
"""
Utilities for bcr_models.
"""

from __future__ import division

import os
import math
import tempfile
import subprocess
import collections

import Bio.PDB
import Bio.SubsMat.MatrixInfo
from Bio.SCOP.Raf import protein_letters_3to1 as pl3to1
import numpy as np
import requests

try:
    from string import maketrans
except ImportError:
    maketrans = str.maketrans

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

#Amino acids
AMINO_ACID_THREE = ('ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS',
      'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR')

AMINO_ACID_ONE = 'ACDEFGHIKLMNPQRSTVWY'

pl1to3 = {o: t for o, t in zip(AMINO_ACID_ONE, AMINO_ACID_THREE)}


#External Programs
HMMALIGN_BIN  = os.environ.get('HMMALIGN_BIN',  'hmmalign')
HMMSEARCH_BIN = os.environ.get('HMMSEARCH_BIN', 'hmmsearch')
SCWRL4_BIN    = os.environ.get('SCWRL4_BIN', 'Scwrl4')


#Nucleoic to peptide conversion
TRANSLATION_TABLE = {'ttt': 'F', 'ttc': 'F', 'tta': 'L', 'ttg': 'L', 'ctt': 'L',
                     'ctc': 'L', 'cta': 'L', 'ctg': 'L', 'att': 'I', 'atc': 'I',
                     'ata': 'I', 'atg': 'M', 'gtt': 'V', 'gtc': 'V', 'gta': 'V',
                     'gtg': 'V', 'tct': 'S', 'tcc': 'S', 'tca': 'S', 'tcg': 'S',
                     'cct': 'P', 'ccc': 'P', 'cca': 'P', 'ccg': 'P', 'act': 'T',
                     'acc': 'T', 'aca': 'T', 'acg': 'T', 'gct': 'A', 'gcc': 'A',
                     'gca': 'A', 'gcg': 'A', 'tat': 'Y', 'tac': 'Y', 'taa': 'X',
                     'tag': 'X', 'cat': 'H', 'cac': 'H', 'caa': 'Q', 'cag': 'Q',
                     'aat': 'N', 'aac': 'N', 'aaa': 'K', 'aag': 'K', 'gat': 'D',
                     'gac': 'D', 'gaa': 'E', 'gag': 'E', 'tgt': 'C', 'tgc': 'C',
                     'tga': 'X', 'tgg': 'W', 'cgt': 'R', 'cgc': 'R', 'cga': 'R',
                     'cgg': 'R', 'agt': 'S', 'agc': 'S', 'aga': 'R', 'agg': 'R',
                     'ggt': 'G', 'ggc': 'G', 'gga': 'G', 'ggg': 'G'}

def fetchpdb(pdbname):
    """Fetch a pdbname from online PDB."""
    req = requests.get('http://www.rcsb.org/pdb/files/{}.pdb'.format(pdbname))
    if req.status_code != 200:
        raise Exception('Error fetching PDB: {!r}'.format(pdbname))
    txt = StringIO(req.text)
    return Bio.PDB.PDBParser(QUIET=True).get_structure(pdbname, txt)[0]


def run_cmd(cmd, input_string=''):
    """Run the cmd with input_string as stdin and return output."""
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                         stderr=subprocess.PIPE, universal_newlines=True)
    out, stderr = p.communicate(input=input_string)
    if p.returncode:
        raise Exception('Cmd {} failed: {}'.format(cmd[0], stderr))
    return out


def pairwise_score(seq1, seq2, subs_matrix=None, norm=False):
    """Score two sequences according to matrix.

    Usage::
        >>> seq1 = 'DIQMTQSPSSLS--SVGD--RVTITCR-SQDIR-NDLGWYQQKPGKA'
        >>> seq2 = 'DIVMTQTPSTLS-ASVGD--RVTLACKASQDI--SYLAWYQQKPGKA'
        >>> blosum, ident = pairwise_score(seq1, seq2)
        >>> round(ident, 7)
        0.7209302
        >>> blosum
        167

    Normalized blosum score::

        >>> blosum, ident = pairwise_score(seq1, seq2, norm=True)
        >>> round(blosum, 7)
        0.7767442

    """
    if len(seq1) != len(seq2):
        raise Exception('Sequence lengths not identical.\n1: {}\n2: {}'.format(seq1, seq2))

    if subs_matrix is None:
        subs_matrix = BLOSUM62

    empty = frozenset('-')
    score = score1 = score2 = ident = match = 0
    for pair in zip(seq1, seq2):
        if norm:
            a = frozenset(pair[0])
            b = frozenset(pair[1])
            #Increase blosum reference scores
            score1 += subs_matrix[a] if a != empty else 0
            score2 += subs_matrix[b] if b != empty else 0

        pair = frozenset(pair)
        #Increase Blosum score
        score += subs_matrix[pair]
        #All insertions are not counted towards ident score
        if pair != empty:
            match += 1
            if len(pair) == 1:
                ident += 1

    if norm:
        return score/math.sqrt(score1*score2), ident/match

    return score, ident/match


def frozen_matrix(subs_matrix):
    """Turn a Biopython subsmatrix into a frozenset subs_matrix."""
    new_matrix = collections.defaultdict(int)
    for entry, score in subs_matrix.items():
        new_matrix[frozenset(entry)] = score
    return new_matrix


BLOSUM62 = frozen_matrix(Bio.SubsMat.MatrixInfo.blosum62)


def resid2biopdb(q):
    """Return a Bio.PDB residue id from a string id.

    >>> resid2biopdb('100A')
    (' ', 100, 'A')
    >>> resid2biopdb('100')
    (' ', 100, ' ')

    """
    try:
        return (' ', int(q), ' ')
    except ValueError:
        return (' ', int(q[:-1]), q[-1])


def convert_residue(res, to, resid=None):
    """Convert residue to another and copy backbone atoms."""
    to    = pl1to3[to] if len(to) == 1 else to
    resid = resid or res.id[:]
    r = Bio.PDB.Residue.Residue(resid, to, '')

    for atid in ('N', 'CA', 'C', 'O'):
        if atid in res:
            r.add(res[atid].copy())

    return r


def align_insertions(seq, start, end, align='c'):
    """Align all insertions between start and end.

    >>> align_insertions('ABC-D--EFG--HI', 1, 13, 'c')
    'ABCD-----EFGHI'
    >>> align_insertions('ABC-D--EFG--HI', 1, 12, 'c')
    'ABCD-----EFGHI'
    >>> align_insertions('ABC-D--EFG--HI', 1, 13, 'l')
    'ABCDEFGH-----I'
    >>> align_insertions('ABC-D--EFG--HI', 1, 13, 'r')
    'A-----BCDEFGHI'

    """
    before = seq[:start]
    after  = seq[end:]
    gap    = seq[start:end]
    no_gap = gap.replace('-', '')
    count  = gap.count('-')

    if align.lower() == 'l':
        return ''.join((before, no_gap, '-' * count, after))
    elif align.lower() == 'r':
        return ''.join((before, '-' * count, no_gap, after))
    elif align.lower() == 'c':
        return ''.join((
            before,
            no_gap[:len(no_gap)//2],
            '-' * count,
            no_gap[len(no_gap)//2:],
            after
        ))

    raise Exception('Unknown alignment: '.format(align))


def save_pdbmodel(model, outfile):
    """Save the modelled structure."""
    new_structure = Bio.PDB.Structure.Structure('BCRSTRUCTURE')

    new_structure.add(model)

    io = Bio.PDB.PDBIO()
    io.set_structure(new_structure)
    io.save(outfile)


class SCWRLError(Exception):
    """Raised when a Scwrl error is detected."""


def run_scwrl(model, method='different'):
    """SCWRLS model"""

    sequence = []
    for chain in model:
        for res in chain:
            letter = pl3to1[res.resname].lower()
            if not getattr(res, 'conserve', False):
                letter = letter.upper()
            sequence.append(letter)

    sequence = ''.join(sequence)

    method = str(method).lower()
    if method == 'all':
        sequence = sequence.upper()
    elif method != 'different':
        raise Exception('Unknown Scwrl4 method: "{}" [different|all]'.format(method))

    temp_pdb_in  = tempfile.NamedTemporaryFile(suffix='.pdb')
    temp_pdb_out = tempfile.NamedTemporaryFile(suffix='.pdb')
    temp_seq     = tempfile.NamedTemporaryFile()
    save_pdbmodel(model, temp_pdb_in.name)
    temp_pdb_in.flush()

    temp_seq.write(sequence.encode('ascii'))
    temp_seq.flush()
    out = run_cmd([SCWRL4_BIN, '-i', temp_pdb_in.name, '-o',
                   temp_pdb_out.name, '-s', temp_seq.name])

    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('X', temp_pdb_out.name)

    try:
        return structure[0]
    except KeyError:
        raise SCWRLError('Scwrl (probably) exited with an error:\n{}'.format(out))


def hmmalign(hmm, sequence, trim=True):
    """Align a sequence to an HMM profile."""
    #Run hmmalign
    if trim:
        cmd = [HMMALIGN_BIN, '--trim', hmm, '-']
    else:
        cmd = [HMMALIGN_BIN, hmm, '-']

    aln = run_cmd(cmd, '>actseq\n' + sequence)

    #Parse output
    seq = []
    for line in aln.splitlines():
        line = line.split()
        if line and line[0] == 'actseq':
            seq.append(line[1])

    return ''.join(seq)


def chain2seq(chain, scwrl=False):
    """Convert a Bio.PDB Chain or list of Residues to a string."""
    if not scwrl:
        return ''.join([pl3to1.get(r.resname, '') if r else '-' for r in chain])
    else:
        seq = []
        for res in chain:
            letter = pl3to1.get(res.resname).lower() if res else '-'
            if not getattr(res, 'conserve', False):
                letter = letter.upper()
            seq.append(letter)

        return ''.join(seq)


def score_model(atoms, outfmt):
    """Score a model based on an atom dict.

    :param dict atoms: Dictionary of atom coordinates. See example.
    :param Bio.PDB.Model ref_pdb: Reference PDB model which is moved to the
                                  superimposed location.
    :param str outfmt:            Output format. rmsd, raw or tm-score.

    Examples::

        >>> import numpy as np
        >>> atoms = {'Framework': (np.array([[51.65, -1.90, 50.07],
        ...                                  [50.40, -1.23, 50.65],
        ...                                  [50.68, -0.04, 51.54],
        ...                                  [50.22, -0.02, 52.85]]),
        ...                        np.array([[51.30, -2.99, 46.54],
        ...                                  [51.09, -1.88, 47.58],
        ...                                  [52.36, -1.20, 48.03],
        ...                                  [52.71, -1.18, 49.38]])),
        ...          'loop': (np.array([[51.65, -1.90, 50.07],
        ...                             [50.22, -0.02, 52.85]]),
        ...                   np.array([[51.30, -2.99, 46.54],
        ...                             [52.71, -1.18, 49.38]]))}
        >>> score = score_model(atoms, 'rmsd')
        >>> score['Framework'] #doctest: +ELLIPSIS
        0.0030426...
        >>> score['loop'] #doctest: +ELLIPSIS
        0.0028796...

    """
    outfmt = outfmt.lower()
    #Run the superimposer
    rot, tran = svd_superimpose(np.array(atoms['Framework'][1]),
        np.array(atoms['Framework'][0]))

    score = {}
    for loop, (moving_atoms, static_atoms) in atoms.items():
        moving_atoms = [np.dot(a, rot) + tran for a in moving_atoms]

        if outfmt == 'rmsd':
            score[loop] = rmsd(moving_atoms, static_atoms)
        elif outfmt == 'tm-score':
            score[loop] = tm_score(moving_atoms, static_atoms)
        else:
            raise bcr.BCRRuntimeError(
                'Unknown distance score: {} (raw, rmsd or tm-score)'.format(outfmt))

    return score


def svd_superimpose(static, moving):
    """Superimpose using singular value decomposition.

        :param ndarray static: List of static coordinate vectors.
        :param ndarray moving: List of moving coordinate vectors.

    Usage::

        >>> import numpy as np
        >>> static = np.array([[51.65, -1.90, 50.07],
        ...                    [50.40, -1.23, 50.65],
        ...                    [50.68, -0.04, 51.54],
        ...                    [50.22, -0.02, 52.85]])
        >>> moving = np.array([[51.30, -2.99, 46.54],
        ...                    [51.09, -1.88, 47.58],
        ...                    [52.36, -1.20, 48.03],
        ...                    [52.71, -1.18, 49.38]])
        >>> rotation, translation = svd_superimpose(static, moving)
        >>> rotation
        array([[ 0.68304983,  0.53664371,  0.49543563],
               [-0.52277295,  0.83293229, -0.18147242],
               [-0.51005037, -0.13504564,  0.84947707]])
        >>> translation
        array([ 38.78608157, -20.65451334, -15.42227366])

    """
    #Center on centroid
    av1 = np.average(moving, axis=0)
    av2 = np.average(static, axis=0)
    moving = moving - av1
    static = static - av2

    #Correlation matrix
    a = np.dot(np.transpose(moving), static)
    u, d, vt = np.linalg.svd(a)
    rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))

    #Check if we have found a reflection
    if np.linalg.det(rot) < 0:
        vt[2] = -vt[2]
        rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    tran = av2 - np.dot(av1, rot)

    return rot, tran


def rmsd(a, b):
    """Calculate RMSD between coordinate set a and coordinate set b.

    Usage::

        >>> import numpy as np
        >>> a = np.array([[51.65, -1.90, 50.07],
        ...               [50.40, -1.23, 50.65],
        ...               [50.68, -0.04, 51.54],
        ...               [50.22, -0.02, 52.85]])
        >>> b = np.array([[51.30, -2.99, 46.54],
        ...               [51.09, -1.88, 47.58],
        ...               [52.36, -1.20, 48.03],
        ...               [52.71, -1.18, 49.38]])
        >>> round(rmsd(a, b), 7)
        3.8784565

    """
    a = np.array(a)
    b = np.array(b)
    return np.sqrt(np.sum((a - b)**2) / a.shape[0])


def tm_score(a, b, d0=None):
    """Calc TM-score."""
    L = a.shape[0]
    if d0 is None:
        d0 = L**0.39 * np.sqrt(0.58 + 0.05*L * np.exp(-L/4.7) - 0.63 * np.exp(-L/37))

    return 1/L * np.sum(1 / (1 + (np.sqrt(np.sum((a - b)**2, axis=1)) / d0)**2))


def encode_sparse(seq):
    """Sparse encode a sequence.

    Usage::

        >>> encode_sparse('GQ') #doctest: +ELLIPSIS
        [1, 0, 0, 0, 0, 0, 0, ..., 0, 0, 0, 1, 0, 0, 0, 0, 0]

    """
    amino_acids = ('G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S',
                   'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H',)
    codings = {AA: [1 if i == j else 0 for j in range(len(amino_acids))]
        for i, AA in enumerate(amino_acids)}
    gap = [0] * len(amino_acids)

    coded = []
    for res in seq.upper():
        coded.extend(codings.get(res, gap))

    return coded


def reverse_complement(x):
    """Reverse complement.

    Checks for an internal reverse_complement() method of the supplied object,
    otherwise performs a string translate. String translate retains case.

        >>> s = "AaTtGgCc"
        >>> reverse_complement(s)
        'gGcCaAtT'

    """
    transtable = maketrans('ATGCRYMKBVDHatgcrymkbvdh',
                           'TACGYRKMVBHDtacgyrkmvbhd')
    return str(x).translate(transtable)[::-1]


def translate_nucleotides(sequence, reading_frame=1, stop=True):
    """Translates DNA sequence of specific reading frame to peptide sequence

    Usage::

        >>> from bcr_models import utils
        >>> translate_nucleotides('atggaaactctcctgggagtgtctttggtgatt', -1)
        'NHQRHSQESFH'
        >>> translate_nucleotides('atggaaactctcctgggagtgtctttggtgatt', -2)
        'ITKDTPRRVS'
        >>> translate_nucleotides('atggaaactctcctgggagtgtctttggtgatttagctg')
        'METLLGVSLVI'
        >>> translate_nucleotides('atggaaactctcctgggagtgtctttggtgatttagctg', stop=False)
        'METLLGVSLVIXL'

    """
    sequence = str(sequence).lower()
    peptide_sequence = ''

    if reading_frame < 0:
        sequence = reverse_complement(sequence)
        reading_frame = abs(reading_frame)

    for i in range((reading_frame - 1), len(sequence), 3):
        codon = sequence[i:(i+3)]
        if len(codon) < 3:
            break
        try:
            residue = TRANSLATION_TABLE[codon]
            if stop and residue == 'X':
                break
            peptide_sequence += residue
        except KeyError:
            raise KeyError('Unknown codon: "{}"'.format(codon))

    return peptide_sequence
