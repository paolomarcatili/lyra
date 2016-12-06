
"""
Database Structures
===================

Ig Databases
------------

All immunoglobolin databases must implement the interface specified by
:class:`TemplateDatabase`.

.. autoclass:: BuiltinTemplateDatabase
    :members:
    :undoc-members:
    :member-order: bysource


.. autoclass:: TemplateCSVDatabase
    :members:
    :undoc-members:
    :member-order: bysource


.. autoclass:: TemplateDatabase
    :members:
    :undoc-members:
    :member-order: bysource


Canonical Structure databases
-----------------------------

Similarly to ig databases, all cs databases must implement the interface
specified by the base class :class:`CsDatabase`.

.. autoclass:: CsDatabase
    :members:
    :undoc-members:
    :member-order: bysource


.. autoclass:: BuiltinCsDatabase
    :members:
    :undoc-members:
    :member-order: bysource


PDB file databases
------------------

.. autoclass:: PDBDirectoryDatabase
    :members:
    :undoc-members:
    :member-order: bysource


.. autoclass:: BuiltinPDBDatabase
    :members:
    :undoc-members:
    :member-order: bysource

"""

import os
import sys
import csv
import glob
import collections

import Bio.PDB

import bcr_models as bcr
from bcr_models.canonical_structures import CsDatabase, CsJSONDatabase, BuiltinCsDatabase


class TemplateDatabase(object):
    """
    Database interface for immunoglobolin sequences.
    """

    def search(self, **search):
        for entry in self._db.values():
            if all([getattr(entry, key) == val for key,val in search.items()]):
                yield entry

    def __getitem__(self, item):
        return self._db[item]

    def __contains__(self, item):
        return item in self._db

    def __iter__(self):
        return iter(self._db.values())

    def realign(self, hmm, chain):
        """Realign the whole database. """
        seqs = []
        for entry in self.search(chain=chain):
            seqs.extend(['>' + entry.pdbname, entry.sequence.replace('-', '')])

        cmd = ['hmmalign', '--trim', hmm, '-']
        aln = bcr.utils.run_cmd(cmd, '\n'.join(seqs))

        aligned_templates = []

        for line in aln.splitlines():
            if not line.strip() or line[0] in ('#', '/'):
                continue

            pdbname, alnseq = line.split()
            self[pdbname, chain].sequence = alnseq

    def dump_csv(self, outfile=None):
        """Dump the database to a file or stdout."""
        if not outfile:
            outfile = sys.stdout

        writer = csv.DictWriter(outfile, self._fields)
        for row in self._db.values():
            writer.writerow(row.asdict())


class TemplateEntry(object):
    """
    Container for a template entry

    Usage::

        >>> sequence = ('DIQMTQSPSSLSASVGDRVTITCRTSQSIS----------SYLNWYQQKPGKAPKLLIYA------'
        ... '--ASSLQSGVPSRFSGSGSG--------TDFTLTISSLQPEDFATYYCQQSYSAP--------RTFGQGTKVEIKRTV')

        >>> entry = TemplateEntry('1DEE', 'L', 'K', '2', '1', '1', sequence, '0')
        >>> entry.cs1, entry.K2, entry['K3']
        (2, 1, 1)
        >>> entry.antigen
        False
        >>> entry
        TemplateEntry(1DEE, K)

    Align::

        >>> import bcr_models as igm
        >>> from bcr_models.tests import IgChain_kappa
        >>> aligned_entry = entry.align(IgChain_kappa())
        >>> aligned_entry
        TemplateEntry(1DEE, K, 80.0%)
        >>> aligned_entry.loop_score
        12

    Blacklist certain loops::

        >>> kappa = IgChain_kappa()
        >>> aligned_entry = entry.align(kappa, blacklist=('K2', ))
        >>> aligned_entry.loop_score
        0

    Get an aligned PDB chain::

        >>> aligned_entry.aligned_seq #doctest: +ELLIPSIS
        'DIQMTQSPSSLSASVGDRVTITCRTSQS----------ISSYLNWYQQKPGK...'
        >>> pdbchain = aligned_entry.get_aligned_pdb(kappa)
        >>> pdbchain[0]
        <Residue ASP het=  resseq=1 icode= >
        >>> igm.utils.chain2seq(pdbchain) #doctest: +ELLIPSIS
        'DIQMTQSPSSLSASVGDRVTITCRTSQS----------ISSYLNWYQQKPGK...'

    With custom sequence::

        >>> aligned_entry.aligned_seq = sequence
        >>> aligned_entry.realign = False
        >>> pdbchain = aligned_entry.get_aligned_pdb(kappa)
        >>> igm.utils.chain2seq(pdbchain) #doctest: +ELLIPSIS
        'DIQMTQSPSSLSASVGDRVTITCRTSQSIS----------SYLNWYQQKPGK...'

    """
    def __init__(self, pdbname, pdbchain, chain, cs1, cs2, cs3, sequence, antigen):
        self.pdbname  = pdbname
        self.pdbchain = pdbchain
        self.chain    = chain
        self.cs1      = int(cs1)
        self.cs2      = int(cs2)
        self.cs3      = int(cs3)
        self.sequence = sequence
        self.antigen  = bool(int(antigen))
        for i in '1', '2', '3':
            setattr(self, self.chain + i, getattr(self, 'cs' + i))
            setattr(self, self.chain + i + '_score', None)

        #Alignment details
        self.aligned_seq    = None
        self.score          = None
        self.blosum_score   = None
        self.identity_score = None
        self.loop_score     = None

        #Realign when getting pdbchain
        self.realign = True

    def __getitem__(self, item):
        """Get an item.

        Raises KeyError::

            >>> entry = TemplateEntry('15C8', 'L', 'K', '6', '1', '1', None, '0')
            >>> entry['nothing']
            Traceback (most recent call last):
                ...
            KeyError: 'nothing'

        """
        try:
            return getattr(self, item)
        except AttributeError:
            raise KeyError(item)

    def __setitem__(self, item, value):
        setattr(self, item, value)

    def copy(self):
        return TemplateEntry(self.pdbname, self.pdbchain, self.chain, self.cs1, self.cs2,
            self.cs3, self.sequence, self.antigen)

    def asdict(self):
        return {'pdbname': self.pdbname, 'pdbchain': self.pdbchain, 'chain': self.chain,
                'cs1': self.cs1, 'cs2': self.cs2, 'cs3': self.cs3,
                'sequence': self.sequence, 'antigen': int(self.antigen)}

    def __repr__(self):
        if self.score:
            return 'TemplateEntry({}, {}, {:.1%})'.format(self.pdbname, self.chain, self.identity_score)

        return 'TemplateEntry({}, {})'.format(self.pdbname, self.chain)

    def __str__(self):
        return self.pdbname

    def csv_line(self):
        return '{},{},{},{:>2},{:>2},{:>2},{},{}'.format(*[self.pdbname,
            self.pdbchain, self.chain, str(self.cs1), str(self.cs2), str(self.cs3),
            self.sequence, int(self.antigen)])

        return ','.join([self.pdbname, self.pdbchain, self.chain, str(self.cs1),
            str(self.cs2), str(self.cs3), self.sequence, str(int(self.antigen))])

    def align(self, target, blacklist):
        """Align/score to target."""
        entry = self.copy()

        #Correct gaps in entry sequence
        entry.aligned_seq = target._gapalign(entry.sequence)

        #Calculate BLOSUM score and %ID score
        entry.score = bcr.utils.pairwise_score(target.aligned_seq, entry.aligned_seq)
        entry.blosum_score, entry.identity_score = entry.score

        #Calculate a loop score based on matching loops
        entry.loop_score = 0
        #Calculate score for each cdr loop
        for csid, csregion in target.cdr_regions.items():
            #Extract loop sequences
            target_cdrseq = target.aligned_seq[slice(*csregion)]
            entry_cdrseq = entry.aligned_seq[slice(*csregion)]

            #Skip if loops are not same length
            if target_cdrseq.count('-') != entry_cdrseq.count('-'):
                entry[csid + '_score'] = None
                continue

            #Score
            entry[csid + '_score'] = bcr.utils.pairwise_score(target_cdrseq, entry_cdrseq)

            #Everytime a loop matches, add the blosum score to the loop score
            if entry[csid] == target.cs[csid] and csid not in blacklist:
                entry.loop_score += entry[csid + '_score'][0]

        return entry

    def get_aligned_pdb(self, target, force_realign=None):
        """Get a list of PDB residues aligned to target (ig_chain).

        :param force_realign: If False, force not realigning, and if True, for
                              realigning. If None, use self.realign as guide.

        """
        #Template entry
        entry = target.template_db[self.pdbname, target.chain_type]
        #Bio.PDB.Chain
        pdb_chain = target.pdb_db.get(self.pdbname)[entry.pdbchain]

        realign = force_realign
        if realign is None:
            realign = self.realign

        aligned_pdbseq = self.aligned_seq
        if not aligned_pdbseq or realign:
            aligned_pdbseq = None

        return target.get_aligned_pdb(pdb_chain, aligned_pdbseq=aligned_pdbseq)


class TemplateCSVDatabase(TemplateDatabase):
    """
    Read an immunoglobolin database from a csv file.

    :param filehandle csvhandle: Open filehandle for database csv file.

    """

    def __init__(self, csvhandle):
        try:
            self._filename = csvhandle.name
        except AttributeError:
            self._filename = 'Unknown'

        self._parse_csv(csvhandle)

    def _parse_csv(self, csvhandle):
        """Parse the csv file."""
        reader = csv.DictReader(csvhandle)
        self._fields = reader.fieldnames
        #Generator to fill in the database
        genr = (((e['pdbname'], e['chain']), TemplateEntry(**e)) for e in reader)
        self._db = collections.OrderedDict(genr)

    def __str__(self):
        return 'TemplateCSVDatabase({})'.format(self._filename)


class BuiltinTemplateDatabase(TemplateCSVDatabase):
  """
  Built in templates.

  """

  def __init__(self):
    self._filename = os.path.join(os.path.dirname(__file__), 'data', 'templates.csv')

    with open(self._filename) as f:
        self._parse_csv(f)


#
# PDB file databases
# ------------------

class PDBDatabase(object):
    """Get PDB files."""


class PDBDirectoryDatabase(PDBDatabase):
    """Get PDB files stored in a directory."""

    def __init__(self, pdbdir):
        self._pdbdir = pdbdir
        self._cache  = {}

    def get(self, pdbname):
        """Return PDB model with pdbname."""
        pdbname = pdbname.upper()

        if pdbname not in self._cache:
            local_path = os.path.join(self._pdbdir, '{}.pdb'.format(pdbname))
            if not os.path.isfile(local_path):
                raise bcr.BCRDatabaseError('Unable to locate {}.'.format(local_path))

            #Parse the PDB
            parser = Bio.PDB.PDBParser(QUIET=True)
            self._cache[pdbname] = parser.get_structure(pdbname, local_path)[0]

        return self._cache[pdbname]

    def __str__(self):
        return 'PDBDirectoryDatabase({})'.format(self._pdbdir)

    def filename(self, pdbname):
        """Get the filename for pdbname."""
        return os.path.join(self._pdbdir, '{}.pdb'.format(pdbname))

    def __iter__(self):
        for filename in os.listdir(self._pdbdir):
            pdbname, ext = os.path.splitext(filename)
            if ext.lower() == '.pdb':
                yield pdbname


class BuiltinPDBDatabase(PDBDirectoryDatabase):
    """Built in PDB files. """

    def __init__(self):
        pdbdir = os.path.join(os.path.dirname(__file__), 'data', 'pdb')
        super(BuiltinPDBDatabase, self).__init__(pdbdir)


#
# HMM "databases"
# ---------------

def builtin_hmms():
    """Get a list of HMM files. """
    return glob.glob(os.path.join(os.path.dirname(__file__), 'data', '*.hmm'))

def builtin_hmm(chain_type):
    """Get the HMM profile that fits chain_type. """
    hmm = {
        'K': 'kappa.hmm',
        'L': 'lambda.hmm',
        'H': 'heavy.hmm',
        'A': 'alpha.hmm',
        'B': 'beta.hmm',
    }.get(chain_type.upper())

    return os.path.join(os.path.dirname(__file__), 'data', hmm)