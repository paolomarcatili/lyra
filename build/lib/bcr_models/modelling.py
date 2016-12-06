# -*- coding: utf-8 -*-
"""
Classes to model ig superfamily molecules.

Ig Complexes
------------

.. autosummary::

   IgComplex
   IgComplex.from_pdb
   IgComplex.hmmsearch
   IgComplex.canonical_structures
   IgComplex.find_templates
   IgComplex.build_structure
   IgComplex.save
   IgComplex.renumber_pdbmodel
   IgComplex.benchmark

How to instantiate:

.. autoclass:: IgComplex
.. automethod:: IgComplex.from_pdb


Modelling Complexes
^^^^^^^^^^^^^^^^^^^

.. automethod:: IgComplex.hmmsearch
.. automethod:: IgComplex.canonical_structures
.. automethod:: IgComplex.find_templates
.. automethod:: IgComplex.build_structure
.. automethod:: IgComplex.save


Convenience Methods
^^^^^^^^^^^^^^^^^^^

.. automethod:: IgComplex.renumber_pdbmodel
.. automethod:: IgComplex.benchmark


Ig Chains
---------

.. autosummary::

    IgChain
    IgChain.from_template
    IgChain.hmmalign
    IgChain.hmmsearch
    IgChain.canonical_structures
    IgChain.find_templates
    IgChain.build_structure
    IgChain.save
    IgChain.show_alignment
    IgChain.show_templates
    IgChain._parse_hmm
    IgChain.remodel_sidechains
    IgChain.set_framework
    IgChain.set_loop
    IgChain.renum
    IgChain.denum
    IgChain.get_aligned_pdbname
    IgChain.get_aligned_pdb
    IgChain.renumber_pdbchain
    IgChain.benchmark

Attribute summary:

How to instantiate:

.. autoclass:: IgChain
.. automethod:: IgChain.from_template


Modelling with IgChain
^^^^^^^^^^^^^^^^^^^^^^

.. automethod:: IgChain.hmmalign
.. automethod:: IgChain.hmmsearch
.. automethod:: IgChain.canonical_structures
.. automethod:: IgChain.find_templates
.. automethod:: IgChain.build_structure
.. automethod:: IgChain.save


Displaying IgChains
^^^^^^^^^^^^^^^^^^^

.. automethod:: IgChain.show_alignment
.. automethod:: IgChain.show_templates


Manual Modelling with IgChain
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automethod:: IgChain._parse_hmm
.. automethod:: IgChain.remodel_sidechains
.. automethod:: IgChain.set_framework
.. automethod:: IgChain.set_loop
.. automethod:: IgChain.renum
.. automethod:: IgChain.denum


Convenience Methods
^^^^^^^^^^^^^^^^^^^

.. automethod:: IgChain.get_aligned_pdbname
.. automethod:: IgChain.get_aligned_pdb
.. automethod:: IgChain.renumber_pdbchain
.. automethod:: IgChain.benchmark


"""

from __future__ import division
from __future__ import print_function

import os
import copy
import logging
import itertools
import re

import Bio.PDB
import numpy as np

import bcr_models as igm
import bcr_models.utils

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO


#Define a log
log = logging.getLogger('LYRA')
log.addHandler(logging.NullHandler())


class IgComplex(object):
    """
    Main class to model a complex of two Ig-like structures.
    """

    def __init__(self, seq1, seq2, template_db=None, pdb_db=None):
        self._chains = [
            IgChain(seq1, template_db, pdb_db) if seq1 else None,
            IgChain(seq2, template_db, pdb_db) if seq2 else None,
        ]

        self.template_db = template_db or igm.db.BuiltinTemplateDatabase()
        self.blacklist   = []
        self.pdbmodel    = None

    @classmethod
    def detect(cls, seq1, seq2, template_db=None, pdb_db=None, hmms=()):
        """Initialize IgComplex and detect Protein/Nucleic acid sequences.

        Autoconverts to the best reading frame::

            >>> from bcr_models.tests import NUCLEO_ALPHA, NUCLEO_BETA
            >>> NUCLEO_BETA #doctest: +ELLIPSIS
            'ATGGGCACCAGGCTCCTCTGCTGGGCAGCCCTGTGCCTCCTGGG...'
            >>> ig_complex = IgComplex.detect(NUCLEO_ALPHA, NUCLEO_BETA)
            >>> ig_complex['A'].aligned_seq #doctest: +ELLIPSIS
            '-QQGEEDP-QALSIQEGENATMNCSYKTS--------INNLQWYRQ...'
            >>> ig_complex['B'].aligned_seq #doctest: +ELLIPSIS
            '-AGVSQTPSNKVTEKGKYVELRCDPIS-------GHTALYWYRQSLGQGP...'

        Peptides too::

            >>> from bcr_models.tests import EX_HEAVY, EX_KAPPA
            >>> ig_complex = IgComplex.detect(EX_HEAVY, EX_KAPPA)
            >>> ig_complex['H'].aligned_seq #doctest: +ELLIPSIS
            'EVQLVESGPGLVQPGKSLRLSCVASGFT-------FSGYGMHW...'
            >>> ig_complex['K'].aligned_seq #doctest: +ELLIPSIS
            'DIVMTQTPSTLSASVGDRVTLTCKASQD-----------ISYLAWYQQK...'

        """
        chain1 = IgChain.detect(seq1, template_db=template_db, pdb_db=pdb_db, hmms=hmms)
        chain2 = IgChain.detect(seq2, template_db=template_db, pdb_db=pdb_db, hmms=hmms)

        ig_complex = cls(chain1.sequence, chain2.sequence,
            template_db=template_db, pdb_db=pdb_db)

        hmms = hmms or igm.db.builtin_hmms()
        ig_complex.hmmsearch(*hmms)

        return ig_complex

    @classmethod
    def from_pdbmodel(cls, pdbmodel, template_db=None, pdb_db=None, hmms=()):
        """Load an IgChain from a Bio.PDB.Model."""
        if len(pdbmodel) >= 2:
            err = 'Expected two (or one) chains in Bio.PDB.Model. (Found {})'
            raise igm.BCRParserError(err.format(len(pdbmodel)))

        ig_complex = cls(None, None, template_db=template_db, pdb_db=pdb_db)
        ig_complex.pdbmodel = Bio.PDB.Model.Model('model')

        for i, pdbchain in enumerate(pdbmodel):
            ig_complex._chains[i] = IgChain.from_pdbchain(pdbchain,
                template_db=template_db, pdb_db=pdb_db, hmms=hmms)
            ig_complex.pdbmodel.add(ig_complex._chains[i].remodelled_chain)

        return ig_complex

    @classmethod
    def from_dirtypdb(cls, pdbmodel, template_db=None, pdb_db=None, hmms=(),
                      l_cutoff=100, score_cutoff=50, interface_cutoff=14.):
        """Load any PDB file and find a complex if possible."""
        hmms         = hmms or igm.db.builtin_hmms()
        l_cutoff     = l_cutoff or float('inf')
        score_cutoff = score_cutoff or float('inf')
        chains       = []
        for pdbchain in pdbmodel:
            try:
                ig_chain = IgChain.from_pdbchain(pdbchain)
            except igm.BCRParserError:
                continue

            if len(ig_chain.aligned_seq.replace('-', '')) < l_cutoff:
                continue

            if ig_chain._hmmsearch_score < score_cutoff:
                continue

            chains.append(ig_chain)

        legal_pairs = {
            frozenset({'A', 'B'}),
            frozenset({'L', 'H'}),
            frozenset({'K', 'H'}),
        }

        #Pick first chain in list and find suitable candidates
        candidates = [c for c in chains
            if frozenset({chains[0].chain_type, c.chain_type}) in legal_pairs]
        chains = [chains[0], ]

        for candidate in candidates:
            if chains[0].min_interface_dist(candidate) < interface_cutoff:
                chains.append(candidate)
                break

        if len(chains) == 2:
            chain_types = set(c.chain_type for c in chains)
            if frozenset(chain_types) not in legal_pairs:
                err = 'Unkown pairing: {!r}'
                raise igm.BCRParserError(err.format(chain_types))

        ig_complex = cls(None, None, template_db=template_db, pdb_db=pdb_db)
        ig_complex.pdbmodel = Bio.PDB.Model.Model('model')
        for i, ig_chain in enumerate(chains):
            ig_complex._chains[i] = ig_chain
            ig_complex.pdbmodel.add(ig_complex._chains[i].remodelled_chain)

        return ig_complex

    def __iter__(self):
        return iter(self._chains)

    def __getitem__(self, q):
        """Return chain based either on integer or chain type query."""
        if isinstance(q, int):
            return self._chains[q]

        for chain in self:
            if chain.chain_type == q:
                return chain

    def hmmsearch(self, *hmms):
        """Find the best alignment from the hmms and align the sequence.

        :param list hmms: Locations of the hmmer profiles to align.
                          When no profiles are supplied, the built-in profiles
                          will be used.

        """
        for chain in self._chains:
            chain.hmmsearch(*hmms)

    def canonical_structures(self, csdb):
        """Calculate canonical structures.

        :param csdb: An instance of a canonical structures database.
        :type  csdb: CsDatabase

        """
        for chain in self._chains:
            chain.canonical_structures(csdb)

    def find_templates(self, single_model=True, force_graft=(), single_id_min=0.8,
                       id_max=None, local_id_max=None, blacklist=(), method='weighted',
                       framework_weight=1):
        """Automatic selection of templates.

        :param int n:        The number of templates to load.
        :param single_model: Prioritise chain pairs from same model.
        :type  single_model: bool
        :param force_graft:  List of canonical structures to force grafting.
        :type  force_graft:  list
        :param float id_min: Never include templates with id lower than this.
        :param float id_max: Discard templates with %ID above this value

        Use ``force_graft`` to enter a list of cdr ids to graft.

        """
        if force_graft:
            raise NotImplementedError

        id_max = id_max if id_max is not None else 1.01
        #local_id_max will follow id_max unless explicitly set.
        local_id_max = local_id_max if local_id_max is not None else id_max

        pairs = []

        #Force blacklist upper-case
        blacklist = [pdbname.upper() for pdbname in blacklist]
        self.blacklist = blacklist

        if single_model:
            aligned_templates = {}
            for chain in self:
                #Calc canonical structures
                if not chain.cs:
                    chain.canonical_structures()

                #Get a list of possible templates
                templates = [t for t in chain.aligned_templates() if
                    #Exclude based on identity and blacklists
                    id_max >= t.identity_score >= single_id_min
                    and t.pdbname not in blacklist]

                #Pick the best based on loop_score
                aligned_templates[chain.chain_type] = {t.pdbname: t for t in templates}

            try:
                primary_templates, secondary_templates = aligned_templates.values()
            except ValueError:
                log.error('Detected two identical chain types in complex: {}.'
                    ''.format(chain.description))
                log.error('Please submit compatible chains.')
                raise igm.BCRModellingError('Single chain complex: {}.'
                    ''.format(chain.description))

            for entry1 in primary_templates.values():
                #There is only one framework template; skip.
                if entry1.pdbname not in secondary_templates:
                    continue

                entry2   = secondary_templates[entry1.pdbname]
                loop_sum = entry1.loop_score     + entry2.loop_score
                bl_sum   = entry1.blosum_score   + entry2.blosum_score

                if method == 'loop_score':
                    pair = (loop_sum, bl_sum, entry1.pdbname)
                elif method == 'weighted':
                    pair = (loop_sum + framework_weight * bl_sum, entry1.pdbname)
                elif method == 'framework':
                    pair = (bl_sum, loop_sum, entry1.pdbname)
                pairs.append(pair)

        for chain in self:
            if pairs:
                fmw_template = aligned_templates[chain.chain_type][max(pairs)[-1]]
                chain.templates['framework'] = fmw_template
                chain._select_cs_templates(id_max=local_id_max, blacklist=blacklist)
            else:
                chain.find_templates(force_graft=force_graft, id_max=id_max,
                    local_id_max=local_id_max, blacklist=blacklist, method=method,
                    framework_weight=framework_weight)

        if not pairs:
            self._select_packing_template(blacklist=blacklist)

        return self

    def _select_packing_template(self, blacklist=()):
        """Select the template used for superimposing/packing two chains."""
        #Make substitution matrix
        submat = igm.utils.BLOSUM62

        scored = {}
        #Search on one of the chains
        for entry0 in self.template_db.search(chain=self[0].chain_type):
            #We only need entries that have both chains
            try:
                entry1 = self.template_db[(entry0.pdbname, self[1].chain_type)]
            except KeyError:
                continue

            score = 0
            #Score first chain
            for chain, entry, in zip(self, (entry0, entry1)):
                if entry.pdbname in blacklist:
                    continue

                for i in chain.interface_residues:
                    #Resi is plural of res
                    resi = frozenset((chain.aligned_seq[i], entry.sequence[i]))
                    if '-' in resi:
                        continue
                    score += submat[resi]

            scored[entry0.pdbname] = score, entry0.pdbname

        try:
            best_packing = max(scored, key=lambda s: scored[s])
        except ValueError:
            raise igm.BCRModellingError('Chains {} and {} are incompatible.'
                ''.format(self[0].description, self[1].description))

        for chain in self:
            if chain.templates['framework'].pdbname != best_packing:
                #Get the id for the packing template
                template_id = (scored[best_packing][1], chain.chain_type)
                #Set the template
                chain.templates['packing'] = self.template_db[template_id]

    def build_structure(self, anchors=2, strip_cdr3_sidechains=True):
        """Build the 3D structure.

        :param int anchors: How many residues to use for superimposing.

        """
        built_chains = []
        for chain in self:
            chain.build_structure(anchors, strip_cdr3_sidechains)
            built_chains.append(chain.remodelled_chain)

        #Create a temporary PDB model for Scwrl4
        self.pdbmodel = Bio.PDB.Model.Model('BCRMODEL')
        for new_chain in built_chains:
            self.pdbmodel.add(new_chain)

        return self

    def remodel_sidechains(self, method='different'):
        #Run Scwrl, replace template sidechains with target sequence sidechains.
        try:
            self.pdbmodel = igm.utils.run_scwrl(self.pdbmodel, method)
        except igm.utils.SCWRLError as e:
            log.error('Scwrl error: ' + str(e))
            for chain in self:
                log.info('ALN: ' + chain.aligned_seq)
                log.info('TPL: ' + chain.template_seq)
            raise

        return self

    def save(self, outfile, with_remarks=True):
        """Save the modelled structure."""
        if with_remarks:
            self._save_remarks(outfile)

        igm.utils.save_pdbmodel(self.pdbmodel, outfile)

    def _save_remarks(self, outfile):
        """Save the modelled structure."""
        remarks = [
            'REMARK     ',
            'REMARK     Output from LYRA: Lymphocyte Receptor Automated Modelling',
            'REMARK     Made at Center for Biological Sequence Analysis (CBS),',
            'REMARK     Technical University of Denmark (DTU)',
            'REMARK     Version: {}'.format(igm.__version__),
        ]

        newrmk = lambda r: remarks.append('REMARK     {}'.format(r)[:80])

        for chain in self:
            newrmk('')
            newrmk('CHAIN {} {}'.format(chain.chain_type, chain.description))
            newrmk('FRAMEWORK TEMPLATE {} {}'.format(chain.chain_type,
                chain.templates.get('framework')))

            #Chain template and aligned sequence
            alnseq = chain.aligned_seq
            tplseq = chain.template_seq
            stkseq = ['|' if a==t and a!='-' else ' ' for a, t in zip(alnseq, tplseq)]
            stkseq = ''.join(stkseq)

            n = 63
            alnseqs = [alnseq[i:i+n] for i in range(0, len(alnseq), n)]
            stkseqs = [stkseq[i:i+n] for i in range(0, len(stkseq), n)]
            tplseqs = [tplseq[i:i+n] for i in range(0, len(tplseq), n)]

            for alnseq, stkseq, tplseq in zip(alnseqs, stkseqs, tplseqs):
                newrmk('TPL {} {}'.format(chain.chain_type, tplseq))
                newrmk('      {}' .format(stkseq))
                newrmk('INP {} {}'.format(chain.chain_type, alnseq))

            for cdr, (start, end) in sorted(chain.cdr_regions.items()):
                newrmk('')
                newrmk('LOOP {} ABSPOS {}-{}'.format(cdr, start, end))
                newrmk('LOOP {} PDBPOS {}-{}'.format(cdr,
                    chain.denum(start), chain.denum(end)))
                newrmk('LOOP {} TEMPLATE {}'.format(cdr, chain.templates.get(cdr)))
                newrmk('LOOP {} TEMPLATE SEQUENCE: {}'.format(cdr,
                    chain.template_seq[start:end].replace('-', '')))
                newrmk('LOOP {} INPUT SEQUENCE:    {}'.format(cdr,
                    chain.aligned_seq[start:end].replace('-', '')))

            newrmk('PACKING TEMPLATE {} {}'.format(chain.chain_type,
                chain.templates.get('packing')))

            for pdbname in self.blacklist:
                newrmk('BLACKLIST {}'.format(pdbname))

            print('\n'.join(remarks), file=outfile)



    def save_pymol(self, outfile):
        """Save the pretty image pymol script."""

        #Make PDB string ready for PyMol
        pdbfile = StringIO()
        self.save(pdbfile)
        pdblines = pdbfile.getvalue().splitlines()
        pdblines = '\\\n'.join(pdblines)
        pdbfile.close()

        #Make all commands for PyMol script
        commands = [
            'cmd.read_pdbstr("""{}""", "LYRA_model")'.format(pdblines),
        ]

        colors = (
            ('0x909090', ('0xFF6633', '0xFF9933', '0xCC0000')),
            ('0x505050', ('0x3399CC', '0x0066FF', '0x9933FF'))
        )

        zit = zip(sorted(self, key=lambda c: c.chain_type), colors)
        for chain, (fmw_color, loop_colors) in zit:
            commands.append('color {}, ///{}'.format(fmw_color, chain.chain_type))

            local_zit = zip(sorted(chain.cdr_regions.items()), loop_colors)
            for (cdr, (start, end)), loop_color in local_zit:
                commands.append('color {}, ///{}/{}-{}'.format(loop_color,
                    chain.chain_type, chain.denum(start), chain.denum(end)))

        commands.append('as cartoon')
        commands.append('set opaque_background, off')
        commands.append('bg_color white')
        commands.append('set ray_trace_mode, 1')

        print('\n'.join(commands), file=outfile)

    def benchmark(self, ref_pdb, outfmt='rmsd', ref_chains=None, excluding=None):
        """Superimpose ref_pdb (in-place) to the built structure framework.

        :param Bio.PDB.Model ref_pdb: Reference PDB model.
        :param dict ref_chain:        Dictionary of reference chains **EXPLAIN**
        :returns: Dictionary of rmsds

        Chain type from PDB must match chain type from IgChains.

        The rmsd is given as a sum of squared distances and number of atoms. To
        calculate the final RMSD use: sqrt(sq_sum / n)

        """
        #Superimpose PDB to model
        atoms = {
            'Total':        ([], []),
            'Framework':    ([], []),
            'Binding site': ([], []),
        }

        if excluding:
            exclname = 'Binding site excl. ' + excluding
            atoms[exclname] = ([], [])

        score = {}
        for chain in self:
            #Get the right reference chain
            ref_chain = chain.chain_type
            if ref_chains:
                ref_chain = ref_chains.get(chain.chain_type, chain.chain_type)

            ref_chain = ref_pdb[ref_chain]

            #Get aligned backbone atoms
            chain_atoms = chain._aligned_atoms(ref_chain, excluding=excluding)
            #Extend list of all atoms
            for region, atom_lists in chain_atoms.items():
                loop = region
                if region in chain.cdr_regions:
                    region = 'Binding site'
                elif region not in ('Framework', 'Total'):
                    continue

                atoms[region][0].extend(atom_lists[0])
                atoms[region][1].extend(atom_lists[1])
                if excluding and region == 'Binding site' and loop != excluding:
                    atoms[exclname][0].extend(atom_lists[0])
                    atoms[exclname][1].extend(atom_lists[1])

            score.update(igm.utils.score_model(chain_atoms, outfmt).items())

        score.update(igm.utils.score_model(atoms, outfmt))
        return score

    def __repr__(self):
        vrs = ', '.join([c.chain_type if c else 'X' for c in self])
        return 'IgComplex({})'.format(vrs)

    @classmethod
    def renumber_pdbmodel(cls, pdbmodel, hmms=None):
        """Renumber a Bio.PDB.Model."""
        new_model = Bio.PDB.Model.Model('')
        for chain in pdbmodel:
            new_chain = IgChain.renumber_pdbchain(chain, hmms)
            new_model.add(new_chain)

        return new_model


class IgChain(object):
    """
    Main class to model a single immunoglobulin-like structure.

    :var sequence:           The original sequence.
    :var template_db:        The template database instance.
    :var pdb_db:             The PDB template database instance.

    :var cs:                 Dictionary of canonical structures.

    :var interface_residues: Residues used to select packing template.
    :var packing_residues:   Residues used to superimpose chain on packing template.

    """

    def __init__(self, sequence, template_db=None, pdb_db=None):
        # -- Set by init()
        self.sequence    = str(sequence)
        self.template_db = template_db or igm.db.BuiltinTemplateDatabase()
        self.pdb_db      = pdb_db      or igm.db.BuiltinPDBDatabase()

        #Set by canonical_structures()
        self.cs                 = {}

        #Set by hmmalign()
        self.aligned_seq        = None
        self.leftaln_seq        = None

        #set by hmmsearch()
        self._hmm               = None
        self._hmmsearch_score   = None
        self.E_value            = None

        #Set by find_templates()
        self.templates          = {}
        self._aligned_templates = None

        #Set by _parse_hmm()
        self.cdr_regions        = {}
        self.chain_type         = None
        self.description        = ''
        self.interface_residues = None
        self.packing_residues   = None
        self.template_seq       = ''
        self._renum_table       = None
        self._numindex          = None
        self._insertion_sites   = []
        self.segi               = None

        #Set by build_chain()
        self.remodelled_chain   = None
        self.aligned_chain      = None
        self.loop_graft_rmsd    = {}

        #Set by get_aligned_pdb()
        self._pdb_templates     = {}

    @classmethod
    def detect(cls, seq, template_db=None, pdb_db=None, hmms=()):
        """Initialize IgComplex and detect Protein/Nucleic acid sequences.

        Autoconverts to the best reading frame.

        Usage::

            >>> from bcr_models.tests import NUCLEO_ALPHA
            >>> NUCLEO_ALPHA #doctest: +ELLIPSIS
            'atggaaactctcctgggagtgtctttggtgattctatggcttca...'
            >>> ig_chain = IgChain.detect(NUCLEO_ALPHA)
            >>> ig_chain.description
            'Alpha chain'
            >>> ig_chain.sequence #doctest: +ELLIPSIS
            'METLLGVSLVILWLQLARVNSQQGEEDPQALSIQEGENATM...'
            >>> ig_chain.aligned_seq #doctest: +ELLIPSIS
            '-QQGEEDP-QALSIQEGENATMNCSYKTS--------INNLQW...'

        """
        hmms = hmms or igm.db.builtin_hmms()

        if re.match(r'^[ATGCU]+$', seq, flags=re.IGNORECASE):
            seq = seq.lower().replace('u', 't')
            chains = {}
            #Make six different IgChains
            for reading_frame in (-1, -2, -3, 1, 2, 3):
                tmp_seq = igm.utils.translate_nucleotides(seq, reading_frame, stop=False)
                ig_chain = cls(tmp_seq, template_db=template_db, pdb_db=pdb_db)
                #Score
                try:
                    ig_chain.hmmsearch(*hmms)
                except igm.BCRParserError:
                    pass
                else:
                    chains[reading_frame] = ig_chain

            try:
                #Pick out the best chain based on HMM score.
                reading_frame, best_chain = max(chains.items(),
                    key=lambda c: c[1]._hmmsearch_score)
            except ValueError:
                raise igm.BCRParserError('Unable to translate "{}" into '
                    'meaningful peptide sequence.'.format(seq))

            log.info('Detected Nucleic Acid sequence for chain {} in reading frame {}'
                ''.format(best_chain.chain_type, reading_frame))

            return best_chain

        elif re.match(r'^\w+$', seq):
            ig_chain = cls(seq, template_db=template_db, pdb_db=pdb_db)
            ig_chain.hmmsearch(*hmms)
            return ig_chain

        raise BCRParserError('Sequence "{}" is invalid'.format())

    @classmethod
    def from_template(cls, pdbname, chain_type, hmm=None, template_db=None, pdb_db=None):
        """Create an IgChain from a PDB Chain."""
        hmm         = hmm         or igm.db.builtin_hmm(chain_type)
        template_db = template_db or igm.db.BuiltinTemplateDatabase()
        pdb_db      = pdb_db      or igm.db.BuiltinPDBDatabase()

        template = template_db[pdbname, chain_type]
        pdbchain = pdb_db.get(pdbname)[template.pdbchain]

        ig = cls(template.sequence.replace('-', ''), template_db, pdb_db)
        ig.hmmalign(hmm)
        ig.cs = {c: template[c] for c in ig.cdr_regions}
        ig.remodelled_chain = pdbchain

        return ig

    @classmethod
    def from_pdbchain(cls, pdbchain, template_db=None, pdb_db=None, hmms=None):
        """Load an IgChain from a Bio.PDB.Chain. """
        hmms     = hmms or igm.db.builtin_hmms()
        sequence = igm.utils.chain2seq(pdbchain)

        ig_chain = cls(sequence, template_db=template_db, pdb_db=pdb_db)
        ig_chain.hmmsearch(*hmms)

        ig_chain.aligned_chain    = []
        ig_chain.remodelled_chain = Bio.PDB.Chain.Chain(ig_chain.chain_type)
        for i, res in enumerate(ig_chain.get_aligned_pdb(pdbchain)):
            if res:
                res    = res.copy()
                res.id = igm.utils.resid2biopdb(ig_chain.denum(i))

                #For some reason this hack is needed
                for at in res:
                    at.disordered_flag = 0

                #Add to the
                ig_chain.remodelled_chain.add(res)

            #Add to aligned chain
            ig_chain.aligned_chain.append(res)

        return ig_chain

    def hmmalign(self, hmm):
        """Align a single HMM to the sequence.

        :param str hmm: Location of hmmer profile to align.

        Usage::

            >>> from bcr_models.tests import KAPPA_HMM, EX_KAPPA
            >>> KAPPA_HMM # doctest: +ELLIPSIS
            '.../kappa.hmm'
            >>> ig = IgChain(EX_KAPPA)
            >>> ig.hmmalign(KAPPA_HMM)
            >>> ig.description
            'Kappa Light chain'
            >>> ig.chain_type
            'K'
            >>> sorted(ig.cdr_regions.items())
            [('K1', [25, 42]), ('K2', [59, 70]), ('K3', [116, 130])]
            >>> ig._numindex #doctest: +ELLIPSIS
            ['1', '2', '3',...'29', '30', '30A', '30B',...'30J', '31',...'95H',...'109']

        """
        self.aligned_seq = igm.utils.hmmalign(hmm, self.sequence, trim=True)

        try:
            self._parse_hmm(hmm)
        except (ValueError, AttributeError) as e:
            err = 'Could not parse DESC record in {}'
            raise igm.BCRParserError(err.format(hmm))

        if not self.aligned_seq.isupper():
            err = '{}|{}'.format(self.chain_type, self.aligned_seq)
            raise igm.BCRInsertionsInAlignmentError(err)

    def _parse_hmm(self, hmm):
        """Parse a hmmer profile extra information.

        :param str hmm: Path to the hmm file.

        """
        self._hmm = hmm
        #Locate the description in the HMM file
        desc     = None
        numlines = None
        intfc    = ''
        pckg     = ''
        with open(hmm) as f:
            for line in f:
                if line[0:4] == 'DESC':
                    desc = line.split(None, 1)[1].strip()
                elif line[0:4] == 'INSS':
                    self._create_renum(line.split(None, 1)[1].strip())
                elif line[0:5] == 'INTFC':
                    intfc, pckg = line.split()[1].split('|')

        self.interface_residues = [self.renum(i) for i in intfc.split(',')]

        #Finish parsing packing residues
        pckg = [i.split('-') for i in pckg.split(',')]
        self.packing_residues = []
        for start, end in pckg:
            rg = (self.renum(start), self.renum(end) + 1)
            self.packing_residues.extend(list(range(*rg)))

        #Parse description record
        self.description, self.chain_type, cdrs = desc.split('|')
        #Parse complementary determining regions
        self.cdr_regions = {}
        for cdr in cdrs.split(';'):
            cdr_id, reg = cdr.split(':')
            self.cdr_regions[cdr_id] = [self.renum(s) for s in reg.split('-', 1)]
            self.cdr_regions[cdr_id][1] += 1

        #Re-align gaps
        self.aligned_seq = self._gapalign(self.aligned_seq)
        self.leftaln_seq = self._gapalign(self.aligned_seq, align='l')
        self._create_segi()

        return self

    def _gapalign(self, seq, insertions=None, align='c'):
        """Redo alignment to center gaps in loop areas.

        Usage::

            >>> from bcr_models.tests import IgChain_heavy
            >>> ig_chain = IgChain_heavy()
            >>> ig_chain._gapalign(ig_chain.aligned_seq, align='c') #doctest: +ELLIPSIS
            'EVQL...KVKFYDP------------------TAPNDYWGQGTLVTVSS'
            >>> ig_chain._gapalign(ig_chain.aligned_seq, align='l') #doctest: +ELLIPSIS
            'EVQL...KVKFYDPTAPND------------------YWGQGTLVTVSS'

        With residues outside alignment::

            >>> extended_seq = 'outside' + ig_chain.aligned_seq
            >>> ig_chain._gapalign(extended_seq, align='c') #doctest: +ELLIPSIS
            'outsideEVQL...KVKFYDP------------------TAPNDYWGQGTLVTVSS'
            >>> ig_chain._gapalign(extended_seq, align='l') #doctest: +ELLIPSIS
            'outsideEVQL...KVKFYDPTAPND------------------YWGQGTLVTVSS'

        """
        offset = 0
        while seq[offset] != '-' and seq[offset].islower():
            offset += 1

        new_seq = seq[offset:]
        if insertions is None:
            insertions = sorted(self.cdr_regions.values())

        for s, e in insertions:
            new_seq = igm.utils.align_insertions(new_seq, s, e, align)

        new_seq = seq[:offset] + new_seq

        assert seq.replace('-', '') == new_seq.replace('-', '')
        return new_seq

    def hmmsearch(self, *hmms):
        """Find the best alignment from the hmms and align the sequence.

        :param list hmms: Locations of the hmmer profiles to align.

        """
        #No point in searching
        if len(hmms) == 1:
            return self.hmmalign(hmms[0])
        #Fetch the default hmm profiles
        elif not hmms:
            raise igm.BCRModellingError('No hmm profiles given to hmmsearch')

        self._hmmsearch_scores = dict()
        for hmm in hmms:
            #Run hmmsearch
            aln_raw = igm.utils.run_cmd([igm.utils.HMMSEARCH_BIN, hmm, '-'],
                '>actseq\n' + self.sequence)
            aln = [line.split() for line in aln_raw.splitlines()]
            #Parse
            score = None
            for i, line in enumerate(aln):
                if line[0:3] == ['E-value', 'score', 'bias'] and aln[i+2]:
                    try:
                        E_value = float(aln[i+2][0])
                        score = float(aln[i+2][1])
                        break
                    except ValueError:
                        E_value = float(aln[i+3][0])
                        score = float(aln[i+3][1])
                        break
            #Error checking
            if score is not None:
                #Register score and E_value
                self._hmmsearch_scores[hmm] = score, E_value

        try:
            #Pick out the winner
            best_hmm = max(self._hmmsearch_scores.items(), key=lambda x: x[1][0])
            self._hmm, (self._hmmsearch_score, self.E_value) = best_hmm
            #Do the final alignment
            self.hmmalign(self._hmm)
        except ValueError:
            raise igm.BCRParserError('Unable to align sequence "{}". '
                'Make sure it is a TCR or BCR sequence.'.format(self.sequence))

    def canonical_structures(self, csdb):
        """Calculate canonical structures.

        :param csdb: An instance of a canonical structures database.
        :type  csdb: CsDatabase

        Usage::

            >>> from pprint import pprint
            >>> from bcr_models.tests import IgChain_kappa
            >>> from bcr_models.tests import TestCsDatabase

            >>> chain = IgChain_kappa()
            >>> chain.canonical_structures(TestCsDatabase())
            >>> pprint(chain.cs)
            {'K1': 1, 'K2': 1, 'K3': 0}

        """
        self.cs = csdb.get(self)

    def aligned_templates(self, realign=False, blacklist=('H3', )):
        """Score all templates.

        :param float id_min: Identity cutoff. Do not allow templates with a
                             total %ID lower than this.

        Return a list of dictionaries of all possible entries with the
        total_score calculated.

        >>> 'Show the dict' # doctest: +SKIP
        'RESULT'

        """
        if not self.aligned_seq:
            raise igm.BCRRuntimeError('Align sequence before selecting templates')

        if self._aligned_templates and not realign:
            return self._aligned_templates

        self._aligned_templates = []
        for entry in self.template_db.search(chain=self.chain_type):
            aligned_entry = entry.align(self, blacklist=blacklist)
            self._aligned_templates.append(aligned_entry)

        return self._aligned_templates

    def find_templates(self, force_graft=(), id_max=None, local_id_max=None,
                       blacklist=(), method='weighted', framework_weight=1):
        """Automatic selection of templates.

        :param force_graft:  List of canonical structures to force grafting.
        :type  force_graft:  list
        :param float id_min: Never include templates with id lower than this.

        Use ``force_graft`` to enter a list of cdr ids to graft.

        """
        if force_graft:
            raise NotImplementedError

        id_max = id_max if id_max is not None else 1.01
        #local_id_max will follow id_max unless explicitly set.
        local_id_max = local_id_max if local_id_max is not None else id_max

        #Choose selection criteria
        if method == 'loop_score':
            criteria = lambda t: (t.loop_score, t.score, t.pdbname)
        elif method == 'weighted':
            fmww = framework_weight
            criteria = lambda t: (t.loop_score + fmww * t.blosum_score, t.pdbname)
        elif method == 'framework':
            criteria = lambda t: (t.blosum_score, t.loop_score, t.pdbname)

        self.templates['framework'] = max([t for t in self.aligned_templates() if
            #Exclude based on idettity and blacklists
            id_max >= t.identity_score and t.pdbname not in blacklist],
            key=criteria)

        self._select_cs_templates(id_max=local_id_max, blacklist=blacklist)

        return self

    def _select_cs_templates(self, id_max=None, blacklist=()):
        """Select templates for CS regions.

        :param highscore: Best template (framework template).
        :type  highscore: IgDb entry (dict)
        :param candidate_templates: List of candidate templates to use for grafting.
        :type  candidate_templates: List of IgDb entries
        :param float id_max:        Maximum %ID (total) of templates to accept.
        :param float id_min_escape: Minimum %ID (total) to accept loop templates
                                    if no templates are found using id_min.
                                    Use None to raise exception if no templates
                                    are found above id_min.

        """
        #Choose CS regions
        for csid in self.cdr_regions:
            self._select_cs_template(csid, id_max=id_max, blacklist=blacklist)

    def _select_cs_template(self, csid, id_max=None, blacklist=()):
        """Select the best loop from candidate_templates.

        :param str csid:            Canonical structure ID to find template for. (Eg L1)
        :param candidate_templates: List of candidate templates:
        :type  candidate_templates: List of IgDb entries.
        :param float id_max:        Maximum %ID (total) of templates to accept.
        :param float id_min_escape: Minimum %ID (total) to accept loop templates
                                    if no templates are found using id_min.
                                    Use None to raise exception if no templates
                                    are found above id_min.

        """
        scorekey = csid + '_score'

        if id_max is None:
            id_max = 1.01

        #Matching CS for the best scorer:
        if (self.cs[csid] and self.templates['framework'][scorekey]  #There is a CS (!= 0)
            and self.templates['framework'][csid] == self.cs[csid]   #CS type match framework
            and self.templates['framework'].pdbname not in blacklist #Blacklisted
            and self.templates['framework'][scorekey][1] <= id_max): #%id is not above id_max
                #.. then choose framework as CS template
                self.templates[csid] = self.templates['framework']
                return

        aligned_templates = [t for t in self.aligned_templates()
            if t[scorekey] and t[scorekey][1] <= id_max and t.pdbname not in blacklist]

        best_loops = sorted(aligned_templates, reverse=True,
            key=lambda t: (t[scorekey], t.score, t.pdbname))

        if not best_loops:
            log.error('Could not find loop for {} length={}'.format(csid, self.cdr_len(csid)))
            log.error('There are no templates in the database that match the '
                'length of the input sequence.')
            log.error('Maybe the sequence is incomplete, the HMM failed or '
                'the loop in question is too exotic for us to model.')
            log.error('Alignment: {}'.format(self.aligned_seq))
            log.error('Loops:     {}'.format(self._get_loop_definition()))
            raise igm.BCRNoTemplateError('Completely unable to find matching loop for {}. '
                ''.format(csid), self, csid)

        self.templates[csid] = None
        #No point looking if cs is 0
        if self.cs[csid]:
            #Pick the best loop with the right canonical structure
            for loop in best_loops:
                #This means mismatch in length is reached
                if not loop[scorekey]:
                    break
                #Matching canonical structure
                if loop[csid] == self.cs[csid]:
                    self.templates[csid] = loop
                    break

        #No matching canonical structures
        if not self.templates[csid]:
            #Pick the best
            self.templates[csid] = best_loops[0]

        assert self.templates[csid].chain == self.chain_type

    def _pack_residuelist(self):
        """Superimpose residuelist onto the packing template structure. """
        moving_reslist = self.aligned_chain

        static = []
        moving = []
        #Get pdb file and align to a residuelist
        static_template = self.templates['packing']
        static_pdbchain = self.template_db[static_template.pdbname, self.chain_type].pdbchain
        static_chain   = self.pdb_db.get(static_template.pdbname)[static_pdbchain]
        static_reslist = self.get_aligned_pdb(static_chain)
        #Extract superimposing atoms.
        for i in self.packing_residues:
            s_res = static_reslist[i]
            m_res = moving_reslist[i]
            if s_res is None or m_res is None:
                continue

            static.append(static_reslist[i]['CA'].get_coord())
            moving.append(moving_reslist[i]['CA'].get_coord())

        #Run the superimposer and aquire rotation matrix and translation vector
        rotran = igm.utils.svd_superimpose(np.array(static), np.array(moving))

        for res in moving_reslist:
            if res is not None:
                res.transform(*rotran)

        self.remodelled_chain = moving_reslist

    def build_structure(self, anchors=2, strip_cdr3_sidechains=True):
        """Build the 3D structure.

        :param int anchor_len: How many residues to use for superimposing.

        """
        if not self.templates.get('framework'):
            raise igm.BCRRuntimeError('No framework structure found.')

        framework = self.templates['framework']
        self.set_framework()
        for loop in self.cdr_regions:
            entry = self.templates.get(loop)

            #Set all loops where there is a template that is not the same as the framework
            if entry and entry.pdbname != framework.pdbname:
                self.set_loop(loop, anchors=anchors)

        if strip_cdr3_sidechains:
            self._strip_cdr_sidechains(self.chain_type + '3')

        if self.templates.get('packing'):
            self._pack_residuelist()

        #Realign the current chain to be more conservative on how to re-align insertions
        realigned_chain = [r for r in self.aligned_chain if r]
        realigned_chain = self.get_aligned_pdb(realigned_chain,
            insertions=self._insertion_sites, align='l')

        self.remodelled_chain = Bio.PDB.Chain.Chain(self.chain_type)
        for i, res in enumerate(realigned_chain):
            if res:
                #Get a Bio.PDB compatible id to renumber the current residue
                resid = igm.utils.resid2biopdb(self.denum(i))
                res.id = resid
                #For some reason this hack is needed
                for at in res:
                    at.disordered_flag = 0
                self.remodelled_chain.add(res)

        return self

    def remodel_sidechains(self, method='different'):
        """Remodel sidechains."""
        #Create a temporary PDB model for Scwrl4
        new_model = Bio.PDB.Model.Model('BCRMODEL')
        new_model.add(new_chain)

        #Run Scwrl, replace template sidechains with target sequence sidechains.
        remodelled_model = igm.utils.run_scwrl(new_model, method)

        self.remodelled_chain = list(remodelled_model)[0]
        self.aligned_chain    = self.get_aligned_pdb(self.remodelled_chain)

    def set_framework(self, pdbname=None):
        """Load a PDB, align it and set is as the current structure.

        :param str pdbname: pdbname to use for framework. Leave as `None` to use
                            the internal template dict.

        Usage::

            >>> from bcr_models.tests import IgChain_kappa
            >>> ig_chain = IgChain_kappa()
            >>> ig_chain.set_framework('1DEE')
            >>> bcr_models.utils.chain2seq(ig_chain.aligned_chain) #doctest: +ELLIPSIS
            'DIVMTQTPSTLSASVGDRVTLTCKASQD-----...CLQQN---------SNWTFGQGTKVDIK---'

        """
        if pdbname:
            framework_chain = self.get_aligned_pdbname(pdbname)
        else:
            framework_chain = self.templates['framework'].get_aligned_pdb(self)
        self.template_seq = igm.utils.chain2seq(framework_chain)

        self.aligned_chain = []
        for r, a in zip(framework_chain, self.aligned_seq):
            r = r if a != '-' else None
            if r:
                if a == igm.utils.pl3to1[r.resname]:
                    r = r.copy()
                    r.detach_parent()
                    r.conserve = True
                else:
                    r = igm.utils.convert_residue(r, a)
                    r.conserve = False

            self.aligned_chain.append(r)

    def set_loop(self, loop, anchors=2, pdbname=None):
        """Load a PDB, align it and graft it onto the framework.

        :param str loop:    Which loop to graft.
        :param int anchors: How many anchor residues to use.
        :param str pdbname: pdbname to use for loop. Leave as `None` to use
                            the internal template dict.

        First a framework must be present::

            >>> from bcr_models.tests import IgChain_kappa
            >>> chain2seq = bcr_models.utils.chain2seq
            >>> ig_chain = IgChain_kappa()
            >>> ig_chain.set_framework('1DEE')
            >>> out = [0, 0, 0]
            >>> out[0] = chain2seq(ig_chain.aligned_chain)
            >>> out[1] = ig_chain.template_seq
            >>> out[2] = chain2seq(ig_chain.aligned_chain, scwrl=True)

        First sequence is the aligned chain converted to the input sequence,
        second sequence is the original template sequence, and the last sequence
        is the input sequence with conserved residues in lower-case::

            >>> print('\\n'.join(out)) #doctest: +ELLIPSIS
            DIVMTQTPSTLSASVGDRVTLTCKASQD-----------ISYLAWYQQKPGKAP...
            DIQMTQSPSSLSASVGDRVTITCRTSQS----------ISSYLNWYQQKPGKAP...
            diVmtqTpsTlsasvgdrvtLtcKAsqD-----------IsylAwyqqkpgkap...

        Then set the loop::

            >>> ig_chain.set_loop('K1', pdbname='1EO8')
            >>> chain2seq(ig_chain.aligned_chain) #doctest: +ELLIPSIS
            'DIVMTQTPSTLSASVGDRVTLTCKASQD-----------ISYLAWYQQKPGKAP...'
            >>> round(ig_chain.loop_graft_rmsd['K1'], 4)
            0.2903
            >>> ig_chain.template_seq #doctest: +ELLIPSIS
            'DIQMTQSPSSLSASVGDRVTITCRTSSD-----------ISYLNWYQQKPGKAP...'
            >>> chain2seq(ig_chain.aligned_chain, scwrl=True) #doctest: +ELLIPSIS
            'diVmtqTpsTlsasvgdrvtLtcKAsQd-----------isylAwyqqkpgkap...'

        """
        pdbname = pdbname or self.templates[loop].pdbname
        aligned_loop = self.get_aligned_pdbname(pdbname)

        start, end = self.cdr_regions[loop]

        #Superimpose anchors
        static = []
        moving = []
        for r in range(start - anchors, start), range(end, end + anchors):
            for i in r:
                for atname in ('CA', 'C', 'N', 'O'):
                    try:
                        static_at = self.aligned_chain[i][atname].get_coord()
                        moving_at = aligned_loop[i][atname].get_coord()
                    except (KeyError, TypeError):
                        pass
                    else:
                        static.append(static_at)
                        moving.append(moving_at)

        static = np.array(static)
        moving = np.array(moving)
        rot, tran = igm.utils.svd_superimpose(static, moving)

        template_seq = list(self.template_seq)

        #Transform and graft residues
        for i in range(start, end):
            res = aligned_loop[i] or None

            #Set template sequence
            template_seq[i] = igm.utils.pl3to1[res.resname] if res else '-'

            if res:
                if self.aligned_seq[i] == igm.utils.pl3to1[res.resname]:
                    res = aligned_loop[i].copy()
                    res.detach_parent()
                    res.conserve = True
                else:
                    res = igm.utils.convert_residue(res, self.aligned_seq[i])
                    res.conserve = False

                #Transform residue
                res.transform(rot, tran)

            self.aligned_chain[i] = res

        rmsd = igm.utils.rmsd(static, np.dot(moving, rot) + tran)
        self.loop_graft_rmsd[loop] = rmsd

        self.template_seq = ''.join(template_seq)

    def _strip_cdr_sidechains(self, loop):
        """Strips sidechains of given cdr region and remodel them with SCWRL4.

        First a framework must be present::

            >>> from bcr_models.tests import IgChain_kappa
            >>> chain2seq = bcr_models.utils.chain2seq
            >>> ig_chain = IgChain_kappa()
            >>> ig_chain.set_framework('1DEE')

        Comparison of before and after the stripping and remodelling::

            >>> chain2seq(ig_chain.aligned_chain, scwrl=True) #doctest: +ELLIPSIS
            'diVmtqTpsTlsasvgdrvtLtcKAsqD-----------IsylAwyqq...'
            >>> ig_chain._strip_cdr_sidechains('K1')
            >>> chain2seq(ig_chain.aligned_chain, scwrl=True) #doctest: +ELLIPSIS
            'diVmtqTpsTlsasvgdrvtLtcKASQD-----------ISYlAwyqq...'

        """
        start, end = self.cdr_regions[loop]

        for i in range(start,end):
            res = self.aligned_chain[i]

            if res:
                res = igm.utils.convert_residue(res, self.aligned_seq[i])
                res.conserve = False
                self.aligned_chain[i] = res

    def get_aligned_pdbname(self, pdbname):
        """Get a list of residues aligned to the current hmm.

        Usage::

            >>> from bcr_models.tests import IgChain_kappa
            >>> reslist = IgChain_kappa().get_aligned_pdbname('1DEE')
            >>> bcr_models.utils.chain2seq(reslist) #doctest: +ELLIPSIS
            'DIQMTQSPSSLSASVGDRVTITCRTSQS-----...CQQSYS--------APRTFGQGTKVEIK---'

        `pdbname` must be present in the template database.

        """
        pdbname = str(pdbname.upper())

        if pdbname not in self._pdb_templates:
            #Get the template entry
            entry = self.template_db[pdbname, self.chain_type]
            #Bio.PDB.Chain
            pdb_chain = self.pdb_db.get(pdbname)[entry.pdbchain]
            #Align
            self._pdb_templates[pdbname] = self.get_aligned_pdb(pdb_chain)

        return self._pdb_templates[pdbname]

    def get_aligned_pdb(self, pdb_chain, insertions=None, align='c', aligned_pdbseq=None):
        """Fetch, align and return a renumbered pdb."""
        #Three to one residue ids
        pdbseq = igm.utils.chain2seq(pdb_chain)

        if not aligned_pdbseq:
            #Align the sequence to the profile
            aligned_pdbseq = igm.utils.hmmalign(self._hmm, pdbseq, trim=False)
            aligned_pdbseq = self._gapalign(aligned_pdbseq,
                insertions=insertions, align=align)

        aligned_chain = []
        pdbchain_iter = iter(list(pdb_chain))

        for resid in aligned_pdbseq:
            #Gap
            if resid == '-':
                aligned_chain.append(None)
            #Outside profile TODO: more explicit?)
            elif resid.lower() == resid:
                next(pdbchain_iter)
                continue
            #Actual residue in profile
            else:
                try:
                    newres = next(pdbchain_iter)
                    #Hetatm
                    if newres.resname not in igm.utils.pl3to1:
                        newres = None
                #Catch exhaustion of the pdb iterator
                except StopIteration:
                    newres = None
                #Add the residue to the aligned list
                aligned_chain.append(newres)

        if len(self.aligned_seq) > len(aligned_chain):
            diff = (len(self.aligned_seq) - len(aligned_chain))
            aligned_chain.extend([None] * diff)
        elif not len(self.aligned_seq) == len(aligned_chain):
            raise igm.BCRRuntimeError('Could not align pdb')

        #Trim from the left
        i = 0
        while self.aligned_seq[i] == '-':
            aligned_chain[i] = None
            i += 1
        #Trim from the right
        j = -1
        while self.aligned_seq[j] == '-':
            aligned_chain[j] = None
            j -= 1

        return aligned_chain

    def save(self, outfile):
        """Save the modelled structure."""
        new_model = Bio.PDB.Model.Model('BCRMODEL')
        new_model.add(self.remodelled_chain)

        bcr_models.utils.save_pdbmodel(new_model, outfile)

    def _get_aln_sticks(self, seq1=None, seq2=None):
        """Get the alignment sticks."""
        seq1 = (seq1 or self.aligned_seq).upper()
        seq2 = (seq2 or self.template_seq or '').upper()
        stx = []
        for a, b in zip(seq1, seq2):
            if a == b and a != '-':
                stx.append('|')
            elif a != b and (a == '-' or b == '-'):
                stx.append('X')
            else:
                stx.append(' ')
        return ''.join(stx)

    def _get_loop_definition(self, use_unicode=False):
        """Get a string of loop definitions."""
        marker_start, marker_end, marker_mdl, = '|', '|', '-'
        if use_unicode:
            marker_start, marker_end, marker_mdl, = u'└', u'┘', u'─'

        aln = [' '] * len(self.aligned_seq)

        #Go through each region
        for cdr, (s, e) in sorted(self.cdr_regions.items(), key=lambda t: -t[1][0]):
            #Mark endpoints
            aln[s] = marker_start
            aln[e-1] = marker_end
            #..and everything in between
            for i in range(s+1, e-1):
                aln[i] = marker_mdl
            #Calculate an offset for the label (cdr id) and fill positions
            o = s + (e - s) //2 - len(cdr)//2
            aln[o - 1] = ' '
            aln[o + len(cdr)] = ' '
            for i, c in enumerate(cdr):
                aln[o + i] = c

        return ''.join(aln)

    def _get_insertions_line(self):
        """Get all the insertion letters in a line, as well as number sticks."""
        letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        ins = [' ' for _ in self.aligned_seq]
        num = [' ' for _ in self.aligned_seq]

        for i in range(len(self.aligned_seq)):
            if self.denum(i)[-1] in letters:
                ins[i] = self.denum(i)[-1]
            elif (i + 1) % 10 == 0 or i == 0:
                ins[i] = '|'
                j = self.denum(i)
                num[i:i+len(j)] = j

        return ''.join(num), ''.join(ins)


    def show_alignment(self, html=False, use_unicode=False, numbering=False,
                       terminal=False, realign=False):
        """Return the alignment as a string.

        Usage::

            >>> from bcr_models.tests import IgChain_kappa
            >>> ig_chain = IgChain_kappa()
            >>> print(ig_chain.show_alignment()) #doctest: +ELLIPSIS
            DIVMTQTPSTLSASVGDRVTLTCKASQD-----------ISYLAWYQQKPG...
                                     |----- K1 ------|         ...

        With added templates::

            >>> ig_chain.set_framework('1DEE')
            >>> ig_chain.set_loop('K1', pdbname='1EO8')
            >>> ig_chain.set_loop('K3', pdbname='1DQL')
            >>> print(ig_chain.show_alignment()) #doctest: +ELLIPSIS
            Template: DIQMTQSPSSLSASVGDRVTITCRTSSD-----------ISYLNWYQQKPGKAPKLL...
                      || ||| || |||||||||| ||  | |           |||| ||||||||||| |...
            Input:    DIVMTQTPSTLSASVGDRVTLTCKASQD-----------ISYLAWYQQKPGKAPKKL...
                                               |----- K1 ------|               ...

        """
        aln = self.aligned_seq
        tpl = self.template_seq
        if realign:
            aln = igm.utils.hmmalign(self._hmm, aln.replace('-', ''), trim=True)
            aln = self._gapalign(aln, insertions=self._insertion_sites, align='l')
            if tpl:
                tpl = igm.utils.hmmalign(self._hmm, tpl.replace('-', ''), trim=True)
                tpl = self._gapalign(tpl, insertions=self._insertion_sites, align='l')

        offset = 0
        if self.template_seq:
            offset = 10
            out = [
                'Template: {}'.format(tpl),
                '          {}'.format(self._get_aln_sticks(tpl, aln)),
                'Input:    {}'.format(aln),
               u'          {}'.format(self._get_loop_definition(use_unicode)),
            ]
        else:
            out = [
                self.aligned_seq,
                self._get_loop_definition(use_unicode)
            ]

        if html or terminal:
            if html:
                begin = '<span class="lyra-aln-loop-{}">'
                end   = '</span>'
            else:
                begin = '\033[32m'
                end   = '\033[0m'
            for cdr, (s, e) in sorted(self.cdr_regions.items(), key=lambda t: -t[1][0]):
                for i, line in enumerate(out):
                    line = list(line)
                    line.insert(e + offset, end)
                    line.insert(s + offset, begin.format(cdr))
                    out[i] = ''.join(line)
            if html and self.template_seq:
                out[1] = '<span class="lyra-aln-sticks">{}</span>'.format(out[1])

        if numbering:
            num, ins = self._get_insertions_line()
            out.insert(0, ' ' * offset + '{}'.format(ins))
            out.insert(0, ' ' * offset + '{}'.format(num))

        return '\n'.join(out)

    def show_templates(self):
        """Show the templates that are selected."""
        out = []
        for region, template in sorted(self.templates.items()):
            if region == 'packing':
                identity == ''
            elif region == 'framework':
                identity = ' ({:.0%})'.format(template.identity_score)
            else:
                try:
                    identity = ' ({:.0%})'.format(template[region + '_score'][1])
                except KeyError:
                    idettity = ''

            out.append('{}: {}{}'.format(region, template.pdbname, identity))

        return ', '.join(out)


    def renum(self, q):
        """Get the python sequence index from a residue number."""
        if not self._renum_table:
            return int(q) - 1
        return self._renum_table.get(q, None)

    def denum(self, i):
        """Get a residue code from a sequence number."""
        if not self._numindex:
            return str(i + 1)
        return self._numindex[i]

    def cdr_len(self, cdr):
        """Get the length of a CDR.

        Simple::

            >>> from bcr_models.tests import IgChain_kappa
            >>> ig_chain = IgChain_kappa()
            >>> print(ig_chain.show_alignment()) #doctest: +ELLIPSIS
            DIVMTQTPSTLSASVGDRVTLTCKASQD-----------ISYLAWYQQKP...
                                     |----- K1 ------|        ...
            >>> ig_chain.cdr_len('K1')
            6

        """
        cdr_seq = self.aligned_seq[slice(*self.cdr_regions[cdr])]
        return len(cdr_seq) - cdr_seq.count('-')

    def min_interface_dist(self, other_ig_chain):
        """Min dist between two CA atoms in interface_residues."""
        m = float('inf')
        for i in self.interface_residues:
            res0 = self.aligned_chain[i]
            if not res0 or 'CA' not in res0:
                continue
            atom0 = res0['CA']

            for j in other_ig_chain.interface_residues:
                res1 = other_ig_chain.aligned_chain[j]
                if not res1 or 'CA' not in res1:
                    continue
                atom1 = res1['CA']

                d = np.sqrt(np.sum((atom0.get_coord() - atom1.get_coord())**2))
                if d < m:
                    m = d
        return m

    def _create_renum(self, renum):
        """Create the renum table."""
        self._insertion_sites = []
        for site in sorted(renum.split(',')):
            sp = site.split('-')
            self._insertion_sites.append((int(sp[0]) - 1, int(sp[1])))

        self._insertion_sites.sort()

        if not self.aligned_seq:
            raise Exception('Align sequence before parsing HMM')

        inscodes = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        self._numindex = []
        i = n = 0
        for j in range(len(self.aligned_seq)):
            for start, end in self._insertion_sites:
                if start <= j < end:
                    self._numindex.append(str(i) + inscodes[n])
                    n += 1
                    break
            else:
                i += 1
                n = 0
                self._numindex.append(str(i))

        self._renum_table = {n: i for i, n in enumerate(self._numindex)}

    def _create_segi(self):
        """Create segment list.

        Called in _parse_hmm::

            >>> from bcr_models.tests import IgChain_kappa
            >>> ig = IgChain_kappa()
            >>> ig.segi #doctest: +ELLIPSIS
            ['F1', 'F1', ... 'K1', ... 'F2', ... 'K2', ... 'F3', ... 'K3', ... 'F4', 'F4']
            >>> len([s for s in ig.segi if s[0] == 'F'])
            102

        """
        self.segi = []
        fmw   = 1 #Current framework region number
        #Loop: keep track of which loops have been iterated
        fmw_l = {r: 1 for r in self.cdr_regions}
        for i in range(len(self.aligned_seq)):
            for loop, (r_start, r_end) in self.cdr_regions.items():
                if r_start <= i < r_end:
                    self.segi.append(loop)
                    #Increase framework region number ig loop is hitherto unseen
                    fmw += fmw_l[loop]
                    fmw_l[loop] -= fmw_l[loop]
                    break
            else:
                self.segi.append('F{}'.format(fmw))

    @property
    def _desc_str(self):
        """Short description compatible with the DESC record hmmer profiles."""
        if not self.description or not self.chain_type or not self.cdr_regions:
            return None
        cdr_regions = sorted(self.cdr_regions.items())
        cdrs = [(a, self.denum(s), self.denum(e-1)) for a, (s, e) in cdr_regions]
        cdrs = ';'.join(['{}:{}-{}'.format(*cdr) for cdr in cdrs])
        return '|'.join([self.description, self.chain_type, cdrs])

    def __str__(self):
        desc = self._desc_str or ''
        cs = ','.join([str(cs) for cdr, cs in sorted(self.cs.items())])
        cs = ', cs={}'.format(cs) if cs else ''
        return 'IgChain({}{})'.format(desc, cs)

    __repr__ = __str__

    def benchmark(self, ref_chain, outfmt='rmsd'):
        """Superimpose ref_pdb (in-place) to the built structure framework.

        :param Bio.PDB.Chain ref_chain: Reference PDB Chain.
        :returns: Dictionary of rmsds

        """
        #ref_chain = ref_chain or self.chain_type
        atoms = self._aligned_atoms(ref_chain)
        return igm.utils.score_model(atoms, outfmt)

    def _aligned_atoms(self, ref_chain, excluding=None,
                       backbone_atoms=('C', 'O', 'N', 'CA')):
        """.."""
        loopname  = '{} loops'.format(self.chain_type)
        total     = 'Total {}'.format(self.chain_type)
        framework = 'Framework {}'.format(self.chain_type)
        atoms = {
            'Total':     ([], []),
            'Framework': ([], []),
            loopname:    ([], []),
            framework: ([], []),
            total: ([], []),
        }

        if excluding and excluding in self.cdr_regions:
            exclname = loopname + ' excl. ' + excluding
            atoms[exclname] = ([], [])
        else:
            excluding = None

        ref_sequence = self.get_aligned_pdb(ref_chain)
        slf_sequence = self.get_aligned_pdb(self.remodelled_chain)

        for i, (r, q) in enumerate(zip(ref_sequence, slf_sequence)):
            if r is None or q is None:
                continue

            moving = []
            static = []

            try:
                for atid in backbone_atoms:
                    moving.append(r[atid].get_coord())
                    static.append(q[atid].get_coord())
            except KeyError:
                #Skip this residue
                continue

            atoms['Total'][0].extend(moving)
            atoms['Total'][1].extend(static)
            atoms[total][0].extend(moving)
            atoms[total][1].extend(static)

            #Check for loops
            for loop, (r_start, r_end) in sorted(self.cdr_regions.items()):
                if r_start <= i < r_end:
                    if loop not in atoms:
                        atoms[loop] = ([], [])
                    #Local loop
                    atoms[loop][0].extend(moving)
                    atoms[loop][1].extend(static)

                    #Collected loop regions
                    atoms[loopname][0].extend(moving)
                    atoms[loopname][1].extend(static)

                    #Loop regions without excluding
                    if excluding and loop != excluding:
                        atoms[exclname][0].extend(moving)
                        atoms[exclname][1].extend(static)

                    break
            #We are in the framework region, append
            else:
                atoms['Framework'][0].extend(moving)
                atoms['Framework'][1].extend(static)
                atoms[framework][0].extend(moving)
                atoms[framework][1].extend(static)

        return atoms

    @classmethod
    def renumber_pdbchain(cls, pdbchain, hmms=None):
        """Renumber a Bio.PDB.Chain."""
        seq  = igm.utils.chain2seq(pdbchain)
        ig   = cls(seq)
        hmms = hmms or igm.db.builtin_hmms()

        ig.hmmsearch(*hmms)
        log.debug('Renumbering using {}'.format(ig._hmm))
        new_chain = Bio.PDB.Chain.Chain(ig.chain_type)
        for i, res in enumerate(ig.get_aligned_pdb(pdbchain)):
            if res:
                res = res.copy()
                resid = igm.utils.resid2biopdb(ig.denum(i))
                res.id = resid

                #For some reason this hack is needed
                for at in res:
                    at.disordered_flag = 0

                new_chain.add(res)

        return new_chain