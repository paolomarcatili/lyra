
"""
DOC
"""

from __future__ import print_function

import os
import glob
import argparse
import logging

import bcr_models as bcr

log = logging.getLogger('igm/main')
log.addHandler(logging.NullHandler())


#Python < 3.2 compat
try:
   FileNotFoundError
except NameError:
   FileNotFoundError = IOError


def get_fasta(filehandle):
    """Aquire sequences from a fasta."""
    seqid = current_seq = None
    sequences = []

    for line in filehandle:
        if line[0] == '>':
            seqid = line[1:].split()[0]
            if current_seq:
                current_seq['seq'] = ''.join(current_seq['seq']).upper()
                sequences.append(current_seq)
            current_seq = {'seq': [], 'seqid': seqid}
        elif seqid:
            current_seq['seq'].append(line.strip())

    current_seq['seq'] = ''.join(current_seq['seq']).upper()
    sequences.append(current_seq)

    return sequences


def entry():
    parser = argparse.ArgumentParser()
    parser.add_argument('seq', help='Light and heavy chain sequence in fasta file')
    parser.add_argument('-o', '--output', help='Output PDB file', default=None)
    parser.add_argument('--no-renumber', action='store_true',help='Renumber the sequences using internal HMM numbering rather that KC numbering scheme')
    parser.add_argument('-q', '--output-noscwrl', default=None,
        help='Output PDB file without running Scwrl')
    parser.add_argument('-p', '--output-pymol', default=None,
        help='Output PyMol script')
    parser.add_argument('--no-scwrl', action='store_true', help='Do not use SCWRL4')
    parser.add_argument('--no-check', action='store_true', help='Do not perform structural check on the model')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('-a','--noalign', action='store_true',help='Use the alignment provided by the user in the input file')
    parser.add_argument('--debug', action='store_true',
        help='Make an effort to keep things working.')
    parser.add_argument('--align', action='store_true', help='Only align the sequences')
    parser.add_argument('--scwrl-method', default='different', choices=('different', 'all'),
        help='Scwrl4 method to remodel sidochains: different or all')
    parser.add_argument('hmms', nargs='*',
        help='HMM profiles to align. Omit to use built-in profiles')
    args = parser.parse_args()

    if args.debug:
        args.verbose = True

    logfmt = '%(asctime)s %(name)-12s: %(levelname)-8s %(message)s'
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, format=logfmt, datefmt='%Y-%m-%d %H:%M')
    else:
        logging.basicConfig(level=logging.INFO, format=logfmt, datefmt='%Y-%m-%d %H:%M')

    if args.hmms:
        hmms = []
        for hmm in args.hmms:
            hmms.extend(glob.glob(hmm))
        if not all([os.path.isfile(h) for h in hmms]) or not hmms:
            print('Unable to locate HMM profiles.')
            exit(1)
    else:
        hmms = bcr.db.builtin_hmms()

    log.debug('hmmalign: '  + bcr.utils.HMMALIGN_BIN)
    log.debug('hmmsearch: ' + bcr.utils.HMMSEARCH_BIN)
    log.debug('Scwrl4: '    + bcr.utils.SCWRL4_BIN)

    log.debug('hmmer profiles:')
    for hmm in hmms:
        log.debug('-- {}'.format(hmm))

    with open(args.seq) as f:
        seqs = get_fasta(f)

    template_db = bcr.db.BuiltinTemplateDatabase()
    log.debug('Template database:')
    log.debug('-- {}'.format(template_db))

    pdb_db = bcr.db.BuiltinPDBDatabase()
    log.debug('PDB database:')
    log.debug('-- {}'.format(pdb_db))

    csdb = bcr.db.BuiltinCsDatabase()
    if args.debug:
        csdb = bcr.canonical_structures.CsDatabase()
    log.debug('Cs database:')
    log.debug('-- {}'.format(csdb))

    ig_complex = bcr.IgComplex.detect(seqs[0]['seq'], seqs[1]['seq'], template_db, pdb_db)
    ig_complex.hmmsearch(*hmms)
    if args.noalign:
       for ig_chain in ig_complex:
          for seqm in seqs:
             if (ig_chain.aligned_seq.replace("-","") in seqm['seq'].replace("-","")):
                log.debug('HMM aligned sequence:\n{}'.format(ig_chain.aligned_seq))
                log.debug('User provided alignment:\n{}'.format(seqm['seq']))
                ig_chain.aligned_seq=seqm['seq']

    ig_complex.canonical_structures(csdb)
    if not args.align:
        if args.debug:
            for ig_chain in ig_complex:
                log.debug('Realign {} with {}'.format(ig_chain.chain_type, ig_chain._hmm))
                ig_chain.template_db.realign(ig_chain._hmm, ig_chain.chain_type)
        ig_complex.find_templates()
        ig_complex.build_structure()
        if (not args.no_check):
           ig_complex.structure_check()
        if args.output_noscwrl:
            with open(args.output_noscwrl, 'w') as f:
                ig_complex.save(f)

        if not args.no_scwrl:
            ig_complex.remodel_sidechains(args.scwrl_method)

    print('>', seqs[0]['seqid'], ' ', ig_complex[0], sep='')
    print()
    print(ig_complex[0].show_alignment(use_unicode=True, numbering=True,
        terminal=True, realign=(not args.no_renumber)))
    print()
    print(ig_complex[0].show_templates())

    print()

    print('>', seqs[1]['seqid'], ' ', ig_complex[1], sep='')
    print()
    print(ig_complex[1].show_alignment(use_unicode=True, numbering=True,
        terminal=True, realign=(not args.no_renumber)))
    print()
    print(ig_complex[1].show_templates())

    if args.output:
        with open(args.output, 'w') as f:
            ig_complex.save(f)

    if args.output_pymol:
        with open(args.output_pymol, 'w') as f:
            ig_complex.save_pymol(f)


if __name__ == '__main__':
    entry()
