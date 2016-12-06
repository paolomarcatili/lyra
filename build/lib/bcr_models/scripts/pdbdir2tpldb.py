#!/usr/bin/env python

"""
Load a directory of PDB files and output a template database.

"""

from __future__ import print_function

import os
import argparse

import bcr_models as igm


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('pdbdir', help='PDB directory to load.')
    parser.add_argument('-s', '--skip-header', help='Do not print header.', action='store_true')
    args = parser.parse_args()

    tpllines = {}

    csdb = igm.db.BuiltinCsDatabase()
    pdbs = igm.db.PDBDirectoryDatabase(args.pdbdir)
    for pdbname in pdbs:
        pdbmodel = pdbs.get(pdbname)
        for pdbchain in pdbmodel:
            ig_chain = igm.IgChain.from_pdbchain(pdbchain)
            ig_chain.canonical_structures(csdb)

            tt = ig_chain.chain_type
            if tt not in tpllines:
                tpllines[tt] = []

            tplline = [pdbname, tt, tt]
            tplline.extend([cstyp for cdr, cstyp in sorted(ig_chain.cs.items())])
            tplline.extend([ig_chain.aligned_seq, 0])

            tpllines[tt].append(tplline)

    if not args.skip_header:
        print('pdbname,pdbchain,chain,cs1,cs2,cs3,sequence,antigen')
    for tt, lines in sorted(tpllines.items()):
        for line in lines:
            print('{},{},{},{:>2},{:>2},{:>2},{},{}'.format(*line))

if __name__ == '__main__':
    main()