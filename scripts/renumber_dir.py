#!/usr/bin/env python

"""
Renumber an entire directory of PDB files.
"""

import os
import logging
import argparse

import Bio.PDB

import bcr_models as igm


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('inp', help='Input dir')
    parser.add_argument('out', help='Output dir')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    logfmt = '%(asctime)s %(name)-15s: %(levelname)-8s %(message)s'
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, format=logfmt, datefmt='%Y-%m-%d %H:%M')
    else:
        logging.basicConfig(level=logging.INFO, format=logfmt, datefmt='%Y-%m-%d %H:%M')


    for f in os.listdir(args.inp):
        if f[-4:].lower() != '.pdb':
            continue

        print('Processing', f)
        full = os.path.join(args.inp, f)
        mdl = Bio.PDB.PDBParser(QUIET=True).get_structure(f, full)[0]

        new_model = igm.IgComplex.renumber_pdbmodel(mdl)

        igm.utils.save_pdbmodel(new_model, os.path.join(args.out, f))