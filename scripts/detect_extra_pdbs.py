#!/usr/bin/env python

"""
Detect PDB files in a dir not present in templates.
"""

from __future__ import print_function

import argparse

import bcr_models as igm


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--missing', help='Find missing PDB files instead.', action='store_true')
    args = parser.parse_args()

    tpls = set([e.pdbname for e in igm.db.BuiltinTemplateDatabase()])
    pdbs = igm.db.BuiltinPDBDatabase()
    pdbset = set(pdbs)

    if args.missing:
        for pdbname in sorted(tpls - pdbset):
            print(pdbname)
    else:
        for pdbname in sorted(pdbset - tpls):
            print(pdbs.filename(pdbname))