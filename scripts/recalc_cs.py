#!/usr/bin/env python

"""
Recalculate canonical structures for all or one
"""

import argparse

import bcr_models as igm

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', help='Template database', default=None)
    parser.add_argument('-c', help='Only calc these chains', default=None, nargs='+')
    args = parser.parse_args()

    if args.t:
        with open(args.t) as f:
            templates = igm.db.TemplateCSVDatabase(f)
    else:
        templates = igm.db.BuiltinTemplateDatabase()

    pdbs = igm.db.BuiltinPDBDatabase()
    csdb = igm.canonical_structures.RandomForestCsDatabase()

    print('pdbname,pdbchain,chain,cs1,cs2,cs3,sequence,antigen')
    for entry in templates:
        if args.c and entry.chain not in args.c:
            print(entry.csv_line())
            continue

        ig_chain = igm.IgChain.from_template(entry.pdbname, entry.chain,
            template_db=templates, pdb_db=pdbs)
        ig_chain.canonical_structures(csdb)

        for cdr, cs in ig_chain.cs.items():
            setattr(entry, 'cs' + cdr[-1], cs)

        #print(entry.cs1)
        print(entry.csv_line())
