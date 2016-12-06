#!/usr/bin/env python

"""
Realign the whole database and output CSV to stdout
"""

import argparse

import bcr_models as igm


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('chain', help='Realign only this chain', nargs='+')
    parser.add_argument('--csv', help='CSV file to parse', default='built-in')
    args = parser.parse_args()

    if args.csv == 'built-in':
        tdb = igm.db.BuiltinTemplateDatabase()
    else:
        with open(args.csv) as f:
            tdb = igm.db.TemplateCSVDatabase(f)

    for chain in args.chain:
        hmm = igm.db.builtin_hmm(chain)
        tdb.realign(hmm=hmm, chain=chain)
    print('pdbname,pdbchain,chain,cs1,cs2,cs3,sequence,antigen')
    tdb.dump_csv()
