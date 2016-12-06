#!/usr/bin/env python

"""
Benchmark a directory of PDB files.

"""

from __future__ import print_function
from __future__ import division

import logging
import argparse

import bcr_models as igm
from bcr_models.scripts import benchmarking

#Define a log
log = logging.getLogger('LYRA/benchdir')


def benchmark_pdbnames(pdbnames, pdb_db, id_max=0.98, local_id_max=None, monitor=None):
    rmsds = {}
    for pdbname in pdbnames:
        log.debug('Trying {}'.format(pdbname))
        try:
            pdbmodel = pdb_db.get(pdbname)
            entry_blacklist = (pdbname, )
            entry_rmsd, igc = benchmarking.benchmark_pdbmodel(pdbmodel, id_max,
                local_id_max, entry_blacklist)
            for region, v in entry_rmsd.items():
                if region == 'K1': region = 'L1'
                if region == 'K2': region = 'L2'
                if region == 'K3': region = 'L3'
                if region == 'K loops': region = 'L loops'
                if region not in rmsds:
                    rmsds[region] = []

                rmsds[region].append(v)
        except (igm.BCRBaseError) as err:
            log.error('{} failed: {}'.format(pdbname, err))
        else:
            if monitor:
                if rmsds[monitor[0]][-1] > monitor[1]:
                    log.info('{} {:.2f}'.format(pdbname, rmsds[monitor[0]][-1]))
            else:
                log.info('{} {:.2f}'.format(pdbname, rmsds.get('Total', [999])[-1]))

    return rmsds


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('pdbdir', help='Directory of pdbfiles to benchmark.')
    parser.add_argument('-v', help='Verbose', action='store_true')
    parser.add_argument('-m', help='Maximum %%ID', default=0.98, type=float)
    parser.add_argument('-M', help='Maximum local %%ID', default=None, type=float)
    parser.add_argument('--plot', help='Plot histogram', default=False, nargs='?', const=True)
    parser.add_argument('--monitor', help='Monitor structure, [structure ..]', default=None, nargs=2)
    args = parser.parse_args()

    logfmt = '%(asctime)s %(name)-15s: %(levelname)-8s %(message)s'
    if args.v:
        logging.basicConfig(level=logging.DEBUG, format=logfmt, datefmt='%Y-%m-%d %H:%M')
    else:
        logging.basicConfig(level=logging.INFO, format=logfmt, datefmt='%Y-%m-%d %H:%M')

    if args.monitor:
        args.monitor = (args.monitor[0], float(args.monitor[1]))

    pdb_db = igm.db.PDBDirectoryDatabase(args.pdbdir)
    rmsds = benchmark_pdbnames(sorted(pdb_db), pdb_db, args.m, args.M, args.monitor)

    benchmarking.print_results(rmsds)

    if args.plot:
        if not args.M:
            args.M = args.m
        elif args.M > 1.:
            args.M = 1.

        if args.plot is True:
            args.plot = ('histogram_bench_m{:0>3.0f}M{:0>3.0f}.pdf'
                ''.format(args.m*100, args.M*100))
        elif args.plot[-4:] != '.pdf':
            args.plot += '.pdf'

        plottitle = r'RMSD ($n = {},\, ID < {:.0f}\%,\, ID_{{local}} < {:.0f}\%$)'.format(
            len(rmsds['Total']), args.m*100, args.M*100)
        benchmarking.plot_results(rmsds, args.plot, plottitle)


if __name__ == '__main__':
    main()