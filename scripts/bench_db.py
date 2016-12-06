#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Benchmark template database
"""

from __future__ import print_function
from __future__ import division

import logging
import argparse
import random

import joblib
import numpy as np
import matplotlib.pyplot as plt

import bcr_models as igm


#Define a log
log = logging.getLogger('igm/benchmark')
#memory = joblib.Memory('/home/micker/Ramdisk', verbose=0)


ORDER = [
    'A1', 'A2', 'A3', 'A loops',
    'B1', 'B2', 'B3', 'B loops',
    'L1', 'L2', 'L3', 'L loops',
    'K1', 'K2', 'K3', 'K loops',
    'H1', 'H2', 'H3', 'H loops', 'H loops excl. H3',
    'Binding site excl. H3',
    'Binding site',
    'Total', 'Total A', 'Total B', 'Total H', 'Total K', 'Total L',
    'Framework', 'Framework A', 'Framework B',
    'Framework H', 'Framework K', 'Framework L',
]


#@memory.cache
def benchmark_db(template_db, chains, id_max=0.98, local_id_max=None,
    blacklist=(), subset=None, monitor=None):
    """Benchmark template_db."""
    pdbnames = dict()
    for entry in template_db:
        if entry.chain not in chains:
            continue
        pdbnames[entry.pdbname] = pdbnames.get(entry.pdbname, 0) + 1

    skipped = sorted([p for p in pdbnames if pdbnames[p] != 2])
    log.info('Skipping:    ' + ', '.join(skipped))
    log.info('Blacklisted: ' + ', '.join(set(pdbnames) & set(blacklist)))

    pdbnames = sorted(set(pdbnames) - set(skipped) - set(blacklist))

    if subset:
        pdbnames = random.sample(pdbnames, subset)

    rmsds = {}
    for pdbname in pdbnames:
        try:
            entry_rmsd = benchmark_entry(pdbname, id_max, local_id_max)
            for region, v in entry_rmsd.items():
                if region == 'K1': region = 'L1'
                if region == 'K2': region = 'L2'
                if region == 'K3': region = 'L3'
                if region == 'K loops': region = 'L loops'
                if region not in rmsds:
                    rmsds[region] = []

                rmsds[region].append(v)
        except (igm.BCRBaseError) as err:
            log.debug('{} failed: {}'.format(pdbname, err))
        else:
            if monitor:
                if rmsds[monitor[0]][-1] > monitor[1]:
                    log.info('{} {:.2f}'.format(pdbname, rmsds.get('Total', [999])[-1]))
            else:
                log.info('{} {:.2f}'.format(pdbname, rmsds.get('Total', [999])[-1]))

    return rmsds


#@memory.cache
def benchmark_entry(pdbname, id_max, local_id_max, get_igc=False):
    """Benchmark a single entry."""
    pdb_db      = igm.db.BuiltinPDBDatabase()
    template_db = igm.db.BuiltinTemplateDatabase()
    hmms        = igm.db.builtin_hmms()
    csdb        = igm.db.BuiltinCsDatabase()
    #csdb        = igm.canonical_structures.RandomForestCsDatabase()

    log.debug('Trying: ' + pdbname)

    ref_pdb = igm.IgComplex.renumber_pdbmodel(pdb_db.get(pdbname))

    #Extract sequences
    seq1, seq2 = [igm.utils.chain2seq(chain) for chain in ref_pdb]

    #Build PDB model
    igc = igm.IgComplex(seq1, seq2, template_db, pdb_db)
    igc.hmmsearch(*hmms)
    igc.canonical_structures(csdb)
    igc.find_templates(id_max=id_max, local_id_max=local_id_max, blacklist=(pdbname, ))
    igc.build_structure()

    excluding = 'H3' if 'H' in [c.chain_type for c in igc] else None

    #Align IgComplex to reference pdb model
    rmsd = igc.benchmark(ref_pdb, excluding=excluding)

    if get_igc:
        return rmsd, igc

    return rmsd


def print_results(rmsds):
    bins = [0, 1, 2, 5, 10]
    print('Region/Loop             ', end='')
    print('        avg',              end='')
    print('   n ',                    end='')
    for i in range(len(bins) - 1):
        print('{:>6} '.format('< {}'.format(bins[i+1])), end='')
    print(' > {}'.format(bins[-1]))

    ordered = sorted(rmsds, key=lambda r: (ORDER.index(r), r) if r in ORDER else (999, r))
    for region in ordered:
        vals = rmsds[region]
        avg  = np.average(vals)
        print('{:<24} '.format(region),   end='')
        print('{:>10.2f}'.format(avg),    end='')
        print('{:>4} '.format(len(vals)), end='')
        hist, edges = np.histogram(vals, bins=bins)

        inhist = 0
        for i, n in enumerate(hist):
            inhist += n
            print('{:>6} '.format(n), end='')

        rest = len(vals) - inhist
        print(' {:>4}'.format(rest))


def plot_results(rmsds, plotname, plottitle):
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 7)

    colors = ['#66cc00', '#ffcc00', '#ff9900', '#990000', '#000000', ]

    x = 0
    bins = [0, 1, 2, 5, 10]
    ordered = sorted(rmsds, key=lambda r: (ORDER.index(r), r) if r in ORDER else (999, r))
    for region in ordered:
        vals = rmsds[region]
        avg = np.average(vals)
        hist, edges = np.histogram(vals, bins=bins)

        inhist = 0
        width = 0.95
        bac = []
        for i, n in enumerate(hist):
            b = ax.bar([x], [n], width, bottom=[inhist], color=colors[i], linewidth=0)
            bac.append(b)
            inhist += n

        rest = len(vals) - inhist
        b = ax.bar([x], [rest], width, bottom=[inhist], color=colors[-1], linewidth=0)
        bac.append(b)
        ax.text(x+width/2, len(vals), u'{:.2f}Å'.format(avg), ha='center', va='bottom')

        x += 1

    ind = np.arange(len(rmsds))
    plt.xticks(ind+width/2., ordered, rotation=30)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment('right')

    ax.set_title(plottitle)
    ax.set_ylabel('n')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.xlim(0, len(rmsds))
    plt.ylim(0, int((len(rmsds['Total']) + 1) * 1.1))

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.legend(bac, (
        'RMSD $\leq 1$',
        'RMSD $\leq 2$',
        'RMSD $\leq 5$',
        'RMSD $\leq 10$',
        'RMSD $> 10$',
    ), loc=3)

    plt.tight_layout()
    plt.savefig(args.plot)


def print_igc(pdbname, igc):
    print(pdbname, igc[0])
    print(igc[0].show_alignment())
    print()
    print(igc[0].show_templates())

    print()

    print(pdbname, igc[1])
    print(igc[1].show_alignment())
    print()
    print(igc[1].show_templates())
    print()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('chains', choices=['H', 'K', 'L', 'A', 'B'], nargs='+',
        help='Which chains to include in the benchmark')
    parser.add_argument('-v', help='Verbose', action='store_true')
    parser.add_argument('-m', help='Maximum %%ID', default=0.98, type=float)
    parser.add_argument('-M', help='Maximum local %%ID', default=None, type=float)
    parser.add_argument('--single', help='Benchmark single pdbid', default=None)
    parser.add_argument('--monitor', help='Monitor structure', default=None)
    parser.add_argument('--save', help='Save [single] output model', default=None)
    parser.add_argument('--plot', help='Plot histogram', default=False, nargs='?', const=True)
    parser.add_argument('--subset', help='Random subset', type=int, default=None)
    parser.add_argument('--random-seed', help='Random seed', default=1337)
    args = parser.parse_args()

    logfmt = '%(asctime)s %(name)-12s: %(levelname)-8s %(message)s'
    if args.v:
        logging.basicConfig(level=logging.DEBUG, format=logfmt, datefmt='%Y-%m-%d %H:%M')
    else:
        logging.basicConfig(level=logging.INFO, format=logfmt, datefmt='%Y-%m-%d %H:%M')

    random.seed(args.random_seed)

                 #TCR
    blacklist = []#'2NW2', '1BWM', '1KTK', '2P1Y', '3ARB', '3TZV']

    template_db = igm.db.BuiltinTemplateDatabase()

    if args.single:
        rmsd, igc = benchmark_entry(args.single, args.m, args.M, get_igc=True)
        print_igc(args.single, igc)

        ordered = sorted(rmsd, key=lambda r: (ORDER.index(r), r) if r in ORDER else (999, r))
        for region in ordered:
            print('{:<22} {:.2f} Å'.format(region + ':', rmsd[region]))

        if args.save:
            raise NotImplementedError
        exit()

    if args.monitor:
        args.monitor = monitor.split(',')
        args.monitor[1] = float(monitor[1])

    #Do the benchmarking
    rmsds = benchmark_db(template_db, args.chains, args.m, args.M,
        blacklist, args.subset, args.monitor)
    print_results(rmsds)

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
        plot_results(rmsds, args.plot, plottitle)

