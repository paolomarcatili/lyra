# -*- coding: utf-8 -*-
"""
Common benchmarking routines.

"""

from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import bcr_models as igm


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


def benchmark_pdbmodel(pdbmodel, id_max, local_id_max, blacklist=()):
    """Benchmark a single entry."""
    ref_pdb = igm.IgComplex.renumber_pdbmodel(pdbmodel)

    #Extract sequences
    seq1, seq2 = [igm.utils.chain2seq(chain) for chain in ref_pdb]

    #Build PDB model
    igc = igm.IgComplex(seq1, seq2)
    igc.hmmsearch(*igm.db.builtin_hmms())
    igc.canonical_structures(igm.db.BuiltinCsDatabase())
    #igc.canonical_structures(igm.db.CsDatabase())
    igc.find_templates(id_max=id_max, local_id_max=local_id_max, blacklist=blacklist, method='framework')
    igc.build_structure()

    excluding = 'H3' if 'H' in [c.chain_type for c in igc] else None

    #Align IgComplex to reference pdb model
    rmsd = igc.benchmark(ref_pdb, excluding=excluding)

    return rmsd, igc


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
    plt.savefig(plotname)