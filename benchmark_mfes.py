#!/usr/bin/env python
# Copyright (C) 2016 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from socket import gethostname
from time import time
from networkx import DiGraph, strongly_connected_components
from utils import print_timestamp, serialize, Stats

TMP_DIR = '/users/ali/benchmarks/' if gethostname() == 'otto' else '/tmp/benchmarks/'

N_SAMPLE = 101

#-------------------------------------------------------------------------------
# Erdos-Renyi

def __main():
    print_timestamp()
    for n, c in gen_n_c():
        _erdos_renyi(n, c)
    print_timestamp()


def gen_n_c():
    for c in (5, 7, 9):
        for n in range(35, 101, 5):
            if n > 75 and c == 9:
                continue
            yield n, c  


def _erdos_renyi(n, c):
    IDX_MEDIAN = int((N_SAMPLE - 1)*0.5)
    stats = []
    for seed, g in deserialize_erdos_renyi_graphs(n, c):
        name = 'erdos renyi'
        params = 'n %d c %d seed %d' % (n, c, seed)
        stat = _solve(g, name, params)
        stats.append(stat)
    times = sorted(s.time for s in stats if s.cost)
    assert len(times) == N_SAMPLE, len(times)
    print('Execution times')
    five_number_stat(times, '%.2f')
    print_group_separator()
    assert times[IDX_MEDIAN] != float('NaN'), times


def deserialize_erdos_renyi_graphs(n, c):
    for seed in gen_erdos_renyi_seed():
        name = 'erdos_renyi_n_{}_c_{}_seed_{}'.format(n, c, seed)
        with open(TMP_DIR + name + '.edges') as f:
            g = DiGraph()
            for line in f:
                u, v = line.split()
                u, v = int(u), int(v)
                g.add_edge(u, v, weight=1, orig_edges=[(u, v)])
        yield seed, g


def gen_erdos_renyi_seed():
    N_SAMPLE = 101
    for seed in range(1, N_SAMPLE + 2 +  1):
        if seed == 61 or seed == 78:
            continue
        yield seed  

#-------------------------------------------------------------------------------
# Tournaments

def _main():
    print_timestamp()
    for n in range(28, 35):
        _tournament(n)
    print_timestamp()

def _tournament(n):
    IDX_MEDIAN = int((N_SAMPLE - 1)*0.5)
    stats = []
    for seed, g in deserialize_tournaments(n):
        name = 'tournament'
        params = 'n %d seed %d' % (n, seed)
        stat = _solve(g, name, params)
        stats.append(stat)
    times = sorted(s.time for s in stats if s.cost)
    assert len(times) == N_SAMPLE, len(times)
    print('Execution times')
    five_number_stat(times, '%.2f')
    print_group_separator()
    assert times[IDX_MEDIAN] != float('NaN'), times


def deserialize_tournaments(n):
    N_SAMPLE = 101
    for seed in range(1, N_SAMPLE +  1):
        name = 'tournament_n_{}_seed_{}'.format(n, seed)
        with open(TMP_DIR + name + '.edges') as f:
            g = DiGraph()
            for line in f:
                u, v = line.split()
                u, v = int(u), int(v)
                g.add_edge(u, v, weight=1, orig_edges=[(u, v)])
        yield seed, g

#-------------------------------------------------------------------------------

def five_number_stat(arr, fmt):
    sz = len(arr)-1
    quantiles = (0, int(sz*0.25), int(sz*0.5), int(sz*0.75), sz)
    for i in quantiles:
        print(fmt % arr[i])
    print() 


def print_group_separator():
    print('===================================================================')    


def _solve(g, name, params):
    from grb_lazy import solve_problem as solve
    print()
    print('***  STARTING A NEW PROBLEM  ***')
    print()
    print('name:', name)
    print('params:', params)
    #
    nontriv_sccs = sum(1 for sc in strongly_connected_components(g) if len(sc)>1)
    assert nontriv_sccs == 1, nontriv_sccs
    assert g.number_of_selfloops() == 0
    #
    stats = Stats(name=name, params=params, is_optimal=True, cost=None, ILP=0, 
                  node=0, iter=0, time=None)
    start = time()
    elims, cost, cycle_matrix = solve(g, stats)
    end = time()
    #
    stats.time = end - start
    #
    stats.cost = cost
    print_stats(g, stats)
    fname = (name + '_' + params).replace(' ', '_')
    serialize(cycle_matrix, TMP_DIR + 'cycle_matrix_' + fname + '.pkl.gz')
    serialize(elims,        TMP_DIR + 'solution_'     + fname + '.pkl.gz')
    #
    return stats


def print_stats(g, stats):
    print()
    print('***  SOLVED  ***')
    print()
    print('name:', stats.name)
    print('params:', stats.params)
    print('n:', g.number_of_nodes())
    print('m:', g.number_of_edges())
    print('cost:', stats.cost)
    print()
    print('ILPs:', stats.ILP)
    print('node:', stats.node)
    print('iter:', stats.iter)
    print('time: {0:0.2f} s'.format(stats.time))
    print('------------------------------------------------------------------')


if __name__ == '__main__':
    _main()
