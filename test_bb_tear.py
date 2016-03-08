# Copyright (C) 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from random import Random
from networkx import Graph
from networkx.algorithms.bipartite import is_bipartite_node_set
from bb_tear import solve_problem, _worst_cases
from order_util import coo_matrix_to_bipartite
from plot_ordering import to_pdf
from py3compat import irange, izip
from test_utils import create_rnd_bipartite, solve_difficult_for_ilp, \
                       solve_test_matrices
from testmatrices import create_difficult_pattern
from utils import print_timestamp, deserialize, marshal_load#, marshal_dump


def log(*args): pass
# log = print


def test_to_hessenberg_none_forbidden(n_eqs, n_vars, seed):
    log('seed =', seed)
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_vars, seed, rng, nonsing=False, c=3)
    solve_problem(g, eqs)
    log()

#def test_proxy(n_eqs, seed):
#    test_nonsing_none_forbidden(n_eqs, seed)
#
#@profile
def test_nonsing_none_forbidden(n_eqs, seed):
    log('seed =', seed)
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_eqs, seed, rng, nonsing=True, c=3)
    #serialize((g, eqs), DATADIR + 'ilp_poor_performance_c=5_r=2997.pkl.gz')
    #plot_bipartite(g, eqs, forbidden)
    solve_problem(g, eqs)
    log()

#-------------------------------------------------------------------------------

def solve_difficult_for_bb():
    filepath = 'data/worst_of_10M_runs/bb_bad_perf_n=10_c=0.5_1.pkl.gz'
    g, eqs, _ = deserialize(filepath)
    rowp, colp = solve_problem(g, eqs)[0:2]
    to_pdf(g, rowp, colp)

#-------------------------------------------------------------------------------

#@profile
def naughty_brute_force():
    for size in irange(3, 6):
        print_timestamp()
        print('Testing (naughty) bipartite graphs of size', size)
        # serialized in test_utils, but optimums can be serialized here
        #opts = [ ]
        #opts = marshal_load('data/all_bips/opt_n'+str(size)+'.bin')
        #all_edgelists = marshal_load('data/all_bips/bip_n'+str(size)+'.bin')
        opts = marshal_load('data/bip_filt/opt_n'+str(size)+'.bin')
        all_edgelists = marshal_load('data/bip_filt/filt_n'+str(size)+'.bin')
        print('Loaded', len(all_edgelists), 'graphs')
        print_timestamp()
        #for edgelist in all_edgelists:
        for i, (edgelist, opt) in enumerate(izip(all_edgelists, opts)):
            assert len(edgelist) % 2 == 0
            g = Graph()
            g.add_edges_from(e for e in izip(edgelist[::2], edgelist[1::2]))
            assert len(g) == 2*size
            g.graph['name'] = str(i)
            _, _, _, tear_set, _ = solve_problem(g, set(irange(size)))
            assert opt == len(tear_set)
            #---
            #solve_problem(g, set(irange(size)))
            #---
            #opt = len(solve_problem(g, set(irange(size)))[3])
            #opts.append(opt)
            #---
            #to_pdf(g, rowp,  colp)
        #assert len(opts) == len(all_edgelists)
        #marshal_dump(opts, '/tmp/opt_n'+str(size)+'.bin')
        print([t[0] for t in _worst_cases])
        #print('Len:', len(_worst_cases))
        _worst_cases.sort(key=sort_patterns)
#         for i, (explored, g, _, rowp, colp, ub) in enumerate(_worst_cases, 1):
#             msg   = 'Index: ' + g.graph['name']
#             fname = '{0:03d}a'.format(i)
#             to_pdf(g, list(irange(size)), irange(size, 2*size), msg, fname)
#             msg   = 'OPT = {}, BT: {}'.format(ub, explored)
#             fname = '{0:03d}b'.format(i)
#             to_pdf(g, rowp, colp, msg, fname)
        _worst_cases[:] = [ ]
        print_timestamp()
        print()


def sort_patterns(tup):
    _, g, _, rowp, colp, ub = tup
    c_index = { name : i for i, name in enumerate(n for n in colp) }
    r_index = { name : i for i, name in enumerate(n for n in rowp) }
    # Last occupied columns rowwise, empty rows allowed
    adj = g.adj
    last_elem  = tuple(max(c_index[c] for c in adj[r]) if adj[r] else -1 for r in rowp)
    n_rows = len(rowp)
    # First occupied row columnwise
    first_elem = tuple(min(r_index[r] for r in adj[c]) if adj[c] else n_rows for c in colp)
    # 
    return ub, last_elem, first_elem, g.number_of_edges()

#-------------------------------------------------------------------------------

def difficult(size):
    print('Solving patterns leading to many ties (backtracking) of size', size)
    msg   = 'Size: ' + str(size)
    fname = '{0:03d}a'.format(size)
    g = create_difficult_pattern(size)
    to_pdf(g, list(irange(size)), irange(size, 2*size), msg, fname)
    #
    solve_problem(g, set(irange(size)))
    #
    explored, g, _, rowp, colp, ub = _worst_cases[0]
    msg   = 'OPT = {}, BT: {}'.format(ub, explored)
    fname = '{0:03d}b'.format(size)
    to_pdf(g, rowp, colp, msg, fname)
    _worst_cases[:] = [ ]
    print_timestamp()


def plot_dm_decomp(g, size):
    rows, cols = [ ], [ ]
    for u, v in g.edges_iter(irange(size)):
        rows.append(u)
        cols.append(v-size)
    values = [1]*len(rows)
    shape = (size, size)
    g_dup, eqs, vrs = coo_matrix_to_bipartite(rows, cols, values, shape)
    assert is_bipartite_node_set(g_dup, eqs)
    assert is_bipartite_node_set(g_dup, vrs)    
    from dm_decomp import blt_with_tearing
    msg = 'Size: {}'.format(size)
    blt_with_tearing(rows, cols, values, shape, [size-1], [0], show=True, name=msg)


################################################################################

if __name__ == '__main__':

    difficult(6)
    difficult(8)
    difficult(10)
    difficult(12)
    difficult(14)
    #difficult(16)
    #quit()
    
    naughty_brute_force()
    solve_difficult_for_bb()
    solve_difficult_for_ilp(solve_problem, print)
    solve_test_matrices(solve_problem, print)
    #quit()
    
    print('Started generative testing...')
    
    import os
    os.environ['HYPOTHESIS_STORAGE_DIRECTORY'] = '/tmp/ht'
    from hypothesis import given, Settings
    from hypothesis.strategies import integers
    
    MAX_VALUE = 12
    MAX_EXAMP = 1000
    
    decor = given(n_eqs  = integers(min_value=1, max_value=MAX_VALUE),
                  n_vars = integers(min_value=0, max_value=MAX_VALUE), 
                  seed   = integers(min_value=0),
                  settings = Settings(max_examples=MAX_EXAMP))
    
    decor(test_to_hessenberg_none_forbidden)()
    
    #-----------------------------------------
    
    SIZE = 12
    MAX_EXAMP = 1000
    
    decor = given(n_eqs  = integers(min_value=SIZE, max_value=SIZE),
                  seed   = integers(min_value=0),
                  settings = Settings(max_examples=MAX_EXAMP, timeout=3600))
    
    print('###  Nonsingular  ###')
    decor(test_nonsing_none_forbidden)()
    
    print('Done!')

