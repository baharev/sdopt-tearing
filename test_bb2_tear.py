# Copyright (C) 2016 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from random import Random
from networkx import Graph
from py3compat import irange, izip
from test_utils import create_rnd_bipartite
from utils import print_timestamp, marshal_load

from bb_tear   import solve_problem as bb_solve
from bb2_tear  import solve_problem as bb2_solve


def log(*args): pass
log = print

def test_to_hessenberg_none_forbidden(n_eqs, n_vars, seed):
    #log('seed =', seed)
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_vars, seed, rng, nonsing=False, c=0.5)
    crosscheck(g, eqs)

def test_nonsing_none_forbidden(n_eqs, seed):
    #log('seed =', seed)
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_eqs, seed, rng, nonsing=True, c=0.5)
    crosscheck(g, eqs)

def crosscheck(g, eqs):
    #log('----------------')
    #log('Solving with bb')
    _, _, _, tear_set, _ = bb_solve(g, eqs)
    bb_cost = len(tear_set)
    #log('----------------')
    #log('Solving with bb2')
    _, _, _, tear_set, _ = bb2_solve(g, eqs)
    bb2_cost = len(tear_set)
    assert bb_cost == bb2_cost, (bb_cost, bb2_cost)
    #log()

#-------------------------------------------------------------------------------

def naughty_brute_force():
    for size in irange(1, 6):
        print_timestamp()
        print('Testing (naughty) bipartite graphs of size', size)
        opts = marshal_load('data/all_bips/opt_n'+str(size)+'.bin')
        all_edgelists = marshal_load('data/all_bips/bip_n'+str(size)+'.bin')
        print('Loaded', len(all_edgelists), 'graphs')
        print_timestamp()
        for i, (edgelist, opt) in enumerate(izip(all_edgelists, opts)):
            g = Graph()
            g.add_edges_from(e for e in izip(edgelist[::2], edgelist[1::2]))
            g.graph['name'] = str(i)
            _, _, _, tear_set, _ = bb2_solve(g, set(irange(size)))
            assert opt == len(tear_set)
            #to_pdf(g, rowp,  colp)
        #print([t[0] for t in _worst_cases])
        #print('Len:', len(_worst_cases))
        #_worst_cases.sort(key=sort_patterns)
#         for i, (explored, g, _, rowp, colp, ub) in enumerate(_worst_cases, 1):
#             msg   = 'Index: ' + g.graph['name']
#             fname = '{0:03d}a'.format(i)
#             to_pdf(g, list(irange(size)), irange(size, 2*size), msg, fname)
#             msg   = 'OPT = {}, BT: {}'.format(ub, explored)
#             fname = '{0:03d}b'.format(i)
#             to_pdf(g, rowp, colp, msg, fname)
        #_worst_cases[:] = [ ]
        print_timestamp()
        print()
        
#-------------------------------------------------------------------------------


if __name__ == '__main__':
    
    naughty_brute_force()
    #quit()
    
    print_timestamp()
    print('Started generative testing...')
    
    import os
    os.environ['HYPOTHESIS_STORAGE_DIRECTORY'] = '/tmp/ht'
    from hypothesis import given, Settings
    from hypothesis.strategies import integers
    
    MAX_VALUE = 6
    MAX_EXAMP = 10000
    
    decor = given(n_eqs  = integers(min_value=1, max_value=MAX_VALUE),
                  n_vars = integers(min_value=0, max_value=MAX_VALUE), 
                  seed   = integers(min_value=0),
                  settings = Settings(max_examples=MAX_EXAMP, timeout=3600))
    
    decor(test_to_hessenberg_none_forbidden)()
    
    #-----------------------------------------
    print_timestamp()
    log('\n###  Testing with non-singular square matrices  ###\n')
    
    SIZE = 6
    MAX_EXAMP = 10000
    
    decor = given(n_eqs  = integers(min_value=SIZE, max_value=SIZE),
                  seed   = integers(min_value=0),
                  settings = Settings(max_examples=MAX_EXAMP, timeout=3600))

    decor(test_nonsing_none_forbidden)()
    
    print_timestamp()
    print('Done!')
