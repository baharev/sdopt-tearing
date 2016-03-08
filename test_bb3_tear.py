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
from bb3_tear  import solve_problem as bb3_solve


def log(*args): pass
log = print

def test_to_hessenberg_none_forbidden(n_eqs, n_vars, seed):
    #log('seed =', seed)
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_vars, seed, rng, nonsing=False, c=2)
    crosscheck(g, eqs)

def test_nonsing_none_forbidden(n_eqs, seed):
    #log('seed =', seed)
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_eqs, seed, rng, nonsing=True, c=2)
    crosscheck(g, eqs)

def crosscheck(g, eqs):
    #log('----------------')
    #log('Solving with bb')
    _, _, _, tear_set, _ = bb_solve(g, eqs)
    bb_cost = len(tear_set)
    #log('----------------')
    #log('Solving with bb3')
    res = bb3_solve(g, eqs)
    bb3_cost = res.ub
    assert bb_cost == bb3_cost, (bb_cost, bb3_cost)
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
            res = bb3_solve(g, set(irange(size)))
            assert opt == res.ub
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
    
    #test_to_hessenberg_none_forbidden(n_eqs=12, n_vars=28, seed=636)
    #quit()
    #test_to_hessenberg_none_forbidden(n_eqs=13, n_vars=29, seed=766)
    #quit()
    
#     test_to_hessenberg_none_forbidden(5, 10,  291)
#     test_to_hessenberg_none_forbidden(6, 11,  344)
#     test_to_hessenberg_none_forbidden(6, 11,  928)
#     test_to_hessenberg_none_forbidden(7,  8, 1631)
#     quit()
    
    naughty_brute_force()
    #quit()
    
    print_timestamp()
    print('Started generative testing...')
    
    import os
    os.environ['HYPOTHESIS_STORAGE_DIRECTORY'] = '/tmp/ht'
    from hypothesis import given, Settings
    from hypothesis.strategies import integers
    
    # fast: size = 6, examp = 10 000
    # ~5 min: size = 30, examp = 1500
    # ~12.5 min: size = 30, examp = 5000
    MAX_VALUE = 30
    MAX_EXAMP = 150
    
    decor = given(n_eqs  = integers(min_value=1, max_value=MAX_VALUE),
                  n_vars = integers(min_value=0, max_value=MAX_VALUE), 
                  seed   = integers(min_value=0),
                  settings = Settings(max_examples=MAX_EXAMP, timeout=3600))
    
    decor(test_to_hessenberg_none_forbidden)()
    
    #-----------------------------------------
    print_timestamp()
    log('\n###  Testing with non-singular square matrices  ###\n')
    
    # ~ 6 min: size 30, examp =  500
    # ~10 min: size 30, examp = 1000 
    SIZE = 30
    MAX_EXAMP = 50
    
    decor = given(n_eqs  = integers(min_value=SIZE, max_value=SIZE),
                  seed   = integers(min_value=0),
                  settings = Settings(max_examples=MAX_EXAMP, timeout=3600))

    decor(test_nonsing_none_forbidden)()
    
    print_timestamp()
    print('Done!')
