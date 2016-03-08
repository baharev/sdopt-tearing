# Copyright (C) 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from random import Random
from networkx import Graph
from ilp_tear import solve_with_pcm
#from plot_ordering import plot_bipartite_no_red_greedy_order as plot_bipartite
from py3compat import irange, izip
from test_utils import create_rnd_bipartite, solve_test_matrices
from utils import print_timestamp, marshal_load

def log(*args): pass
log = print


def solve_problem(g, eqs, forbidden=set()):
    return solve_with_pcm(g, eqs, forbidden, use_min_degree=True)

# FIXME - Has allowed edge in the TORN graph
#       - Profile
#       - Refactor, move things to the appropriate modules
#-------------------------------------------------------------------------------

def test_to_hessenberg_none_forbidden(n_eqs, n_vars, seed):
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_vars, seed, rng)
    solve_problem(g, eqs)


def test_to_hessenberg_some_forbidden(n_eqs, n_vars, seed):
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_vars, seed, rng)
    edges = g.edges(eqs)
    log('Edges:', edges)
    forbidden = set(rng.choice(edges) for _ in range(len(edges)//2) )
    log('Forbidden:', list(forbidden))
    #plot_bipartite(g, eqs, forbidden, 'input')  # see commented out import!
    solve_problem(g, eqs, forbidden)


def test_to_hessenberg_all_forbidden(n_eqs, n_vars, seed):
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_vars, seed, rng)
    forbidden = set(g.edges_iter(eqs))
    solve_problem(g, eqs, forbidden)

#-------------------------------------------------------------------------------

def test_nonsing_none_forbidden(n_eqs, seed):
    print('seed =', seed)
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_eqs, seed, rng, nonsing=True)
    #serialize((g, eqs), DATADIR + 'ilp_poor_performance_c=5_r=2997.pkl.gz')
    #plot_bipartite(g, eqs, forbidden)  # see commented out import!
    solve_problem(g, eqs)


def test_nonsing_some_forbidden(n_eqs, seed):
    print('seed =', seed)
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_eqs, seed, rng, nonsing=True)
    edges = g.edges(eqs)
    log('Edges:', edges)
    forbidden = set(rng.choice(edges) for _ in range(len(edges)//2) )
    log('Forbidden:', list(forbidden))
    #plot_bipartite(g, eqs, forbidden)  # see commented out import!
    solve_problem(g, eqs, forbidden)


def test_nonsing_all_forbidden(n_eqs, seed):
    print('\nseed =', seed)
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_eqs, seed, rng, nonsing=True)
    forbidden = set(g.edges_iter(eqs))
    solve_problem(g, eqs, forbidden)

################################################################################
# FIXME - This code is now triplicated, move it to test_utils
#       - Use pickle instead so that I can test with Python 3 too
#       - Split this code into two: one for serialization, and another for testing

def naughty_brute_force():
    for size in irange(1, 7):
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
            _, _, _, tear_set, _ = solve_problem(g, set(irange(size)))
            assert opt == len(tear_set)
        print('Done with size', size)
        print_timestamp()
        print()


################################################################################

if __name__ == '__main__':
    
    from gurobipy import setParam
    setParam('LogFile', '/tmp/gurobi.log')
    setParam('LogToConsole', 0)

    #naughty_brute_force()
    #quit()

    #solve_difficult_for_ilp(solve_problem, log)
    solve_test_matrices(solve_problem, log)
    #quit()
    
    print('Started generative testing...')
    
    import os
    os.environ['HYPOTHESIS_STORAGE_DIRECTORY'] = '/tmp/ht'
    from hypothesis import given, Settings
    from hypothesis.strategies import integers
    
    MAX_VALUE = 4
    MAX_EXAMP = 100
    
    decor = given(n_eqs  = integers(min_value=1, max_value=MAX_VALUE),
                  n_vars = integers(min_value=0, max_value=MAX_VALUE), 
                  seed   = integers(min_value=0),
                  settings = Settings(max_examples=MAX_EXAMP))
    
    decor(test_to_hessenberg_none_forbidden)()
    decor(test_to_hessenberg_all_forbidden)()
    decor(test_to_hessenberg_some_forbidden)()
    
    #-----------------------------------------
    
    SIZE = 4
    MAX_EXAMP = 100
    
    decor = given(n_eqs  = integers(min_value=SIZE, max_value=SIZE),
                  seed   = integers(min_value=0),
                  settings = Settings(max_examples=MAX_EXAMP, timeout=3600))
    
    print('###  None forbidden  ###')
    decor(test_nonsing_none_forbidden)()
    print('###  All forbidden  ###')
    decor(test_nonsing_all_forbidden)()
    print('###  Some forbidden  ###')
    decor(test_nonsing_some_forbidden)()
    
    print('Done!')
