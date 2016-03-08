# Copyright (C) 2016 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from random import Random
from test_utils import create_rnd_bipartite

from bb_tear       import solve_problem as bb_solve
from test_ilp_tear import solve_problem as ilp_solve


def log(*args): pass
log = print


def test_to_hessenberg_none_forbidden(n_eqs, n_vars, seed):
    log('seed =', seed)
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_vars, seed, rng, nonsing=False, c=2)
    crosscheck(g, eqs)

def test_nonsing_none_forbidden(n_eqs, seed):
    log('seed =', seed)
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_eqs, seed, rng, nonsing=True, c=2)
    crosscheck(g, eqs)

def crosscheck(g, eqs):
    log('----------------')
    log('Solving with B&B')
    _, _, _, tear_set, _ = bb_solve(g, eqs)
    bb_cost = len(tear_set)
    log('----------------')
    log('Solving with ILP')
    _rowp, _colp, _match, tear_set, _sink_set = ilp_solve(g, eqs)
    ilp_cost = len(tear_set)
    assert bb_cost == ilp_cost, (bb_cost, ilp_cost)
    log()

#===============================================================================

if __name__ == '__main__':

    from gurobipy import setParam
    setParam('LogFile', '/tmp/gurobi.log')
    setParam('LogToConsole', 0)

    print('Started generative testing...')
    
    import os
    os.environ['HYPOTHESIS_STORAGE_DIRECTORY'] = '/tmp/ht'
    from hypothesis import given, Settings
    from hypothesis.strategies import integers
    
    MAX_VALUE = 6
    MAX_EXAMP = 100
    
    decor = given(n_eqs  = integers(min_value=1, max_value=MAX_VALUE),
                  n_vars = integers(min_value=0, max_value=MAX_VALUE), 
                  seed   = integers(min_value=0),
                  settings = Settings(max_examples=MAX_EXAMP))
    
    decor(test_to_hessenberg_none_forbidden)()
    
    #-----------------------------------------
    log('\n###  Testing with non-singular square matrices  ###\n')
    
    # For cross checking: size = 6, examp = 1000, c = 2
    SIZE = 6
    MAX_EXAMP = 100
    
    decor = given(n_eqs  = integers(min_value=SIZE, max_value=SIZE),
                  seed   = integers(min_value=0),
                  settings = Settings(max_examples=MAX_EXAMP, timeout=3600))

    decor(test_nonsing_none_forbidden)()
    
    print('Done!')
