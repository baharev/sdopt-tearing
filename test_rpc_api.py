# Copyright (C) 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from random import Random
#from dm_decomp import dbg_show_coomat
from order_util import get_row_weights
from py3compat import irange
from rpc_api import hessenberg, fine_dulmage_mendelsohn, tearing_hand_guided
from test_utils import create_coomat
from utils import pairwise

def log(*args): pass
#log = print

def test_rpc(n_rows, n_cols, seed):
    rng = Random(seed)
    #
    g, rows, cols, values = create_coomat(n_rows, n_cols, rng)
    #
    #print('Input:')
    #dbg_show_coomat(rows, cols, values, (n_rows, n_cols))
    #
    #---------------------------------------------------------------------------
    result = hessenberg(rows, cols, values, n_rows, n_cols, 'IGNORE')
    check(result, n_rows, n_cols, values)
    #
    #---------------------------------------------------------------------------
    result = hessenberg(rows, cols, values, n_rows, n_cols, 'MIN_FIRST')
    compare = lambda c1, w1, c2, w2: (c1,w1) <= (c2,w2)
    check_hessenberg_tie_breaking(g, n_rows, n_cols, values, result, compare)
    #
    #---------------------------------------------------------------------------
    result = hessenberg(rows, cols, values, n_rows, n_cols, 'MAX_FIRST')
    compare = lambda c1, w1, c2, w2: (c1,w2) <= (c2,w1)
    check_hessenberg_tie_breaking(g, n_rows, n_cols, values, result, compare)
    #
    #---------------------------------------------------------------------------
    result = fine_dulmage_mendelsohn(rows, cols, values, n_rows, n_cols, upper=False, minimize=True)
    check(result, n_rows, n_cols, values)
    #
    #---------------------------------------------------------------------------
    torn_rows, torn_cols = [ ], [ ]
    #
    if n_rows:
        torn_rows = rng.sample(irange(n_rows), rng.randint(0, n_rows-1))
    #
    if n_cols:
        torn_cols = rng.sample(irange(n_cols), rng.randint(0, n_cols-1))
    #
    result = tearing_hand_guided(rows, cols, values, n_rows, n_cols, torn_rows, torn_cols)
    check(result, n_rows, n_cols, values)


def check(result, n_rows, n_cols, values):
    msg = result.get('error_msg')
    if msg is not None:
        assert n_rows==0 or n_cols==0 or not values, msg 


def check_hessenberg_tie_breaking(g, n_rows, n_cols, values, result, compare):
    if result.get('error_msg'):
        check(result, n_rows, n_cols, values)
        return
    #
    rowp, colp = result['rowp'], result['colp']
    rperm = [-1]*len(rowp)
    for i, r in enumerate(rowp):
        rperm[r] = i
    # Last occupied columns rowwise, empty rows allowed
    n_rows = len(rperm)
    last_elem = [max(colp[c-n_rows] for c in g[r]) if g[r] else -1 for r in rperm]
    log('Last elems:', last_elem)
    row_weights = get_row_weights(g, n_rows)
    perm_weights = [row_weights[r] for r in rperm]
    # Checking if we really break ties if the last elements equal
    for i, ((c1,c2), (w1,w2)) in enumerate(zip(pairwise(last_elem), pairwise(perm_weights))):
        log('i:', i)
        log('Last elems:', c1, c2)
        log('Weights:', w1, w2)
        #assert (c1, w1) <= (c2, w2)
        assert compare(c1, w1, c2, w2)


if __name__ == '__main__':
    
    print('Started generative testing...')
    
    import os
    os.environ['HYPOTHESIS_STORAGE_DIRECTORY'] = '/tmp/ht'
    from hypothesis import given, Settings
    from hypothesis.strategies import integers
    
    MAX_VALUE = 5
    MAX_EXAMP = 1000
    
    decor = given(n_rows = integers(min_value=0, max_value=MAX_VALUE),
                  n_cols = integers(min_value=0, max_value=MAX_VALUE), 
                  seed   = integers(min_value=0),
                  settings = Settings(max_examples=MAX_EXAMP))
    
    decor(test_rpc)()
    print('Done!')
