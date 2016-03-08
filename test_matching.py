# Copyright (C) 2016 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from random import Random
from test_utils import create_rnd_bipartite
from utils import print_timestamp

from matching import maxmatch_len, nx_maxmatch_len

def log(*args): pass
log = print

def test_with_random_matrices(n_eqs, n_vars, seed):
    #log('seed =', seed)
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_vars, seed, rng, nonsing=False, c=1)
    nrows, ncols, nzeros = len(eqs), len(g) - len(eqs), g.number_of_edges()
    log('rows: {}, cols: {}, nzeros: {}'.format(nrows, ncols, nzeros))
    max_len = maxmatch_len(g, eqs)
    nx_len  = nx_maxmatch_len(g)
    assert max_len == nx_len, (max_len, nx_len)
    log('underdetermined by:', min(nrows, ncols) - max_len)

#-------------------------------------------------------------------------------


if __name__ == '__main__':
    
    print_timestamp()
    print('Started generative testing...')
    
    import os
    os.environ['HYPOTHESIS_STORAGE_DIRECTORY'] = '/tmp/ht'
    from hypothesis import given, Settings
    from hypothesis.strategies import integers
    
    MAX_VALUE = 50
    MAX_EXAMP = 10
    
    decor = given(n_eqs  = integers(min_value=1, max_value=MAX_VALUE),
                  n_vars = integers(min_value=0, max_value=MAX_VALUE), 
                  seed   = integers(min_value=0),
                  settings = Settings(max_examples=MAX_EXAMP, timeout=3600))
    
    decor(test_with_random_matrices)()

    print_timestamp()
    print('Done!')
