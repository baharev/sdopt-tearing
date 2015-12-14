# Copyright (C) 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from random import Random
from dm_decomp import dm_decomp
from test_utils import create_rnd_bipartite

def test_DM(n_eqs, n_vars, seed):
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_vars, seed, rng)
    dm_decomp(g, eqs)

if __name__ == '__main__':

    print('Started generative testing...')

    import os
    os.environ['HYPOTHESIS_STORAGE_DIRECTORY'] = '/tmp/ht'
    from hypothesis import given, Settings
    from hypothesis.strategies import integers
    
    MAX_VALUE = 12
    MAX_EXAMP = 10000
    
    DM_decor = given(n_eqs  = integers(min_value=1, max_value=MAX_VALUE),
                     n_vars = integers(min_value=0, max_value=MAX_VALUE), 
                     seed   = integers(min_value=0),
                     settings = Settings(max_examples=MAX_EXAMP))
    
    DM_decor(test_DM)()
    
    print('Done!')
