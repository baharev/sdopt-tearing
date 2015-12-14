# Copyright (C) 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from random import Random
from networkx.algorithms.bipartite import is_bipartite_node_set
from heap_md import to_hessenberg_form, to_spiked_form
from order_util import to_bipart_w_weights, check_nonincreasing_envelope
from test_utils import create_rnd_bipartite, create_diagonal_matrix, \
                       raw_rnd_bipartite
from utils import pairwise

def log(*args): pass
#log = print

#-------------------------------------------------------------------------------

def test_to_hessenberg_none_forbidden(n_eqs, n_vars, seed):
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_vars, seed, rng)
    to_hessenberg_form(g, eqs)

def test_to_hessenberg_some_forbidden(n_eqs, n_vars, seed):
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_vars, seed, rng)
    edges = g.edges(eqs)
    log('Edges:', edges)
    forbidden = set(rng.choice(edges) for _ in range(len(edges)//2) )
    log('Forbidden:', list(forbidden))
    to_hessenberg_form(g, eqs, forbidden)

def test_to_hessenberg_all_forbidden(n_eqs, n_vars, seed):
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n_eqs, n_vars, seed, rng)
    to_hessenberg_form(g, eqs, set(g.edges_iter(eqs)))

#-------------------------------------------------------------------------------

def test_to_spiked_none_forbidden(n, seed):
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n, n, seed, rng)
    to_spiked_form(g, eqs)

def test_to_spiked_some_forbidden(n, seed):
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n, n, seed, rng)
    edges = g.edges(eqs)
    log('Edges:', edges)
    forbidden = set( rng.choice(edges) for _ in range(len(edges)//2) )
    log('Forbidden:', list(forbidden))
    to_spiked_form(g, eqs, forbidden)

def test_to_spiked_all_forbidden(n, seed):
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n, n, seed, rng)
    to_spiked_form(g, eqs, set(g.edges_iter(eqs)))

def test_nonsingular_spiked(n, seed):
    rng = Random(seed)
    g, eqs = create_rnd_bipartite(n, n, seed, rng, nonsing=True)
    edges = g.edges(eqs)
    log('Edges:', edges)
    forbidden = set( rng.choice(edges) for _ in range(len(edges)//2) )
    log('Forbidden:', list(forbidden))
    singular = to_spiked_form(g, eqs, forbidden)[0]
    assert not singular

#-------------------------------------------------------------------------------

def test_to_hessenberg_none_forbidden_diagonal(n, seed):
    rng = Random(seed)
    g, eqs = create_diagonal_matrix(n, seed, rng)
    rowp = to_hessenberg_form(g, eqs)[0]
    assert rowp == sorted(rowp)

def test_to_hessenberg_all_forbidden_diagonal(n, seed):
    rng = Random(seed)
    g, eqs = create_diagonal_matrix(n, seed, rng)
    rowp = to_hessenberg_form(g, eqs, set(g.edges_iter(eqs)))[0]
    assert rowp == sorted(rowp)

#-------------------------------------------------------------------------------

def test_to_spiked_none_forbidden_diagonal(n, seed):
    rng = Random(seed)
    g, eqs = create_diagonal_matrix(n, seed, rng)
    singular, rowp = to_spiked_form(g, eqs)[:2]
    assert not singular
    assert rowp == sorted(rowp)

def test_to_spiked_all_forbidden_diagonal(n, seed):
    rng = Random(seed)
    g, eqs = create_diagonal_matrix(n, seed, rng)
    singular, rowp = to_spiked_form(g, eqs, set(g.edges_iter(eqs)))[:2]
    assert not singular
    assert rowp == sorted(rowp)

#-------------------------------------------------------------------------------

#def test_proxy(n_eqs, n_vars, seed):
#    test_weighted_bipart(n_eqs, n_vars, seed)

#@profile
def test_weighted_bipart(n_eqs, n_vars, seed):
    log('---------------------------------------------------------------------')
    bip = raw_rnd_bipartite(n_eqs, n_vars, seed)
    assert is_bipartite_node_set(bip, range(n_eqs))
    #
    rng = Random(seed)
    cols_rowwise = [ list(bip[i]) for i in range(n_eqs) ]
    vals_rowwise = [ rnd_weights(rng, len(cols)) for cols in cols_rowwise ]
    row_weights  = [ sum(vals, 0.0) for vals in vals_rowwise ]
    #
    g, eqs, mapping, _ = to_bipart_w_weights(cols_rowwise, vals_rowwise)
    rowp, colp, _, sinks, row_matches, _ = to_hessenberg_form(g, eqs)
    #
    rowp  = [ mapping[r] for r in rowp ]
    sinks = [ mapping[r] for r in sinks ]
    row_matches = [ mapping[r] for r in row_matches ]
    #
    log('bip rowp:', rowp)
    log('eqs =', n_eqs, ' vars =', n_vars, ' edges =', bip.number_of_edges())
    #
    check_nonincreasing_envelope(bip, rowp, colp)
    check_nondecreasing_row_weights(bip, rowp, colp, row_weights)
    log('---------------------------------------------------------------------')

def rnd_weights(rng, length):
    return [rng.random() for _ in range(length)]

def check_nondecreasing_row_weights(bip, rowp, colp, row_weights):
    # Last occupied columns rowwise, empty rows allowed
    c_index = { name : i for i, name in enumerate(n for n in colp) }
    last_elem = [max(c_index[c] for c in bip[r]) if bip[r] else -1 for r in rowp]
    log('Last elems:', last_elem)
    perm_weights = [row_weights[r] for r in rowp]
    # Checking if we really break ties if the last elements equal
    for i, ((c1,c2), (w1,w2)) in enumerate(zip(pairwise(last_elem), pairwise(perm_weights))):
        log('i:', i)
        log('Last elems:', c1, c2)
        log('Weights:', w1, w2)
        assert (c1, w1) <= (c2, w2)

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    
    print('Started generative testing...')
    
    import os
    os.environ['HYPOTHESIS_STORAGE_DIRECTORY'] = '/tmp/ht'
    from hypothesis import given, Settings
    from hypothesis.strategies import integers
    
    MAX_VALUE = 8
    MAX_EXAMP = 500
    
    hessenberg_decor = given(n_eqs  = integers(min_value=1, max_value=MAX_VALUE),
                             n_vars = integers(min_value=0, max_value=MAX_VALUE), 
                             seed   = integers(min_value=0),
                             settings = Settings(max_examples=MAX_EXAMP))
    
    #hessenberg_decor(test_proxy)()
    #quit()
    hessenberg_decor(test_weighted_bipart)()
    #
    hessenberg_decor(test_to_hessenberg_none_forbidden)()
    hessenberg_decor(test_to_hessenberg_some_forbidden)()
    hessenberg_decor(test_to_hessenberg_all_forbidden)()
    
    #-----
    
    hessenberg_decor = given(n    = integers(min_value=1, max_value=MAX_VALUE),
                             seed = integers(min_value=0),
                             settings = Settings(max_examples=MAX_EXAMP))
    
    # Checking whether the row order is retained (hashing could mess it up)
    hessenberg_decor(test_to_hessenberg_none_forbidden_diagonal)()
    hessenberg_decor(test_to_hessenberg_all_forbidden_diagonal)()
    
    #-----
   
    spiked_decor = given(n    = integers(min_value=1, max_value=MAX_VALUE),
                         seed = integers(min_value=0),
                         settings = Settings(max_examples=MAX_EXAMP))
    
    spiked_decor(test_nonsingular_spiked)()
    spiked_decor(test_to_spiked_none_forbidden)()
    spiked_decor(test_to_spiked_some_forbidden)()
    spiked_decor(test_to_spiked_all_forbidden)()
    
    # Checking whether the row order is retained (hashing could mess it up)    
    spiked_decor(test_to_spiked_none_forbidden_diagonal)()
    spiked_decor(test_to_spiked_all_forbidden_diagonal)()
    
    print('Done!')
