# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>

# Heap-based minimum-degree ordering with NO lookahead.
# 
# See also min_degree.py which uses lookahead, and simple_md.py which is a 
# hacked version of min_degree.py that still uses repeated linear scans to find 
# the minimum degree nodes but does not do lookahead. 
from __future__ import print_function
from py3compat import cPickle_loads, cPickle_dumps, cPickle_HIGHEST_PROTOCOL
from itertools import chain
from networkx import max_weight_matching
from pqueue import PriorityQueue as heapdict
from py3compat import irange
from order_util import colp_to_spiked_form, get_hessenberg_order, check_spiked_form,\
                       coo_matrix_to_bipartite, partial_relabel, argsort, \
                       get_inverse_perm, get_row_weights


def hessenberg(rows, cols, values, n_rows, n_cols, tie_breaking):
    'Tie breaking options: MIN_FIRST, MAX_FIRST, IGNORE'
    assert tie_breaking in ('IGNORE', 'MIN_FIRST', 'MAX_FIRST'), tie_breaking
    # The col IDs in cols are shifted by n_rows, must undo later
    g, eqs, _ = coo_matrix_to_bipartite(rows, cols, values, (n_rows, n_cols))
    if tie_breaking != 'IGNORE':
        # Relabel the rows such that they are ordered by weight
        row_weights = get_row_weights(g, n_rows)
        reverse = True if tie_breaking == 'MAX_FIRST' else False
        row_pos = argsort(row_weights, reverse)
        mapping = {n: i for i, n in enumerate(row_pos)}
        #
        eqs = set(mapping[eq] for eq in eqs)
        g = partial_relabel(g, mapping)
    #
    rperm, cperm, _, _, _, _ = to_hessenberg_form(g, eqs)
    # Finally, shift the colp such that it is a permutation of 0 .. n_cols-1
    cperm = [c-n_rows for c in cperm]
    #
    if tie_breaking != 'IGNORE':
        rperm = [row_pos[r] for r in rperm]
    #
    rowp, colp = get_inverse_perm(rperm, cperm)
    assert sorted(rowp) == list(irange(n_rows))
    assert sorted(colp) == list(irange(n_cols))  
    return rowp, colp

################################################################################
#
# TODO Hereafter, rowp and colp here seems to be consistently used for 
# rperm and cperm, the permuted row and col identifiers.
#
################################################################################

def to_spiked_form(g, eqs, forbidden=None):
    '''Returns the tuple of: bool singular, [row permutation], 
    [column permutation], [spike variables], [residual equations]. The spikes 
    and the residuals are ordered according to the permutation.'''
    # Check singularity, apparently only the permutation to spiked form needs it
    assert 2*len(eqs) == len(g),  'Not a square matrix!'
    matches = max_weight_matching(g)
    if len(matches) != 2*len(eqs):
        return (True, [], [], [], [])
    if forbidden is None:
        forbidden = set()
    rowp, colp_hess, matches, tear_set, sink_set = min_degree(g, eqs, forbidden)
    colp = colp_to_spiked_form(rowp, colp_hess, matches, tear_set, sink_set)
    check_spiked_form(g, rowp, colp, tear_set)
    #from plot_ordering import plot_hessenberg, plot_bipartite
    #plot_hessenberg(g, rowp, colp_hess, [], '')
    #plot_bipartite(g, forbidden, rowp, colp)
    tears = [c for c in colp if c in tear_set]
    sinks = [r for r in rowp if r in sink_set]
    return (False, rowp, colp, tears, sinks)


def to_hessenberg_form(g, eqs, forbidden=None):
    '''Returns the tuple of: [row permutation], [column permutation], 
    [guessed variables], [residual equations], [row matches], [col matches]. 
    Everything is ordered according to the permutation.'''
    rowp, colp, matches, tear_set, sink_set = min_degree(g, eqs, forbidden)
    tears = [c for c in colp if c in tear_set]
    sinks = [r for r in rowp if r in sink_set]
    row_matches = [r for r in rowp if r in matches]
    col_matches = [c for c in colp if c in matches]
    return (rowp, colp, tears, sinks, row_matches, col_matches)


def min_degree(g_orig, eqs, forbidden=None):
    '''Returns: tuple([row permutation], [column permutation], 
    {eq: var and var: eq matches}, set(tear vars), set(residual equations)).'''
    # Duplicated in bb_tear.initial_solution with none forbidden
    assert eqs
    if forbidden is None:
        forbidden = set()
    if not isinstance(eqs, (set, dict)):
        eqs = set(eqs)  # Make sure that `n in eqs` will be O(1).
    g_allowed, g = setup_graphs(g_orig, forbidden)
    eq_tot = create_heap(g_allowed, g, eqs)
    rowp, matches = [ ], { }
    while eq_tot:
        (_cost, _tot, _eq), eq = eq_tot.popitem()
        #assert _eq == eq, (_eq, eq)
        #print('Eq:', eq)
        rowp.append(eq)
        
        if g_allowed[eq]:
            var = sorted(g_allowed[eq])[0] # or [-1] for last
            assert eq  not in matches
            assert var not in matches
            matches[eq]  = var
            matches[var] = eq
            #print('Var:', var)
        
        vrs = sorted(g[eq])
        
        eqs_update = set(chain.from_iterable(g[v] for v in vrs))
        eqs_update.discard(eq)
        
        g_allowed.remove_node(eq)
        g.remove_node(eq)
    
        g_allowed.remove_nodes_from(vrs)
        g.remove_nodes_from(vrs)
    
        for e in sorted(eqs_update): # keep in sync with create_heap
            tot = len(g[e])
            cost = tot-1 if g_allowed[e] else tot
            eq_tot[e]  = (cost, tot, e)
    
    assert len(rowp) == len(eqs)
    # The row permutation determines the column permutation, let's get it!
    # get_hessenberg_order also asserts non-increasing envelope, among others
    colp = get_hessenberg_order(g_orig, eqs, rowp)
    sink_set = { n for n in rowp if n not in matches }
    tear_set = { n for n in colp if n not in matches }
    #
    #print('Number of tears:', len(tear_set))
    #print('Row permutation:', rowp)
    #print('Col permutation:', colp)
    #
    return rowp, colp, matches, tear_set, sink_set


def setup_graphs(g_orig, forbidden):
    # g is a copy of g_orig; g_allowed contains only the allowed edges of g_orig
    g_pkl = cPickle_dumps(g_orig, cPickle_HIGHEST_PROTOCOL)
    g = cPickle_loads(g_pkl)
    g_allowed = cPickle_loads(g_pkl)
    adj = g_allowed.adj
    for u, v in forbidden:
        del adj[u][v]
        del adj[v][u] # assumes no self loops    
    return g_allowed, g


def create_heap(g_allowed, g, eqs):
    eq_tot  = heapdict()
    for e in sorted(eqs):
        tot = len(g[e])
        cost = tot-1 if g_allowed[e] else tot
        eq_tot[e]  = (cost, tot, e)
    return eq_tot


def run_tests():
    from test_tearing import gen_testproblems
    for g, eqs, forbidden in gen_testproblems():
        _, _, tears, sinks, _, _ = to_hessenberg_form(g, eqs, forbidden)
        print('Tears:', tears)
        print('Residuals:', sinks)


if __name__=='__main__':
    run_tests()
