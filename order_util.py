# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from six import itervalues
from networkx import Graph, NetworkXUnfeasible, relabel_nodes
from networkx.algorithms.bipartite import is_bipartite_node_set
from py3compat import irange
from utils import izip, pairwise

#-------------------------------------------------------------------------------
# The input is assumed to be correct and is not checked!

def to_bipartite(rows, cols_rowwise):
    'Returns: (g, eqs). Assumes disjoint row and column identifier sets.'
    edge_list = ((r, c) for r, cols in izip(rows, cols_rowwise) for c in cols)
    g = Graph(edge_list)
    g.add_nodes_from(rows) # Empty rows are allowed (but not empty columns)
    eqs = set(rows)
    assert is_bipartite_node_set(g, eqs)
    return g, eqs

def to_bipartite_w_forbidden(rows, cols_rowwise, vals_rowwise):
    '''Returns: (g, eqs, forbidden). Assumes disjoint row and column 
    identifier sets.'''
    g, forbidden = Graph(), set() 
    g.add_nodes_from(rows) # Empty rows are allowed (but not empty columns)
    for r, cols, vals in izip(rows, cols_rowwise, vals_rowwise):
        for c, v in izip(cols, vals):
            g.add_edge(r, c)
            if v > 1:
                forbidden.add((r, c))
    assert is_bipartite_node_set(g, rows)
    return g, set(rows), forbidden

def to_bipart_w_weights(cols_rowwise, vals_rowwise):
    '''Returns the tuple of: g, eqs, mapping (a list) to undo the row 
    permutation by weight, and the row weights in the same order as in the 
    input. This function does not receive the row identifiers but makes up
    new ones: 0, 1, ..., n_rows-1.'''
    n_rows = len(cols_rowwise)
    rows = list(irange(n_rows))
    row_weights = [sum(vals, 0.0) for vals in vals_rowwise]
    row_pos = argsort(row_weights)
    #print('Row weights: ', row_weights)
    #print('Row position:', row_pos)
    g = Graph()
    g.add_nodes_from(rows) # Empty rows are allowed (but not empty columns)
    # Apply row permutation row_pos
    edges = ((i,c) for i, r in enumerate(row_pos) for c in cols_rowwise[r])
    g.add_edges_from(edges)
    assert is_bipartite_node_set(g, rows) # Same ID for both a col and a row?    
    return g, set(rows), row_pos, row_weights

def argsort(seq, reverse=False):
    return sorted(range(len(seq)), key=seq.__getitem__, reverse=reverse)

#-------------------------------------------------------------------------------

def get_row_weights(g, n_rows):
    return [sum(d['weight'] for d in itervalues(g[r])) for r in irange(n_rows)]

def partial_relabel(g, mapping):
    return relabel_nodes(g, mapping, copy=True)

#-------------------------------------------------------------------------------

def bipartite_from_empty_matrix(shape):
    n_rows, n_cols = shape
    assert n_rows >= 0 and n_cols >= 0
    r_nodes = list(irange(n_rows))
    cols = [ ]    
    g = Graph()
    g.add_nodes_from(r_nodes)
    g.add_nodes_from(irange(n_rows, n_rows+n_cols))
    assert len(g) == n_rows + n_cols
    return g, set(r_nodes), cols


def coo_matrix_to_bipartite(rows, cols, values, shape):
    check_coordinate_format(rows, cols, values, shape)
    # See _check_coordinate_format in rpc_api too
    n_rows, n_cols = shape
    r_nodes = list(irange(n_rows))
    # relabel the columns, the caller must undo it later
    cols = [c+n_rows for c in cols]    
    g = Graph()
    g.add_nodes_from(r_nodes)
    g.add_nodes_from(irange(n_rows, n_rows+n_cols))
    assert len(g) == n_rows + n_cols
    g.add_edges_from(izip(rows, cols, ({'weight': int(v)} for v in values)))
    return g, set(r_nodes), cols


def check_coordinate_format(rows, cols, values, shape):
    n_rows, n_cols = shape
    assert n_rows >= 1 and n_cols >= 1, 'At least one row and one column are expected'
    assert rows and cols and values, 'There should be at least one entry'
    assert len(rows) == len(cols) and len(rows) == len(values), 'Array length mismatch'    
    check_if_indices_are_in_range(rows, shape[0], 'row')
    check_if_indices_are_in_range(cols, shape[1], 'column')    


def check_if_indices_are_in_range(idx, n_elems, elem_type, empty_allowed=False):
    if empty_allowed and not idx:
        return
    min_idx = min(idx)
    assert min_idx >= 0, 'The smallest {} index is {}'.format(elem_type, min_idx)
    max_idx = max(idx)
    assert max_idx < n_elems, 'The largest {} index is {}'.format(elem_type, max_idx)


def get_inverse_perm(rperm, cperm):
    rowp = [-1]*len(rperm)
    colp = [-1]*len(cperm)
    for i, r in enumerate(rperm):
        rowp[r] = i
    for i, c in enumerate(cperm):
        colp[c] = i
    assert rowp.count(-1) == 0
    assert colp.count(-1) == 0
    return rowp, colp

#-------------------------------------------------------------------------------
# Helper functions for the bipartite case, without any block structure

def get_row_col_perm(eqs, dag, tears, sinks, order):
    '''Returns the row and the column identifiers in permuted order for the 
    spiked form, given by the bipartite dag.'''
    gen_cols = (n for n in order if n not in eqs)
    indexof = { name : i for i, name in enumerate(gen_cols) }
    tear_stack = list(tears)
    tear_stack.sort(key=lambda v: indexof[v], reverse=True)
    resids = set(sinks)
    rows = [ n for n in order if n in eqs ]
    cols = [ ]
    for n in rows: # we walk along the envelope
        if n not in resids:
            # eq must have exactly one out edge in a valid elimination order
            (var,) = dag[n]
        else:
            var = tear_stack.pop() # Throws if structurally singular
        cols.append(var)
    assert not tear_stack
    assert len(rows) == len(eqs)
    assert len(cols) == len(dag) - len(eqs)
    return (rows, cols)

def colp_to_spiked_form(rowp, colp_hess, matches, tear_set, sink_set):
    '''Returns the column permutation vector to spiked form, given by the 
    bipartite matching and row permutation vector rowp.'''
    assert len(rowp) == len(colp_hess)
    tear_stack = list(n for n in reversed(colp_hess) if n in tear_set)
    colp = [ ]
    for n in rowp:
        if n not in sink_set:
            var = matches[n]
            assert var not in tear_set
        else:
            var = tear_stack.pop() # Throws if structurally singular
        colp.append(var)
    assert not tear_stack
    assert len(colp) == len(colp_hess)
    return colp

def check_spiked_form(g_orig, rowp, colp, tear_set):
    assert len(rowp) == len(colp)
    r_index = { name : i for i, name in enumerate(n for n in rowp) }
    first_elem = [ min(r_index[r] for r in g_orig[c]) for c in colp ]
    for i, c in enumerate(colp):
        first_nonzero = first_elem[i]
        if c in tear_set:
            assert first_nonzero <= i # on or above the diagonal
        else:
            assert first_nonzero == i # on the diagonal

def get_hessenberg_order(g, eqs, rowp, matches):
    colp = build_colp(g, rowp, matches)
    # append all isolated columns at the back
    isolated_cols = sorted(n for n in g if n not in eqs and len(g[n])==0)
    colp.extend(isolated_cols)
    assert len(rowp) == len(eqs)
    assert len(colp) == len(g) - len(eqs)
    check_nonincreasing_envelope(g, rowp, colp)
    return colp

def build_colp(g, rowp, matches):
    colp, seen_cols = [ ], set()
    adj = g.adj
    for r in rowp:
        cols = set(adj[r]) - seen_cols
        if cols:
            to_append = sorted(cols)
            colp.extend(to_append)
            seen_cols.update(to_append)
    return colp

def check_nonincreasing_envelope(g, rowp, colp):
    c_index = { name : i for i, name in enumerate(n for n in colp) }
    r_index = { name : i for i, name in enumerate(n for n in rowp) }
    # Last occupied columns rowwise, empty rows allowed
    adj = g.adj
    last_elem  = [max(c_index[c] for c in adj[r]) if adj[r] else -1 for r in rowp]
    c_viol = non_monotone_indices(last_elem)
    #if c_viol:
    #    from plot_ordering import plot_bipartite
    #    plot_bipartite(g, set(), rowp, colp)
    assert not c_viol, 'Non-monotone last elements in rows:  {}'.format(c_viol)
    # First occupied rows columnwise
    n_rows = len(rowp)
    first_elem = [min(r_index[r] for r in adj[c]) if adj[c] else n_rows for c in colp]
    r_viol = non_monotone_indices(first_elem)
    assert not r_viol, 'Non-monotone first elements in cols: {}'.format(r_viol)

def non_monotone_indices(lst):
    return [i for i, (u, v) in enumerate(pairwise(lst)) if u > v]

#-------------------------------------------------------------------------------
# Compare with dag_util in SDOPT

def deterministic_topological_sort(dag, nbunch=None):
    # This function is stolen from networkx.algorithms.dag.topological_sort.
    # Made the returned order deterministic by pre-sorting the nodes by node ID.
    assert dag.is_directed()
    seen = set()
    order = []
    explored = set()
    if nbunch is None:
        nbunch = sorted(dag.nodes_iter())                           # <-- SORTED
    for v in nbunch:     # process all vertices in dag
        if v in explored:
            continue
        fringe = [v]   # nodes yet to look at
        while fringe:
            w = fringe[-1]  # depth first search
            if w in explored: # already looked down this branch
                fringe.pop()
                continue
            seen.add(w)     # mark as seen
            # Check successors for cycles and for new nodes
            new_nodes = []
            for n in sorted(dag[w]):                                # <-- SORTED
                if n not in explored:
                    if n in seen: #CYCLE !!
                        raise NetworkXUnfeasible("Graph contains a cycle.")
                    new_nodes.append(n)
            if new_nodes:   # Add new_nodes to fringe
                fringe.extend(new_nodes)
            else:           # No new nodes so w is fully explored
                explored.add(w)
                order.append(w)
                fringe.pop()    # done considering this node
    return list(reversed(order))
