# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from itertools import chain, product
import six
import networkx as nx
from networkx.algorithms.bipartite import projected_graph
from py3compat import izip, irange

from order_util import deterministic_topological_sort, coo_matrix_to_bipartite,\
                       get_inverse_perm, check_coordinate_format,\
                       bipartite_from_empty_matrix, check_if_indices_are_in_range

from test_utils import to_bipartite_from_test_string

def test_rowwise():
    print('DM decomposition of sparse matrices, stored row-wise')
    for name, mat_str in sorted(six.iteritems(TEST_CASES)):
        print(name)
        g, eqs = to_bipartite_from_test_string(mat_str)
        dm_decomp(g, eqs, upper=False)


def test_coo():
    print('Sparse matrices in coordinate form, DM decomposition')
    problems = sorted(six.iteritems(COORDINATE_TEST_CASES))
    for upper, (name, mat_str) in product((True, False), problems):
        print('---------------------------------------------------------------')
        print(name)
        rows, cols, values, shape = to_sparse_mat(mat_str)
        dm_coordinate(rows, cols, values, shape, upper, minimize=upper, show=True, name=name)
    print('Done!')


def test_blt_with_tearing():
    print('BLT with tearing')
    problems = sorted(six.iteritems(COORDINATE_TEST_CASES))
    for name, mat_str in problems:
        print('---------------------------------------------------------------')
        print(name)
        rows, cols, values, shape = to_sparse_mat(mat_str)
        blt_with_tearing(rows, cols, values, shape, [0,1], [2,3], show=True, name=name)
    print('Done!')

#-------------------------------------------------------------------------------

def _dbg_show_input(rows_orig, cols_orig, values_orig, shape, torn_rows, torn_cols):
    print('-------------------------------------------------------------------')
    print('Input matrix')
    dbg_show_coomat(rows_orig, cols_orig, ['x']*len(values_orig), shape)
    print('Torn rows:', torn_rows)
    print('Torn cols:', torn_cols)
    print()

def _dbg_show_submatrix(rows, cols, values, shape_sub):
    print('After removing the torn rows and columns:')
    dbg_show_coomat(rows, cols, ['x']*len(values), shape_sub)


def blt_with_tearing(rows_orig, cols_orig, values_orig, shape, torn_rows, torn_cols, \
                     show=False, name=''):
    #_dbg_show_input(rows_orig, cols_orig, values_orig, shape, torn_rows, torn_cols)
    n_rows, n_cols = shape
    check_if_indices_are_in_range(torn_rows, n_rows, 'torn row', empty_allowed=True)
    check_if_indices_are_in_range(torn_cols, n_cols, 'torn column', empty_allowed=True)
    # Filter out the not torn rows and columns (the active submatrix)
    rows_to_remove, cols_to_remove = set(torn_rows), set(torn_cols)
    rows, cols, values = [ ], [ ], [ ]
    for r, c, v in izip(rows_orig, cols_orig, values_orig):
        if r not in rows_to_remove and c not in cols_to_remove:
            rows.append(r)
            cols.append(c)
            values.append(v)
    # Indices in the submatrix, also the inverse maps for undoing relabel()
    row_idx  = [r for r in irange(n_rows) if r not in rows_to_remove]
    col_idx  = [c for c in irange(n_cols) if c not in cols_to_remove]
    # Relabel the rows and cols in the active submatrix: 
    # The smallest index is 0 and the indices are continuous.
    rows, cols = relabel(rows, cols, row_idx, col_idx)
    shape_sub = (len(row_idx), len(col_idx))
    #_dbg_show_submatrix(rows, cols, values, shape_sub)
    # Do the BLT decomposition on the relabeled submatrix 
    tup = dm_coordinate(rows, cols, values, shape_sub, upper=False, minimize=True)
    _, _, _,  R0, C0, rdiag, cdiag, rowpart, colpart, _ = tup
    # Undo relabel()
    R0 = [row_idx[r] for r in R0]
    C0 = [col_idx[c] for c in C0]
    rdiag = [row_idx[r] for r in rdiag]
    cdiag = [col_idx[c] for c in cdiag]
    rowpart = undo_relabel(rowpart, row_idx)
    colpart = undo_relabel(colpart, col_idx)
    #
    rperm = get_IDs_in_permuted_order(rdiag, R0, torn_rows)
    cperm = get_IDs_in_permuted_order(cdiag, C0, torn_cols)
    rowp, colp = get_inverse_perm(rperm, cperm)
    #
    gray_rows = chain(R0, torn_rows)
    gray_cols = chain(C0, torn_cols)
    colors = get_color_groups(rows_orig, cols_orig, gray_rows, gray_cols, rdiag, cdiag)
    assert sorted(rowp) == list(irange(n_rows))
    assert sorted(colp) == list(irange(n_cols))
    assert len(colors) == len(values_orig)
    #
    if show:
        pretty_print(R0, C0, rowpart, colpart)
        show_plot(name, rows_orig, cols_orig, colors, rowp, colp, rowpart, colpart)
    return rowp, colp,  R0, C0, rdiag, cdiag, rowpart, colpart, colors


def relabel(rows, cols, row_idx, col_idx):
    # This function should be moved to ordering_util?
    rmap = {r: i for i, r in enumerate(row_idx)}
    cmap = {c: i for i, c in enumerate(col_idx)}
    rows = [rmap[r] for r in rows]
    cols = [cmap[c] for c in cols]
    return rows, cols


def undo_relabel(part, idx_map):
    if isinstance(part, list):
        return [ undo_relabel(child, idx_map) for child in part ]
    else:
        return idx_map[part]


def get_IDs_in_permuted_order(diag, unmatched, torn):
    perm = list(diag)
    perm.extend(unmatched)
    perm.extend(torn)
    return perm

#-------------------------------------------------------------------------------

def dm_coordinate(rows, cols, values, shape, upper, minimize=False, show=False, name=''):
    # rows, cols, values: a sparse matrix in coordinate format
    # shape: (n_rows, n_cols)
    # The corresponding bipartite graph is relabeled, the col IDs start at 
    # n_rows instead of 0.  
    # rowpart, colpart: row and col partition, a tree given as nested lists,
    #      the unmatched entries are inserted just before returning
    # rperm, cperm: the indices in permuted order
    # rowp, colp: permutation vectors, putting the input in permuted order 
    # rdiag, cdiag: the diagonal entries in the SCCs, marked red on the plot
    # colors: color groups, 1: black, 2: red, 3: gray
    #---------------------------------------------------------------------------
    if values: # non-empty submatrix
        # Re-weight to get a minimum or maximum weight diagonal
        values = re_weight(values, minimize)
        assert all(v==int(v) and v >= 1 for v in values)
        # The col IDs in cols are shifted by n_rows, must undo later
        g, eqs, cols = coo_matrix_to_bipartite(rows, cols, values, shape)
    else:
        g, eqs, cols = bipartite_from_empty_matrix(shape)
    R0, C0, rowpart, colpart, rdiag, cdiag = dm_decomp(g, eqs, upper=upper) 
    rperm, cperm = get_permIDs(R0, C0, rdiag, cdiag, upper)
    # Undo shift
    n_rows = shape[0]
    C0 = [c-n_rows for c in C0]
    cperm = [c-n_rows for c in cperm]
    cdiag = [c-n_rows for c in cdiag]
    cols =  [c-n_rows for c in cols ]
    #pretty_indent(colpart)
    colpart = undo_shift(colpart, n_rows)
    #pretty_indent(colpart)
    # Build the permutation vectors that permute the *input* into permuted form
    rowp, colp = get_inverse_perm(rperm, cperm)
    assert len(rowp) == n_rows
    assert len(colp) == shape[1]
    colors = get_color_groups(rows, cols, R0, C0, rdiag, cdiag)
    if show:
        pretty_print(R0, C0, rowpart, colpart)
        show_plot(name, rows, cols, colors, rowp, colp, rowpart, colpart)
    # Careful: the g has shifted IDs for the cols!
    return g, rowp, colp,  R0, C0, rdiag, cdiag, rowpart, colpart, colors


def sccs_on_plot(rowpart, colpart, rowp, colp):
    # SCC squares: list of the tuple (r index, c index, size)
    scc_squares = [ ]
    for rpart, cpart in izip(rowpart, colpart):
        for rscc, cscc in izip(rpart, cpart):
            r, c = rowp[rscc[0]], colp[cscc[0]]
            assert len(rscc) == len(cscc)
            scc_squares.append((r, c, len(rscc)))
    return scc_squares


def show_plot(name, rows, cols, colors, rowp, colp, rowpart, colpart):
    from plot_ordering import plot_dm
    scc_squares = sccs_on_plot(rowpart, colpart, rowp, colp)
    plot_dm(name, rows, cols, rowp, colp, colors, scc_squares)


def re_weight(values, minimize):
    if minimize:
        big_M = max(values) + 1
        return [big_M-v for v in values]
    else:
        shift = 1 - min(values)
        return [v + shift for v in values]


def undo_shift(part, n_rows):
    if isinstance(part, list):
        return [ undo_shift(child, n_rows) for child in part ]
    else:
        return part-n_rows


def get_color_groups(rows, cols, gen_gray_rows, gen_gray_cols, rdiag, cdiag):
    # Put the diagonal entries of the SCCs into color group 2 (red)
    diags = set(izip(rdiag, cdiag))
    # Put the unmatched / torn entries in color group 3 (gray), 
    # and all others in 1 (black)
    gray_rows = set(gen_gray_rows)
    gray_cols = set(gen_gray_cols)
    def color_group(r_c):
        if r_c in diags:
            return 2
        if r_c[0] in gray_rows or r_c[1] in gray_cols:
            return 3
        return 1
    #
    return [color_group(r_c) for r_c in izip(rows, cols)]    


def get_permIDs(R0, C0, rdiag, cdiag, upper):
    # Row and column IDs in permuted order
    if upper:
        rperm = list(rdiag)
        rperm.extend(R0)
        cperm = list(C0)
        cperm.extend(cdiag)
    else:
        rperm = list(R0)
        rperm.extend(rdiag)
        cperm = list(cdiag)
        cperm.extend(C0)
    return rperm, cperm

#-------------------------------------------------------------------------------

def dm_decomp(g, eqs, upper=False):
    coarse_tup = coarse_DM(g, eqs)
    return fine_DM(*coarse_tup, upper=upper)


def coarse_DM(g, eqs):
    #from matching import max_weight_matching
    mate = nx.max_weight_matching(g, maxcardinality=True)
    bipart = orient(g, eqs, mate)
    # unmatched rows and columns
    R0 = sorted(n for n in g if n     in eqs and n not in mate)
    C0 = sorted(n for n in g if n not in eqs and n not in mate)
    # See Timothy A. Davis: Direct Methods for Sparse Linear Systems
    #                       7.4 Dulmage-Mendelsohn decomposition, p. 122.
    C3 = sorted(reachable_in_other_nodeset(bipart.pred, R0))
    R1 = sorted(reachable_in_other_nodeset(bipart.succ, C0))
    C1 = sorted(mate[r] for r in R1)
    R3 = sorted(mate[c] for c in C3)
    R2complement = set(chain(R0, R1, R3))
    C2complement = set(chain(C0, C1, C3))    
    R2 = sorted(n for n in bipart if n     in eqs and n not in R2complement)
    C2 = sorted(n for n in bipart if n not in eqs and n not in C2complement)
    check_coarse_decomp(g, eqs, mate, R0, R1, R2, R3, C0, C1, C2, C3)
    #print('[A11 A12] = ', (R1, (C0, C1)))
    #print('[A23] = ', (R2, C2))
    #print('[A34; A44] = ', ((R3, R0), C3))
    return mate, bipart, R0, R1, R2, R3, C0, C1, C2, C3


def fine_DM(mate, bipart, R0, R1, R2, R3, C0, C1, C2, C3, upper=False):
    rowpart, colpart = [ ], [ ]
    if upper:
        coarse_blocks_in_order = ((R1,C1), (R2,C2), (R3,C3))
    else:
        coarse_blocks_in_order = ((R3,C3), (R2,C2), (R1,C1))
    #
    for rowidx, colidx in coarse_blocks_in_order:
        rpart, cpart = partition_by_sccs(bipart, mate, rowidx, colidx, upper)
        rowpart.append(rpart)
        colpart.append(cpart)
    rdiag, cdiag = get_diagonal(rowpart, colpart, mate)
    #pretty_print(R0, C0, rowpart, colpart)
    return R0, C0, rowpart, colpart, rdiag, cdiag


def orient(g, eqs, mate):
    eq_var_matches = {e: v for e, v in six.iteritems(mate) if e in eqs}
    assert len(eq_var_matches) == len(mate)//2
    bipart = nx.DiGraph()
    bipart.add_nodes_from(g)
    # orient all edges eq <- var (all tears)
    bipart.add_edges_from((var, eq) for eq, var in g.edges_iter(eqs))
    # reverse only the matched edges
    for eq, var in six.iteritems(eq_var_matches):
        bipart.remove_edge(var, eq)
        bipart.add_edge(eq, var)
    assert len(bipart) == len(g)
    assert bipart.number_of_edges() == g.number_of_edges()
    return bipart


def reachable_in_other_nodeset(g_nbrs, sources):
    # Returns: nodes reachable from sources in the other node set of g. 
    # Only for bipartite graphs!
    visited = set(sources)
    reachable = set()
    level = 0
    nextlevel = set(chain.from_iterable(g_nbrs[s] for s in sources))    
    while nextlevel:
        thislevel = nextlevel
        nextlevel = set(chain.from_iterable(g_nbrs[n] for n in thislevel))
        level += 1
        if level % 2:
            reachable.update(thislevel)
        visited.update(thislevel)
        nextlevel.difference_update(visited)
    return reachable


def check_coarse_decomp(g, eqs, mate, R0, R1, R2, R3, C0, C1, C2, C3):
    all_rows = set(chain(R0, R1, R2, R3))
    all_cols = set(chain(C0, C1, C2, C3))
    assert eqs == all_rows
    assert not all_rows & all_cols
    assert len(g) == len(all_rows) + len(all_cols)
    assert len(all_rows) == len(R0) + len(R1) + len(R2) + len(R3)
    assert len(all_cols) == len(C0) + len(C1) + len(C2) + len(C3)
    assert len(R1) == len(C1)
    assert len(R2) == len(C2)
    assert len(R3) == len(C3)
    C2set = set(C2)
    R2set = set(R2)
    assert all(mate[r] in C2set for r in R2)
    assert all(mate[c] in R2set for c in C2)
    assert len(R1) + len(R2) + len(R3) == len(mate)//2
    assert len(C1) + len(C2) + len(C3) == len(mate)//2


def partition_by_sccs(bipart, mate, rowindices, colindices, upper):
    if not rowindices and not colindices:
        return [ ], [ ]  # otherwise condensation crashes
    sub_g = bipart.subgraph(chain(rowindices, colindices))
    eq_sccs = nx.condensation(projected_graph(sub_g, rowindices))
    precedence_order = deterministic_topological_sort(eq_sccs)
    rowpart, colpart = [ ], [ ]
    if upper:
        precedence_order = reversed(precedence_order) 
    for scc in precedence_order:
        equations = sorted(eq_sccs.node[scc]['members'])
        variables = [mate[eq] for eq in equations]    
        rowpart.append(equations)
        colpart.append(variables)
    return rowpart, colpart


def get_diagonal(rowpart, colpart, mate):
    rdiag, cdiag = [ ], [ ]
    for rpart, cpart in izip(rowpart, colpart):
        rdiag.extend(chain.from_iterable(rpart))
        cdiag.extend(chain.from_iterable(cpart))

    match_len = len(mate)//2
    assert len(rdiag) == match_len
    assert len(cdiag) == match_len
    for i in irange(match_len):
        assert rdiag[i] == mate[cdiag[i]]
        assert cdiag[i] == mate[rdiag[i]]
        
    return rdiag, cdiag

#-------------------------------------------------------------------------------

def dbg_show_coomat(rows, cols, values, shape, rowp=None, colp=None):
    
    n_rows, n_cols = shape
    print('Shape: {}x{}'.format(n_rows, n_cols))
    
    row_padding  = 2 if n_rows > 10 else 1
    #col_padding  = 2 if n_cols > 10 else 1
    row_fmt  =  '%2s: |' if n_rows > 10 else  '%s: |'
    col_fmt  = ' %2s'    if n_cols > 10 else ' %s'
    zero_padding = '  ' if n_cols > 10 else ' '
    
    if rowp is None:
        rowp = list(irange(n_rows))
    
    if colp is None:
        colp = list(irange(n_cols))
    
    dense_mat = [[zero_padding]*n_cols for _ in irange(n_rows)]
    for r, c, v in izip(rows, cols, values):
        dense_mat[rowp[r]][colp[c]] = v
  
    # print the col identifiers on the top
    print(' '* (row_padding + len(': |')), end='')
    for c in irange(n_cols):
        print(col_fmt % c, end='')
    print()
    
    for i in irange(n_rows):
        print(row_fmt % i, end='')
        for c in dense_mat[i]:
            print(col_fmt % c, end='')
        print(' |')
    print()


def pretty_print(R0, C0, rpartition, cpartition):
    print('Unmatched rows:', R0)
    print('Unmatched cols:', C0)
    for i, (rowparts, colparts) in enumerate(izip(rpartition, cpartition), start=1):
        print('Diagonal block', i)
        for j, (rpart, cpart) in enumerate(izip(rowparts, colparts), start=1):
            print('  SCC', j, '[', end='')
            for r, c in izip(rpart, cpart):
                print('', r, '->', c, '', end='')
            print(']')
        print()
    #
    print('rpartition:', rpartition)
    print('cpartition:', cpartition)
    #---------------------
    #print_node(rpartition, 0)
    #print()
    #---------------------
    #pretty_indent(rpartition)


def print_node(node, level):
    print(' {}[ '.format(level), end='')
    for child in node:
        if isinstance(child, list):
            print_node(child, level+1)
        else:
            print(child, '', end='')
    print(']{} '.format(level), end='')


PRE_ORDER, IN_ORDER, POST_ORDER = 'PRE_ORDER', 'IN_ORDER', 'POST_ORDER' 

def depth_first_traversal(node, level):
    if isinstance(node, list):
        yield PRE_ORDER, level, None
        for child in node:
            for move, depth, value in depth_first_traversal(child, level+1):
                yield move, depth, value
        yield POST_ORDER, level, None
    else:
        yield IN_ORDER, level, node

def pretty_indent(tree):
    for event, level, value in depth_first_traversal(tree, level=0):
        if event is PRE_ORDER:
            print(' {}[ '.format(level), end='')
        elif event is POST_ORDER:
            print(']{} '.format(level), end='')
        else:
            print(value, '', end='')
    print()
    for event, level, value in depth_first_traversal(tree, level=0):
        padding = ' '*level*2 
        if event is PRE_ORDER:
            print(padding, '{}[ '.format(level), sep='')
        elif event is POST_ORDER:
            print(padding, ']{} '.format(level), sep='')
        else:
            print(padding, value, sep='')
    print()

#-------------------------------------------------------------------------------
# See test_utils.to_bipartite_from_test_string to create the corresponding graph

TEST_CASES = {
              
    'test_1' : (    'e1 e2 e3 e4',
                  '''v1 v2 v3 v4
                           v3 v4
                              v4
                              v4''' ),
    
    'test_2' : (   'e1  e2  e3  e4',
                 '''x1  x2  x3  x4
                    x1  x2  x3  x4
                            x3  x4
                            x3  x4'''),
    
    'test_3' : (  'e01  e02  e03  e04  e05  e06  e07  e08  e09  e10  e11  e12  e13  e14',
                '''x01  x02  x03  x04  x05  x06
                   x01  x02  x03  x04  x05  x06
                   x01  x02            x05  x06
                   x01  x02            x05  x06
                                                 x07
                                                      x08  x09  x10
                                                      x08  x09  x10
                                                      x08  x09  x10
                                                                     x11  x12  x13  x14
                                                                     x11  x12  x13  x14
                                                                               x13  x14
                                                                               x13  x14
                                                                     x11  x12  x13  x14
                                                                     x11  x12  x13  x14''')
}

#-------------------------------------------------------------------------------

def to_sparse_mat(mat_str):
    rows = mat_str[0].split()
    cols = mat_str[1].split()
    rows = [int(r) for r in rows]
    cols = [int(c) for c in cols]
    if len(mat_str) == 2:
        vals = [1]*len(rows)
    else:
        vals = mat_str[2].split()
        vals = [int(v) for v in vals]
    shape = ( len(set(rows)), len(set(cols)) )
    msg = check_coordinate_format(rows, cols, vals, shape)
    assert not msg, msg    
    #---
    ## A hack to apply random permutation
    #from random import Random
    #rng = Random(3) # 7
    #mapping = list(irange(shape[0]))
    #rng.shuffle(mapping)
    #print(mapping)
    #rows = [mapping[r] for r in rows]
    #rng.shuffle(mapping)
    #cols = [mapping[c] for c in cols]
    #_dbg_show_input(rows, cols, vals, shape, [], []) 
    #print(rows)
    #print(cols)
    #---
    return rows, cols, vals, shape

COORDINATE_TEST_CASES = {
                         
    'test_1': ('''0 0 0 0
                      1 1
                        2
                        3''',
                         
               '''0 1 2 3
                      2 3
                        3
                        3'''),
 
    'test_2': ('''0 0 0 0
                  1 1 1 1
                      2 2
                      3 3''',
                         
               '''0 1 2 3
                  0 1 2 3
                      2 3
                      2 3'''),
                          
    'test_3' : ('''0  0  0  0  0  0
                   1  1  1  1  1  1
                   2  2        2  2
                   3  3        3  3
                                     4
                                        5  5  5
                                        6  6  6
                                        7  7  7
                                                 8  8  8  8
                                                 9  9  9  9
                                                      10 10
                                                      11 11
                                                12 12 12 12
                                                13 13 13 13''',
                                                 
                '''0  1  2  3  4  5
                   0  1  2  3  4  5
                   0  1        4  5
                   0  1        4  5
                                     6
                                        7  8  9
                                        7  8  9
                                        7  8  9
                                                10 11 12 13
                                                10 11 12 13
                                                      12 13
                                                      12 13
                                                10 11 12 13
                                                10 11 12 13'''),
 
    'test_4' : ('''0 0 0 0
                   1 1 1 1
                   2 2 2 2''',
                    
                '''0 1 2 3
                   0 1 2 3
                   0 1 2 3''',
                    
                '''9 1 1 6 
                   3 8 1 2
                   5 4 1 7'''),
                          
     'test_5' : ('''0  0  0
                       1  1
                             2
                                3  3  3
                                4  4  4
                                5  5  5
                                          6  6
                                          7  7
                                          8  8''',
                     
                 '''0  1  2
                       1  2
                             3
                                4  5  6
                                4  5  6
                                4  5  6
                                          7  8
                                          7  8
                                          7  8'''),
                         
     'test_6' : ('''0  0  0        0
                       1  1  1               1
                             2
                                3     3
                                4  4      4
                                5  5  5
                                          6  6
                                          7  7
                                          8  8''',
                    
                 '''0  1  2        5
                       1  2  3               8
                             3
                                4     6
                                4  5      7  
                                4  5  6
                                          7  8
                                          7  8
                                          7  8'''),
    
    'test_7' : ('''0  0
                    1  1  1
                       2  2  2
                          3  3  3
                             4  4  4
                                5  5''',
                    
                 '''0  1
                    0  1  2
                       1  2  3
                          2  3  4
                             3  4  5
                                4  5''',
                    
                 '''10   1
                     1  10   1
                         1  10   1
                             1  10   1
                                 1  10   1 
                                     1  10'''),

    'test_8' : ('''0  0
                    1  1  1
                       2  2  2
                          3  3  3
                             4  4  4
                                5  5  5
                                   6  6''',
                    
                 '''0  1
                    0  1  2
                       1  2  3
                          2  3  4
                             3  4  5
                                4  5  6
                                   5  6''',
                    
                 '''10   1
                     1  10   1
                         1  10   1
                             1  10   1
                                 1  10   1 
                                     1  10   1
                                         1  10'''),

    'test_9' : ('''0   0   0                            
                   1   1   1   1                        
                   2   2   2   2   2                    
                       3   3   3   3   3                
                           4   4   4   4   4            
                               5   5   5   5   5        
                                   6   6   6   6   6    
                                       7   7   7   7   7
                                           8   8   8   8
                                               9   9   9''',
                '''0   1   2                            
                   0   1   2   3                        
                   0   1   2   3   4                    
                       1   2   3   4   5                
                           2   3   4   5   6            
                               3   4   5   6   7        
                                   4   5   6   7   8    
                                       5   6   7   8   9
                                           6   7   8   9
                                               7   8   9''',
                '''10   1   1                            
                    1  10   1   1                        
                    1   1  10   1   1                    
                        1   1  10   1   1                
                            1   1  10   1   1            
                                1   1  10   1   1        
                                    1   1  10   1   1    
                                        1   1  10   1   1
                                            1   1  10   1
                                                1   1  10'''),

    'test_10' : (   '''0   0   0                                                                                                            
                       1   1   1   1                                                                                                        
                       2   2   2   2   2                                                                                                    
                           3   3   3   3   3                                                                                                
                               4   4   4   4   4                                                                                            
                                   5   5   5   5   5                                                                                        
                                       6   6   6   6   6                                                                                    
                                           7   7   7   7   7                                                                                
                                               8   8   8   8   8                                                                            
                                                   9   9   9   9   9                                                                        
                                                      10  10  10  10  10                                                                    
                                                          11  11  11  11  11                                                                
                                                              12  12  12  12  12                                                            
                                                                  13  13  13  13  13                                                        
                                                                      14  14  14  14  14                                                    
                                                                          15  15  15  15  15                                                
                                                                              16  16  16  16  16                                            
                                                                                  17  17  17  17  17                                        
                                                                                      18  18  18  18  18                                    
                                                                                          19  19  19  19  19                                
                                                                                              20  20  20  20  20                            
                                                                                                  21  21  21  21  21                        
                                                                                                      22  22  22  22  22                    
                                                                                                          23  23  23  23  23                
                                                                                                              24  24  24  24  24            
                                                                                                                  25  25  25  25  25        
                                                                                                                      26  26  26  26  26    
                                                                                                                          27  27  27  27  27
                                                                                                                              28  28  28  28
                                                                                                                                  29  29  29''',
                    '''0   1   2                                                                                                            
                       0   1   2   3                                                                                                        
                       0   1   2   3   4                                                                                                    
                           1   2   3   4   5                                                                                                
                               2   3   4   5   6                                                                                            
                                   3   4   5   6   7                                                                                        
                                       4   5   6   7   8                                                                                    
                                           5   6   7   8   9                                                                                
                                               6   7   8   9  10                                                                            
                                                   7   8   9  10  11                                                                        
                                                       8   9  10  11  12                                                                    
                                                           9  10  11  12  13                                                                
                                                              10  11  12  13  14                                                            
                                                                  11  12  13  14  15                                                        
                                                                      12  13  14  15  16                                                    
                                                                          13  14  15  16  17                                                
                                                                              14  15  16  17  18                                            
                                                                                  15  16  17  18  19                                        
                                                                                      16  17  18  19  20                                    
                                                                                          17  18  19  20  21                                
                                                                                              18  19  20  21  22                            
                                                                                                  19  20  21  22  23                        
                                                                                                      20  21  22  23  24                    
                                                                                                          21  22  23  24  25                
                                                                                                              22  23  24  25  26            
                                                                                                                  23  24  25  26  27        
                                                                                                                      24  25  26  27  28    
                                                                                                                          25  26  27  28  29
                                                                                                                              26  27  28  29
                                                                                                                                  27  28  29''',
                   '''10   1   1                                                                                                            
                       1  10   1   1                                                                                                        
                       1   1  10   1   1                                                                                                    
                           1   1  10   1   1                                                                                                
                               1   1  10   1   1                                                                                            
                                   1   1  10   1   1                                                                                        
                                       1   1  10   1   1                                                                                    
                                           1   1  10   1   1                                                                                
                                               1   1  10   1   1                                                                            
                                                   1   1  10   1   1                                                                        
                                                       1   1  10   1   1                                                                    
                                                           1   1  10   1   1                                                                
                                                               1   1  10   1   1                                                            
                                                                   1   1  10   1   1                                                        
                                                                       1   1  10   1   1                                                    
                                                                           1   1  10   1   1                                                
                                                                               1   1  10   1   1                                            
                                                                                   1   1  10   1   1                                        
                                                                                       1   1  10   1   1                                    
                                                                                           1   1  10   1   1                                
                                                                                               1   1  10   1   1                            
                                                                                                   1   1  10   1   1                        
                                                                                                       1   1  10   1   1                    
                                                                                                           1   1  10   1   1                
                                                                                                               1   1  10   1   1            
                                                                                                                   1   1  10   1   1        
                                                                                                                       1   1  10   1   1    
                                                                                                                           1   1  10   1   1
                                                                                                                               1   1  10   1
                                                                                                                                   1   1  10'''),

      'test_11' : ( '''0   0   0                                                                                                                
                       1   1   1   1                                                                                                            
                       2   2   2   2   2                                                                                                        
                           3   3   3   3   3                                                                                                    
                               4   4   4   4   4                                                                                                
                                   5   5   5   5   5                                                                                            
                                       6   6   6   6   6                                                                                        
                                           7   7   7   7   7                                                                                    
                                               8   8   8   8   8                                                                                
                                                   9   9   9   9   9                                                                            
                                                      10  10  10  10  10                                                                        
                                                          11  11  11  11  11                                                                    
                                                              12  12  12  12  12                                                                
                                                                  13  13  13  13  13                                                            
                                                                      14  14  14  14  14                                                        
                                                                          15  15  15  15  15                                                    
                                                                              16  16  16  16  16                                                
                                                                                  17  17  17  17  17                                            
                                                                                      18  18  18  18  18                                        
                                                                                          19  19  19  19  19                                    
                                                                                              20  20  20  20  20                                
                                                                                                  21  21  21  21  21                            
                                                                                                      22  22  22  22  22                        
                                                                                                          23  23  23  23  23                    
                                                                                                              24  24  24  24  24                
                                                                                                                  25  25  25  25  25            
                                                                                                                      26  26  26  26  26        
                                                                                                                          27  27  27  27  27    
                                                                                                                              28  28  28  28  28
                                                                                                                                  29  29  29  29
                                                                                                                                      30  30  30''',
                    
                    '''0   1   2                                                                                                                
                       0   1   2   3                                                                                                            
                       0   1   2   3   4                                                                                                        
                           1   2   3   4   5                                                                                                    
                               2   3   4   5   6                                                                                                
                                   3   4   5   6   7                                                                                            
                                       4   5   6   7   8                                                                                        
                                           5   6   7   8   9                                                                                    
                                               6   7   8   9  10                                                                                
                                                   7   8   9  10  11                                                                            
                                                       8   9  10  11  12                                                                        
                                                           9  10  11  12  13                                                                    
                                                              10  11  12  13  14                                                                
                                                                  11  12  13  14  15                                                            
                                                                      12  13  14  15  16                                                        
                                                                          13  14  15  16  17                                                    
                                                                              14  15  16  17  18                                                
                                                                                  15  16  17  18  19                                            
                                                                                      16  17  18  19  20                                        
                                                                                          17  18  19  20  21                                    
                                                                                              18  19  20  21  22                                
                                                                                                  19  20  21  22  23                            
                                                                                                      20  21  22  23  24                        
                                                                                                          21  22  23  24  25                    
                                                                                                              22  23  24  25  26                
                                                                                                                  23  24  25  26  27            
                                                                                                                      24  25  26  27  28        
                                                                                                                          25  26  27  28  29    
                                                                                                                              26  27  28  29  30
                                                                                                                                  27  28  29  30
                                                                                                                                      28  29  30''',
                    
                    
                    '''10   1   1                                                                                                                
                       1  10   1   1                                                                                                            
                       1   1  10   1   1                                                                                                        
                           1   1  10   1   1                                                                                                    
                               1   1  10   1   1                                                                                                
                                   1   1  10   1   1                                                                                            
                                       1   1  10   1   1                                                                                        
                                           1   1  10   1   1                                                                                    
                                               1   1  10   1   1                                                                                
                                                   1   1  10   1   1                                                                            
                                                       1   1  10   1   1                                                                        
                                                           1   1  10   1   1                                                                    
                                                               1   1  10   1   1                                                                
                                                                   1   1  10   1   1                                                            
                                                                       1   1  10   1   1                                                        
                                                                           1   1  10   1   1                                                    
                                                                               1   1  10   1   1                                                
                                                                                   1   1  10   1   1                                            
                                                                                       1   1  10   1   1                                        
                                                                                           1   1  10   1   1                                    
                                                                                               1   1  10   1   1                                
                                                                                                   1   1  10   1   1                            
                                                                                                       1   1  10   1   1                        
                                                                                                           1   1  10   1   1                    
                                                                                                               1   1  10   1   1                
                                                                                                                   1   1  10   1   1            
                                                                                                                       1   1  10   1   1        
                                                                                                                           1   1  10   1   1    
                                                                                                                               1   1  10   1   1
                                                                                                                                   1   1  10   1
                                                                                                                                       1   1  10'''),
                    
       'test_12': ('''0                          0
                         1                       1
                            2                    2
                               3                 3
                                  4              4
                                     5           5
                                        6        6
                                           7     7
                                              8  8
                      9  9  9  9  9  9  9  9  9  9''',
                    
                   '''0                          9
                         1                       9
                            2                    9
                               3                 9
                                  4              9
                                     5           9
                                        6        9
                                           7     9
                                              8  9
                      0  1  2  3  4  5  6  7  8  9''',
                    
                   '''1                          1
                         1                       1
                            1                    1
                               1                 1
                                  1              1
                                     1           1
                                        1        1
                                           1     1
                                              1  1
                     10  1  1  1  1  1  1  1  1  1'''),

    'test_14' :   ('''0                          0
                         1                       1
                            2                    2
                               3                 3
                                  4              4
                                     5           5
                                        6        6
                                           7     7
                                              8  8
                      9  9  9  9  9  9  9  9  9  9''',
                    
                   '''0                          9
                         1                       9
                            2                    9
                               3                 9
                                  4              9
                                     5           9
                                        6        9
                                           7     9
                                              8  9
                      0  1  2  3  4  5  6  7  8  9''',
                    
                   '''1                          1
                         1                       1
                            1                    1
                               1                 1
                                  1              1
                                     1           1
                                        1        1
                                           1     1
                                              1  1
                      1  1  1  1  1  1  1  1  1 10'''),
}


def run_tests():
    test_rowwise()
    test_coo()
    test_blt_with_tearing()    

if __name__ == '__main__':
    run_tests()
