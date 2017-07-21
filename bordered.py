# Copyright (C) 2017 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from itertools import chain
from heap_md import min_degree
from order_util import coo_matrix_to_bipartite, get_inverse_perm


def to_bordered_form(rows, cols, values, shape):
    # values must be integers
    # cols are shifted by n_rows, must undo later
    g, eqs, cols = coo_matrix_to_bipartite(rows, cols, values, shape)
    forbidden = {(i, j) for i, j, v in zip(rows, cols, values) if v <= 0}
    assert all(g.has_edge(*e) for e in forbidden)
    rowp, colp, size = to_bordered(g, eqs, forbidden)
    # undo the col shift
    n_rows = shape[0]
    colp = [c-n_rows for c in colp]
    # find the inverse permutation
    rowp, colp = get_inverse_perm(rowp, colp)
    return rowp, colp, size

def to_bordered(g_orig, eqs, forbidden):
    # TODO Matching preference based on weights, then original order?
    rowp, colp, matches, tear_set, sink_set = min_degree(g_orig, eqs, forbidden)
    # get the ordering on the lower triangular part (lt_*)
    lt_row_order = [n for n in rowp if n in matches]
    lt_col_order = [n for n in colp if n in matches]
    assert len(lt_row_order) == len(lt_col_order)
    # TODO Order the border in non-decreasing order, the sinks as in rowp?
    #      Column border: break ties according to original order?
    #      Figure out a better node order: the order in the coo matrix?
    node_order = {n: i for i, n in enumerate(chain(rowp, colp))}
    def key(n):
        return node_order[n]
    row_border = sorted(sink_set, key=key)
    col_border = sorted(tear_set, key=key)
    # The row / col permutation that we are returning has the lower triangular
    # part first, and then the row / col border.
    res_rowp, res_colp = lt_row_order + row_border, lt_col_order + col_border
    assert sorted(rowp) == sorted(res_rowp)
    assert sorted(colp) == sorted(res_colp)
    return res_rowp, res_colp, len(lt_row_order) # size of the LT part

def run_tests():
    from test_tearing import gen_testproblems
    for g, eqs, forbidden in gen_testproblems():
        to_bordered(g, eqs, forbidden)

if __name__=='__main__':
    run_tests()
