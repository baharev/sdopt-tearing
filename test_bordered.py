# Copyright (C) 2017 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from random import Random
from bordered import to_bordered_form
from dm_decomp import COORDINATE_TEST_CASES, to_sparse_mat, dbg_show_coomat
from order_util import get_inverse_perm
from test_utils import create_coomat


def run():
    # test_5 and test_3 seem to be buggy
    #test_cases = {'test_5': COORDINATE_TEST_CASES['test_5']}
    test_cases = COORDINATE_TEST_CASES
    for name, mat_str in test_cases.items():
        print(name)
        rows, cols, values, shape = to_sparse_mat(mat_str)
        dbg_show_coomat(rows, cols, values, shape)
        rowp, colp, size = to_bordered_form(rows, cols, values, shape)
        dbg_show_coomat(rows, cols, values, shape, rowp, colp)
        rowp, colp = get_inverse_perm(rowp, colp)
        print('rows:', rowp)
        print('cols:', colp)
        print('LT size:', size)
        check_leading_lt_part(rows, cols, values, rowp, colp, size)
        print('---------------------------------------------------------------')

def run_property_based_tests():
    print('Started generative testing...')
    import os
    os.environ['HYPOTHESIS_STORAGE_DIRECTORY'] = '/tmp/ht'
    from hypothesis import settings
    with settings(max_examples=10, timeout=0):
        _run_tests()
    print('Done!')

def _run_tests():
    from hypothesis import given
    from hypothesis.strategies import integers
    
    MAX_VALUE = 8
    
    decor = given(n_rows = integers(min_value=1, max_value=MAX_VALUE),
                  n_cols = integers(min_value=1, max_value=MAX_VALUE), 
                  seed   = integers(min_value=0))

    decor(test_random_coo)()

def test_random_coo(n_rows, n_cols, seed):
    print('n_rows:', n_rows, 'n_cols:', n_cols, 'seed:', seed)
    rng = Random(seed)
    _g, rows, cols, values = create_coomat(n_rows, n_cols, rng)
    values = [v if rng.random() > 0.5 else -v for v in values]
    shape = (n_rows, n_cols)
    dbg_show_coomat(rows, cols, values, shape)
    rowp, colp, size = to_bordered_form(rows, cols, values, shape)
    dbg_show_coomat(rows, cols, values, shape, rowp, colp)
    rowp, colp = get_inverse_perm(rowp, colp)
    print('rows:', rowp)
    print('cols:', colp)
    print('LT size:', size)
    check_leading_lt_part(rows, cols, values, rowp, colp, size)
    print('---------------------------------------------------------------')
    
#-------------------------------------------------------------------------------

def check_leading_lt_part(rows, cols, values, rowp, colp, size):
    prows, pcols = permute(rows, cols, rowp, colp)
    rows, cols, values = get_lt_part(prows, pcols, values, size)
    print()
    print('Lower triangular part')
    dbg_show_coomat(rows, cols, values, (size, size))
    check_lt_with_pos_diagonal(rows, cols, values)

def permute(rows, cols, rowp, colp):
    rmap = {r: i for i, r in enumerate(rowp)}
    cmap = {c: i for i, c in enumerate(colp)}
    prows = [rmap[r] for r in rows]
    pcols = [cmap[c] for c in cols]
    return prows, pcols

def get_lt_part(prows, pcols, orig_values, size):
    rows, cols, values = [], [], []
    for r, c, v in zip(prows, pcols, orig_values):
        if r < size and c < size:
            rows.append(r)
            cols.append(c)
            values.append(v)
    return rows, cols, values

def check_lt_with_pos_diagonal(rows, cols, values):
    for r, c, v in zip(rows, cols, values):
        assert r >= c, (r, c)
        assert r != c or v > 0, (r, c, v)

if __name__ == '__main__':
    #run()
    run_property_based_tests()

