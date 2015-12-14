# Copyright (C) 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from random import shuffle
from itertools import chain
from dm_decomp import dm_coordinate, blt_with_tearing #, dbg_show_coomat
import heap_md as md
from order_util import get_inverse_perm, argsort,\
                       check_coordinate_format, check_if_indices_are_in_range
from py3compat import irange
from utils import duplicates


def fine_dulmage_mendelsohn(rows, cols, values, n_rows, n_cols, upper, minimize):
    shape = (n_rows, n_cols)
    error_msg = _check_coordinate_format(rows, cols, values, shape)
    if error_msg is not None:
        return {'error_msg': error_msg}
    tup = dm_coordinate(rows, cols, values, shape, upper, minimize=minimize)
    _, rowp, colp,  R0, C0, _, _, rowpart, colpart, colors = tup
    # add unmatched to the partition, per request
    if upper:
        rowpart.append(R0)
        colpart.insert(0, C0)
    else:
        rowpart.insert(0, R0)
        colpart.append(C0)
    #pretty_indent(rowpart)
    #pretty_indent(colpart)
    _check_partition(rowp, colp, rowpart, colpart)
    return _pack(rowp, colp, rowpart, colpart, colors)


def tearing_hand_guided(rows, cols, values, n_rows, n_cols, torn_rows, torn_cols):
    shape = (n_rows, n_cols)
    error_msg = _check_coordinate_format(rows, cols, values, shape) or \
                _check_torn_indices(shape, torn_rows, torn_cols)
    if error_msg is not None:
        return {'error_msg': error_msg}
    tup = blt_with_tearing(rows, cols, values, shape, torn_rows, torn_cols)
    rowp, colp, R0, C0, _, _, rowpart, colpart, colors = tup
    # Similarly to fine_dulmage_mendelsohn, we add the unmatched but also the 
    # torn rows/cols as "unmatched".
    rowpart.append(list(chain(R0, torn_rows)))
    colpart.append(list(chain(C0, torn_cols)))
    #print('Returning:')
    #dbg_show_coomat(rows, cols, colors, shape, rowp, colp)
    #print('Row partition:', rowpart)
    #print('Col partition:', colpart)
    _check_partition(rowp, colp, rowpart, colpart)
    return _pack(rowp, colp, rowpart, colpart, colors)


def _check_partition(rowp, colp, rowpart, colpart):
    dbg_rperm = list(_flatten(rowpart))
    dbg_cperm = list(_flatten(colpart))
    dbg_rowp, dbg_colp = get_inverse_perm(dbg_rperm, dbg_cperm)
    assert dbg_rowp == rowp
    assert dbg_colp == colp


def _flatten(partition):
    if isinstance(partition, list):
        for elem in partition:
            for leaf in _flatten(elem):
                yield leaf
    else:
        yield partition # leaf


def hessenberg(rows, cols, values, n_rows, n_cols, tie_breaking):
    shape = (n_rows, n_cols)
    error_msg = _check_coordinate_format(rows, cols, values, shape)
    if error_msg is not None:
        return {'error_msg': error_msg}
    rowp, colp = md.hessenberg(rows, cols, values, n_rows, n_cols, tie_breaking)
    return _pack(rowp, colp)


def lexicographical(row_ids, col_ids, sort_rows, sort_cols):
    dups = duplicates(row_ids) + duplicates(col_ids)
    if dups:
        return {'error_msg': 'Duplicate identifiers: {}'.format(dups) }
    rpart = argsort(row_ids) if sort_rows else list(irange(len(row_ids)))
    cpart = argsort(col_ids) if sort_cols else list(irange(len(col_ids)))
    rowp, colp = get_inverse_perm(rpart, cpart)
    return _pack(rowp, colp)


def random(n_rows, n_cols):
    rowp, colp = list(irange(n_rows)), list(irange(n_cols))
    shuffle(rowp)
    shuffle(colp)
    return _pack(rowp, colp)


def _check_coordinate_format(rows, cols, values, shape):
    try:
        check_coordinate_format(rows, cols, values, shape)
    except AssertionError as e:
        return str(e)


def _check_torn_indices(shape, torn_rows, torn_cols):
    try:
        check_if_indices_are_in_range(torn_rows, shape[0], 'torn row', empty_allowed=True)
        check_if_indices_are_in_range(torn_cols, shape[1], 'torn column', empty_allowed=True)
    except AssertionError as e:
        return str(e)


def _pack(rowp, colp, rpart=[], cpart=[], colors=[]):
    return { 'rowp': rowp, 'colp': colp, 'rpart': rpart, 'cpart': cpart, 
             'color_groups': colors }
