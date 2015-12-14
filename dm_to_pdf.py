# Copyright (C) 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
#from glob import glob
from cjson import decode as loads
import six
from dm_decomp import sccs_on_plot, blt_with_tearing #dm_coordinate
from order_util import check_coordinate_format, coo_matrix_to_bipartite
from plot_ordering import plot_dm

TOO_SMALL = 1
TOO_BIG   = 100000

# input: decompress data/exploring/testenv_all_merged.zip into /tmp 
#        to get /tmp/testenv/*.json
# output: the PDFs will be created in /tmp/pics

TESTENV_DIR = '/tmp/testenv/'
EXT = '.json'


class TestsetStatistics:
    def __init__(self):
        self.too_small = [ ]
        self.too_big   = [ ]
        self.prob_info = { }


class ProblemStatistics:
    # A Data Transfer Object (DTO)
    def __init__(self, size, n_nzeros, h_index, c_factor, d, n_sccs, max_scc, 
                 max_scc_percentage):
        self.size = size # (n_rows, n_cols)
        self.n_nzeros = n_nzeros
        self.h_index = h_index
        self.c_factor = c_factor
        self.d = d        
        self.n_sccs = n_sccs
        self.max_scc = max_scc
        self.max_scc_percentage = max_scc_percentage

    def __str__(self):
        # Not shown: nz={}, 
        fmt = '{}x{}, h={}, c={:.1f}, d={}, SCCs={} (max={}, {:.1f}%)'
        return fmt.format( self.size[0], self.size[1], self.h_index, 
                           self.c_factor, self.d, self.n_sccs, self.max_scc, 
                           self.max_scc_percentage )

    def to_csv_str(self):
        # rows, cols, h, d
        fmt = '{};{};{};{}'
        return fmt.format(self.size[0], self.size[1], self.h_index, self.d)


def test_set_image_export():
    #test_problems = sorted(glob(TESTENV_DIR+'*'+EXT))
    test_problems = [TESTENV_DIR + name + EXT for name in sorted(HAND_GUIDED_TEARING)]
    ext_lenght = len(EXT)
    stats = TestsetStatistics()
    for problem in test_problems:
        i = problem.rfind('/') + 1
        name = problem[i:-ext_lenght]
        print(name)
        with open(problem) as json_file:
            data = loads(json_file.read().rstrip())
        rows, cols, values = data['rows'], data['cols'], data['weights']
        rowIDs, colIDs = data['rowIDs'], data['colIDs']
        process_problem(name, rowIDs, colIDs, rows, cols, values, stats)
        #just_basic_stats(name, rowIDs, colIDs, rows, cols, values, stats)
    print_skipped(stats)
    print('name;n_rows;n_cols;h;d')
    for name, info in sorted(six.iteritems(stats.prob_info)):
        line = ';'.join((name, info.to_csv_str())) 
        print(line) 
    print('Done!')


def just_basic_stats(name, rowIDs, colIDs, rows, cols, values, stats):
    n_rows, n_cols = len(rowIDs), len(colIDs)
    shape = (n_rows, n_cols)
    error_msg = check_coordinate_format(rows, cols, values, shape)
    assert error_msg is None, error_msg
    # The col IDs in cols are shifted by n_rows, must undo later
    g, eqs, _ = coo_matrix_to_bipartite(rows, cols, values, shape)
    entries_per_row = sorted((len(g[eq]) for eq in eqs), reverse=True)
    h_row, _ = max(enumerate(entries_per_row, start=1), key=lambda (i, e): min(i, e))
    entries_per_col = sorted((len(g[n]) for n in g if n not in eqs) , reverse=True)
    h_col, _ = max(enumerate(entries_per_col, start=1), key=lambda (i, e): min(i, e))
    h_index = min(h_row, h_col)
    # the max of the min row/col count puts a lower bound on optimum in tearing
    min_count_row = min(eqs, key=lambda eq: len(g[eq]))
    min_count_col = min((n for n in g if n not in eqs), key=lambda n: len(g[n])) 
    min_rowcount = len(g[min_count_row])
    min_colcount = len(g[min_count_col])
    d = min(min_rowcount, min_colcount) # Careful: here min, in an SCC max!
    stats.prob_info[name] = ProblemStatistics(shape, -1, h_index, -1, d, -1, -1, -1)


HAND_GUIDED_TEARING = {

'bratu_020_A' : (['eq_20_'], 
                 ['x_1_']),

'extr22_A' : (['M_eq_1_22_', 'M_eq_2_22_'], 
              ['v_1_1_', 'v_2_1_', 'v_3_1_']),

'extr30_A' : (['M_eq_1_30_', 'M_eq_2_30_'],
              ['v_1_1_', 'v_2_1_', 'v_3_1_']),

'JacobsenChemProcSim1D_A' : (['VLE_eq_7_'], 
                             ['zeta']),

'JacobsenChemProcSim2D_A' : (['VLE_eq_8_', 'spec_V_N'], 
                             ['v_1_1_', 'v_2_1_']),

'JacobsenTorn2D_A' : (['component_balance_0_', 'heat_balance_0_'], 
                      ['L_8_', 'x_8_']),

'mss20_A' : (['M_eq_19_1_', 'M_eq_19_2_'], 
             ['x_1_1_', 'x_2_1_', 'x_3_1_']),

'mss20Bloat_A' : (['M_eq_19_1_', 'M_eq_19_2_'], 
                  ['x_1_1_', 'x_2_1_', 'x_3_1_']),

'mss20Torn_A' : (['M_eq_19_1_', 'M_eq_19_2_'], 
                 ['x_1_1_', 'x_2_1_', 'x_3_1_']),

'mss20heatBalance_A' : (['M_eq_1_20_', 'M_eq_2_20_', 'M_eq_3_20_','sum_y_minus_sum_x_20_', 'sum_z_20_'], 
                        ['x_1_1_', 'x_2_1_', 'x_3_1_', 'T_1_', 'lambda_1_', 'V_1_']),

'mss20heatBalanceBloat_A' : (['M_eq_1_20_', 'M_eq_2_20_', 'M_eq_3_20_','sum_y_minus_sum_x_20_', 'sum_z_20_'], 
                             ['x_1_1_', 'x_2_1_', 'x_3_1_', 'T_1_', 'lambda_1_', 'V_1_']),

'ReactorLuyben_A' : (['total_condenser'],
                     ['TR']),

'TunnelDiodesSum20_A' : (['def_s'],
                         ['s']),

}

def process_problem(name, rowIDs, colIDs, rows, cols, values, stats):
    n_rows, n_cols = len(rowIDs), len(colIDs)
    if n_rows < TOO_SMALL or n_cols < TOO_SMALL:
        stats.too_small.append(name)
        return
    if n_cols > TOO_BIG or n_rows > TOO_BIG:
        stats.too_big.append(name)
        return
    shape = (n_rows, n_cols)
    error_msg = check_coordinate_format(rows, cols, values, shape)
    assert error_msg is None, error_msg
    #---
    rname_to_idx = {r: i for i, r in enumerate(rowIDs)}
    cname_to_idx = {c: i for i, c in enumerate(colIDs)}
    torn_rows, torn_cols = HAND_GUIDED_TEARING[name]
    torn_rows = [rname_to_idx[r] for r in torn_rows]
    torn_cols = [cname_to_idx[c] for c in torn_cols]
    tup = blt_with_tearing(rows, cols, values, shape, torn_rows, torn_cols)
    rowp, colp,  R0, C0, _, _, rowpart, colpart, colors = tup
    g, _, _ = coo_matrix_to_bipartite(rows, cols, values, shape)
    #---
    #tup = dm_coordinate(rows, cols, values, shape, upper=False)
    ## Careful: g has shifted col IDs
    #g, rowp, colp,  R0, C0, _, _, rowpart, colpart, colors = tup
    #---
    sccs = sccs_on_plot(rowpart, colpart, rowp, colp)
    prob_stat = problem_statistics(g, rows, n_rows, n_cols, R0, C0, sccs)
    stats.prob_info[name] = prob_stat
    # Saving plot as PDF
    msg = str(prob_stat)
    plot_dm(name, rows, cols, rowp, colp, colors, sccs, show=False, msg=msg)


def problem_statistics(g, rows, n_rows, n_cols, R0, C0, sccs):
    # building a new graph without the unmatched rows and columns
    # re-shift the col IDs as g has them with shifted IDs
    R0, C0 = set(R0), set(c+n_rows for c in C0)
    g = g.subgraph(n for n in g if n not in R0 and n not in C0)
    eqs = set(r for r in rows if r not in R0) # each matched row *exactly* once
    # h index: function f in decreasing order, h(f) = max_i min(f(i), i)
    entries_per_row = sorted((len(g[eq]) for eq in eqs), reverse=True)
    h_index, _ = max(enumerate(entries_per_row, start=1), key=lambda (i, e): min(i, e))
    n_entries = g.number_of_edges()
    assert n_entries <= len(rows)
    c_factor = n_entries / float(n_rows)
    # the max of the min row/col count puts a lower bound on optimum in tearing
    min_count_row = min(eqs, key=lambda eq: len(g[eq]))
    min_count_col = min((n for n in g if n not in eqs), key=lambda n: len(g[n])) 
    min_rowcount = len(g[min_count_row])
    min_colcount = len(g[min_count_col])
    #if min_rowcount!=min_colcount:
    #    pass
    d = max(min_rowcount, min_colcount)
    n_sccs = len(sccs)
    _, _, max_scc_size = max(sccs, key=lambda t: t[2]) # sccs: [(r,c,size)]
    prob_size = min(n_rows, n_cols)
    max_scc_percentage = (max_scc_size/float(prob_size)) * 100.0
    return ProblemStatistics( (n_rows, n_cols), n_entries, h_index, c_factor, \
                              d, n_sccs, max_scc_size, max_scc_percentage )


def print_skipped(stats):
    print('-------------------------------------------------------------------')
    if stats.too_small:
        print('Skipped because it was too small (<{})'.format(TOO_SMALL))
        print([name for name in sorted(stats.too_small)])
    if stats.too_big:
        print('Skipped because it was too big (>{})'.format(TOO_BIG))
        print([name for name in sorted(stats.too_big)])


if __name__ == '__main__':
    test_set_image_export()
