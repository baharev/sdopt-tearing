#!/usr/bin/env python
# Copyright (C) 2016 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
#
# Derived from: sdopt/datagen/generate_gjh.py
from __future__ import print_function, division
from collections import namedtuple
from functools import partial
from itertools import groupby
from multiprocessing import Pool
from os import listdir, mkdir, makedirs, remove #, _exit
from os.path import join, isfile, isdir
from random import Random
from shutil import copy, rmtree
from six import iteritems
from socket import gethostname
from string import Template
from subprocess import call
from sys import stderr
from traceback import format_exc
from networkx import Graph, relabel_nodes
from networkx.algorithms.bipartite import is_bipartite_node_set
from matching import maxmatch_len
from plot_ordering import to_pdf
from py3compat import ifilter, imap, irange, izip
from utils import create_dir_w_parents, serialize, deserialize, pairwise, \
                  print_timestamp

from bb4_tear import *

warning = partial(print, file=stderr)

#===============================================================================

def clean(directory):
    rmtree(directory, ignore_errors=True)
    mkdir(directory)
    
def generate_gjh(presolve):
    clean(GJH_DIR)
    clean(TMPDIR)
    modfiles = sorted(f for f in listdir(COCONUT_MOD_DIR) if f.endswith('.mod'))
    for modfile in modfiles:
        #if modfile.startswith(('lib1_', 'lib2_')):
        #    continue
        #if modfile.startswith('lib2_'):
        #    break
        print(modfile)
        try_generation(modfile, presolve)
    #return
    copy_output()

def try_generation(modfile, presolve):
    # With seed 0 we first try the initial point that is in the mode, 
    # but see also the string template below.
    SEEDS = [0, 1, 31, 42, 13, 26, 87, 59, 64, 77, 95]    
    for seed in SEEDS:
        #
        create_ampl_input(modfile, seed, presolve)
        #
        error = run_ampl(modfile, AMPL_COMMAND)
        if error == ERR_LICENSE:
            # The error has already been logged by AMPL
            return
        elif error == ERR_SEGFAULT:
            warning('###  AMPL SEGFAULT  ###')
            return
        assert not error, 'Unknown error code: {}'.format(error)
        #
        error = run_gjh(modfile)
        if not error:
            return

def create_ampl_input(modfile, seed, presolve):
    with open(join(TMPDIR, modfile), 'w') as dst:
        dst.write(content_until_solve_cmd(modfile))
        presolve_str = PRESOLVE if presolve else NO_PRESOLVE
        script = TEMPLATE.substitute(problem_name=modfile[:-4], seed=seed,
                                     presolve_str=presolve_str) 
        dst.write(script)

def content_until_solve_cmd(modfile):
    with open(join(COCONUT_MOD_DIR, modfile)) as src:
        return ''.join(relevant_lines(src))

def relevant_lines(f):
    for lineno, line in enumerate(f, 1):
        # hs057.mod had display inside a multiline comment!
        if line.startswith('solve;') or line.startswith('display '):
            line = line.rstrip()
            print('Parsing stopped at line {}: {}'.format(lineno, line))
            break
        yield line

def run_ampl(modfile, CMD):
    logfile_name = join(TMPDIR, modfile[:-4]+'.log')
    with open(logfile_name, 'w') as logfile:
        return call([CMD, modfile], cwd=TMPDIR, stdout=logfile)

def run_gjh(modfile):
    basename = modfile[:-4]
    logfile_name = join(TMPDIR, basename + '.log')
    nl_file_name = join(TMPDIR, basename + '.nl')
    with open(logfile_name, 'a') as logfile:
        return call([GJH_COMMAND, nl_file_name], cwd=TMPDIR, stdout=logfile)

def copy_output(): 
    to_cp = sorted(f for f in listdir(TMPDIR) if f.endswith(('.gjh', '.log')))
    for fname in to_cp:
        copy(join(TMPDIR, fname), GJH_DIR)


PRESOLVE = '''
option presolve 100;
option presolve_eps 1.0e-10;
option var_bounds 2;
'''

NO_PRESOLVE = '''
option presolve 0;
'''

TEMPLATE = Template('''##########################################################

option show_stats 1;

$presolve_str

option substout 0;
option nl_comments 0;

option randseed $seed;

# We need j_ instead of j as ampl apparently does not understand scopes

if $seed != 0 then {

for {j_ in 1.._snvars} {
  let _svar[j_] := Uniform(max(-10000, _svar[j_].lb), min(10000,_svar[j_].ub));
  if _svar[j_] == 0 then {
    let _svar[j_] := Uniform(-1, 1);
  }
}

}

display _svar;

param n_ineq_cons_;
let n_ineq_cons_ := 0;
for {i_ in 1.._sncons} {
  if _scon[i_].lbs != _scon[i_].ubs then {
    let n_ineq_cons_ := n_ineq_cons_ + 1;
  }
}


# Force model generation, otherwise the log becomes messed up
print "ignore this line", _snzcons;

print "@@@ Problem statistics: rows  cols  nonzeros  ineq_cons";
print _sncons, _snvars, _snzcons, n_ineq_cons_;
print "";

option solver gjh;
#option gjh_options 'sparse';
write g$problem_name;

''')

#===============================================================================
# nineq: number of inequality constraints
# nentries: the number of nonzeros in the system of the first-order optimality 
#           conditions (well, sort of, assuming equality constraints only)

ProblemStats = namedtuple('ProblemStats', 'nrows  ncols  nzeros  nineq  nentries')


def create_problem_stats_with_segment_check():
    names = get_names_with_check()
    to_eval = {}
    for problem in sorted(names):
        print(problem, end='')
        stats = problem_statistics(join(GJH_DIR, problem+'.log'))
        print('\t', stats.nrows, stats.ncols, stats.nzeros)
        J = checked_J_segment(problem, stats.nrows, stats.ncols, stats.nzeros)
        H = checked_H_segment(problem, stats.ncols)
        entries = len(H) + 2*len(J) + stats.ncols
        stats = stats._replace(nentries=entries)
        to_eval[problem] = stats._asdict()
    create_dir_w_parents(DATABANK_DIR)
    #
    serialize(to_eval, join(DATABANK_DIR, PROBLEM_STATS))
    #
    nmods = len([1 for f in listdir(COCONUT_MOD_DIR) if f.endswith('.mod')])
    fmt = 'All: {}, some output: {}, can be evaluated: {}'
    print(fmt.format(nmods, len(names), len(to_eval)))
    probs_by_nzeros = sorted(iteritems(to_eval), key=lambda x: x[1]['nzeros'], reverse=True)
    fmt = '{}\trows {}\tcols {}\tnnz {}\tineq {}'
    for name, stats in probs_by_nzeros:
        stats = ProblemStats(**stats)
        print(fmt.format(name, stats.nrows, stats.ncols, stats.nzeros, stats.nineq))

def get_names_with_check():
    modfiles = problem_name_set(COCONUT_MOD_DIR, '.mod')
    logfiles = problem_name_set(GJH_DIR, '.log')
    gjhfiles = problem_name_set(GJH_DIR, '.gjh')
    missing_logpair = gjhfiles - logfiles
    assert not missing_logpair, sorted(missing_logpair)
    missing_gjh = logfiles - gjhfiles
    if missing_gjh: # most likely evaluation failed
        warning('gjh is missing for:', sorted(missing_gjh))
    missing = logfiles ^ modfiles 
    assert not missing, sorted(missing)
    return gjhfiles

def problem_name_set(directory, extension):
    ext_len = len(extension)
    return {f[:-ext_len] for f in listdir(directory) if f.endswith(extension)}

def problem_statistics(filename):
    with open(filename) as f:
        for line in f:
            if line.startswith('@@@ Problem statistics: rows  cols  nonzeros'):
                break
        else:
            raise AssertionError('Problem statistics not found in {}'.format(filename))
        line = next(f)
        rows, cols, nzeros, ineq = line.split() # we do not know nentries yet
        return ProblemStats(int(rows), int(cols), int(nzeros), int(ineq), None)

def checked_J_segment(problem, nrows, ncols, nzeros):
    sparse_J = sparse_mat_tuples(join(GJH_DIR, problem + '.gjh'), 'J')
    assert sparse_J or not nrows, problem
    #print(sparse_J)
    assert_eq(len(sparse_J), nzeros)
    row_ids = {i for i, _ in sparse_J}
    col_ids = {j for _, j in sparse_J}
    assert len(row_ids) <= nrows
    assert len(col_ids) <= ncols
    return sparse_J

def checked_H_segment(problem, ncols):
    sparse_H = sparse_mat_tuples(join(GJH_DIR, problem + '.gjh'), 'H')
    #print(sparse_H)
    check_indices_and_symmetry(sparse_H, ncols)
    return sparse_H

def assert_eq(a, b):
    assert a == b, '{} != {}'.format(a, b)

def assert_le(a, b):
    assert a <= b, '{} > {}'.format(a, b)

def check_indices_and_symmetry(mat, size):
    # Zero-based indices assumed, row indices go from 0...size, col indices too.
    # The `mat` is an edgelist: (i, j) tuples
    idx_rowwise = {i for i, _ in mat}
    idx_colwise = {j for _, j in mat}
    assert idx_rowwise == idx_colwise, sorted(idx_rowwise ^ idx_colwise)
    if mat:
        min_idx = min(idx_rowwise)
        max_idx = max(idx_rowwise)
        assert 0 <= min_idx and min_idx < size, (min_idx, size)
        assert 0 <= max_idx and max_idx < size, (max_idx, size)
    transpose = {(j,i) for i, j in mat}
    assert transpose == set(mat)     

#-------------------------------------------------------------------------------
# Parse the sparse matrix format of AMPL

def sparse_mat_tuples(filename, segment):
    sparse_mat = [ ]
    with open(filename, 'r') as f:
        for line in gen_lines(f, segment):
            if line[0] == '[':
                # [2,*] -> i = 2
                i = int(line[1:].split(',', 1)[0])
                continue
            else:
                # line: "   3    51.72" and we want: (2, 3)
                j = int(line.split(None, 1)[0]) # j = 3 in our example
                sparse_mat.append((i-1, j-1)) # zero based
    return sparse_mat

def gen_lines(f, segment):
    assert segment == 'J' or segment == 'H', segment
    # discard until segment, including the param J / H declaration
    declaration = 'param %s :=' % segment
    found = False
    for line in f:
        if line.startswith(declaration):
            found = True
            break
    if found:
        # yield the segment (if any) until the closing ; is reached
        for line in f:
            line = line.rstrip() 
            if line == ';':
                break
            yield line

#===============================================================================

def get_probs_w_stats():
    stats = deserialize(join(DATABANK_DIR, PROBLEM_STATS))
    return [(name, ProblemStats(**(stats[name]))) for name in sorted(stats)]

def delete_all_not_hidden_files(directory):
    to_del = [f for f in listdir(directory) if not f.startswith('.')]
    for f in to_del:
        remove(join(directory, f))
    print('Deleted', len(to_del), 'files in', directory)

def serialize_first_order_cond():
    delete_all_not_hidden_files(DATABANK_DIR)
    counter = 0
    for name, stats in get_probs_w_stats():
        nrows, ncols, nzeros, nineq = stats.nrows, stats.ncols, stats.nzeros, stats.nineq 
        print('---------------------------------------------------------------')
        print('Problem:', name)
        print('rows: {}, cols: {}, nzeros: {}'.format(nrows, ncols, nzeros))
        size = nrows + ncols
        if size < 8:
            print('Too small matrix\n')
            continue
        if nrows*ncols == nzeros:
            print('Full Jacobian, skipping problem\n')
            continue
        if nineq:
            print('Problem has inequality constraints, skipping it\n')
            continue
        sparse_H = sparse_mat_tuples(join(GJH_DIR, name + '.gjh'), 'H')
        if not sparse_H:
            print('Hessian empty, skipping problem\n')
            continue
        H_nzeros = len(sparse_H)
        if H_nzeros == ncols*ncols:
            print('Full Hessian, skipping problem\n')
            continue           
        #if 8 <= size: # and H_nzeros + nzeros <= 10000:
        sparse_J = sparse_mat_tuples(join(GJH_DIR, name + '.gjh'), 'J')
        #
        edgelist = first_order_cond(stats, sparse_J, sparse_H)
        #
        if size != maxmatch_len(*to_sym_bipartite_graph(size, edgelist)):
            print('Structurally singular matrix')
            continue
        #
        serialize((size, edgelist), join(DATABANK_DIR, name + '.pkl.gz'))
        #
        counter += 1
        print()
    print('==================================================================')
    print('Wrote', counter, 'problems')

def first_order_cond(stats, sparse_J, sparse_H):
    nrows, ncols, nzeros = stats.nrows, stats.ncols, stats.nzeros
    size = nrows + ncols
    assert_eq(len(sparse_J), nzeros)
    # Add H
    edgelist = [(i, j+size) for i, j in sparse_H]
    # Add J below H
    edgelist.extend((i+ncols, j+size) for i, j in sparse_J)
    # Add J^T next to H
    shift = size + ncols
    edgelist.extend((j, i+shift) for i, j in sparse_J)
    # Add the identity matrix to the bottom right corner
    gen_i = irange(ncols, ncols+nrows)
    gen_j = irange(shift, shift+nrows)
    edgelist.extend(izip(gen_i, gen_j)) 
    # Debugging:
    tups = [(i, j-size) for i, j in edgelist]
    check_indices_and_symmetry(tups, size)
    return edgelist

#-------------------------------------------------------------------------------

def serialize_constraint_jacobian():
    delete_all_not_hidden_files(DATABANK_DIR)
    counter = 0
    for name, stats in get_probs_w_stats():
        nrows, ncols, nzeros, nineq = stats.nrows, stats.ncols, stats.nzeros, stats.nineq
        print('---------------------------------------------------------------')
        print('Problem:', name)
        print('rows: {}, cols: {}, nzeros: {}'.format(nrows, ncols, nzeros))
        #if nrows != ncols:
        #    print('Rectangular problem\n')
        #    continue
        if nrows < 8 or ncols < 8:
            print('Too small matrix\n')
            continue
        if nrows*ncols == nzeros:
            print('Full Jacobian, skipping problem\n')
            continue
        if nineq:
            print('Problem has inequality constraints, skipping it\n')
            continue
        #
        sparse_J = sparse_mat_tuples(join(GJH_DIR, name + '.gjh'), 'J')
        assert_eq(len(sparse_J), nzeros)
        jac = [(i, j+nrows) for i, j in sparse_J]
        #
        if min(nrows, ncols) != maxmatch_len(*to_bipartite(nrows, ncols, jac)):
            warning('Structurally singular matrix: {}'.format(name))
            continue
        #
        serialize(jac, join(DATABANK_DIR, name + '.pkl.gz'))
        #
        counter += 1
        print()
    print('==================================================================')
    print('Wrote', counter, 'problems')

#===============================================================================

SUCCESS, TIMED_OUT, CRASHED = 1, 2, 3

def solve_problems(dir_name, deserialize_g_eqs, solve_func):
    # setup directory
    directory = join(SOLUTION_DIR, dir_name)
    if isdir(directory):
        print(directory, 'already exists, exiting...')
        return
    makedirs(directory)
    # prepare tasks
    tasks = [ ]
    for name, stats in get_probs_w_stats():
        if not isfile(join(DATABANK_DIR, name + '.pkl.gz')):
            continue
        tasks.append((directory, name, stats, deserialize_g_eqs))
    # run the tasks in a process pool
    pool = Pool(processes=4)
    #tasks = tasks[:10]
    res = list(pool.imap_unordered(solve_func, tasks))
    #res = [solve_func(tasks[0])]
    #
    print('-------------------------------------------------------------')
    print('Success:  ', res.count(SUCCESS))
    print('Timed out:', res.count(TIMED_OUT))
    print('Crashed:  ', res.count(CRASHED))

# solve_func is either solve_rnd_instance or solve_instance. With the 
# solve_rnd_instance we solve the original input, the reversed original input,
# and 10 random permutations of the input, and we bogusly always report success.
#-------------------------------------------------------------------------------

def solve_rnd_instance(_tup):
    try:
        directory, name, stats, deserialize_g_eqs = _tup
        g_orig, eqs = deserialize_g_eqs(name, stats)
        solve_permutations(g_orig, eqs, name, directory)
    except:
        print('\n'.join(('Crashed:', name, format_exc())))
        return CRASHED        
    return SUCCESS # This is a lie

def solve_permutations(g, eqs, name_orig, directory):
    for name, mapping in iteritems(create_mappings(name_orig, len(eqs))):
        solve_internal(g, eqs, mapping, name, directory)

def create_mappings(name, nrows):
    # Careful: get_rnd_g_eqs and base_of_rnd knows these naming conventions!
    eq_ids = list(irange(nrows))
    mappings = { }
    # Input, mapping is identity
    mappings[name+'-IN'] = {e: e for e in eq_ids}
    # Reversed input 
    mappings[name+'-RE'] = {e: re for e, re in izip(eq_ids, reversed(eq_ids))}
    # Random permutations
    rng = Random(31)
    for i in irange(0, 10):
        rng.shuffle(eq_ids)    
        mappings[name+'-{:02}'.format(i)] = {k: eq_ids[k] for k in irange(nrows)}
    return mappings    

def solve_internal(g_orig, eqs, mapping, name, directory):
    # get_rnd_g_eqs duplicates some of this code!
    g = relabel_nodes(g_orig, mapping)
    assert is_bipartite_node_set(g, eqs)
    res = solve_problem(g, eqs) #, log=print)
    serialize(res._asdict(), join(directory, name + '.result.pkl.gz'))
    serialize(mapping,       join(directory, name + '.mapping.pkl.gz'))

#-------------------------------------------------------------------------------

def solve_instance(_tup):
    directory, name, stats, deserialize_g_eqs = _tup
    g, eqs = deserialize_g_eqs(name, stats)
    try:
        res = solve_problem(g, eqs) #, log=print)
    except:
        print('\n'.join(('Crashed:', name, format_exc())))
        return CRASHED
    serialize(res._asdict(), join(directory, name + '.pkl.gz'))
    return SUCCESS if res.optimal else TIMED_OUT

#-------------------------------------------------------------------------------
# Stuff for deserialization

def get_g_eqs(name, stats):
    edgelist = deserialize(join(DATABANK_DIR, name + '.pkl.gz'))
    assert_eq(len(edgelist), stats.nzeros)
    return to_bipartite(stats.nrows, stats.ncols, edgelist)

def get_rnd_g_eqs(name, stats, dir_name=None):
    # Duplication with get_rnd_first_order_cond_g_eqs
    assert dir_name is not None
    edgelist = deserialize(join(DATABANK_DIR, base_of_rnd(name) + '.pkl.gz'))
    assert_eq(len(edgelist), stats.nzeros)
    g, eqs = to_bipartite(stats.nrows, stats.ncols, edgelist)
    mapping = deserialize(join(SOLUTION_DIR, dir_name, name + '.mapping.pkl.gz'))
    return relabel_nodes(g, mapping), eqs

def base_of_rnd(name):
    # Find the original graph; name-XX -> name
    return name.split('-', 1)[0] 
 
def get_first_order_cond_g_eqs(name, stats):
    size, edgelist = deserialize(join(DATABANK_DIR, name + '.pkl.gz'))
    assert size == stats.nrows + stats.ncols
    return to_sym_bipartite_graph(size, edgelist)

def get_rnd_first_order_cond_g_eqs(name, stats, dir_name=None):
    # Duplication with get_rnd_g_eqs
    assert dir_name is not None
    size, edgelist = deserialize(join(DATABANK_DIR, base_of_rnd(name) + '.pkl.gz'))
    assert size == stats.nrows + stats.ncols
    g, eqs = to_sym_bipartite_graph(size, edgelist)
    mapping = deserialize(join(SOLUTION_DIR, dir_name, name + '.mapping.pkl.gz'))
    return relabel_nodes(g, mapping), eqs

def to_sym_bipartite_graph(size, edgelist):
    return to_bipartite(size, size, edgelist)

def to_bipartite(nrows, ncols, edgelist):
    g = Graph()
    g.add_nodes_from(range(nrows+ncols))
    g.add_edges_from(edgelist)
    eqs = set(range(nrows))
    assert is_bipartite_node_set(g, eqs)
    return g, eqs

#===============================================================================

def read_solutions(dir_name):
    #
    directory = join(SOLUTION_DIR, dir_name)
    #
    optimal, timedout = 0, 0
    for name, stats in get_probs_w_stats():
        if not isfile(join(directory, name + '.pkl.gz')):
            continue
        print('------------------------------------------------------------')
        print('Problem:', name)
        nrows, ncols, nzeros = stats.nrows, stats.ncols, stats.nzeros
        print('rows: {}, cols: {}, nzeros: {}'.format(nrows, ncols, nzeros))
        res = deserialize(join(directory, name + '.pkl.gz'))
        print('Explored:', res['explored'])
        print('Cost: ', res['ub'])
        print('Time: {0:0.1f}'.format(res['time']))
        if res['optimal']:
            if res.get('heur'):
                assert res['explored'] == 0
                print('Solved by:', res['heur'])
            print('Optimal')
            optimal += 1
        else:
            print('Gap:', res['gap'])
            timedout += 1
    print('=============================================================')
    print('Success:  ', optimal)
    print('Timed out:', timedout)
    
#===============================================================================

def read_rnd_solutions(dir_name):
    #
    directory = join(SOLUTION_DIR, dir_name)
    prob_stats = {name: stats for name, stats in get_probs_w_stats()}
    # names, with the '.result.pkl.gz' cropped off
    names = get_names_of_rnd(directory)
    optimal, timedout = 0, 0
    for name in names:
        stats = prob_stats[base_of_rnd(name)]
        print('------------------------------------------------------------')
        print('Problem:', name)
        nrows, ncols, nzeros = stats.nrows, stats.ncols, stats.nzeros
        print('rows: {}, cols: {}, nzeros: {}'.format(nrows, ncols, nzeros))
        res = deserialize(join(directory, name + '.result.pkl.gz'))
        print('Explored:', res['explored'])
        print('Cost: ', res['ub'])
        print('Time: {0:0.1f}'.format(res['time']))
        if res['optimal']:
            if res.get('heur'):
                assert res['explored'] == 0
                print('Solved by:', res['heur'])
            print('Optimal')
            optimal += 1
        else:
            print('Gap:', res['gap'])
            timedout += 1
    print('=============================================================')
    print('Success:  ', optimal)
    print('Timed out:', timedout)

#-------------------------------------------------------------------------------
# This is for the randomized results.

def show_perf_stats(dir_name):
    print_timestamp()
    #
    pool = Pool(processes=4)
    res = list(pool.imap_unordered(stats_of, groups(dir_name)))
    # elements of res: (basename, delta_ub, delta_gap)
    #
    print('Processed {} problems'.format(len(res)))
    ub_delta_sorted  = sorted(res, key=lambda t: t[1], reverse=True)
    gap_delta_sorted = sorted(res, key=lambda t: t[2], reverse=True)
    #
    N = 20
    print('Top', N,'UB deltas:')
    for i in range(N):
        name, diff = ub_delta_sorted[i][0], ub_delta_sorted[i][1]
        print('{}: {}, '.format(name, diff))
    print()
    print('Top', N,'gap deltas:')
    for i in range(N):
        name, diff = gap_delta_sorted[i][0], gap_delta_sorted[i][2]
        print('{}: {}, '.format(name, diff))
    #
    print_timestamp()   

def groups(dir_name):
    prob_stats = {name: stats for name, stats in get_probs_w_stats()}
    directory = join(SOLUTION_DIR, dir_name)
    # names, with the '.result.pkl.gz' cropped off
    names = get_names_of_rnd(directory)
    for basename, group in groupby(names, key=base_of_rnd):
        stats = prob_stats[basename]
        if stats.nentries > 10000:
            continue
        yield directory, basename, sorted(group)

def stats_of(directory_basename_group):
    directory, basename, group = directory_basename_group
    def deser(name):
        return deserialize(join(directory, name + '.result.pkl.gz'))
    #stats = prob_stats[basename]
    #print('------------------------------------------------------------')
    #print('Problem:', basename)
    #nrows, ncols, nzeros = stats.nrows, stats.ncols, stats.nzeros
    #print('rows: {}, cols: {}, nzeros: {}'.format(nrows, ncols, nzeros))
    res  = list(imap(deser, group))
    ubs  = list(imap(lambda r: r['ub'], res))
    lbs  = list(imap(lambda r: r['ub'] - r['gap'], res))
    gaps = list(imap(lambda r: r['gap'], res)) 
    min_ubs = min(ubs)
    max_lbs = max(lbs)
    assert min_ubs >= max_lbs, (basename, min_ubs, max_lbs)
    delta_ubs  = max(ubs)  - min_ubs
    delta_gaps = max(gaps) - min(gaps)    
    return (basename, delta_ubs, delta_gaps)

def get_names_of_rnd(directory):
    # names, with the '.result.pkl.gz' cropped off
    ext = '.result.pkl.gz'
    ext_len = len(ext)
    return sorted(f[:-ext_len] for f in listdir(directory) if f.endswith(ext))

#-------------------------------------------------------------------------------
# This is for the randomized results.
#
# size_func: number of rows, number of nonzeros in the constraint Jacobian, or
#     KKT conditions size: number of rows + cols, nonzero entries in this matrix
# skip_func: to select only the square problems, filtered_groups

ANY_FAILED, ALL_SUCCEEDED = 0, 1 

def show_size_distribution(dir_name, size_func, skip_func):
    print_timestamp()
    prob_stats = {name: stats for name, stats in get_probs_w_stats()}
    skip = partial(skip_func, prob_stats)
    pool = Pool(processes=4)
    res = list(pool.imap_unordered(outcome, filtered_groups(dir_name, skip)))
    any_failed  = [name for name, opt in res if opt==ANY_FAILED]
    all_success = [name for name, opt in res if opt==ALL_SUCCEEDED] 
    print('Processed {} problems'.format(len(res)))
    #print('\n'.join(sorted(name for name, _opt in res)))
    print(len(all_success), 'succeeded in all cases')
    print(len(any_failed),  'timed out at least once')
    size = partial(size_func, prob_stats)  
    print_distrib('Success:', all_success, prob_stats, size)
    print_distrib('Failed:',  any_failed,  prob_stats, size)
    gen_input = (name for _dir, name, _group in filtered_groups(dir_name, skip))
    print_distrib('Input:', gen_input, prob_stats, size)
    print_timestamp()

def size_KKT(prob_stats, name):
    stats = prob_stats[name]
    return (stats.nrows+stats.ncols, stats.nentries)

def size_jac(prob_stats, name):
    stats = prob_stats[name]
    return (stats.nrows, stats.nzeros)

def no_skip(_1, _2):
    return False

def skip_non_square(prob_stats, name):
    stats = prob_stats[name]
    return stats.nrows != stats.ncols

def filtered_groups(dir_name, skip):
    directory = join(SOLUTION_DIR, dir_name)
    # names, with the '.result.pkl.gz' cropped off
    names = get_names_of_rnd(directory)
    for basename, group in groupby(names, key=base_of_rnd):
        if skip(basename):
            continue
        yield directory, basename, sorted(group)

def outcome(directory_basename_group):
    directory, basename, group = directory_basename_group
    def deser(name):
        return deserialize(join(directory, name + '.result.pkl.gz'))
    has_failed = any(ifilter(lambda res: not res['optimal'], imap(deser, group)))
    if has_failed:
        return basename, ANY_FAILED
    return basename, ALL_SUCCEEDED

def print_distrib(msg, iterable, prob_stats, size):
    print()
    print(msg)
    arr = sorted(iterable, key=size)
    print_quantiles(arr, prob_stats, size)

def print_quantiles(arr, prob_stats, size):
    print('Distribution (0%, 25%, 50%, 75%, 100%)')
    assert arr, 'Expects an array with at least one element'
    sz = len(arr)-1
    quantiles = (0, int(sz*0.25), int(sz*0.5), int(sz*0.75), sz) 
    for i in quantiles:
        print('   ', prob_stats[arr[i]], '(%d)' % (i+1))
        print('   ', size(arr[i]))

#===============================================================================
# Currently unused. It compares two directories and asserts that results in the
# second and not worse than those in the first. It should be true if we run the 
# same algorithm for significantly longer time (10s vs. 30s). It was useful for
# debugging mainly the reported results (sometimes bogus results were reported). 

def cross_check():
    counter = 0
    dir1, dir2 = join(SOLUTION_DIR, '10s'), join(SOLUTION_DIR, '30s')
    for name, stats in get_probs_w_stats():
        f1, f2 = join(dir1, name + '.pkl.gz'), join(dir2, name + '.pkl.gz')
        if not isfile(f1):
            continue
        assert isfile(f2)
        res1, res2 = deserialize(f1), deserialize(f2)
        if not res2['optimal']:
            continue
        if res2.get('heur'):
            assert res2['explored'] == 0
        print('------------------------------------------------------------')
        print('Problem:', name)
        nrows, ncols, nzeros = stats.nrows, stats.ncols, stats.nzeros
        print('rows: {}, cols: {}, nzeros: {}'.format(nrows, ncols, nzeros))
        print('Cost: ', res2['ub'])
        print('Improved:', res1['ub'] - res2['ub'])
        assert res1['ub'] >= res2['ub']
        lb1 = res1['ub'] - res1['gap']
        lb2 = res2['ub'] - res2['gap']
        if lb1 > lb2:
            warning('Bogus lower bounds: {} > {}'.format(lb1, lb2))
        assert lb1 <= lb2
        counter += 1
    print('=============================================================')
    print(counter, 'plots written')
    print('Done!')

# Similar to cross_check but on a sequence of directories.

def check_monotonicity():
    #
    selected = [ ]
    directory = join(SOLUTION_DIR, 'fix_0.5')
    for name, _ in get_probs_w_stats():
        fname = join(directory, name + '.pkl.gz')
        if isfile(fname):
            selected.append(name)
    assert len(selected) == 273
    dirs = ['fix_0.0', 'fix_0.1', 'fix_0.25', 'fix_0.5', 'fix_1.0', 'fix_2.0', 'fix_4.0', 'fix_8.0', 'fix_16.0', 'fix_32.0']
    for name in selected:
        print('------------------------------------------------------------')
        print('Problem:', name)
        lbs, ubs = [ ], [ ]
        for d in dirs:
            res = deserialize(join(SOLUTION_DIR, d, name + '.pkl.gz'))
            lbs.append(res['ub']-res['gap'])
            ubs.append(res['ub'])
        print('LB:', lbs)
        print('UB:', ubs)
        if not all(a <= b for a, b in pairwise(lbs)):
            print('Non-monotone lbs', file=stderr)
            assert False
        if not  all(a >= b for a, b in pairwise(ubs)):
            print('Non-monotone ubs', file=stderr)
            assert False
        if max(lbs) > min(ubs):
            print('Disjoint bounds',  file=stderr)
            assert False

#===============================================================================

SKIPPED, PLOTTED = 0, 1

def solutions_to_pdf(dir_name, deser, select_func):
    clean(PICS_TMP_DIR)
    #
    pool = Pool(processes=4)
    res = list(pool.imap_unordered(plot, select_func(dir_name, deser)))
    #res = [ plot(next(select_func(dir_name, deser))) ]    
    #
    print('=============================================================')
    plotted = res.count(PLOTTED)
    skipped = res.count(SKIPPED)
    print(plotted, 'plots written')
    print(skipped, 'skipped')
    print('Concatenating all PDF files into a single file')
    call('pdftk *.pdf cat output 000.pdf', shell=True, cwd=PICS_TMP_DIR)
    print('Done!')

def select_plots(dir_name, deser):
    for name, stats in get_probs_w_stats():
        fname = join(SOLUTION_DIR, dir_name, name + '.pkl.gz')
        if not isfile(fname):
            continue
        res = deserialize(fname)
        #if not res['optimal']:
        #    continue
        yield (name, stats, res, deser)

def select_rnd_plots(dir_name, deser):
    directory = join(SOLUTION_DIR, dir_name)
    prob_stats = {name: stats for name, stats in get_probs_w_stats()}
    # names, with the '.result.pkl.gz' cropped off
    names = get_names_of_rnd(directory)
    for name in names:
        stats = prob_stats[base_of_rnd(name)]
        fname = join(SOLUTION_DIR, dir_name, name + '.result.pkl.gz')
        res = deserialize(fname)
        #if not res['optimal']:
        #    continue
        yield (name, stats, res, deser)

def plot(_args):
    name, stats, res, deser = _args
    MAX_NZEROS = 500
    #
    if stats.nzeros > MAX_NZEROS:
        warning('Skipping {}, too many non-zeros to plot!'.format(name))
        return SKIPPED
    #
    g, eqs = deser(name, stats)
    n_edges = g.number_of_edges() # Not necessarily the same as stats.nzeros
    #
    if n_edges > MAX_NZEROS:
        warning('Skipping {}, too many non-zeros to plot!'.format(name))
        return SKIPPED
    #
    left_lbl  = left_label(g, eqs, n_edges, name, stats)
    rowp, colp = list(irange(len(eqs))), list(irange(len(eqs), len(g)))
    to_pdf(g, rowp, colp, msg=left_lbl,  fname=name+'_a', path=PICS_TMP_DIR)
    #
    right_lbl = right_label(g, eqs, n_edges, name, stats, res)
    rowp, colp = res['rowp'], res['colp']
    to_pdf(g, rowp, colp, msg=right_lbl, fname=name+'_b', path=PICS_TMP_DIR)
    #
    return PLOTTED

def left_label(g, eqs, n_edges, name, stats):
    rows, cols = len(eqs), len(g)-len(eqs)
    lbl = '{}  {}x{} ({})'.format(name, rows, cols, n_edges) 
    ncon, nvar = stats.nrows, stats.ncols
    if rows != ncon and cols != nvar:
        lbl = lbl + ' con: {}, var: {}'.format(ncon, nvar)
    return lbl

def right_label(g, eqs, n_edges, name, stats, res):
    cost = res['ub']
    cost_str = 'cost: {}'.format(cost)
    diff_str = ''
    if len(g)-len(eqs) >= len(eqs): # ncols >= nrows or True for KKT conditions
        diff = cost - (len(g)-2*len(eqs)) # nrows - ncols, or 0 for KKT cond.
        diff_str = 'diff: {}'.format(diff)
    explored = 'explored: {}'.format(res['explored'])
    time_str = 'time: {:0.2f}s'.format(res['time'])
    heur_str = ''
    if res.get('heur'):
        assert res['explored'] == 0
        heur_str = res['heur']
    gap_str = ''
    if not res['optimal']:
        gap_str = 'gap: {}'.format(res['gap'])
    lbl = (s for s in (name, cost_str, diff_str, explored, time_str, heur_str, gap_str) 
              if s)
    return ', '.join(lbl)

#===============================================================================
# Currently unused stuff that was useful for debugging

def solve_one_specific_problem():
    #name = 'aug3d'
    #nrows, ncols, nzeros = 1000, 3873, 6546
    #
    name = 'sosqp2'
    nrows, ncols, nzeros = 10001, 20000, 40000
    g, eqs = get_g_eqs(name, nrows, ncols, nzeros)
    solve_problem(g, eqs, log=print)

def plot_dulmage_mendelsohn():
    from test_bb_tear import plot_dm_decomp # This function is in the wrong module
    for name, stats in get_probs_w_stats():
        fname = join(DATABANK_DIR, name + '.pkl.gz')
        if not isfile(fname):
            continue
        g, eqs = get_first_order_cond_g_eqs(name, stats)
        if g.number_of_edges() > 1000:
            warning('Too many non-zeros: {}'.format(name))
            continue
        plot_dm_decomp(g, len(eqs))

#===============================================================================

TMPDIR  = '/tmp/gjh/'
COCONUT_MOD_DIR = 'data/coconut/prefixed_modfiles'
GJH_DIR = 'data/coconut/prefixed_presolve_gjh/'

#DATABANK_DIR = 'data/coconut/jac_databank/'
#SOLUTION_DIR = 'data/coconut/jac_rnd_solutions/'
DATABANK_DIR = 'data/coconut/first_order_opt_databank/'
SOLUTION_DIR = 'data/coconut/first_order_opt_solutions/'


PROBLEM_STATS = '.stats.pkl.gz' # lives in DATABANK_DIR
PICS_TMP_DIR = '/users/ali/tmp/' if gethostname() == 'otto' else '/tmp/pics/'

# We give up on the problem if AMPL returns a specific code
ERR_LICENSE  = 2
ERR_SEGFAULT = -11

AMPL_COMMAND = '/home/ali/ampl_old/ampl'
GJH_COMMAND = 'gjh'

def run():
    ###  Works in TMPDIR, writes (populates) GJH_DIR
    #generate_gjh(presolve=True)
    #---------------------------------------
    ###  Creates the PROBLEM_STATS file in the DATABANK_DIR but nothing else
    #create_problem_stats_with_segment_check()
    #---------------------------------------
    ###  Creates the edgelists in DATABANK_DIR based on PROBLEM_STATS
    #serialize_constraint_jacobian()  # Deletes all non-hidden files
    #deser = get_g_eqs
    #-
    #serialize_first_order_cond() # Deletes ditto
    #deser = get_first_order_cond_g_eqs
    #---------------------------------------
    TIMEOUT = str(TIME_LIMIT)
    dirname = TIMEOUT + 's'
    #---------------------------------------
    ###  Writes the SOLUTION_DIR. All problems in PROBLEM_STATS solved with 
    #    either solve_instance or solve_rnd_instance. The latter also writes the
    #    random mapping (relabelling of the eqs) into the solution directory.
    #solve_func = solve_instance
    #solve_func = solve_rnd_instance
    #solve_problems(dirname, deser, solve_func)
    #---------------------------------------
    #read_solutions(dirname)
    #
    #read_rnd_solutions(dirname)
    #show_perf_stats(dirname)
    show_size_distribution(dirname, size_KKT, no_skip) 
    #show_size_distribution(dirname, size_jac, no_skip) # no_skip, skip_non_square 
    #---------------------------------------
    ###  Writes the PICS_TMP_DIR
    #    Non-random:
    #select_func = select_plots
    #    Random: needs a new deser since we have to relabel the equations too 
    #    and the mapping for that is in the dir_name directory.
    #    We also have to reassing deser: a curried version of either 
    #    get_rnd_g_eqs or get_rnd_first_order_cond_g_eqs.
    #deser = partial(get_rnd_first_order_cond_g_eqs, dir_name=dirname)
    #deser = partial(get_rnd_g_eqs, dir_name=dirname)
    #select_func = select_rnd_plots
    #
    #solutions_to_pdf(dirname, deser, select_func)
    #---------------------------------------
    #cross_check()
    #check_monotonicity()
    #solve_one_specific_problem()

if __name__=='__main__':
    run()
