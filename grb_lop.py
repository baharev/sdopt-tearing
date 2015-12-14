# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from itertools import combinations, product
from time import time
from benchmarks import gen_benchmark_digraphs
from grb_simplifier import iteratively_remove_runs_and_bypasses
from grb_pcm import get_orig_edges_map
from mfes import noncopy_split_to_nontrivial_sccs
from py3compat import irange
from utils import double_check, solve_ilp, info_short

################################################################################
#
# The code is admittedly ugly: This is a port of the CPLEX OPL script to Python. 
# The linear ordering problem (LOP) formulation is ineffective; the goal of this
# script is to demonstrate this. 
#
################################################################################


def main():
    for g_input in gen_benchmark_digraphs():
        # Giving up on Problem 10 after 2 hours:
        # 0     0   11.00000    0  402   12.00000   11.00000  8.33%     - 7981s
        if g_input.graph['name'] == 'Problem 10 (opt=12)':
            continue
        solve_problem(g_input)


def solve_problem(g_orig):
    'Returns: [torn edges], cost.'
    elims, cost = [ ], 0
    for sc in noncopy_split_to_nontrivial_sccs(g_orig.copy()): # <- Copy passed!
        partial_elims, partial_cost = solve_triangle_inequalities(sc)
        elims.extend(partial_elims)
        cost += partial_cost
    double_check(g_orig, cost, elims, is_labeled=True)
    print('Input graph')
    info_short(g_orig)
    return elims, cost 


def solve_triangle_inequalities(g):
    iteratively_remove_runs_and_bypasses(g)
    # The above simplification removes edges, but the d['orig_edges'] still 
    # refers to these edges. As a result, the double_check calls in this module
    # and in mfes will fail in the cost calculation, when trying to access 
    # those not present edges to figure out their edge weights. So we have to 
    # replace d['orig_edges'] appropriately and undo this at the end of this 
    # function. Additionally, as runs are removed, new edges are introduced
    # that are not in the original graph. It is much simpler to just relabel 
    # the graph and then undo it on the elimination order. 
    origedges_map = get_orig_edges_map(g)
    #
    model, y = build_lp(g)
    start = time()
    success = solve_ilp(model)
    end = time()
    print('Overall solution time: {0:0.1f} s'.format(end-start))
    assert success, 'Solver failures are not handled at the moment...'
    cost = int(round(model.getObjective().getValue()))
    elim_order = recover_order(model, y, len(g))
    elims = get_torn_edges(g, elim_order)
    # Done. Now undo the  d['orig_edges'] mess.
    final_elims = [ ]
    for edge in elims:
        final_elims.extend(origedges_map[edge])
    double_check(g, cost, elims, is_labeled=True)
    return final_elims, cost


def build_lp(g):
    from gurobipy import LinExpr, GRB, Model #, setParam
    # We introduce an integer index per node ID. Sometimes we will use this 
    # index, and sometimes the node ID; this makes the code a bit messy. 
    i_nodeid = {i: n for i, n in enumerate(g)}
    n_nodes  = len(i_nodeid)
    model = Model()
    # The lower half of the n^2 binary variables, the rest: y_{j,i} = 1-y_{i,j}
    y = { (i, j) : model.addVar(vtype=GRB.BINARY) 
                     for i, j in combinations(irange(n_nodes), 2) }
    model.update()
    # The triangle inequalities
    for i, j, k in combinations(irange(n_nodes), 3):
        lhs = LinExpr([( 1,y[(i,j)]), ( 1,y[(j,k)]), (-1,y[(i,k)])])
        model.addConstr(lhs, GRB.LESS_EQUAL, 1.0)
        lhs = LinExpr([(-1,y[(i,j)]), (-1,y[(j,k)]), ( 1,y[(i,k)])])
        model.addConstr(lhs, GRB.LESS_EQUAL, 0.0)
    # The objective
    shift = 0
    c = [[0]*n_nodes for _ in irange(n_nodes)]
    indices = ((i, j) for i, j in product(irange(n_nodes), repeat=2) 
                       if g.has_edge(i_nodeid[j], i_nodeid[i]))
    for j, k in indices:
        w = g[i_nodeid[k]][i_nodeid[j]]['weight']
        if k < j:
            c[k][j] += w
        elif k > j:
            c[j][k] -= w
            shift   += w
    obj = [(c[i][j], y[(i,j)]) for i, j in combinations(irange(n_nodes), 2)]
    model.setObjective(LinExpr(obj), GRB.MINIMIZE)
    model.setAttr('ObjCon', shift)
    model.update()
    return model, y


def recover_order(model, y, n_nodes):
    # Port of the OPL script: Admittedly ugly and inefficient
    # order[j in 1..n] = sum(i in 1..j-1) x[i,j] + sum(k in j+1..n) (1-x[j,k])+1;
    order   = [-1]*n_nodes
    for j in irange(n_nodes): # counting the number of nodes preceeding j 
        s = 0
        for i in irange(0, j):
            s += int(round(y[(i,j)].X))
        for i in irange(j+1, n_nodes):
            s += int(round((1-y[(j,i)].X)))
        order[j] = s
    # elim_order[k in 1..n] = sum(i in 1..n) i*(order[i]==k);
    print(order)
    elim_order = [-1]*n_nodes
    for k in irange(n_nodes): # the permutation realizing `order`
        elim_order[k] = sum(i if order[i]==k else 0 for i in irange(n_nodes))
    print(elim_order)
    return elim_order


def get_torn_edges(g, elim_order):
    # Torn edges are the back edges (pointing backwards)
    i_nodeid = {i: n for i, n in enumerate(g)}
    keys = {i_nodeid[n]: i for i, n in enumerate(elim_order)}
    n_edges = 0
    edges = [ ]
    for pos, i in enumerate(elim_order):
        n = i_nodeid[i]
        for nbr in g[n]:
            if keys[nbr] > pos: # apparently the elim_order should be reversed?
                n_edges +=1
                edges.append((n, nbr))
    print(n_edges)
    return edges


if __name__ == '__main__':
    main()
