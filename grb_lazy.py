# Copyright (C) 2016 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from traceback import print_exc
from sys import stderr
from gurobipy import GRB, LinExpr, Model
from six import iteritems
from networkx import is_directed_acyclic_graph, shortest_path
from grb_simplifier import iteratively_remove_runs_and_bypasses
from mfes import noncopy_split_to_nontrivial_sccs
from utils import info_short as info
from utils import double_check, to_cycle, TimeLimit
from grb_pcm import get_orig_edges_map, feasible_solution, initial_loop_set, \
                    missed_edges

def log(*args, **kwargs): pass
#    print('+ ', *args, **kwargs) if args or kwargs else print()

# Stats only contains counters that will be mutated in-place if passed as the 
# keyword argument stats.

def solve_problem(g_orig, stats, feasible_sol=None):
    'Returns: [torn edges], cost.'
    elims, cost = [ ], 0
    g2 = g_orig.copy()
    # Remove self-loops
    for u, v, d in g_orig.selfloop_edges(data=True):
        g2.remove_edge(u, v)
        elims.append((u,v))
        cost += d['weight']
    cycle_matrix = []
    # Make each SCC acyclic
    for sc in noncopy_split_to_nontrivial_sccs(g2):
        partial_elims, partial_cost, loops = solve_with_pcm(sc, stats, feasible_sol)
        elims.extend(partial_elims)
        cost += partial_cost
        cycle_matrix.extend(loops)
    double_check(g_orig, cost, elims, is_labeled=True, log=log)
    log('Input graph')
    info(g_orig, log=log)
    return elims, cost, cycle_matrix

def solve_with_pcm(g, stats, feasible_sol):
    g_orig = g.copy()
    
    iteratively_remove_runs_and_bypasses(g) # d['orig_edges'] mess,
    origedges_map = get_orig_edges_map(g)   # see explanation in grb_pcm.py
    #
    # FIXME It assumes that we only have a single SCC
    elims, ub = feasible_solution(g) if feasible_sol is None else feasible_sol
    # A shortest path around each edge, see grb_pcm.py for possible improvements
    loops = initial_loop_set(g)
    print('Initial cycle matrix size:', len(loops))
    # Build the model, set the elims as initial solution
    m, vrs = build_ilp(g, loops, elims)
    # Put into the model dict everything we need in the callback
    m._vrs, m._g, m._loops, m._elims, m._ub = vrs, g, loops, elims, ub
    #
    solve(m, stats)
    #
    loops, elims, ub = m._loops, m._elims, m._ub 
    print('Final cycle matrix size:', len(loops))    
    #simplify(g, loops, elims, ub)
    
    #run_IIS(g, loops, elims, ub)
    # Done. Now undo the  d['orig_edges'] mess.
    final_elims = [ ]
    for edge in elims:
        final_elims.extend( origedges_map[edge] )
    
    double_check(g_orig, ub, final_elims, is_labeled=True, log=print)
    return final_elims, ub, loops

def solve(m, stats):
    m.params.LazyConstraints = 1
    m.optimize(callback)
    stats.ILP += 1
    stats.node += int(m.nodeCount)
    stats.iter += int(m.iterCount)
    #
    objective = int(round(m.getObjective().getValue()))
    status = m.status
    if status == GRB.status.TIME_LIMIT:
        print('Objective when time limit reached:', objective)
        raise TimeLimit()
    assert status == GRB.status.OPTIMAL, status
    assert objective == m._ub

def callback(m, where):
    if where == GRB.Callback.MIPSOL:
        try:
            extend_cm(m)
        except:
            print_exc(file=stderr)
            import os
            os._exit(1)

def extend_cm(m):   
    g, ub, ub2 = m._g, m._ub, int(round(m.cbGet(GRB.Callback.MIPSOL_OBJ))) 
    #if ub2 >= ub:
    #    return
    log('Callback with obj: ', ub2)
    #
    relax_elims, new_feas_elims, missed_loops = get_solution(m)
    if not missed_loops:
        if ub2 < ub:
            log('Relaxation became feasible and improved UB {} -> {}'.format(ub, ub2))
            m._elims, m._ub = relax_elims, ub2
            double_check(g, ub2, relax_elims, is_labeled=True, log=log) 
        return
    #
    new_ub = sum(g[u][v]['weight'] for u,v in new_feas_elims)
    if new_ub < ub:
        log('Improved UB: {} -> {}'.format(ub, new_ub))
        double_check(g, new_ub, new_feas_elims, is_labeled=True, log=log)
        m._ub, m._elims = new_ub, new_feas_elims  
    #
    extend_the_cycle_matrix(m, missed_loops)

def get_solution(m):
    # Returns: relaxed eliminations, new feasible elimination, missed loops.
    relax_elims = [e for e, y in iteritems(m._vrs) if m.cbGetSolution(y) > 0.5]
    #
    g_ruins = m._g.copy()
    for e in relax_elims:
        g_ruins.remove_edge(*e)
    #
    if is_directed_acyclic_graph(g_ruins):
        return relax_elims, relax_elims, set()
    #
    missed = missed_edges(g_ruins)
    assert missed
    new_feas_elims = relax_elims + missed
    missed_loops = {to_cycle(shortest_path(g_ruins,v,u)) for u,v in missed}
    #missed_loops |= extra_loops(g_ruins, m._g)
    return relax_elims, new_feas_elims, missed_loops

def extend_the_cycle_matrix(m, missed_loops):
    loops = m._loops
    #missed_loops -= loops
    vrs = m._vrs
    for edgelist in missed_loops:
        y = [vrs[edge] for edge in edgelist]
        sum_y = LinExpr([1.0]*len(y), y)
        m.cbLazy(sum_y >= 1.0)
    loops |= missed_loops  # The naive way: no subset selection
    log('Added', len(missed_loops), 'new rows to the cycle matrix')

def build_ilp(g, edgelist_per_cycle, elims):
    m = Model()
    binary = GRB.BINARY
    edges = g.edges(data=True)
    # add all variables first, keep vrs= {edge: gurobi var} for the constraints 
    vrs = { (u,v) : m.addVar(vtype=binary,obj=d['weight']) for u,v,d in edges }
    m.update()
    for edgelist in edgelist_per_cycle:
        # add the set covering constraint of each simple cycle
        y = [vrs[edge] for edge in edgelist]
        sum_y = LinExpr([1.0]*len(y), y)
        m.addConstr(sum_y >= 1.0)
    elim_set = set(elims)
    # set the greedy heuristic's solution as initial point
    for e, y in iteritems(vrs):
        y.start = 1.0 if e in elim_set else 0.0
    return m, vrs
