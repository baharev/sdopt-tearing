# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import division, print_function
from contextlib import contextmanager
from time import time
from six import iteritems, itervalues
from gurobipy import LinExpr, GRB, Model, setParam
import networkx as nx
from benchmarks import gen_benchmark_digraphs, digraph_to_undirected_bipartite,\
                       gen_digraphs_as_rectangular_bipartite
from equations import info_on_bipartite_graph, read_bipartite_graph
from heap_md import min_degree
from mfes import run_mfes_heuristic
from order_util import deterministic_topological_sort, permute_to_hessenberg
from test_tearing import gen_testproblems
from utils import edges_of_cycle, rotate_min_to_first_pos, solve_ilp
from plot_ordering import plot_hessenberg

__all__ = [ 'solve_problem' ]

# Forked from grb_tear
#
# Here (and also in grb_tear), unlike in other modules, a loop is just a simple
# path of nodes. Elsewhere, the loops are represented as sequence of edges.

def log(*args, **kwargs):  pass
log = print

def main():
    setParam('LogFile', '/tmp/gurobi.log')
    #setParam('OutputFlag', 0)
    setParam('LogToConsole', 0)
    #
    solve_digraphs_as_rectangular_bipartite()
    real_main(use_min_degree=False)
    real_main(use_min_degree=True)

def solve_digraphs_as_rectangular_bipartite():
    for g, eqs in gen_digraphs_as_rectangular_bipartite():
        forbidden = set()
        solve_problem(g, eqs, forbidden)

def real_main(use_min_degree):
    #pname = 'JacobsenShortestSimpBounds'
    pname = 'JacobsenILOSimpBounds'
    g, eqs, forbidden = read_bipartite_graph(pname)
    #feasible_solution = read_feasible_solution(pname, g, eqs, forbidden)
    #
    solve_problem(g, eqs, forbidden, use_min_degree)
    #
    for g, eqs, forbidden in gen_testproblems():
        info_on_bipartite_graph(g, eqs, forbidden, log=log)
        solve_problem(g, eqs, forbidden, use_min_degree)
    #
    # A bipartite graph, triggering poor performance
    dig = next( gen_benchmark_digraphs() )
    for n in dig.nodes():
        if n > 20:
            dig.remove_node(n)
    g, eqs, forbidden = digraph_to_undirected_bipartite(dig)
    solve_problem(g, eqs, forbidden)

def solve_problem(g, eqs, forbidden, use_min_degree=True): 
    # Returns rowp, colp, matches, tear_set, sink_set in Hessenberg form.
    start = time()
    ret = solve_with_pcm(g, eqs, forbidden, use_min_degree)
    end = time()
    log('Overall solution time: {0:0.1f} s'.format(end-start))
    return ret

def solve_with_pcm(g, eqs, forbidden, use_min_degree):
    eqs = set(eqs)
    # Invariant: match belongs to the best known feasible solution (dag), and 
    # match is (hopefully) a near optimal solution to the relaxation
    dag, tears, ub, match, loops = \
                     create_feasible_solution(g, eqs, forbidden, use_min_degree)
    log('UB <=', ub)
    # Build a shortest path loop around each allowed edge of the tear variables
    loops |= path_around_tears(dag, eqs, forbidden, tears)
    #
    candids, dag, match, ub = step(g, eqs, forbidden, loops, dag, match, ub)
    while candids:
        loops |= candids  # Replace with a greedy heuristic?
        # Try to put a shortest path loop around each erroneously matched edge,
        # and try to improve the lower bound and/or the feasible solution (upper
        # bound).
        candids, dag, match, ub = step(g, eqs, forbidden, loops, dag, match, ub)
    #
    variables = sorted(n for n in dag if n not in eqs)
    # all tears are sources, no other variable is source
    sources = sorted(n for n,indeg in dag.in_degree_iter(variables) if indeg==0)
    rowp, colp, matches, tear_set, sink_set = \
                                 permute_to_hessenberg(g, eqs, forbidden, sources)
    assert ub==len(tear_set), (sorted(sources), sorted(tear_set))
    log()
    log('***  Optimal solution found  ***')
    log('Number of tear variables:    ', len(tear_set))
    log('Number of residual equations:', len(sink_set))
    return rowp, colp, matches, tear_set, sink_set

def create_feasible_solution(g, eqs, forbidden, use_min_degree):
    if use_min_degree:
        method_name = 'minimum-degree ordering'
    else:
        method_name = 'small loop subset selection'
    log('Feasible solution will be generated with', method_name)
    start = time()
    #
    if use_min_degree:
        # TODO Suboptimal: match recomputed; type of tears is either set or list
        rowp, colp, matches, tears, sinks = min_degree(g, eqs, forbidden)
        dag = matching_to_dag(g, eqs, forbidden, rowp, colp, matches, tears, sinks)
        ub = len(tears)
        loops = set()
    else:
        dag, tears, ub, loops = feasible_sol_from_loop_subset(g, eqs, forbidden)
    #
    end = time()
    log('Generating a feasible solution took: {0:0.2f} s'.format(end-start))
    match = get_matching(dag, eqs)
    return dag, tears, ub, match, loops  

def feasible_sol_from_loop_subset(g, eqs, forbidden):
    loops = initial_subset(g, eqs, forbidden)
    match, lb = solve_relaxation(g, eqs, forbidden, loops)
    dig = orient_wrt_matching(g, eqs, match, lb)
    tears = mfes_heuristic(dig, eqs)
    ub = len(tears) + lb
    digraph_to_dag(dig, eqs, tears, lb)
    return dig, tears, ub, loops

def get_matching(dag, eqs):
    # The heap-based min-degree has the matching, so it is somewhat wasteful...
    matching = [ ]
    non_sink_eqs = (eq for eq in sorted(eqs) if dag.succ[eq])
    for eq in non_sink_eqs:
        (var,) = dag.succ[eq]  # exactly one out edge in a valid elimination 
        matching.append( (eq,var) )
    return matching

def path_around_tears(dag, eqs, forbidden, tears):
    # Build a shortest path loop around each allowed edge of the tear variables.
    # A loop is a candidate iff *all* of its edges would go into an allowed 
    # direction.
    assert nx.is_directed_acyclic_graph(dag)
    loops = set()
    for var in tears:
        edges = [(eq,var) for eq in dag[var] if (eq,var) not in forbidden]
        for eq,var in edges:
            # Removing var -> eq is necessary in all cases
            with removed_edge(dag, var, eq):
                loops.update( loops_around_edge(dag, eqs, eq, var) )
    return loops

@contextmanager
def removed_edge(dag, u, v):
    dag.remove_edge(u,v)
    try:
        yield
    finally:
        dag.add_edge(u,v)

def loops_around_edge(dag, eqs, eq, var):
    assert not dag.has_edge(eq, var), 'var is supposed to be a tear'
    loops = set()
    # First, remove none of the neighbors of eq or var:
    add_a_loop(dag, eqs, eq, var, loops)
    # TODO Unclear which loop generation strategy is "the best". This one tries 
    # to put more than one loop around the edges by systematically removing
    # edges incident to eq or var or both.
    #
    ## Only neighbor of eq is removed:
    #for eq_nbr in dag.predecessors(eq):
    #    with removed_edge(dag, eq_nbr, eq):
    #        add_a_loop(dag, eqs, eq, var, loops)
    ## Try removing neighbors of var or both var and eq
    #for var_nbr in dag.successors(var):
    #    with removed_edge(dag, var, var_nbr):
    #        # Only neighbor of var is removed:
    #        add_a_loop(dag, eqs, eq, var, loops)
    #        # Both neighbor of eq and var removed:
    #        for eq_nbr in dag.predecessors(eq):
    #            with removed_edge(dag, eq_nbr, eq):
    #                add_a_loop(dag, eqs, eq, var, loops)
    #
    return loops

def add_a_loop(dag, eqs, eq, var, loops):
    # create a small cycle around edge eq -> var with a shortest path var ~> eq
    try:
        simple_path = nx.shortest_path(dag,var,eq)
        loops.add( to_normalized_path(simple_path, eqs) )
    except nx.NetworkXNoPath:
        pass # That's OK

def to_normalized_path(simple_path, eqs):
    # Only for bipartite graphs. Rotate the equation with the smallest id into 
    # the first position, then traverse the loop, starting with that neighbor 
    # variable that has the smaller id.
    equations = [ (n,i) for i,n in enumerate(simple_path) if n in eqs ]
    idx = min(equations)[1]
    left, right = simple_path[idx-1], simple_path[(idx+1) % len(simple_path)]
    step = -1 if left < right else 1
    path = simple_path[idx::step] + simple_path[:idx:step]
    #--- TODO Remove when done, only for debugging
    assert sorted(simple_path)==sorted(path)
    inarg = rotate_min_to_first_pos(simple_path)
    out_1 = rotate_min_to_first_pos(path)
    out_2 = rotate_min_to_first_pos(list(reversed(path)))
    assert inarg==out_1 or inarg==out_2
    #---
    return tuple(path)

def step(g, eqs, forbidden, loops, feas_dag, match, ub):
    # loops give a partial cycle matrix, match a feasible solution
    relax_matching, lb = solve_relaxation(g, eqs, forbidden, loops, match)
    assert lb <= ub
    if lb == ub:
        log('The best feasible solution is proved to be optimal, LB = UB =', ub)
        return [ ], feas_dag, match, ub  # feas_dag is now proved to be optimal
    #
    # Check whether the relax_matching gives an acyclic orientation of the 
    # original undirected graph
    dig = orient_wrt_matching(g, eqs, relax_matching, lb)
    if nx.is_directed_acyclic_graph(dig):
        log('The relaxation became feasible, meaning UB = LB =', lb)
        return [ ], dig, relax_matching, lb
    #
    log('LB >= {}  (UB <= {})'.format(lb, ub))
    # Neither lb==ub nor acyclic; adding the missed loops
    tears = mfes_heuristic(dig, eqs)
    assert tears
    digraph_to_dag(dig, eqs, tears, lb) # checks: new_ub == len(tears) + lb
    candidates = path_around_tears(dig, eqs, forbidden, tears)
    candidates -= loops
    assert candidates, 'We are stuck in an infinite loop'
    #plot_relax_solution(g, dig, eqs, forbidden, candidates)
    #
    # Has the UB improved?
    new_ub = len(tears) + lb  # already double-checked by digraph_to_dag 
    if new_ub < ub:
        log('Improved UB: {} -> {}'.format(ub, new_ub))
        feas_dag, match, ub = dig, get_matching(dig, eqs), new_ub
    return candidates, feas_dag, match, ub

def plot_relax_solution(g, dag, eqs, forbidden, candidates):
    # Also marks the missed loops (candidates) in red.
    tears, _, order = get_solution(dag, eqs, forbidden)
    rowp = [ r for r in order if r in eqs ]
    rindex = { name : i for i, name in enumerate(rowp) }
    #
    colp = sorted(tears, key=lambda c: min(rindex[r] for r in g[c]))
    seen = set(colp)
    colp.extend( c for c in order if c not in eqs and c not in seen )
    # Or, if we want the spiked form, uncomment the lines below:
    #from order_util import get_row_col_perm
    #rowp, colp = get_row_col_perm(eqs, dag, tears, sinks, order)
    #
    cindex = { name : i for i, name in enumerate(colp) }
    mark_red = [ ]
    for loop in candidates:
        for a, b in edges_of_cycle(loop):
            r, c = (a, b) if a in eqs else (b, a)
            mark_red.append( (rindex[r], cindex[c]) )
    #
    msg, partitions = '', [ ]
    plot_hessenberg(g, rowp, colp, partitions, msg, mark_red)

def solve_relaxation(g, eqs, forbidden, loops, match=None):
    log()
    log('The cycle matrix has', len(loops), 'rows')
    m, y = build_ilp(g, eqs, forbidden, loops, match)
    #dump(m, '.lp', increment=True)
    success = solve_ilp(m)
    #dump(m, '.sol')
    assert success, 'Solver failures are not handled'
    objective = int(round(m.getObjective().getValue()))
    n_vars, n_matched_vars = len(g)-len(eqs), objective  
    lb = n_vars - n_matched_vars # unmatched vars are tear variables
    matches = \
           sorted(edge for edge, var in iteritems(y) if int(round(var.x))==1)
    #log('matches:', matches)
    return matches, lb

def dump(m, extension, increment=False):
    m.update()
    if not hasattr(dump, 'counter'): # emulating a static variable
        dump.counter = 0
    if increment:
        dump.counter += 1
    m.write('/tmp/relaxation_'+str(dump.counter)+extension)

def orient_wrt_matching(g_orig, eqs, relax_matching, lb):
    dig = nx.DiGraph()
    dig.add_nodes_from(g_orig)
    matched_edges = { eq_var for eq_var in relax_matching }
    for eq_var in g_orig.edges_iter(eqs):
        u, v = eq_var if eq_var in matched_edges else (eq_var[1], eq_var[0])
        dig.add_edge(u, v)
    variables = sorted(n for n in g_orig if n not in eqs)
    # all tears are sources, no other variable is source
    sources = sorted(n for n,indeg in dig.in_degree_iter(variables) if indeg==0)
    assert len(variables) == len(sources) + len(relax_matching)
    assert lb == len(sources)
    return dig

def matching_to_dag(g_orig, eqs, forbidden, rowp, colp, matches, tears, sinks):
    matched_edges = set(edge for edge in iteritems(matches) if edge[0] in eqs)
    len_matches = len(matched_edges)
    assert not (matched_edges & forbidden)
    
    dag = nx.DiGraph()
    dag.add_nodes_from(rowp) # Empty (isolated) equations are allowed
    #dag.add_nodes_from(variables)
    for eq_var in g_orig.edges_iter(rowp):
        u, v = eq_var if eq_var in matched_edges else (eq_var[1], eq_var[0])
        dag.add_edge(u, v)
        matched_edges.discard(eq_var)
    
    assert not matched_edges
    # FIXME Comparing str and int breaks on Py 3
    has_all_nodes = sorted(dag, key=str) == sorted(g_orig, key=str)
    assert has_all_nodes # Isolated (degree zero) var nodes?
    assert nx.is_directed_acyclic_graph(dag)

    # Check whether the matching is sane
    assert len_matches == len(eqs) - len(sinks)
    assert len_matches == len(g_orig) - len(eqs) - len(tears)    
    
    more_than_one_outedge = [ eq for eq in rowp if len(dag.succ[eq]) > 1 ] 
    assert not more_than_one_outedge, more_than_one_outedge
    
    more_than_one_inedge = [var for var in colp if len(dag.pred[var]) > 1]
    assert not more_than_one_inedge, more_than_one_inedge
    
    return dag

def mfes_heuristic(dig, eqs):
    objective, elims = run_mfes_heuristic(dig, try_one_cut=True)
    ##--- TODO Hack here to run with grb_pcm instead of the heuristic
    #from grb_pcm import solve_problem as mfes_rigorous
    #g = nx.DiGraph()
    #for u, v in dig.edges_iter():
    #    g.add_edge(u, v, { 'weight' : 1, 'orig_edges' : [ (u,v) ] })
    #elims, objective = mfes_rigorous(g)
    ##---
    tears = sorted( v if v not in eqs else u for u,v in elims )
    assert all(v not in eqs for v in tears)
    assert objective == len(tears)
    log('It is still necessary to guess', objective, 'variables')
    if objective <= 5:
        log(tears)
    return tears

def digraph_to_dag(dig, eqs, tears, lb):
    for v in tears:
        eq_nbrs = list(dig.pred[v])
        for e in eq_nbrs:
            dig.remove_edge(e, v)
            dig.add_edge(v, e)
    # all tears are sources, no other variable is source
    variables = sorted(n for n in dig if n not in eqs)
    sources = sorted(n for n,indeg in dig.in_degree_iter(variables) if indeg==0)
    # double-check the upper bound:
    assert len(sources)==len(tears)+lb

def get_solution(dag, eqs, forbidden):
    # Somewhat duplicate of get_final_order in min_degree
    equations = sorted(eqs)
    variables = sorted(n for n in dag if n not in eqs)
    # all tears are sources, no other variable is source
    sources = sorted(n for n,indeg in dag.in_degree_iter(variables) if indeg==0)
    # all residuals are sinks, no other equation is sink
    sinks   = sorted(n for n, out in dag.out_degree_iter(equations) if out==0)
    # do a topological sort;  
    # nbunch tries to break ties by ordering the equations alphabetically.
    nbunch = list(reversed(equations))
    nbunch.extend(reversed(variables))
    order = deterministic_topological_sort(dag, nbunch)
    # Sanity check: each equation is matched with a single allowed variable
    resids = set(sinks)
    non_sink_eqs = (eq for eq in order if eq in eqs and eq not in resids)
    for eq in non_sink_eqs:
        # eq must have exactly one allowed out edge in a valid elimination
        (var,) = dag.succ[eq] 
        assert (eq,var) not in forbidden, (eq,var)
    return sources, sinks, order

#-------------------------------------------------------------------------------

def build_ilp(g, eqs, forbidden, simple_cycles, match):
    m = Model()
    m.setAttr('ModelSense', GRB.MAXIMIZE)
    edges = sorted(g.edges(eqs))
    # y: is edge e in the matching? The worth of each allowed match is 1. 
    obj_coeff = { e : (1 if e not in forbidden else 0) for e in edges }
    y = { e : m.addVar(vtype=GRB.BINARY, obj=obj_coeff[e]) for e in edges }
    m.update()
    # An equation is matched at most once
    for eq in eqs:
        add_matched_at_most_once_con(m, y, g.edges(eq))
    # A variable is matched at most once
    set_eqs = set(eqs)
    vrs = sorted(v for v in g if v not in set_eqs)
    for v in vrs:
        add_matched_at_most_once_con(m, y, list((eq,v) for v,eq in g.edges(v)))
    # 
    break_each_cycle_at_least_once(eqs, m, y, simple_cycles)
    #
    set_lower_bound_if_any(g, eqs, m, y)
    #
    # Set feasible solution from the matching, if any
    if match is not None:
        match = set(match)
        for edge, var in iteritems(y):
            var.setAttr('Start', 1 if edge in match else 0)
    return m, y

def add_matched_at_most_once_con(m, y, edges):
    y_in_con = [ y[e] for e in edges ]
    lhs = LinExpr([1.0]*len(y_in_con), y_in_con)
    m.addConstr(lhs, GRB.LESS_EQUAL, 1.0)

def break_each_cycle_at_least_once(eqs, m, y, simple_cycles):
    equations = set(eqs)
    for cyc in simple_cycles:
        y_in_con = get_all_y_in_cycle(cyc, y, equations)
        n_ys = len(y_in_con)
        assert n_ys % 2 == 0, n_ys
        matched_max = n_ys//2 - 1
        lhs = LinExpr([1.0]*n_ys, y_in_con)
        m.addConstr(lhs, GRB.LESS_EQUAL, matched_max) 

def get_all_y_in_cycle(cyc, y, equations):
    y_in_con = [ ]
    for u,v in edges_of_cycle(cyc):
        # y has (eq,var) as keys: Swap is necessary if u,v is in (var,eq) order
        e = (u,v) if u in equations else (v,u)
        y_in_con.append(y[e])
    return y_in_con

def set_lower_bound_if_any(g, eqs, m, y):
    lb = g.graph.get('trivial_lb')
    if lb is None:
        return
    lhs = LinExpr([1]*len(y), list(itervalues(y)))
    # See solve_relaxation: objective is to maximize the sum of matched 
    # variables, the lb is on the unmatched variables: n_vars - lb is an upper 
    # bound on the variables that can be matched at most
    n_vars = len(g)-len(eqs) 
    m.addConstr(lhs, GRB.LESS_EQUAL, n_vars - lb)

#-------------------------------------------------------------------------------
# TODO Cleanup this code if it turns out that we actually need it. As it stands,
# it is inferior to the min-degree ordering.

def initial_subset(g, eqs, forbidden):
    initial_loops = build_small_loops(g, eqs, forbidden)
    # Remove those small loops that have forbidden edges anyway? Although
    # it is very likely that the presolve in Gurobi will throw them away anyway
    # so maybe it would be wasted developer time, not sure...
    loops = loops_from_path_to_edge_repr(initial_loops, eqs)
    log('Initially we have:', len(loops), 'loops')
    # Subset selection with edge representation
    start = time()
    subset = select_subset_of_loops(g, eqs, loops, max_coverage=2)
    end = time()
    log('Number of selected loops:', len(subset))
    log('Selection computed in {0:0.1f} s'.format(end-start))
    # Convert the subset back to simple path representation
    return convert_loops_from_edge_to_path_repr(subset)

def build_small_loops(g, eqs, forbidden):
    # Build small loops around each edge
    start = time()
    initial_loops = set()
    for eq,var in g.edges_iter(eqs):
        if (eq,var) in forbidden:
            continue
        g.remove_edge(eq,var)
        try:
            simple_path = nx.shortest_path(g,var,eq)
            initial_loops.add( to_normalized_path(simple_path, eqs) )
        except nx.NetworkXNoPath:
            pass # That's OK   
        g.add_edge(eq,var)
    end = time()
    log('Small loops computed in {0:0.1f} s'.format(end-start))
    return initial_loops

def loops_from_path_to_edge_repr(initial_loops, eqs):
    # Convert simple path representation to edge representation 
    start = time()
    loops = [ ]
    # The similar utils function is not good for bipartite setting
    for l in initial_loops:
        loops.append( tuple( (u,v) if u in eqs else (v,u) 
                              for u,v in edges_of_cycle(l) ) )
    end = time()
    log('Edge representation computed in {0:0.1f} s'.format(end-start))
    return loops

def convert_loops_from_edge_to_path_repr(loops_edge_repr):
    # Convert the loops back to simple path representation
    loops_path_repr = set()
    for loop in loops_edge_repr:
        spath = [ ]
        for u, v in loop[::2]:
            spath.append(u)
            spath.append(v)
        loops_path_repr.add( tuple(spath)  )
        # build_small_loops already normalized path
        # assert tuple(spath) == to_normalized_path(spath, eqs)
    return loops_path_repr

# TODO Mostly a duplicate of select_subset_of_loops in grb_pcm but adjusted to 
# bipartite graphs
def select_subset_of_loops(g, eqs, small_loops, max_coverage=1):
    # Select a maximum set of loops such that each edge participates in at most
    # max_coverage loops. The max_coverage=1 means the loops are
    # independent, they do not share edges. The max_coverage=2 means the 
    # selected loops share at most 1 edge, in other words, loops share edges at
    # most pairwise.
    m = Model()
    m.setAttr('ModelSense', GRB.MAXIMIZE)
    # * Bipartite, apart from that, duplicate of grb_pcm select_subset_of_loops*
    edges = g.edges(eqs)
    loop_vars = {loop:m.addVar(vtype=GRB.BINARY,obj=1) for loop in small_loops}
    edge_vars = { edge : m.addVar(vtype=GRB.BINARY) for edge in edges }
    m.update()
    # An edge can participate in at most max_coverage loops
    # TODO Takes very long to build the MILP; most likely e in loop which 
    # searches in linear time in the loop tuple.
    start = time()
    for e in edges:
        in_loops = [var for loop,var in iteritems(loop_vars) if e in loop]
        if in_loops:
            lhs = LinExpr([1]*len(in_loops), in_loops)
            m.addConstr(lhs, GRB.LESS_EQUAL, max_coverage)
    end = time()
    log('Building the ILP model took {0:0.1f} s'.format(end-start))
    # If a loop is chosen, all of its edges are chosen
    for cycle in small_loops:
        for edge in cycle:
            m.addConstr(edge_vars[edge], GRB.GREATER_EQUAL, loop_vars[cycle])
    success = solve_ilp(m)
    assert success, 'Solver failures are not handled'
    objective = int(round(m.getObjective().getValue()))
    loop_subset = {l for l,var in iteritems(loop_vars) if int(round(var.x))==1}
    assert len(loop_subset)==objective
    #log()
    #log('Number of rows in the cycle matrix:', objective)
    return loop_subset

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
