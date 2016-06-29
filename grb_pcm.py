# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from itertools import chain
from time import time
import six
import networkx as nx
from benchmarks import gen_benchmark_digraphs
from grb_set_cover import solve_cm
from grb_simplifier import iteratively_remove_runs_and_bypasses
from mfes import noncopy_split_to_nontrivial_sccs, run_mfes_heuristic
from utils import info_short as info, deserialize, DATADIR
from utils import double_check, to_cycle

def log(*args, **kwargs):  pass
#    print('+ ', *args, **kwargs) if args or kwargs else print()

# Stats only contains counters that will be mutated in-place if passed as the 
# keyword argument stats.

def solve_problem(g_orig, stats=None):
    'Returns: [torn edges], cost.'
    elims, cost = [ ], 0
    for sc in noncopy_split_to_nontrivial_sccs(g_orig.copy()): # <- Copy passed!
        partial_elims, partial_cost = solve_with_pcm(sc, stats)
        elims.extend(partial_elims)
        cost += partial_cost
    double_check(g_orig, cost, elims, is_labeled=True, log=log)
    log('Input graph')
    info(g_orig, log=log)
    return elims, cost 

def solve_with_pcm(g, stats=None):
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
    elims, ub = feasible_solution(g) 
    # Build a shortest path loop around each edge
    # TODO Further improvements: (1) Try to build more than one loop around 
    # those edges that the greedy MFES heuristic would favor. (2) Use the 
    # select_subset_of_loops to add the independent loops. (3) Better 
    # simplification than the iteratively_remove_runs_and_bypasses.
    loops = initial_loop_set(g)
    # The inefficient greedy heuristic, see the comments at its implementation
    #loops = append_loops_greedily(g, set(), loops, ub)
    missed, ruins, elims, ub = step(g, loops, elims, ub, stats)
    while missed:
        # Put a shortest path loop around each missed edge, and try to improve
        # the lower bound and/or the feasible solution (upper bound). 
        candidates = { to_cycle(nx.shortest_path(ruins,v,u)) for u,v in missed }
        # The inefficient greedy heuristic:
        #loops = append_loops_greedily(g, loops, candidates, ub)
        loops = loops | candidates  # <- The naive way: no subset selection
        #
        missed, ruins, elims, ub = step(g, loops, elims, ub, stats)
    #run_IIS(g, loops, ub)
    # Done. Now undo the  d['orig_edges'] mess.
    final_elims = [ ]
    for edge in elims:
        final_elims.extend( origedges_map[edge] )
    return final_elims, ub

def get_orig_edges_map(g):
    # The followings is only needed so that double_check succeeds in this module
    # and in the mfes heuristic
    edge_map = { }
    for u,v,d in g.edges_iter(data=True):
        edge_map[(u,v)] = d['orig_edges']
        d['orig_edges'] = [(u,v)]
    return edge_map

def feasible_solution(g):
    # Only an upper bound estimate is needed, heuristic suffices
    objective, elims = run_mfes_heuristic(g, try_one_cut=True, is_labeled=True)
    #error_msg, elims, obj, _ = rigorous_mfes(g, CUTOFF)
    #assert not error_msg, error_msg
    return  elims, objective

def initial_loop_set(g):
    # Build an initial loop set: we create small cycles around each edge 
    # u -> v with a shortest path v ~> u (we are supposed to be in an SCC) 
    small_loops = set()
    #for cyc in nx.simple_cycles(g):  # <- Uncomment if all loops are needed 
    #    small_loops.add(to_cycle(cyc))
    for u,v in g.edges_iter():
        path = nx.shortest_path(g, source=v, target=u)
        small_loops.add( to_cycle(path) )
    log(len(small_loops), 'small loops')
    return small_loops

def step(g, loops, feas_elims, ub, stats):
    # Solve relaxation with the cycle matrix of loops. This solution 
    # MUST BE rigorous.
    relax_elims, lb = solve_relaxation(g, loops, stats)
    assert lb <= ub
    if lb == ub:
        log('***  Optimal solution found  ***')
        return None, None, feas_elims, ub
    #
    g_ruins = g.copy()
    for e in relax_elims:
        g_ruins.remove_edge(*e)
    #
    if nx.is_directed_acyclic_graph(g_ruins):
        log('***  Relaxation became feasible  ***')
        return None, None, relax_elims, lb,
    #
    log()
    log('Remaining graph')
    info(g_ruins, log=log)
    #
    # Get the missed edges and the ruins of the relaxation: We need new 
    # candidate loops. The missed_edges does not have to be rigorous. We may 
    # also improve the currently best feasible solution. 
    missed = missed_edges(g_ruins)
    assert missed
    new_feas_elims = relax_elims + missed
    new_ub = sum(g[u][v]['weight'] for u,v in new_feas_elims)
    assert new_ub > lb
    if new_ub < ub:
        log('Improved UB: {} -> {}'.format(ub, new_ub))
        ub, feas_elims = new_ub, new_feas_elims 
        double_check(g, ub, feas_elims, is_labeled=True, log=log)
    return missed, g_ruins, feas_elims, ub

def solve_relaxation(g, loops, stats):
    relax_elims, lb = solve_cm(g, loops, stats)
    assert relax_elims is not None, 'Solver failures are not handled'
    log()
    log('LB >=', lb)
    log('Cycle matrix size:', len(loops))
    return relax_elims, lb

def missed_edges(ruins):
    # We have a choice here: If there are not too many simple cycles, we can 
    # solve the remaining graph rigorously. Otherwise, we can only call the 
    # heuristic.
    # missed: additional edges that we had to take out to make ruins acyclic
    #_, missed, cost, _ = rigorous_mfes(ruins, CUTOFF) # <- Hack
    cost, missed = run_mfes_heuristic(ruins, try_one_cut=True, is_labeled=True)
    log('Cost:', cost)
    return missed

#===============================================================================
#
# Below are attempts aiming at greedy loop selections. They did not prove to be
# efficient / fast enough.
#
#-------------------------------------------------------------------------------
# This block is here to test the iteration in step (when we add loops around 
# tears). The need_more_loops can be hacked accordingly. This greedy heuristic 
# is otherwise obsolete and is subject to removal.
#
# Loop subset selection seemed like a good idea. However, as it is implemented 
# in Python, it runs slower than just running Gurobi and leaving up to its 
# presolve phase to throw out the unnecessary loops.
#
# In particular the coverage computation is slow (25-30% of the time, the 
# coverage = [ len(selected & bipart[e].viewkeys()) for e in edges ] line), and 
# the need_more_loops call (65-70% of the time), but not because of Gurobi. 
# Profiling shows that just querying the objective value is slower than solving
# the ILP.

def append_loops_greedily(g, already_selected, candidate_loops, ub_current):
    # candidate_loops were created by putting loops around the tears. None of 
    # these loops can be in the already_selected (it is also asserted below).
    # Adding just one candidate loop is sufficient to make progress. So even a 
    # poor heuristic would not fail: If need_more_loops always returns False,
    # the algorithm will still reach convergence.  
    #
    # The loops are first mapped to integers; selected and candidates are sets 
    # of these integers. The node ids of the loops in the bipartite graph are 
    # these integers too, use nodeid_loop to map back the integer to the loop.
    selected, candidates, bipart, nodeid_loop = \
                      setup_loop_selection(g, already_selected, candidate_loops)
    ub_not_reached = True
    while candidates and ub_not_reached:
        scores = [ ]
        for loop in candidates:
            assert loop not in selected
            edges = bipart[loop] # edges participating in loop
            # coverage = [ len(selected & bipart[e].viewkeys()) for e in edges ]
            coverage = [len(selected & set(bipart[e])) for e in edges] # Py3
            coverage.sort(reverse=True)
            itr = chain(coverage, [0, 0])
            # score: (most shared, second most shared, number of edges, id)
            scores.append( (next(itr), next(itr), len(edges), loop) )
            # TODO edge weights are ignored: use (shared / edge weight)?
        most_shared, second_most_shared, n_edges, loop = min(scores)
        log('Selected score:', (most_shared, second_most_shared, n_edges))
        selected.add(loop)
        candidates.remove(loop)
        loop_array = [ nodeid_loop[i] for i in selected ]
        ub_not_reached = need_more_loops(g, loop_array, ub_current)
    return { nodeid_loop[i] for i in selected }

def setup_loop_selection(g, already_selected, candidate_loops):
    # map the loops to integers, build a bipartite graph of loops -- edges
    all_loops = already_selected | candidate_loops
    candidate_loops = candidate_loops - already_selected
    nodeid_loop = { i : loop for i, loop in enumerate(all_loops) }
    loop_nodeid = { loop : i for i, loop in enumerate(all_loops) }
    bipart = nx.Graph()
    for nodeid, loop in six.iteritems(nodeid_loop):
        for u,v in loop:
            bipart.add_edge((u,v), nodeid, {'weight': g[u][v]['weight']})
    candidates = { loop_nodeid[l] for l in candidate_loops }
    selected   = { loop_nodeid[l] for l in already_selected }
    return selected, candidates, bipart, nodeid_loop

def need_more_loops(g, selected, ub_current):
    # We could run the select most shared edge or the mfes heuristic first, and
    # if that already indicates that we need more loops, we could return early.
    # In fact, rigorous solution is not needed here at all.
    # Gurobi already implements some greedy heuristic, so not calling Gurobi
    # would be wasted developer time.
    _, lb = solve_cm(g, selected) # <- unnecessarily rigorous
    assert lb is not None, 'Solver failures are not handled'
    log('LB >=', lb)
    return False # TODO Hacked here to test the iteration
    return lb < ub_current

#===============================================================================
# For the paper to demonstrate that the proposed method works on the notoriously 
# difficult problem, where the set cover and the LOP formulation fails.

def test_on_Jaconbsen_with_50_stages():
#     # Just to document how the test DAG was created
#     from equations import read_bipartite_graph
#     g, eqs, forbidden = read_bipartite_graph('JacobsenILOSimpBounds')
#     for eq, var in g.edges_iter(eqs):
#         g[eq][var]['weight'] = 1 if (eq, var) in forbidden else 10
#     mate = nx.max_weight_matching(g, maxcardinality=True)
#     # Orient according to the matching, and also label
#     dig = nx.DiGraph()
#     for eq, var in g.edges_iter(eqs):
#         if mate[eq]==var:
#             dig.add_edge(eq, var, weight=1, orig_edges=[(eq,var)])
#         else: 
#             dig.add_edge(var, eq, weight=1, orig_edges=[(var,eq)])
#     assert not nx.is_directed_acyclic_graph(dig)
#     from utils import serialize
#     serialize(dig, 'data/JacobsenILOSimpBounds_as_DAG.pkl.gz')
    dig = deserialize(DATADIR+'JacobsenILOSimpBounds_as_DAG.pkl.gz')
#     # Uncomment to prove that this graph has more than 10M simple cycles: 
#     cutoff = 10000000
#     from itertools import islice
#     n_cycles = sum(1 for _ in islice(nx.simple_cycles(dig), cutoff+1))
#     if n_cycles == cutoff+1:
#         print('More than', cutoff, 'simple cycles, giving up...')
#     else:
#         print('There are', n_cycles, 'simple cycles in total')
    _, cost = solve_problem(dig) 
    print('Cost with ILP:', cost) # 107 with this matching; optimal tearing 53
    cost, _ = run_mfes_heuristic(dig, try_one_cut=True, is_labeled=True)    
    print('Cost with heuristic:', cost) # 160

def run_IIS(g, loops, ub):
    from gurobipy import GRB, LinExpr, setParam
    from grb_set_cover import build_ilp
    #setParam("OutputFlag", 1)
    setParam("IISMethod", 0)
    #success, loops = get_all_cycles(g, cutoff=13747) # if all loops are needed
    #assert success, 'Too many cycles'
    loops = sorted(loops) # loops needs to be a list, and sorting for stability
    m, vrs = build_ilp(g, loops)
    a, y = [ ], [ ]
    for u,v,d in g.edges_iter(data=True):
        a.append(d['weight']), y.append(vrs[u,v])
    lhs = LinExpr(a, y)
    m.addConstr(lhs, GRB.LESS_EQUAL, ub-1)
    #m.update()
    #m.write('JacobsenInfeas.lp')
    m.computeIIS()
    keep = [i for i,cn in enumerate(m.getConstrs()) if cn.getAttr('IISCONSTR')]
    assert keep[-1] == len(loops) # The last one is the dummy constraint
    keep.pop()
    #print(keep)
    cyc_subset = [loops[i] for i in keep]
    #m, _ = build_ilp(g, cyc_subset)
    #m.update()
    #m.write('JacobsenMinimal.lp')
    _, obj = solve_cm(g, cyc_subset)
    print('Obj:', obj)    

#-------------------------------------------------------------------------------

def feas_relax(g, loops, ub):
    from gurobipy import Column, GRB, LinExpr
    from grb_set_cover import build_ilp
    #
    loops = sorted(loops) # loops needs to be a list, and sorting for stability
    m, vrs = build_ilp(g, loops)
    a, y = [ ], [ ]
    for u,v,d in g.edges_iter(data=True):
        a.append(d['weight']), y.append(vrs[u,v])
    lhs = LinExpr(a, y)
    m.addConstr(lhs, GRB.EQUAL, ub)
    m.update()
    m.setObjective(0.0)
    # add slack variables, except the last artificial constraint
    for c in m.getConstrs()[:-1]:
        m.addVar(obj=-1.0, lb=0.0, ub=1.0, column=Column([-1], [c]))
    m.optimize()

#-------------------------------------------------------------------------------

def main():
    #test_on_Jaconbsen_with_50_stages()
    #return
    start = time()
    for g_input in gen_benchmark_digraphs():
        solve_problem(g_input)
    test_on_Jaconbsen_with_50_stages()
    end = time()
    print('Overall solution time: {0:0.1f} s'.format(end-start))


if __name__ == '__main__':
    main()
