# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from copy import deepcopy
from gurobipy import setParam
import networkx as nx
from benchmarks import gen_benchmark_digraphs
from grb_set_cover import rigorous_mfes
from mfes import info_after_cleanup
from plot import dummy as plot
from py3compat import irange
from utils import info, double_check, split_to_nontrivial_sccs

class NoProgress(Exception):
    pass

def main():
    setParam("LogFile", "/tmp/gurobi.log")
    for g in gen_benchmark_digraphs():
        #if g.graph['name'] == 'Problem 10 (opt=12)':
        #if g.graph['name'] != 'Problem 8 (opt=5)':
        #if g.graph['name'] != 'Subproblem 8 (opt=3)':
        #    continue
        STATISICS.clear()
        # Problem 10
        #-----------
        # Sparse cut
        #g.remove_edge(38, 36)
        #g.remove_edge(39, 34)
        #-----------
        # Peeling off the DAG of the bypasses
        #dag_connecting = [ (11,4), (17,16), (1,20), (26,25), (12,29), (35,34),
        #                   (21,38), (44,43), (30,48), (54,53), (39,57) ]
        #for e in dag_connecting:
        #    g.remove_edge(*e)
        #-----------
        # Breaking the small bypasses
        #bypass_breakers = [ (17,11), (26,20), (35,29), (44,38), (54,48) ] 
        #for e in bypass_breakers:
        #    g.remove_edge(*e)
        #-----------
        # Problem 8
        #g.remove_edge(4,6)
        #g.remove_edge(24,16)
        #g.remove_edge(17,22)
        #g.remove_edge(21,20)
        #g.remove_edge(2,6)
        #g.remove_edge(24,26)
        #iteratively_remove_runs_and_bypasses(g)
        #g.remove_edge(2,17)
        #g.remove_edge(8,27)
        #g.remove_edge(27,7)
        #g.remove_edge(8,14)  # These latter two are insane but still OK
        #g.remove_edge(18,21)
        #plot(g, prog='sfdp')
        #-----
        #from utils import deserialize
        #g = deserialize('data/JacobsenILOSimpBounds_as_DAG.pkl.gz')
        #orig_input = deserialize('data/JacobsenILOSimpBounds_as_DAG.pkl.gz')
        # And toggle comment on the # <- lines !!!
        #-----
        orig_input = deepcopy(g) # <-
        simplify(g)
        print('Info on the original input:')
        print(orig_input.graph['name']) # <-
        print('Nodes:', orig_input.number_of_nodes())
        print('Edges:', orig_input.number_of_edges())
        print('Loops:', len(list(nx.simple_cycles(orig_input)))) # <-
        if 'max cycles' in STATISICS:
            print('Loop budget:', STATISICS['max cycles'])
            print('Cutoff:', STATISICS['max cutoff'])
        plot(orig_input, prog='sfdp')
        #break # <-

MAX_CYCLES = 187
BFS_CUTOFF = 5
STATISICS = { }

def update_cyc_max(n_cyc):
    STATISICS['max cycles'] = max(n_cyc, STATISICS.get('max cycles', n_cyc))

def simplify(g):
    print('Nodes:', g.number_of_nodes())
    print('Edges:', g.number_of_edges())
    cycles = list(nx.simple_cycles(g)) # <-
    print('Loops:', len(cycles))       # <- 
    orig_input = deepcopy(g)
    sccs = split_to_nontrivial_sccs(g)
    # Clean-up the new smaller SCCs without splitting them
    for sc in sccs:
        iteratively_remove_runs_and_bypasses(sc)
    #
    running_cost, elims = 0, [ ]
    info_after_cleanup(sccs, running_cost)
    #
    clean_sccs = [ ]
    for scc in sccs:
        plot(scc, prog='sfdp')
        new_sccs, running_cost = try_neighborhood(scc, running_cost, elims)
        clean_sccs.extend( new_sccs )
    # check what remains
    if clean_sccs:
        info_after_cleanup(clean_sccs, running_cost)
        print(sorted(n for n in clean_sccs[0]))
        dbg_dump_as_edgelist(clean_sccs[0])
    else:
        print('The simplifier eliminated the whole input at cost', running_cost)
        double_check(orig_input, running_cost, elims)
        print('Chosen:', elims)

def try_neighborhood(sc, running_cost, elims):
    sccs_to_process = [ sc ]
    rejected = { }
    print('\n*** Looking for sub-SCCs that are safe to eliminate ***\n')
    for cutoff in irange(1, BFS_CUTOFF+1):
        dirty = sccs_to_process
        sccs_to_process = [ ]
        while dirty:
            scc = dirty.pop()
            progressed, cost, new_sccs, rejected_relaxations = \
                                               try_each_node(scc, elims, cutoff)
            if progressed:
                running_cost += cost
                dirty.extend(new_sccs)
            else:
                rejected.update( rejected_relaxations )
                sccs_to_process.append(scc)
            print('-----------------------------------------------------------')
    return sccs_to_process, running_cost

def try_each_node(scc, elims, cutoff):
    rejected = { }
    for n in sorted(scc):
        print('------------------')
        print()
        try:
            safe_elims, cost = try_node(scc, n, cutoff, rejected)
        except NoProgress as e:
            print(e, '\nGiving up on node', n, '\n')
            continue
        # Record the eliminations
        assert safe_elims
        cost = 0
        edges_to_tear = [ (u,v,scc[u][v]) for u,v in safe_elims ]
        for u, v, d in edges_to_tear:
            elims.extend(d['orig_edges'])
            cost += d['weight']
            scc.remove_edge(u, v)
        # Split to SCCs
        new_sccs = split_to_nontrivial_sccs(scc)
        # Clean-up the new smaller SCCs without splitting them
        for new_sc in new_sccs:
            iteratively_remove_runs_and_bypasses(new_sc)
        return True, cost, new_sccs, rejected
    assert scc
    return False, None, None, rejected

def try_node(scc_orig, source, cutoff, rejected):
    STATISICS['max cutoff'] = cutoff
    print('Testing neighbors of node {}, cutoff = {}'.format(source, cutoff))
    neighborhood, _ = get_bfs_neighborhood(scc_orig, source, cutoff)
    #neighborhood = bfs_w_local_search(scc_orig, source, cutoff)
    #neighborhood = simple_path_neighborhood(scc_orig, source, cutoff)
    return simplify_neighborhood(scc_orig, neighborhood, rejected)

def simplify_neighborhood(scc_orig, neighborhood, rejected):
    g_sub = deepcopy(scc_orig.subgraph(neighborhood))
    info(g_sub)
    nontriv_sccs = split_to_nontrivial_sccs(g_sub)
    n_sccs = len(nontriv_sccs)
    if n_sccs == 0: # TODO Try to simplify it if it is single src-target DAG?
        raise NoProgress('Nothing to do, only nontrivial SCCs')
    elif n_sccs != 1:
        raise NoProgress('More than 1 nontrivial SCC, giving up')
    scc = nontriv_sccs[0]
    print('Attempting to simplify the following SCC')
    info(scc)
    return ilp_simplify(scc, scc_orig, rejected)

def ilp_simplify(scc, scc_orig, rejected):
    elims, obj_scc_alone, n_cyc_scc = _solve(scc)
    if len(scc)==len(scc_orig) and \
       scc.number_of_edges()==scc_orig.number_of_edges():
        print('Simplifier running on the entire SCC')
        update_cyc_max(n_cyc_scc)
        return elims, obj_scc_alone
    has_inedge, has_outedge = get_data_for_relaxation(scc, scc_orig)
    relaxed_scc = build_relaxation(scc, has_inedge, has_outedge)
    elims_relax, obj_relax, n_cyc_relax = _solve(relaxed_scc)
    print('Relaxation:', obj_relax)
    print('SCC alone: ', obj_scc_alone)
    if obj_scc_alone >= obj_relax:
        print('Safe to simplify SCC:')
        print(sorted(n for n in scc))
        update_cyc_max(n_cyc_relax)
        return elims_relax, obj_relax
    else:
        ##if obj_relax == obj_scc_alone+1: # border line rejection
        #print('Rejected by:', obj_relax-obj_scc_alone)
        #min_new_loops = max(len(has_inedge), len(has_outedge))
        #print('Min loops:', min_new_loops) 
        #scc_nodes = frozenset(scc)
        #rejected[scc_nodes] = (frozenset(has_inedge),frozenset(has_outedge))
        raise NoProgress('Not safe to simplify')

def get_data_for_relaxation(scc, scc_orig):
    # Returns the ( [nodes in scc having inedge], [nodes in scc with outedge] ).
    has_inedge  = [ t for t in scc for s in scc_orig.pred[t] if s not in scc ]
    has_outedge = [ s for s in scc for t in scc_orig.succ[s] if t not in scc ]
    if not has_inedge or not has_outedge: # sort of "isolated" SCC
        print('SCC can be eliminated on its own')
        raise AssertionError() # TODO Check how this can happen and remove if OK
    assert has_inedge and has_outedge
    return has_inedge, has_outedge

def build_relaxation(scc, has_inedge, has_outedge):
    io_set = set(has_inedge)
    io_set.intersection_update(has_outedge)
    if io_set:
        # n_1 -> (node in scc) -> n_2, where n_1 and n_2 are outside: No way to 
        # break this inside the scc. Requires that the node in scc has at least 
        # 2 in- and 2 out edges 
        raise NoProgress('Giving up (T -> n -> S corner case)')    
    relaxed_scc = deepcopy(scc)
    aux_node = 'auxiliary node'
    assert not relaxed_scc.has_node(aux_node)
    sum_of_w = sum(d['weight'] for _,_,d in relaxed_scc.edges_iter(data=True))
    M = sum_of_w + 1
    for u in has_outedge:
        relaxed_scc.add_edge(u, aux_node, weight=M)
    for v in has_inedge:
        relaxed_scc.add_edge(aux_node, v, weight=M)
    return relaxed_scc

def _solve(g):
    error_msg, elims, obj, ncyc = rigorous_mfes(g, MAX_CYCLES)
    if error_msg:
        raise NoProgress(error_msg)
    return elims, obj, ncyc

#-------------------------------------------------------------------------------

def iteratively_remove_runs_and_bypasses(g):
    while True:
        n_nodes = len(g)
        remove_runs_and_bypasses(g)        
        if n_nodes == len(g):
            return

def remove_runs_and_bypasses(g):
    # Code triplication! See cleanup_siso_nodes in mfes.
    siso_nodes = [ n for n in g if len(g.pred[n])==1 and len(g.succ[n])==1 ]
    for n in siso_nodes:
        (pred,), (succ,) = g.pred[n], g.succ[n]
        cost, (u,v) = min( (g[pred][n]['weight'], (pred,n)),
                           (g[n][succ]['weight'], (n,succ)) )
        edge_d = g[u][v]
        if pred==n and succ==n:
            pass # self-loop
        elif pred == succ:
            pass # 2-loop 
        elif g.has_edge(pred, succ): # A bypass, it would create multiple edges
            #print('Increasing weight of {} -> {} by {}'.format(pred,succ,cost))
            d = g.edge[pred][succ]
            d['orig_edges'].extend( edge_d['orig_edges'] )
            d['weight'] += cost
            g.remove_node(n)
        else: # either a 3-loop pred -> n -> succ -> pred, or n is a junk node
            g.add_edge(pred, succ, deepcopy(edge_d)) # deepcopy needed?
            #print('Adding edge: {} {}'.format((pred,succ), edge_d))
            g.remove_node(n)

#-------------------------------------------------------------------------------

def simple_path_neighborhood(g, source, cutoff):
    # Works a lot worse than BFS
    neighborhood = set()
    for target in g.pred[source]:
        # These will give the loops
        paths = list( nx.all_simple_paths(g, source, target, cutoff) )
        for p in paths:
            neighborhood.update( p )
        # These will give the bypasses
        paths = list( nx.all_simple_paths(g, target, source, cutoff) )
        for p in paths:
            neighborhood.update( p )
    return neighborhood

def get_bfs_neighborhood(g, source, cutoff):
    assert cutoff >= 0
    visited = set()
    level = 0
    visited.add(source)
    nextlevel =  set(g.pred[source])
    nextlevel.update(g.succ[source])
    nextlevel.discard(source)
    while nextlevel and level < cutoff:
        thislevel = nextlevel
        nextlevel = set()
        level += 1
        visited.update(thislevel)
        for n in thislevel:
            nextlevel.update(g.pred[n])
            nextlevel.update(g.succ[n])
        nextlevel.difference_update(visited)
    return visited, nextlevel

def bfs_w_local_search(g, source, cutoff):
    selected, nextlevel = get_bfs_neighborhood(g, source, cutoff)
    for n in nextlevel:
        #nbrs = g.pred[n].viewkeys() | g.succ[n].viewkeys()
        nbrs = set(g.pred[n]) | set(g.succ[n]) # Py3 compatibility
        n_in  = len(nbrs & selected)
        n_out = len(nbrs - selected)
        if n_in > n_out:
        #if n_out == 0:
            selected.add(n)
    return selected

def improve_w_local_search(scc, neighborhood):
    selected = set(neighborhood)
    nextlevel = set()
    for n in selected:
        #nextlevel.update(scc.pred[n].viewkeys() - selected)
        #nextlevel.update(scc.succ[n].viewkeys() - selected)
        nextlevel.update(set(scc.pred[n]) - selected) # Py3 compatibility
        nextlevel.update(set(scc.succ[n]) - selected)
    #
    for n in nextlevel:
        #nbrs = scc.pred[n].viewkeys() | scc.succ[n].viewkeys()
        nbrs = set(scc.pred[n]) | set(scc.succ[n])  # Py3 compatibility
        n_in  = len(nbrs & selected)
        n_out = len(nbrs - selected)
        if n_in > n_out:
        #if n_in > n_out + 1:
        #if n_out == 0:
            selected.add(n)
    return selected  

#-------------------------------------------------------------------------------

def dbg_dump_as_edgelist(g):
    for n in sorted(n for n in g):
        print(n, end='  ')
        for nbr in sorted(g[n]):
            print(nbr, end='  ')
        print()

if __name__ == '__main__':
    main()
