# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from copy import deepcopy
import networkx as nx
from utils import double_check
from utils import info_short as info
# The info from utils will try to enumerate all simple cycles but this may not 
# be tracktable. Use info_short to avoid this problem.

__all__ = [ 'run_mfes_heuristic' ]

dbg = print

#log = print
def log(*args, **kwargs): pass

def run_mfes_heuristic(g_input, try_one_cut=False, is_labeled=False):
    # Set the edge attributes 'weight' and 'orig_edges' if necessary
    g, copy_g = label_edges(g_input, is_labeled)
    #
    greedy_choice = no_lookahead if try_one_cut else with_lookahead
    #
    running_cost, elims = 0, [ ]
    sccs, running_cost = iterate_cleanup(g, running_cost, elims, copy_g)
    info_after_cleanup(sccs, running_cost)
    while True:
        sccs_to_process = [ ]
        for sc in sccs:
            log('-----------------------------------------------------------')
            #distributions(sc)
            new_state = greedy_choice(sc, running_cost, elims)
            new_sccs, elims, running_cost = new_state
            sccs_to_process.extend(new_sccs)
        sccs = sccs_to_process
        if len(sccs)==0:
            # Make sure that what we return is at least consistent
            double_check(g_input, running_cost, elims, is_labeled, log=log)
            return running_cost, elims

#-------------------------------------------------------------------------------

def label_edges(g_input, is_labeled):
    if not is_labeled:
        dig = nx.DiGraph()
        for u, v in g_input.edges_iter():
            dig.add_edge(u, v, { 'weight' : 1, 'orig_edges' : [ (u,v) ] })
        # No need to copy g later: the caller of mfes doesn't know about dig
        g, copy_g = dig, False 
    else:
        # The caller has already labeled the graph but we still need to copy g
        g, copy_g = g_input, True 
    return g, copy_g    

def info_after_cleanup(sccs, running_cost):
    if sccs:
        log('------------')
        log('After invoking the simplifier', len(sccs),'SCCs remain')
        log('Info:\n')
        for sc in sccs:
            info(sc, log=log)
        log('Cost:', running_cost)
    else:
        log('\nSimplifier eliminated the input, cost:', running_cost)

#-------------------------------------------------------------------------------

def no_lookahead(sc, running_cost, elims):
    # See with_lookahead for documentation.
    inc, outc = cheapest_nodes(sc)
    edges_to_break = select_one_edgelist(sc, inc, outc)
    for u, v, d in edges_to_break:
        elims.extend(d['orig_edges'])
        running_cost += d['weight']
        sc.remove_edge(u, v)
    sccs, running_cost = iterate_cleanup(sc, running_cost, elims, copy=False)
    return sccs, elims, running_cost 

def select_one_edgelist(sc, inc, outc):
    # We look for nodes that are highly asymmetric: cheapest to cut on one side 
    # (e.g. all in-edges) but very expensive on the other side (e.g. 
    # all out-edges). With some luck, we can find neighbors, connected by a 
    # single edge, where the tail node has high in-degree, and the head node has
    # high out-degree. The hope is that it breaks many simple cycles.
    if inc:
        min_indeg    = min_in_degree(sc, inc)
        nodes_minin  = nodes_with_indeg(sc, inc, min_indeg)
    if not outc:
        # Sub-optimal: We could pick the max in- + out-cardinality, or max in- 
        # and max out weight. 
        return sc.in_edges(nodes_minin[0], data=True)
    
    if outc:
        min_outdeg   = min_out_degree(sc, outc)
        nodes_minout = nodes_with_outdeg(sc, outc, min_outdeg)
    if not inc:
        return sc.out_edges(nodes_minout[0], data=True) # Sub-optimal, see above
    
    assert min_indeg == min_outdeg
    
    # There might be neighbor nodes, connected with a single edge
    if min_indeg==1: # and min_outdeg == 1 as asserted above
        nbrs = { }
        for v in nodes_minin:
            (u,) = sc.pred[v]
            nbrs[u] = v
        candidates = [ ]
        for u in nodes_minout:
            if u in nbrs:
                u, v = u, nbrs[u] 
                assert sc.has_edge(u,v), (u,v)
                assert u in outc
                assert v in inc
                score = len(sc.pred[u]) + len(sc.succ[v])
                candidates.append( (score,(u,v)) )
        if candidates:
            score, (u,v) = max(candidates, key=lambda t: t[0])
            log('Tearing {} -> {} with score {}'.format(u,v,score))
            return [ (u,v,sc[u][v]) ]
    
    # No luck: Just pick the max fan out (most asymmetric node).
    max_outw   = max_out_weight(sc, inc)
    max_inw    = max_in_weight(sc, outc)
    if max_outw > max_inw:
        n = a_node_with_out_weight(sc, inc, max_outw)
        return sc.in_edges(n, data=True)
    else:
        n = a_node_with_in_weight(sc, outc, max_inw)
        return sc.out_edges(n, data=True)
    
def min_in_degree(g, nbunch):
    return min(indeg  for _, indeg  in g.in_degree_iter(nbunch))

def max_in_degree(g, nbunch):
    return max(indeg  for _, indeg  in g.in_degree_iter(nbunch))

def min_out_degree(g, nbunch):
    return min(outdeg for _, outdeg in g.out_degree_iter(nbunch))

def max_out_degree(g, nbunch):
    return max(outdeg for _, outdeg in g.out_degree_iter(nbunch))

def nodes_with_indeg(g, nbunch, deg):
    return [ n for n, indeg  in g.in_degree_iter(nbunch)  if indeg==deg  ]

def nodes_with_outdeg(g, nbunch, deg):
    return [ n for n, outdeg in g.out_degree_iter(nbunch) if outdeg==deg ]

def max_in_weight(g, nbunch):
    return max(in_w  for _,in_w  in g.in_degree_iter(nbunch, weight='weight'))

def max_out_weight(g, nbunch):
    return max(out_w for _,out_w in g.out_degree_iter(nbunch, weight='weight'))

def a_node_with_in_weight(g, nbunch, w):
    gen = (n for n,in_w in g.in_degree_iter(nbunch,weight='weight') if in_w==w)
    n = next(gen, None)
    assert n is not None
    return n

def a_node_with_out_weight(g, nbunch, w):
    gen = (n for n,outw in g.out_degree_iter(nbunch,weight='weight') if outw==w)
    n = next(gen, None)
    assert n is not None
    return n

#-------------------------------------------------------------------------------
# Implementation boundary for the greedy choice: only with_lookahead is 
# referenced outside of this code block.

def with_lookahead(sc, running_cost, elims):
    'Returns: (sccs, elims, running_cost) for the best greedy choice.'
    # We look for nodes that are cheap to turn into a source or a sink: Either 
    # all in- or all out-edges of a node are cut (these are the in- and out-edge
    # lists, per node). There is a brute force lookahead step: We try all 
    # minimum cost cuts (breaking all edges in the in- or out-edge list of the 
    # node), and pick the one causing the biggest damage to the SCC. 
    inc, outc = cheapest_nodes(sc)
    log('Cheapest in: ', inc)
    log('Cheapest out:', outc)
    inedge_lists  = [ sc.in_edges( n, data=True) for n in inc  ]
    outedge_lists = [ sc.out_edges(n, data=True) for n in outc ]        
    seen = set() 
    incuts  = try_edge_lists(sc,  inedge_lists, running_cost, elims, seen)
    outcuts = try_edge_lists(sc, outedge_lists, running_cost, elims, seen)
    # The returned incuts and outcuts have been sorted by try_edge_lists
    #show_best_ones(incuts,  'Best in cuts: ')
    #show_best_ones(outcuts, 'Best out cuts:')
    # Some statistics
    nedges = sc.number_of_edges()
    percentage = len(seen)/float(nedges) * 100.0
    msg = '{:0.1f} % of the edges have been tried ({} edges out of {})'
    log(msg.format(percentage, len(seen), nedges))
    return pick_best(incuts, outcuts) 
 
def cheapest_nodes(sc):
    min_in  = min_in_weight(sc)
    min_out = min_out_weight(sc)
    assert min_in!=0 and min_out!=0, (min_in, min_out) # Not a strong component!
    cheapest_in  = nodes_with_in_w( sc, min_in ) if min_in  <= min_out else [ ]
    cheapest_out = nodes_with_out_w(sc, min_out) if min_out <= min_in  else [ ]
    return cheapest_in, cheapest_out

def min_in_weight(g):
    return min(in_w  for _,in_w  in g.in_degree_iter( weight='weight'))

def min_out_weight(g):
    return min(out_w for _,out_w in g.out_degree_iter(weight='weight'))

def nodes_with_in_w(g, w):
    return [ n for n,in_w in  g.in_degree_iter(weight='weight')  if in_w==w  ]

def nodes_with_out_w(g, w):
    return [ n for n,out_w in g.out_degree_iter(weight='weight') if out_w==w ]

def try_edge_lists(g, edge_lists, running_cost, elims, seen):
    # Since the inedge_lists and outedge_lists overlap, it is beneficial to save
    # the already tried ones in seen. The implementation is somewhat inefficient
    # as we always push the edges in inedge_lists to seen: The seen only helps
    # when we are processing the outedge_lists. 
    nodecuts = [ ]
    for edges in edge_lists:
        edges_sorted = tuple(sorted((u,v) for u,v,_ in edges))
        if edges_sorted not in seen:
            seen.add(edges_sorted)
            nodecuts.append(try_edges(g, edges, running_cost, elims))
    nodecuts.sort(key=worseness)
    return nodecuts

def try_edges(g_orig, edges_to_break, running_cost, elims_orig):
    #dbg('---')
    g, elims = deepcopy(g_orig), list(elims_orig)
    for u, v, d in edges_to_break:
        elims.extend(d['orig_edges'])
        running_cost += d['weight']
        g.remove_edge(u, v)
        #dbg(u, '->', v)
    sccs, running_cost = iterate_cleanup(g, running_cost, elims, copy=False)
    max_size = size_of(max(sccs, key=size_of)) if sccs else (0,0)
    # TODO Move this data clump into a struct like class?
    return sccs, elims, running_cost, max_size        

def size_of(g):
    return (g.number_of_nodes(), g.number_of_edges())

def worseness(nodecut):
    _, _, running_cost, max_size = nodecut
    return (max_size, running_cost)
    
def show_best_ones(nodecuts, msg, cutoff=5):
    log('\n=====\n', msg, sep='')
    for nodecut in nodecuts[:cutoff]:
        show_nodecut(nodecut)

def show_nodecut(nodecut):
    sccs, _, running_cost, max_size = nodecut
    for sc in sccs: 
        info(sc)
    log('Max. size:', max_size)
    log('Cost:', running_cost)
    log('-----')

def pick_best(incuts, outcuts):
    # Assumes that incuts and outcuts have been sorted
    if len(incuts)==0:
        new_state = outcuts[0] 
    elif len(outcuts)==0:
        new_state =  incuts[0]
    else:
        inbest, outbest = incuts[0], outcuts[0]
        new_state = inbest if worseness(inbest)<worseness(outbest) else outbest
    log('\n***  Picked best node cut  ***')
    show_nodecut(new_state)
    sccs, elims, running_cost, _ = new_state # cut off max_size
    return sccs, elims, running_cost

#-------------------------------------------------------------------------------
# Implementation boundary, only iterate_cleanup is referenced outside this block

def iterate_cleanup(g_orig, running_cost, elims, copy=True):
    # See noncopy_split_to_nontrivial_sccs as to why the g_orig is copied by 
    # default. The elims argument is mutated: eliminations are appended to it.
    dirty = [ deepcopy(g_orig) ] if copy else [ g_orig ]
    clean_sccs = [ ]
    while dirty:
        part = dirty.pop()
        nnodes, nedges = part.number_of_nodes(), part.number_of_edges() 
        new_parts, running_cost = clean_up(part, running_cost, elims)
        if len(new_parts)==0:
            continue
        p = new_parts[0]
        # we are done with part if there was no progress
        if len(new_parts)==1 and nnodes==len(p) and nedges==p.number_of_edges():
            clean_sccs.append(p) # p == part
        else:
            dirty.extend(new_parts)
    return clean_sccs, running_cost

def clean_up(g, running_cost, elims):
    running_cost = remove_self_loops(g, running_cost, elims)
    clean_sccs = [ ]
    for scc in noncopy_split_to_nontrivial_sccs(g):
        running_cost = cleanup_siso_nodes(scc, running_cost, elims)
        if len(scc) > 0:
            clean_sccs.extend( noncopy_split_to_nontrivial_sccs(scc) )
    return clean_sccs, running_cost

def remove_self_loops(g, running_cost, elims):
    for n, n, d in g.selfloop_edges(data=True):
        log('Self-loop at node {}, cost: {}'.format(n, d['weight']))
        elims.extend(d['orig_edges'])
        running_cost += d['weight']
        g.remove_edge(n,n)
    return running_cost

# TODO noncopy_split_to_nontrivial_sccs is still the slowest part of the 
# heuristic, according to the profiling. However, in clean_up it is often (but
# not always) unnecessary to call it because iterate_cleanup has already split 
# the graph to clean SCCs in the previous iteration.

def noncopy_split_to_nontrivial_sccs(g):
    # The edge dictionaries are NOT copied. It is the caller's responsibility 
    # that these attributes of g are not used after the call to this function.
    sccs = list(nx.strongly_connected_components(g))
    return [ g.subgraph(sc) for sc in sccs if len(sc) > 1 ]

# TODO Two clean-up steps should be written: one that can never invalidate SISO
# nodes and another that breaks loops (so it can invalidate SISO nodes). A more 
# sophisticated loop breaking strategy should be written and tested on a hand
# crafted graph; this would only be useful in deriving lower bounds on the 
# minimum feedback edge set.  
#
# Code triplication! The grb_relaxation and grb_simplifier modules use a very 
# similar cleanup: They require the SCC not to be split.

def cleanup_siso_nodes(g, running_cost, elims):
    # g: a single nontrivial SCC with no self-loops
    siso_nodes=[n for n in sorted(g) if len(g.pred[n])==1 and len(g.succ[n])==1]
    for n in siso_nodes:
        n_pred, n_succ = len(g.pred[n]), len(g.succ[n])
        if (n_pred==0 and n_succ==1) or (n_pred==1 and n_succ==0):
            continue # Breaking isolated 3-loops can invalidate SISOs, see above
        (pred,), (succ,) = g.pred[n], g.succ[n]
        cost, edge = min( (g[pred][n]['weight'],(pred,n)),
                          (g[n][succ]['weight'],(n,succ)) )
        edge_dict = g.get_edge_data(*edge)
        orig_edges = edge_dict['orig_edges']
        assert pred!=n and succ!=n, (pred, n, succ) # self-loop
        if pred == succ: # Would create a self-loop, breaking it instead 
            log('Breaking 2-loop at {} - {}, cost: {}'.format(pred, n, cost))
            elims.extend(orig_edges)
            running_cost += cost
            remove_node(g, n)
            if len(g)==1 and g.number_of_edges()==0:
                g.clear()
                break
        elif g.has_edge(pred, succ): # A bypass, it would create multiple edges
            log('Increasing weight of {} -> {} by {}'.format(pred,succ,cost))
            d = g.edge[pred][succ]
            d['orig_edges'].extend(orig_edges)
            d['weight'] += cost
            remove_node(g, n)
        elif g.has_edge(succ, pred): # 3-loop: pred -> n -> succ -> pred
            args=(g, pred, n, succ, cost, edge, orig_edges, running_cost, elims)
            running_cost = handle_3loop(*args)
            if len(g)==0:
                break
        else: # junk node n
            #dbg('New edge: {} -> {}; d: {}'.format(pred,succ,edge_dict))
            g.add_edge(pred, succ, deepcopy(edge_dict)) # Is deepcopy needed?
            remove_node(g, n)
        # TODO Remove when done:
        #assert_sane(g, running_cost, elims, n)
    #assert_sane(g, running_cost, elims, n) # there are breaks in the above loop
    return running_cost

def assert_sane(g, running_cost, elims, n=''):
    # Only if all the initial edge weights were 1!
    for u,v,d in g.edges_iter(data=True):
        assert d['weight']==len(d['orig_edges']), (n,u,v,d) # iff initial w=1 !!
    assert running_cost == len(elims), (n, running_cost, len(elims))    

def handle_3loop(g, pred, n, succ, cost, edge, orig_edges, running_cost, elims):
    # 3-loop: pred -> n -> succ -> pred
    d = g[succ][pred]
    w = d['weight']
    if len(g)==3: # The whole SCC is just this 3-loop
        log('Removing last 3-loop: {} - {} - {}'.format(pred, n, succ))
        assert g.number_of_edges()==3 and g.number_of_selfloops()==0 
        if cost < w:
            elims.extend(orig_edges)
            running_cost += cost
        else:
            elims.extend(d['orig_edges'])
            running_cost += w
        g.clear() # We are done!
    # Print cost, and edge chosen?
    elif w <= cost and (in_card(g, pred)==1 or out_card(g, succ)==1):
        log('Breaking 3-loop  (I): {} - {} - {}'.format(pred, n, succ))
        elims.extend(d['orig_edges'])
        running_cost += w
        g.remove_edge(succ, pred)
    elif cost <= w and (in_card(g, succ)==1 or out_card(g, pred)==1):
        log('Breaking 3-loop (II): {} - {} - {}'.format(pred, n, succ))
        elims.extend(orig_edges)
        running_cost += cost
        g.remove_edge(*edge)
    else: # Leave it for the heuristic.
        log('Unchanged 3-loop: {} - {} - {}'.format(pred, n, succ))
    return running_cost

def remove_node(g, n):
    #log('Removing node:', n)
    g.remove_node(n)    

def in_card(g, n):
    return len(g.pred[n])

def out_card(g, n):
    return len(g.succ[n])

#-------------------------------------------------------------------------------

def main():
    from benchmarks import gen_benchmark_digraphs
    for g in gen_benchmark_digraphs():
        info(g)
        run_mfes_heuristic(g, try_one_cut=False, is_labeled=True)

if __name__ == '__main__':
    main()
