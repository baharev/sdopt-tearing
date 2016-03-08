# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from copy import deepcopy
import six
import networkx as nx
from blt_decomp import Dulmage_Mendelsohn
from plot import dummy as plot

__all__ = [ 'tearing' ]

def dbg_log(_): pass # logging decisions made by the heuristic
# Dulmage-Mendelsohn has its own logging for the matching results
#log = print
def log(*args, **kwargs): pass

# TODO Prioritize equations and variables
# Order the equations as they were in the input; the variables should be
# prioritized too (the user made the var explicit, safe transformations, unsafe
# transformations, forbidden). 
# Ordering according to the position in some other list:
#   position = { name : i for i, name enumerate(equations) }
# and the key function to sort is lambda name: position[name].

def tearing(g, eqs, forbidden=set()):
    # Returns: list of (tears list, (eq,var) elimination list, residuals list)
    torn_blocks = [ ] 
    diagonal_blocks = Dulmage_Mendelsohn(g, eqs)        
    log('Tearing the diagonal blocks of the BLT form')
    for eqs_in_block, vars_in_block in diagonal_blocks:
        equations, variables = sorted(eqs_in_block), sorted(vars_in_block)
        log('-----------------------------------')
        log('Equations:', equations)        
        log('Variables:', variables)
        subgraph = g.subgraph(equations+variables)
        tears_elim_resids = tear_scc(subgraph, equations, variables, forbidden)
        torn_blocks.append( tears_elim_resids )
    return torn_blocks

def tear_scc(graph, equations, variables, forbidden):
    state = deepcopy_input(graph, equations, variables, forbidden)
    # state: g, equations, variables, elim, bwd_elim
    elim_no_choice(*state)
    while True:
        new_states, scores = tear_iteratively(state)
        if scores: # There was at least one black edge, find best var and repeat
            index, score_var = max(enumerate(scores), key=lambda t: t[1][0])
            log('Best score: {} (var: {})'.format(*score_var))
            state = new_states[index]
        else: # No black edge left, done
            break
    # In the future: pass on the relevant part of the input and final states
    dbg_diagnostics(state)
    return to_elimination(state)

def elim_no_choice(g, equations, variables, elim, bwd_elim):
    while True:
        made_progress = False
        # variables with only forbidden black edges left -> tear them
        for var in itr_nodes_with_only_forbidden_black_edges(g, variables):
            dbg_log('{}  blue (tear, only forbidden edges left)'.format(var))
            paint_edges_BLUE(g, var)
            made_progress = True
        # Equations with only forbidden black edges left? Probably taken care of
        # from the variables' side.
        #
        # equations with only one black allowed edge left -> eliminate (eq,var)
        for eq, var in itr_one_black_edge_which_is_allowed(g, equations):
            dbg_log('{}  blue'.format(var))
            dbg_log('{} -> {} red (eq with 1 black edge)'.format(eq,var))
            elim.append((eq,var))
            paint_edges_BLUE(g, var)
            g.edge[eq][var]['color'] = RED
            made_progress = True
        # variables with only one black allowed edge left -> eliminate (eq,var)
        for var, eq in itr_one_black_edge_which_is_allowed(g, variables):
            dbg_log('{}  blue'.format(eq))
            dbg_log('{} -> {} red (var with 1 black edge)'.format(eq,var))    
            bwd_elim.append((eq,var))
            paint_edges_BLUE(g, eq)
            g.edge[eq][var]['color'] = RED
            made_progress = True
        if not made_progress:
            return

def tear_iteratively(state):
    new_states, scores = [ ], [ ]
    g_orig, variables_orig, = state[0], state[2]
    # var has at least one black edge that elim_no_choice could not 
    # eliminate and tearing is necessary; trying each var on a copy of state    
    for var, _ in itr_black_edges(g_orig, variables_orig):
        state_copy = deepcopy(state)
        g = state_copy[0]
        # Score calculation: initialization
        n_red_before = n_red_edges(g)
        not_allowed_black_edges = n_not_allowed_black_edges(g, var)
        # Tear var and see what can be done now: elim_no_choice
        paint_edges_BLUE(g, var)
        dbg_log('{}  probing tear, blue'.format(var))         
        elim_no_choice(*state_copy)
        # Score calculation logic
        made_causal = n_red_edges(g) - n_red_before
        score = made_causal + not_allowed_black_edges # MC3 score
        # Save the result of tearing in var
        scores.append( (score, var) )
        new_states.append( state_copy )
    return new_states, scores

def to_elimination(final_state):
    'Returns: ([tear variables], [elimination (eq,var)], [residual equations])'
    # There is some code duplication with dbg_diagnostics
    g, equations, variables, elim, bwd_elim = final_state
    tears     = [ var for var in itr_nodes_with_only_blue_edges(g, variables) ]
    residuals = [ eq  for eq  in itr_nodes_with_only_blue_edges(g, equations) ]
    elimination = list(elim)
    elimination.extend( reversed(bwd_elim) )
    return tears, elimination, residuals

def dbg_diagnostics(final_state):
    g, equations, variables, elim, bwd_elim = final_state
    
    tears     = [ var for var in itr_nodes_with_only_blue_edges(g, variables) ]
    residuals = [ eq  for eq  in itr_nodes_with_only_blue_edges(g, equations) ]
    
    # Diagnostic information
    log('\nElimination')
    if tears:
        log('Tears:', tears)
    for eq, var in elim:
        log(eq, '->', var)
    if bwd_elim:
        log('---')
    for eq, var in reversed(bwd_elim):
        log(eq, '->', var)
    if residuals:
        log('Residuals:', residuals)
    
    # Direct edges of g in such a way that we get a DAG
    dag = nx.DiGraph()
    for eq, var, d in g.edges_iter(equations, data=True):
        color = d['color']
        assert color==RED or color==BLUE, 'Missed edge: {} - {}'.format(eq,var) 
        if color==RED:    # eq -> var
            dag.add_edge(eq, var)
            assert g[eq][var]['allowed'], (eq,var)
        else: # color==BLUE eq <- var
            dag.add_edge(var, eq)
    
    # Corner case: a single sink equation
    if len(dag)==0 and len(g)==1:
        eq0 = equations[0]
        assert eq0 in g, 'Expected a sink equation, got: {}'.format(eq0)
        dag.add_node(eq0)
    
    # elimination should contain all nodes in topologically sorted order
    elimination = list(tears)
    for eq_var in elim:
        elimination.extend(eq_var)
    for eq_var in reversed(bwd_elim):
        elimination.extend(eq_var)
    elimination.extend(residuals)
    
    # checking the triviality first: elimination and dag must contain all nodes 
    all_nodes = sorted(g, key=str) # FIXME Comparing str and int breaks on Py 3
    assert sorted(elimination, key=str) == all_nodes
    assert sorted(dag, key=str) == all_nodes#,(sorted(n for n in dag),all_nodes)
    
    # checking whether elimination is a valid topological sort
    topsort_pos = { n : i for i, n in enumerate(elimination) }
    for u, v in dag.edges_iter():
        assert  topsort_pos[u] < topsort_pos[v]

    # all tears are sources, no other variable is source
    sources = sorted(n for n,indeg in dag.in_degree_iter(variables) if indeg==0)
    assert sorted(tears)==sources, sources
    
    # all residuals are sinks, no other equation is sink
    sinks = sorted(n for n,out in dag.out_degree_iter(equations) if out==0)
    assert sorted(residuals)==sinks, sinks
    
    log('\nTest passed, valid ordering!')

#-------------------------------------------------------------------------------
# Auxiliary stuff, only interesting on the implementation level

def deepcopy_input(graph, equations, variables, forbidden):
    g = deepcopy(graph)
    initialize_edges(g, equations, forbidden)
    return (g,list(equations),list(variables),[],[]) # last two: elim, bwd_elim

BLACK = 0
RED   = 1
BLUE  = 2

def initialize_edges(g, equations, forbidden):
    for eq in equations:
        for var, d in six.iteritems(g.edge[eq]):
            d['color']   = BLACK
            d['allowed'] = (eq,var) not in forbidden

def itr_one_black_edge_which_is_allowed(g, nodes):
    for n in nodes:
        black_nbr = ((nbr,d['allowed']) for nbr,d in six.iteritems(g.edge[n])
                                         if d['color']==BLACK)
        first_black, allowed = next(black_nbr,(None,None))
        if allowed and not next(black_nbr, None):
            yield n, first_black

def itr_black_edges(g, nodes):
    # yields iff n has at least one black edge 
    for n in nodes:
        nbr_allow = [(nbr,d['allowed']) for nbr,d in six.iteritems(g.edge[n])
                                         if d['color']==BLACK]
        if nbr_allow:
            yield n, nbr_allow

def itr_nodes_with_only_forbidden_black_edges(g, nodes):
    for n, nbr_allow in itr_black_edges(g, nodes):
        # n has at least one black edge 
        first_allowed = next((nbr for nbr, allow in nbr_allow if allow), None)
        if first_allowed is None: # none of the black edges is allowed
            yield n

def itr_nodes_with_only_blue_edges(g, nodes):
    # It also considers isolated nodes as nodes with blue edges. That's fine, 
    # as not doing so would cause problems with isolated nodes (sink equations).
    # The name of this generator should be nodes with no red or black edge?
    for n in nodes:
        notblue = (v for v,d in six.iteritems(g.edge[n]) if d['color']!=BLUE)
        if next(notblue, None) is None:
            yield n

def n_red_edges(g):
    return sum(1 for _,_,d in g.edges_iter(data=True) if d['color']==RED)

def n_not_allowed_black_edges(g, n):
    return sum( 1 for d in six.itervalues(g.edge[n]) 
                   if d['color']==BLACK and not d['allowed'] )

def paint_edges_BLUE(g, node):
    for nbr, d in six.iteritems(g.edge[node]):
        assert d['color'] != RED or plot(g), (node, nbr)
        d['color'] = BLUE

#-------------------------------------------------------------------------------

def run_tests():
    from test_tearing import gen_testproblems
    for g, eqs, forbidden in gen_testproblems():
        print()
        for block in tearing(g, eqs, forbidden):
            tears, _, residuals = block
            # Print only nontrivial strong components
            if tears:
                print('Tears:', tears)
            if residuals:
                print('Residuals:', residuals)
            if tears or residuals:
                print()
    modelica_tearing()

#-------------------------------------------------------------------------------
# Stuff finally removed from demo.py. Subject to deletion.

# Classic tearing, as in Modelica tools. Perform bipartite matching, then 
# find the strongly connected components (SCCs); in short, do a block lower 
# triangular (BLT) decomposition. Then, apply a variant of Cellier's greedy
# tearing heuristic to break all algebraic loops in each SCC. 
#result = model.run_classic_tearing()
#show(result)

def modelica_tearing():
    from model_fmux import BlockTearingResult #, ModelWithInletsAndOutlets
    from total_ordering import to_spiked_form
    #model = ModelWithInletsAndOutlets('demo.xml.gz')
    from utils import DATADIR, deserialize #, serialize
    #serialize(model, 'demo.pkl.gz')
    from os.path import isfile
    fname = DATADIR + 'demo.pkl.gz'
    if not isfile(fname):
        return
    model = deserialize(fname)
    equations, g, eqs, forbidden = _create_bipartite_repr(model)
    torn_blocks = tearing(g, eqs, forbidden)
    alltears, allresids, blocks = to_spiked_form(equations, torn_blocks)
    bounds = model.bounds.copy()
    res = BlockTearingResult(alltears, allresids, blocks, bounds)
    res.plot()
    res.generate_ampl_code()
    res.generate_python_code()

def _create_bipartite_repr(m):
    from equations import gen_nonnesting_eqs, to_bipartite_graph
    # The ordering algorithms may write the equations, so do a deepcopy
    equations = deepcopy(m.equations)
    equations = list(gen_nonnesting_eqs(equations))
    g, eqs, forbidden =  to_bipartite_graph(equations)
    return equations, g, eqs, forbidden

#-------------------------------------------------------------------------------

if __name__=='__main__':
    run_tests()
