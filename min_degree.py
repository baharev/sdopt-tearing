# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from copy import deepcopy
import six
from networkx import DiGraph
from order_util import deterministic_topological_sort
from tearing import BLACK, RED, BLUE, initialize_edges, paint_edges_BLUE, \
                    itr_nodes_with_only_forbidden_black_edges, \
                    itr_one_black_edge_which_is_allowed, itr_black_edges, \
                    itr_nodes_with_only_blue_edges

def min_degree(g_orig, eqs, forbidden=set()):
    'Returns: (dag, [tear vars], [residual equations], [topsort of dag]).'
    if not eqs:
        return None, [ ], [ ], [ ]
    # g is a deep copy of g_orig, equations and variables are sorted lists
    g, equations, variables = initialization(g_orig, eqs, forbidden)
    # Run the heuristic
    remove_no_choice(g, equations, variables)
    _, eq_eligvars = eqs_with_minimum_black_edges(g, equations, variables)
    while eq_eligvars:
        # eligvars: variables eligible for elimination
        # there are ties; all matches are examined and the best one is picked
        g, eq_eligvars = try_all_matches(g, eq_eligvars, equations, variables)
    # The heuristic is done, now turn the colored graph g into an ordering.
    return get_final_order(g, forbidden, equations, variables)

def initialization(g_orig, eqs, forbidden):
    if not isinstance(eqs, (set, dict)):
        eqs = set(eqs) 
    g = deepcopy(g_orig)
    equations = sorted(eqs)
    variables = sorted(n for n in g if n not in eqs)
    initialize_edges(g, equations, forbidden)
    for _,_,d in g.edges_iter(data=True):
        d.pop('weight', None) # could interfere with the degree calculation
    return g, equations, variables

def remove_no_choice(g, equations, variables):
    while True:
        made_progress = False
        # variables with only forbidden black edges left -> tear them
        for var in itr_nodes_with_only_forbidden_black_edges(g, variables):
            paint_edges_BLUE(g, var)
            made_progress = True
        # equations with only forbidden black edges left -> make them sinks
        for eq in itr_nodes_with_only_forbidden_black_edges(g, equations):
            paint_edges_BLUE(g, eq)
            made_progress = True
        # equations with only one black allowed edge left -> eliminate (eq,var)
        for eq, var in itr_one_black_edge_which_is_allowed(g, equations):
            paint_edges_BLUE(g, var)
            g.edge[eq][var]['color'] = RED
            made_progress = True
        # variables with only one black allowed edge left -> eliminate (eq,var)
        for var, eq in itr_one_black_edge_which_is_allowed(g, variables):
            paint_edges_BLUE(g, eq)
            g.edge[eq][var]['color'] = RED
            made_progress = True
        if not made_progress:
            return

def eqs_with_minimum_black_edges(g, equations, variables):
    # Returns: (minimum black edge count, [(eq, [allowed black neighbors])]).
    # _evars stands for eligible variables, eligible for elimination
    eq_cnt_evars = [ ]
    for eq, nbr_allowed in itr_black_edges(g, equations):
        eligible_vars = sorted( nbr for nbr, allowed in nbr_allowed if allowed )
        assert eligible_vars, 'No choice equations should have been processed'
        eq_cnt_evars.append( (eq, len(nbr_allowed), eligible_vars) )
    # In the first pass, get the minimum black edge count
    if eq_cnt_evars:
        minelem = min(eq_cnt_evars, key=lambda t: t[1])
        mincnt = minelem[1]
    else:
        mincnt = 0
    # In the second pass, collect the group of min count equations (partition).
    # The returned list (group) will be empty if mincnt == 0.
    return mincnt, [ (eq,nbrs) for eq,cnt,nbrs in eq_cnt_evars if cnt==mincnt ]

def try_all_matches(g_orig, eq_nbrs, equations, variables):
    matches = [ (eq,var) for eq, nbrs in eq_nbrs for var in nbrs ]
    assert len(matches) >= 2, 'No choice cases should have been processed'
    candidates = [ ]
    for eq, var in matches:
        # Try the eq -> var match on the copy of the original graph
        g = deepcopy(g_orig)
        # Make eq -> var RED, all other eq <- var BLUE, so potentially tear vars
        match(g, eq, var)
        remove_no_choice(g, equations, variables)
        # What will be necessary to try in the *next* round?
        mincnt, eq_evars = eqs_with_minimum_black_edges(g, equations, variables)
        # cost: necessary guesses to continue + torn variables
        guesses = max(mincnt-1, 0)
        n_torn_vars = count_nodes_with_only_blue_edges(g, variables)
        cost = guesses + n_torn_vars
        # better progress = less equations remaining 
        eqs_remaining = count_nodes_with_black_edges(g, equations)
        # favor more future choices
        choices = sum(len(evars) for _, evars in eq_evars)
        score = (cost, eqs_remaining, guesses, -choices)
        candidates.append( (score, g, eq_evars) )
    (score, g, eq_evars) = min(candidates, key=lambda t: t[0]) 
    return g, eq_evars

def match(g, eq, var):
    paint_edges_BLUE(g, eq)
    paint_edges_BLUE(g, var)
    g.edge[eq][var]['color'] = RED # eq -> var
    # make all other variables connected to eq input variables (eq <- var) 
    # and potentially tear variables
    for nbr in g.edge[eq]:
        for d in six.itervalues(g.edge[nbr]):
            if d['color'] == BLACK:
                d['color'] = BLUE

def get_final_order(g, forbidden, equations, variables):
    'Returns: (dag, [tear vars], [residual equations], [topsort of dag]).'
    # DAG from g: red edges: eq -> var; blue edges: eq <- var; no black edges
    dag = to_dag(g, equations, forbidden)
    # all tears are sources, no other variable is source
    sources = sorted(n for n,indeg in dag.in_degree_iter(variables) if indeg==0)
    # all residuals are sinks, no other equation is sink
    sinks   = sorted(n for n, out in dag.out_degree_iter(equations) if out==0)
    # do a topological sort;  
    # nbunch tries to break ties by ordering the equations alphabetically.
    nbunch = list(reversed(equations))
    nbunch.extend(reversed(variables))
    order = deterministic_topological_sort(dag, nbunch)    
    return dag, sources, sinks, order

def to_dag(g, equations, forbidden):
    dag = DiGraph()
    dag.add_nodes_from(n for n in g)
    for eq, var, d in g.edges_iter(equations, data=True):
        color = d['color']
        assert color==RED or color==BLUE, 'Missed edge: {} - {}'.format(eq,var) 
        if color==RED:    # eq -> var
            dag.add_edge(eq, var)
            assert g[eq][var]['allowed'], (eq,var)
        else: # color==BLUE eq <- var
            dag.add_edge(var, eq)
    return dag

#-------------------------------------------------------------------------------
# Auxiliary stuff, only interesting on the implementation level

def count_nodes_with_only_blue_edges(g, nodes):
    return sum(1 for _ in itr_nodes_with_only_blue_edges(g, nodes))

def nodes_with_black_edges(g, nodes):
    for n in nodes:
        for d in six.itervalues(g.edge[n]):
            if d['color']==BLACK:
                yield n
                break

def count_nodes_with_black_edges(g, nodes):
    return sum(1 for _ in nodes_with_black_edges(g, nodes))

#-------------------------------------------------------------------------------

def run_tests():
    from test_tearing import gen_testproblems
    for g, eqs, forbidden in gen_testproblems():
        _, sources, sinks, _ = min_degree(g, eqs, forbidden)
        print('Tears:', sources)
        print('Residuals:', sinks)

if __name__=='__main__':
    run_tests()
