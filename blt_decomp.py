# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
import networkx as nx
from networkx.algorithms.bipartite import is_bipartite_node_set, projected_graph
from plot import dummy as plot

__all__ = [ 'Dulmage_Mendelsohn' ]

#log = print
def log(*args, **kwargs): pass

# When the real plot is used from plot, the order of the plots is as follows:
# (1) The test problem, (2) the bipartite graph after matching, 
# (3) the condensed graph (the SCCs of the equations), and (4) the SCC subgraph
# with all nodes that are connected to this SCC (under/over-determined blocks).

def Dulmage_Mendelsohn(g, eqs):
    '''The input graph g is assumed to be a bipartite graph with no isolated 
    nodes. Returns the diagonal blocks as a list of (equations, variables).'''
    assert_all_in_graph(g, eqs)
    assert_no_isolates(g)
    assert eqs, 'At least one equation is expected'
    assert is_bipartite_node_set(g, eqs)
    # Maximum matching
    mate = nx.max_weight_matching(g, maxcardinality=True)
    matches = sorted((k,mate[k]) for k in mate if k in eqs)
    log('Matches:')
    for eq, var in matches:
        log(eq, var)
    # Direct the edges of g according to the matching
    bipart = to_digraph(g, matches)
    plot(bipart)
    # Find the strongly connected components (SCCs) of the equations
    eq_sccs = nx.condensation( projected_graph(bipart, eqs) )
    plot(eq_sccs)
    # Q: With proper implementation, shouldn't the SCCs be already top. sorted?
    precedence_order = nx.topological_sort(eq_sccs)
    # Collect the diagonal blocks as a list of (equations, variables)
    diagonal_blocks = [ ]
    seen = set()
    for scc in precedence_order:
        equations = eq_sccs.node[scc]['members']
        variables = {n for eq in equations for n in g.edge[eq] if n not in seen}
        seen.update(variables)        
        diagonal_blocks.append( (equations,list(variables)) )
    return diagonal_blocks

def to_digraph(g, eq_var_matches):
    bipart = nx.DiGraph()
    for eq, var in eq_var_matches:
        # eq -> var
        bipart.add_edge(eq, var)
        # eq <- var_k (dependencies of eq other than var)
        deps = [n for n in g.edge[eq] if n!=var]
        for var_k in deps:
            bipart.add_edge(var_k, eq)
        # eq_k <- var (other equations involving var)
        substitute = [(var, n) for n in g.edge[var] if n!=eq]
        bipart.add_edges_from(substitute)
    # If there are no isolates in g, the above loop should have added all 
    # unmatched equations or variables as well
    missed_nodes = [n for n in g if n not in bipart]
    assert not missed_nodes, 'Failed to copy nodes: {}'.format(missed_nodes)
    return bipart

def assert_all_in_graph(g, eqs):
    missing = [ eq for eq in eqs if eq not in g ]
    assert not missing, 'equations not in the input graph: {}'.format(missing)

def assert_no_isolates(g):
    isolated_nodes = nx.isolates(g)
    assert not isolated_nodes, isolated_nodes

#-------------------------------------------------------------------------------

def run_tests():
    from test_tearing import gen_testproblems
    for g, eqs, _ in gen_testproblems():
        plot(g)
        diagonal_blocks = Dulmage_Mendelsohn(g, eqs)
        print('Equation SCCs in topologically sorted order:')
        for equations, variables in diagonal_blocks:
            equations, variables = sorted(equations), sorted(variables)
            print('-----------------------------------')
            print('Equations:    ', equations)        
            print('New variables:', variables)
            plot( g.subgraph(equations+variables) )

if __name__ == '__main__':
    run_tests()
