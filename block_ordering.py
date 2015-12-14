# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
#from networkx import simple_cycles
from plot import plot
#from tearing_opl import dump_opl_input

__all__ = [ 'run_elimination' ]

# The fact that the input graph is a digraph is not used anywhere; it only makes
# the implementation slightly more complicated. We may want to try tie breaking
# heuristics in the future and then we may want to use the edge directions.

def info(g):
    #print('Elementary circuits of the digraph:\n%s' % list(simple_cycles(g)) )
    print('Node weights as (node id, weight):\n%s' %
          [(n,d['weight']) for n, d in g.nodes(data=True)])    
    print('Edge weights as (src node, dst node, weight):\n%s' %
          [(u,v,d['weight']) for u,v,d in g.edges(data=True)])
    print('-------------------------------------------------------------------')

def run_elimination(g, hints=set()):
    '''Returns blocks_in_order which is a list of lists: A list of unit 
    equations is followed by the connections eliminated with the unit (a 
    sequence of lists). Use gen_unit_conn to iterate over blocks_in_order.'''
    info(g)
    assert g.number_of_selfloops()==0
    #dump_opl_input(g)
    #hints = { 'cascade.stages[7].mixer' } # TODO Hack here to overwrite hints
    blocks_in_order = [ ]
    candidates = set()
    total_cost = 0
    while g:
        if candidates:
            # TODO Try to preserve the order so that we don't jump around too 
            # much; hashing unfortunately messes up the order.
            n = select_next_node(candidates, hints)
        else:
            minweight_nodes = find_minweight_nodes(g)
            n = select_next_node(minweight_nodes, hints)
        eliminated_blocks, new_candidates, cost = eliminate(g, n)
        blocks_in_order += eliminated_blocks 
        candidates.update( new_candidates )
        total_cost += cost
    print('===================================================================')
    print('Total cost:', total_cost)
    print('Done!\n\n')
    return blocks_in_order

def find_minweight_nodes(g):
    INT_MAX = int(2**31-1)
    minweight, nodes = INT_MAX, set()
    for n, w in gen_nodes_with_weights(g):
        if w < minweight:
            minweight = w
            nodes.clear()
        if w == minweight:
            nodes.add(n)
    print('Nodes %s with minimum weight %d' % (sorted(nodes), minweight))
    ASSERT(minweight!=INT_MAX and nodes, 'Something strange, no candidate', g)
    return nodes

def gen_nodes_with_weights(g):
    for n, d in g.nodes_iter(data=True):
        w = d['weight']
#         pred = g.pred[n] 
#         succ = g.succ[n]
#         if len(pred)==1 and len(succ)==1: # A SISO node
#             in_weight  = g[next(iter(pred))][n]['weight']
#             out_weight = g[n][next(iter(succ))]['weight']
#             if w==in_weight and w==out_weight:
#                 continue  # Skip trivial SISO nodes # Fails with simple loops!
        yield n, w

# A helper function to identify boring 2-degree nodes; len of the returned list
# would be 2 in that case.
# def get_edgeweights(g, n):
#     inedgeweights  = [ d['weight'] for d in six.itervalues( g.pred[n] ) ]
#     outedgeweights = [ d['weight'] for d in six.itervalues( g.succ[n] ) ]
#     return inedgeweights + outedgeweights

def select_next_node(candidates, hints):
    if len(candidates)==1: # There is nothing to select
        return candidates.pop()
    # TODO Figure out a tie breaking rule!
    favored_candidates = sorted(candidates & hints)
    if favored_candidates:
        print('Favored candidates:', favored_candidates)
        if len(favored_candidates) > 1:
            print('*** Ties, breaking it arbitrarily! (Loc 1)')
        n = favored_candidates[-1]
        print('Selected:', n)
        candidates.remove(n)
        return n
    print('*** Ties, breaking it arbitrarily! (Loc 2)') # len(candidates) > 1
    return candidates.pop()

def eliminate(g, n):
    # TODO Try to reduce the duplication
    cost = g.node[n]['weight']
    print('===================================================================')
    print('Eliminating node: {}, cost: {}'.format(n, cost))
    node_eqs = g.node[n].get('eqs')
    eliminated_blocks = [ node_eqs  ] if node_eqs else [ ]
    new_candidates = set()
    for k in g.pred[n]:
        print('Deleting edge {} -> {}'.format(k,n))
        edge = g[k][n]
        eliminated_blocks.append( edge.get('eqs') )
        update_neighbor(g, k, edge, new_candidates)
    for k in g.succ[n]:
        print('Deleting edge {} -> {}'.format(n,k))
        edge = g[n][k]
        eliminated_blocks.append( edge.get('eqs') )        
        update_neighbor(g, k, edge, new_candidates)
    g.remove_node(n) # <- This deletes all the edges as well
    print('New candidates:', sorted(new_candidates))
    return eliminated_blocks, new_candidates, max(0, cost)

def update_neighbor(g, k, edge, new_candidates):
    node = g.node[k]
    node['weight'] -= edge['weight']
    if node['weight']  <= 0:
        new_candidates.add(k)        

def ASSERT(cond, msg, g):
    if not cond:
        print(msg)
        plot(g, prog='sfdp')
        raise AssertionError()

if __name__ == '__main__':
    from test_ordering import run_tests
    run_tests()
