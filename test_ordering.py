# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from networkx import DiGraph
from benchmarks import gen_benchmark_digraphs
from block_ordering import run_elimination
from plot import dummy as plot
from utils import info

__all__ = [ 'run_tests' ]

def run_tests():
    for g, hints in gen_testgraphs():
        plot(g, prog='sfdp')
        run_elimination(g, hints)
    # Comparing my block ordering to tearing (comparing apples to oranges).
    for g in gen_benchmark_digraphs():
        if g.number_of_selfloops():
            print('Skipping it, as this digraph has self-loops\n')
            continue
        plot(g, prog='sfdp')
        info(g)
        run_elimination(g)

def gen_testgraphs():
    # Returns a test digraph and a set of nodes that might be preferred when 
    # breaking ties (e.g. sources and sinks in the (imaginary) original graph).
    testgraphs = sorted(g for g in dir(TestGraphs) if not g.startswith('_'))
    obj = TestGraphs()
    for graph_name in testgraphs:
        print("Creating test graph '%s'" % graph_name)
        yield getattr(obj, graph_name)()

#-------------------------------------------------------------------------------

class TestGraphs:
#     def connected_cycles(self):
#         g = nx.DiGraph()
#         g.add_cycle([0, 1, 2, 3])  
#         g.add_cycle([4, 5, 6, 7])
#         g.add_edge(1,4)
#         return g
#     
#     def three_cycles(self):
#         g = nx.DiGraph()
#         g.add_cycle([0, 1, 2, 3])
#         g.add_cycle([1, 2, 7, 4])    
#         g.add_cycle([4, 5, 6, 7])
#         return g
#     
#     def cascade(self):
#         g = nx.DiGraph()
#         g.add_cycle([0, 1, 2, 3])
#         g.add_cycle([1, 2, 7, 4])    
#         g.add_cycle([4, 5, 6, 7])
#         g.add_cycle([5, 6, 10, 11])
#         g.add_cycle([8, 9, 10, 11])
#         g.add_edge(3,12)
#         g.add_edge(9,13)
#         return g    
#     
#     def nested(self):
#         g = nx.DiGraph()
#         g.add_cycle([0, 1, 2])    
#         g.add_cycle([3, 0, 1, 4, 5])
#         return g 
#     
#     def grid(self):
#         g = nx.DiGraph()
#         g.add_cycle([0, 1, 2, 3])
#         g.add_cycle([2, 3, 4, 5])
#         g.add_cycle([5, 2, 6, 7])
#         g.add_cycle([1, 2, 6, 8])
#         return g
#     
#     def four_branch_loop(self):
#         g = nx.DiGraph()
#         g.add_path([0, 1, 2, 3,  4,  5])
#         g.add_path([6, 7, 8, 9, 10, 11])    
#         g.add_cycle([2, 3, 8, 9])
#         return g

    def process(self):
        g = DiGraph()
        g.add_cycle([1, 2, 3])
        g.add_cycle([0, 1, 2, 4])
        g.add_cycle([5, 6, 7])
        g.add_cycle([0, 1, 2, 4, 5, 6, 8])
        hints = { 0, 8 }
        set_edge_and_node_weights(g, extra_guesses_at=[2, 4, 6, 8])
        return g, hints        

    def cascade(self):
        g = DiGraph()
        g.add_cycle([0, 1, 2, 3])
        g.add_cycle([2, 3, 4, 5])    
        g.add_cycle([4, 5, 6, 7])
        set_edge_and_node_weights(g, extra_guesses_at=[1, 3, 5, 7])
        return g, { 1 }

    def loop(self):
        # This only has trivial SISO nodes and no hints 
        g = DiGraph()
        g.add_cycle([0, 1])
        g.node[0]['weight'] = 1
        g.node[1]['weight'] = 1        
        g[0][1]['weight'] = 1
        g[1][0]['weight'] = 1
        return g, set()
 
#-------------------------------------------------------------------------------

def set_edge_and_node_weights(g, extra_guesses_at):
    C = 3
    for n in g:
        g.node[n]['weight'] = C*g.in_degree(n)
    for n in extra_guesses_at:
        g.node[n]['weight'] += 1            
    for _, _, d in g.edges_iter(data=True):
        d['weight'] = C

if __name__ == '__main__':
    run_tests()
