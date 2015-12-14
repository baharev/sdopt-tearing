# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from networkx import Graph
from plot_ordering import plot_bipartite
from test_tearing import grouper
from order_util import get_row_col_perm

__all__ = [ 'create_test_problem']

def main():
    # lazy import, otherwise gurobipy creates gurobi.log
    from ilp_tear import solve_problem as ilp_tearing
    for g, eqs, forbidden in gen_testproblems():
        dag, tears, sinks, order = ilp_tearing(g, eqs, forbidden)
        # Get the spiked form (row and column permutation) from the ordering:
        row_perm, col_perm = get_row_col_perm(eqs, dag, tears, sinks, order)
        plot_bipartite(g, forbidden, row_perm, col_perm)

def gen_testproblems():
    # yields: (undirected bipartite graph, equations, forbidden edges)
    print('===============================================================')    
    for problem_name, edge_list in TESTPROBLEMS.iteritems():
        print("Creating test problem '%s'" % problem_name)
        yield edgelist_to_graph_eqs_forbidden(edge_list) 
        print('===============================================================')

def edgelist_to_graph_eqs_forbidden(edge_list):
    W = 'weight'
    edges = [ (e,v,{W:1}) for e,v in grouper(edge_list.split(),2) ]
    forbidden = set()
    g = Graph(edges)
    eqs = [ n for n in g if n[0]=='e' ]
    return g, eqs, forbidden


TESTPROBLEMS = {

'b 3, s 1, N 2, d 2' : '''
e0  x2
e0  x0
e1  x0
e1  x1
e2  x1
e2  x2
e3  x5
e3  x3
e4  x3
e4  x4
e5  x4
e5  x5
e0  x6
e0  x7
e1  x6
e1  x7
e2  x6
e2  x7
e6  x3
e6  x4
e6  x5
e7  x3
e7  x4
e7  x5
e3  x0
e3  x1
e3  x2
e4  x0
e4  x1
e4  x2
e5  x0
e5  x1
e5  x2
''',

'b 3, s 1, N 3, d 1' : '''
e0  x2
e0  x0
e1  x0
e1  x1
e2  x1
e2  x2
e3  x5
e3  x3
e4  x3
e4  x4
e5  x4
e5  x5
e6  x8
e6  x6
e7  x6
e7  x7
e8  x7
e8  x8
e0  x9
e1  x9
e2  x9
e9  x6
e9  x7
e9  x8
e3  x0
e3  x1
e3  x2
e4  x0
e4  x1
e4  x2
e5  x0
e5  x1
e5  x2
e6  x3
e6  x4
e6  x5
e7  x3
e7  x4
e7  x5
e8  x3
e8  x4
e8  x5
''',

}


def create_test_problem(opt):
    b = opt['block size'] 
    s = opt['spike width']
    N = opt['number of blocks']
    d = opt['border width']
    spp = [ ]
    
    block_pattern = [ ]
    for i in xrange(b):
        block_pattern.extend((i, j if j>=0 else j+b) for j in xrange(i-s, i+1))
    
    # Diagonal blocks  
    for shift in xrange(0, b*N, b):
        blk = map(lambda t: (t[0]+shift,t[1]+shift), block_pattern)
        spp.extend(blk)
    
    # Border
    spp.extend((i,j) for i in xrange(b) for j in xrange(b*N, b*N+d))
    
    # Final equations
    pattern = list(xrange(b*(N-1), b*N))
    spp.extend((i,j) for i in xrange(b*N, b*N+d) for j in pattern)
    
    # Dense connecting blocks
    dense_pattern = [(i,j) for i in xrange(b) for j in xrange(b)]
    for i_shift in xrange(b, b*N, b):
        j_shift = i_shift - b
        blk = map(lambda t: (t[0]+i_shift,t[1]+j_shift), dense_pattern)
        spp.extend(blk)
    
    g = Graph( ('e%02d'%i, 'x%02d'%j) for i, j in spp )
    eqs = [ n for n in g if n[0]=='e' ]
       
    row_perm = sorted(eqs)
    col_perm = sorted(n for n in g if n not in eqs)
    plot_bipartite(g, set(), row_perm, col_perm)

    return g, eqs


if __name__ == '__main__':
    main()
