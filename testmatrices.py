# Copyright (C) 2016 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from itertools import product, repeat
from networkx import Graph
from networkx.algorithms.bipartite import is_bipartite_node_set
from py3compat import irange, izip


def create_difficult_pattern(size):
    '''The eq ids go from 0..size-1, the column ids from size..2*size-1.
    A pathological pattern, resulting in many ties:
    | x x         | 
    |     x x     |
    |         x x |
    | x x x x     |
    |     x x x x |
    | x x     x x |  '''
    assert size % 2 == 0, size 
    rows, cols = list(irange(size)), list(irange(size, 2*size))
    g = Graph()
    half_size = size//2
    # build upper half
    for i in irange(half_size):
        g.add_edges_from(((i, size + 2*i), (i, size+ 2*i+1)))
    # build lower half
    for i in irange(half_size, size):
        k = 2*(i - half_size)
        vrs = [ size + v % size for v in irange(k, k+4) ]
        g.add_edges_from(izip(repeat(i),vrs))
    assert is_bipartite_node_set(g, rows)
    assert is_bipartite_node_set(g, cols)
    #to_pdf(g, rows, cols, '', str(size))
    #plot_dm_decomp(g, size)
    return g


def create_block_pattern(n_blocks):
    '''The eq ids go from 0..size-1, the column ids from size..2*size-1.
    A pathological pattern, resulting in many ties:
    | x x x         | 
    | x x x         |
    | x x x x x     |
    |     x x x     |
    |     x x x x x |
    |         x x x | 
    |         x x x | '''
    size = 2*n_blocks + 1
    g = Graph()
    g.add_nodes_from(irange(2*size))
    for i in irange(0, size-2, 2):
        eqs = [i, i+1, i+2]
        j = i + size
        vrs = [j, j+1, j+2]
        g.add_edges_from(product(eqs, vrs))
    #print('Nodes:', g.nodes())
    #print('Edges:', sorted(g.edges()))
    assert is_bipartite_node_set(g, set(irange(size)))
    return g, size
