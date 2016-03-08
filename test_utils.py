# Copyright (C) 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>

from __future__ import print_function
from itertools import chain
import string
import six
import networkx as nx
from networkx.algorithms.bipartite import is_bipartite_node_set
from utils import duplicates, deserialize, DATADIR, marshal_dump
from py3compat import irange, izip

def log(*args): pass
#log = print

def create_diagonal_matrix(n, seed, rng):
    n_nodes = 2*n
    g = nx.Graph()
    eqs = list(irange(n))
    vrs = list(irange(n, n_nodes))
    g.add_nodes_from(eqs, bipartite=0)
    g.add_nodes_from(vrs, bipartite=1)
    g.add_edges_from(zip(eqs, vrs))
    return _finalize(g, n, n, rng)

def create_rnd_bipartite(n_eqs, n_vars, seed, rng, nonsing=False, c=1.0):
    # Node identifiers: random strings
    g = raw_rnd_bipartite(n_eqs, n_vars, seed, c)
    # Fill up the diagonal if requested, makes the matrix nonsingular
    if nonsing:
        eqs = (n for n, d in g.nodes_iter(data=True) if d['bipartite']==0)
        vrs = (n for n, d in g.nodes_iter(data=True) if d['bipartite']==1)
        g.add_edges_from(zip(eqs, vrs))
    return _finalize(g, n_eqs, n_vars, rng)

def raw_rnd_bipartite(n_eqs, n_vars, seed, c=1.0):
    # Node identifiers: 0, 1, ..., n_eqs + n_vars -1
    n_nodes = n_eqs + n_vars
    p = n_nodes/(n_eqs*n_vars + 1.0e-17)*c
    p = min(p, 1.0 - 1.0e-6) 
    return nx.bipartite_random_graph(n_eqs, n_vars, p, seed)

def _finalize(g, n_eqs, n_vars, rng):
    log('eqs =', n_eqs, ' vars =', n_vars, ' edges =', g.number_of_edges())
    # Relabel the nodes (issue with hashing)
    names = _get_rnd_names(n_eqs+n_vars, rng)
    mapping = { node : name for node, name in zip(g, names) }
    nx.relabel_nodes(g, mapping, copy=False)
    # Getting the bipartite node set
    eqs = set(n for n, d in g.nodes_iter(data=True) if d['bipartite']==0)
    assert len(eqs) == n_eqs
    # Fixing isolated var nodes, isolated equations are fine
    isolated_vars = (n for n in g if n not in eqs and not g[n])
    eq_list = list(eqs)
    for n in isolated_vars:
        eq = rng.choice(eq_list)
        log('Isolated node:', n)
        g.add_edge(eq, n)
    log('Nodes:', list(g))
    log('Eqs:', list(eqs))
    log()
    return g, eqs

def _get_rnd_names(n_nodes, rng):
    names = set()
    while len(names) < n_nodes:
        size = rng.randint(1, 6)
        names.add(_str_generator(size, rng))
    return names       

ALPHABET = string.ascii_letters + string.digits

def _str_generator(size, rng):
    return ''.join(rng.choice(ALPHABET) for _ in range(size))

#-------------------------------------------------------------------------------

def create_coomat(n_rows, n_cols, rng):
    g = raw_rnd_bipartite(n_rows, n_cols, rng.randint(0, 2**32))
    #
    rows, cols = [ ], [ ]
    for r in irange(n_rows):
        for c in g[r]:
            rows.append(r)
            cols.append(c-n_rows)
    #
    n_nonzeros = g.number_of_edges()
    values = [rng.randint(1, 9) for _ in irange(n_nonzeros)]
    for r, c, v in izip(rows, cols, values):
        g[r][c+n_rows]['weight'] = v
    #
    return g, rows, cols, values

################################################################################
# FIXME Clean up this mess

TEST_MATRICES = {
                 
    # The first 3 matrices do not trigger poor performance in ILP, only the full
    # graph does, that yields (n choose 2)^2 rows, which is 4356 for n=12.

    'test_1' : (   'e01  e02  e03  e04  e05  e06  e07  e08  e09  e10  e11  e12',
                 '''A                       G                   
                        B                       H               
                            C                       I           
                                D                       J       
                                    E                       K   
                                        F                       L
                    A                       G                   
                        B                       H               
                            C                       I           
                                D                       J       
                                    E                       K   
                                        F                       L''',
                    6),

    'test_2' : (   'e01  e02  e03  e04  e05  e06  e07  e08  e09  e10  e11  e12',
                 '''A                        G    H                
                        B                        H    I            
                            C                        I    J        
                                D                        J    K    
                                    E                        K    L
                    A                    F                        L
                    A    B                    G                    
                        B    C                    H                
                            C    D                    I            
                                D    E                    J        
                                    E    F                    K    
                                        F    G                    L''',
                    2),

    'test_3' : (   'e01  e02  e03  e04  e05  e06  e07  e08  e09  e10  e11',
                 '''A   B   C                               
                    A   B   C                               
                    A   B   C   D   E                       
                            C   D   E                       
                            C   D   E   F   G               
                                    E   F   G               
                                    E   F   G   H   I       
                                            G   H   I       
                                            G   H   I   J   K
                                                    I   J   K
                                                    I   J   K''',
                    6),

#     'test_4' : (   'e01  e02  e03  e04  e05  e06  e07  e08  e09  e10  e11  e12',
#                  '''A   B   C   D   E   F   G   H   I   J   K   L
#                     A   B   C   D   E   F   G   H   I   J   K   L
#                     A   B   C   D   E   F   G   H   I   J   K   L
#                     A   B   C   D   E   F   G   H   I   J   K   L
#                     A   B   C   D   E   F   G   H   I   J   K   L
#                     A   B   C   D   E   F   G   H   I   J   K   L
#                     A   B   C   D   E   F   G   H   I   J   K   L
#                     A   B   C   D   E   F   G   H   I   J   K   L
#                     A   B   C   D   E   F   G   H   I   J   K   L
#                     A   B   C   D   E   F   G   H   I   J   K   L
#                     A   B   C   D   E   F   G   H   I   J   K   L
#                     A   B   C   D   E   F   G   H   I   J   K   L'''),

    'test_5' : (   'e01  e02  e03  e04  e05  e06  e07  e08  e09  e10  e11 e12 e13 e14 e15 e16 e17 e18 e19 e20 e21 e22 e23 e24 e25',
                 '''A   B   C                                                                                       
                    A   B   C                                                                                       
                    A   B   C   D   E                                                                               
                            C   D   E                                                                               
                            C   D   E   F   G                                                                       
                                    E   F   G                                                                       
                                    E   F   G   H   I                                                               
                                            G   H   I                                                               
                                            G   H   I   J   K                                                       
                                                    I   J   K                                                       
                                                    I   J   K   L   M                                               
                                                            K   L   M                                               
                                                            K   L   M   N   O                                       
                                                                    M   N   O                                       
                                                                    M   N   O   P   Q                               
                                                                            O   P   Q                               
                                                                            O   P   Q   R   S                       
                                                                                    Q   R   S                       
                                                                                    Q   R   S   T   U               
                                                                                            S   T   U               
                                                                                            S   T   U   V   W       
                                                                                                    U   V   W       
                                                                                                    U   V   W   X   Y
                                                                                                            W   X   Y
                                                                                                            W   X   Y''',
                   13),                    

}

def to_bipartite_from_test_string(mat_str):
    # Unaware of the optional opt in the tuple (dm_decomp does not have opt)
    rows = mat_str[0].split()
    cols_rowwise = [line.split() for line in mat_str[1].splitlines()]
    # check rows for typos
    eqs = set(rows)
    assert len(eqs) == len(rows), (sorted(eqs), sorted(rows))
    assert len(rows) == len(cols_rowwise)
    # check cols for typos
    all_cols = set(chain.from_iterable(cols for cols in cols_rowwise))
    both_row_and_col = sorted( eqs & all_cols ) 
    assert not both_row_and_col, both_row_and_col
    # check cols for duplicates
    for r, cols in izip(rows, cols_rowwise):
        dups = duplicates(cols)
        assert not dups, 'Duplicate column IDs {} in row {}'.format(dups, r)
    #print(rows)
    #print(cols_rowwise)
    g = nx.Graph()
    g.add_nodes_from(rows)
    g.add_nodes_from(all_cols)
    g.add_edges_from((r,c) for r, cols in izip(rows, cols_rowwise) for c in cols)
    assert is_bipartite_node_set(g, eqs)
    return g, eqs

def solve_test_matrices(solve_function, log, skip={}):
    for name, (rows, cols_rowwise, opt) in sorted(six.iteritems(TEST_MATRICES)):
        if name in skip:
            continue
        log(name)
        g, eqs = to_bipartite_from_test_string((rows, cols_rowwise))
        _, _, _, tear_set, _ = solve_function(g, eqs)
        cost = len(tear_set)
        assert opt == cost, (opt, cost)
        log() 

DIFFICULT_FOR_ILP = [
    # 12x12 nonsingular square systems. c rougly: nonzeros per row, r: rows
    # As for c, see raw_rnd_bipartite in test_utils, computation of p.
    # This graph was the worst out of 500 random tests:
    ('ilp_poor_performance_c=2_r=509.pkl.gz',  5),
    # This graph was the worst out of 100 random tests:
    ('ilp_poor_performance_c=3_r=1337.pkl.gz', 7),
    # This graph was the worst out of 100 random tests:
    ('ilp_poor_performance_c=5_r=2997.pkl.gz', 9),
] 

def solve_difficult_for_ilp(solve_function, log):
    for filename, opt in DIFFICULT_FOR_ILP:
        log(filename)
        g, eqs = deserialize(DATADIR + filename)
        _, _, _, tear_set, _ = solve_function(g, eqs)
        cost = len(tear_set)
        assert opt == cost, (opt, cost)
        log()

#-------------------------------------------------------------------------------
# FIXME Document better what's going on

def parse_nauty_edgelist(size):
#     all_edgelists = [ ]
#     for edgelist in edgelists(size):
#         #print(edgelist)
#         assert len(edgelist) % 2 == 0
#         g = nx.Graph()
#         g.add_edges_from(e for e in izip(edgelist[::2], edgelist[1::2]))
#         assert len(g) == 2*size
#         assert is_bipartite_node_set(g, irange(size))
#         #print(len(g), g.number_of_edges())
#         all_edgelists.append(edgelist)
    # For size > 6:
    all_edgelists = list(edgelists(size))
    print('There were', len(all_edgelists), 'graphs')   
    marshal_dump(all_edgelists, '/tmp/filt_n'+str(size)+'.bin')


def edgelists(size):
    with open('/tmp/filt_n'+str(size)+'.txt') as f:
        next(f) # ugly, we ignore the first empty line
        edgelist, state = [ ], 'graph_line'
        for line in f:
            line = line.rstrip()
            if not line:
                yield edgelist
                edgelist = [ ]
                state = 'graph_line'
                continue
            if state == 'graph_line':
                #print(line)
                state = 'info_line'
                continue
            if state == 'info_line':
                #print(line)
                state = 'edge_list_lines'
                continue
            if state == 'edge_list_lines':
                edgelist.extend(int(s) for s in line.split())
        yield edgelist    


def serialize_all_bip():
    # Generate input:
    # genbg -d1 -q 3 3 | showg -e >/tmp/bip_n3.txt
    # Or assuming preprocessing:
    # genbg -d2 -D6 -c -q 7 7 | showg -e >/tmp/bip_n7.txt
    # Then transform into networkx graphs
    #parse_nauty_edgelist(1)
    #parse_nauty_edgelist(2)
    parse_nauty_edgelist(3)
    parse_nauty_edgelist(4)
    parse_nauty_edgelist(5)
    parse_nauty_edgelist(6)

# if __name__ == '__main__':
#     serialize_all_bip()
