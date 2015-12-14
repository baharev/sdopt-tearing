# Copyright (C) 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>

import string
import networkx as nx
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

def create_rnd_bipartite(n_eqs, n_vars, seed, rng, nonsing=False):
    # Node identifiers: random strings
    g = raw_rnd_bipartite(n_eqs, n_vars, seed)
    # Fill up the diagonal if requested, makes the matrix nonsingular
    if nonsing:
        eqs = (n for n, d in g.nodes_iter(data=True) if d['bipartite']==0)
        vrs = (n for n, d in g.nodes_iter(data=True) if d['bipartite']==1)
        g.add_edges_from(zip(eqs, vrs))
    return _finalize(g, n_eqs, n_vars, rng)

def raw_rnd_bipartite(n_eqs, n_vars, seed):
    # Node identifiers: 0, 1, ..., n_eqs + n_vars -1
    n_nodes = n_eqs + n_vars
    p = n_nodes/(n_eqs*n_vars + 1.0e-17)
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

ALPHABET = string.letters + string.digits

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
