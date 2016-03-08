# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
import six
import networkx as nx
from networkx.algorithms import bipartite
from block_ordering import run_elimination
from equations import dbg_ordered_blocks, gen_unit_conn, is_connection
from flatten import DATADIR
from min_degree import min_degree
from utils import deserialize

__all__ = [ 'total_order', 'total_ordering_to_dag' ]

#-------------------------------------------------------------------------------

class Block:
    def __init__(self, sources, conn_triples):
        self.eqs = [ ]
        self.vars = [ ] # variable names only
        self.conn_triples = conn_triples
        self.tears = sources
        self.resids = [ ]

#-------------------------------------------------------------------------------

# TODO Figure out an abstraction for the spiked form because tear_blocks, 
# get_row_col_perm, to_spiked_form and get_spiked_form_rowwise all do 
# essentially the same.

#-------------------------------------------------------------------------------
# Implementation boundary: only total_order is referenced outside of this block

def total_order(process_graph):
    '''Orders the blocks, then the equations within the blocks. Returns a tuple
    of ([alltears], [allresids], [blocks]) where the type of blocks is Block.
    The returned blocks hold references to the equations, which in turn were 
    passed in with the process_graph. Side-effect: sets eq.solved_for on the 
    equations  which interferes with later orderings if not reset with 
    reset_solved_for.'''
    blocks_in_order = run_elimination(process_graph)
    dbg_ordered_blocks(blocks_in_order)
    # Order the equations within each block
    alltears_allresids_blocks = tear_blocks(blocks_in_order)
    return alltears_allresids_blocks

def tear_blocks(blocks_in_order, log=print):
    seen, blocks, tear_stack, alltears, allresids = set(), [ ], [ ], [ ], [ ]
    for _, unit_block, conn_triples in gen_unit_conn(blocks_in_order):
        # eqs is { eq_id : eq }, as min_degree requires node IDs for building a
        # graph
        g, eqs, sink_eqs, forbidden = create_bipartite_graph(unit_block, seen)
        dag, sources, sinks, order = min_degree(g, eqs, forbidden)
        # Record the tear variables first
        blk = Block(sources, conn_triples)
        tear_stack.extend(blk.tears)
        alltears.extend(blk.tears)
        # Elimination along the upper envelope, in topologically sorted order
        elim_envelope(dag, eqs, order, sinks, tear_stack, blk)
        # Add the isolated nodes (sink equations)
        add_sink_eqs(sink_eqs, tear_stack, blk)
        # Finish up the bookkeeping 
        allresids.extend(blk.resids)
        blocks.append(blk)
        for eq in unit_block:
            seen.update(eq.names)
        for y, x, _ in conn_triples:
            seen.update( (y,x) )
    assert len(tear_stack)==0, tear_stack
    return alltears, allresids, blocks

def elim_envelope(dag, eqs, order, sinks, tear_stack, blk):
    eqs_in_order = ((n,eqs[n]) for n in order if n in eqs)
    resids = set(sinks)
    for eq_id, eq in eqs_in_order:
        if eq_id not in resids:
            # eq must have exactly one out edge in a valid elimination order
            (var,) = dag[eq_id]
            assert var in eq.elims, (var, eq.names)
            eq.solved_for = var
        else:
            var = tear_stack.pop() # Throws if structurally singular
            blk.resids.append(eq)
        blk.eqs.append(eq)
        blk.vars.append(var)

def add_sink_eqs(sink_eqs, tear_stack, blk):
    # The sink_eqs are isolated nodes, with no in-edges in this block
    # Should I move this above the tears?
    for _, eq in sorted(six.iteritems(sink_eqs), key=lambda t: t[0]):
        var = tear_stack.pop() # Throws if structurally singular
        blk.resids.append(eq)
        blk.eqs.append(eq)
        blk.vars.append(var)

def create_bipartite_graph(unit_block, seen):
    g = nx.Graph()
    eqs, sink_eqs, forbidden = {}, {}, [] # forbidden will be converted to set
    for eq in unit_block:
        eq_id = eq.id
        assert eq_id is not None, eq.names
        assert eq_id not in eqs and eq_id not in sink_eqs, (eq.id, eq.names)
        edges, notelim = get_edgelist_forbidden(eq_id, eq.names, eq.elims, seen)
        if edges:
            eqs[eq_id] = eq
            g.add_edges_from(edges)
        else:
            sink_eqs[eq_id] = eq
        forbidden += notelim
    return g, eqs, sink_eqs, set(forbidden)

def get_edgelist_forbidden(eq_id, var_names, allowed, seen_vars):
    edge_list, forbidden = [ ], [ ]
    unseen_vars = ( v for v in var_names if v not in seen_vars )
    for var in unseen_vars:
        if var not in allowed:
            forbidden.append((eq_id, var))
            w = -1
        else:
            w = 1
        edge_list.append( (eq_id,var,{'weight':w}) )
    return edge_list, forbidden

#-------------------------------------------------------------------------------

def to_spiked_form(equations, torn_blocks):
    # The torn_blocks come from BLT + Cellier's heuristic with lookahead
    id_eq = { eq.id : eq for eq in equations }
    alltears, allresids, blocks = [ ], [ ], [ ]
    tear_stack = [ ]
    for block in torn_blocks:
        tears, elimination, residuals = block
        # Record the tear variables first
        blk = Block(tears, conn_triples=[])
        tear_stack.extend(tears)
        alltears.extend(tears)
        # Elimination along the upper envelope, in topologically sorted order
        for eq_id, var in elimination:
            eq = id_eq[eq_id]
            set_fake_elimination_if_necessary(eq, var)
            blk.eqs.append(eq)
            blk.vars.append(var)
            assert var in eq.elims, (var, eq.names)
            eq.solved_for = var
        # Add the isolated nodes (sink equations)
        for r in residuals:
            eq = id_eq[r]
            var = tear_stack.pop() # Throws if structurally singular
            blk.eqs.append(eq)
            blk.vars.append(var)
            blk.resids.append(eq)
        #
        allresids.extend(blk.resids)
        blocks.append(blk)
    return alltears, allresids, blocks

#-------------------------------------------------------------------------------

def reset_solved_for(equations):
    # Computing an ordering has the side-effect that it sets eq.solved_for,
    # we undo that here
    for eq in equations:
        eq.solved_for = None

#-------------------------------------------------------------------------------
# A helper function for the bipartite case, without any block structure

def to_one_block(row_perm, col_perm, equations, match, tears, sinks):
    '''Returns the whole system as one big block, the type of the returned
    object is Block. The returned block holds references to the equations. 
    Side-effect: sets eq.solved_for on the equations which interferes with later
    orderings if not reset with reset_solved_for.'''
    # Hack: Instead of writing a custom codegen for the bipartite case, we 
    # simply put everything into one big block and re-use the codegen written 
    # for the case with blocks. 
    blk = Block(tears, conn_triples=[])
    id_eq = { eq.id : eq for eq in equations }
    blk.eqs = [ id_eq[r] for r in row_perm ]
    blk.vars = col_perm
    blk.resids = [ id_eq[r] for r in sinks ]
    resids = set(sinks)
    for n in row_perm:
        if n not in resids:
            var = match[n]
            eq = id_eq[n]
            set_fake_elimination_if_necessary(eq, var)
            eq.solved_for = var
    return blk

#-------------------------------------------------------------------------------

def set_fake_elimination_if_necessary(eq, var):
    if is_connection(eq):
        # Connection equations have no eq.elims and are treated 
        # differently in codegen; so we do yet another hack to be able 
        # to treat them like unit equations
        x, y = eq.names
        eq.elims[var] = var + ' = ' + (x if x != var else y)    

#-------------------------------------------------------------------------------
# These are helper functions, written when the ILP-based tearing was 
# implemented. The purpose of these functions is to pass a good quality feasible
# solution to the ILP-based tearing as initial point.

def total_ordering_to_dag(blocks):
    # blocks: as returned from tear_blocks (called by total_order). Returns: the
    # tuple of (eq-var bipartite DAG corresponding to the total ordering, the 
    # equations node set of this graph, the forbidden edges, the tear variables
    # (sources) in the DAG).
    dag, eqs, forbidden = nx.DiGraph(), set(), set() 
    for blk in blocks:
        #
        for eq in blk.eqs:
            eqs.add(eq.id)
            for name in eq.names:
                u, v = (eq.id, name) if name == eq.solved_for else (name, eq.id)
                dag.add_edge(u, v)
                if name not in eq.elims:
                    forbidden.add( (eq.id, name) )
        #
        for y, x, eq_id in blk.conn_triples:
            eqs.add(eq_id)
            dag.add_edge(x, eq_id)
            dag.add_edge(eq_id, y)
    #
    tears = sorted(n for n in dag if n not in eqs and len(dag.pred[n])==0)
    #
    info_on_total_order_dag(dag, eqs, forbidden, tears)
    #from utils import serialize
    #serialize((dag, eqs, forbidden, tears), 'JacobsenILOSimpBounds_TO.pkl.gz')
    return dag, eqs, forbidden, tears

def info_on_total_order_dag(dag, eqs, forbidden, tears):
    print()
    print('Ordered DAG (bipartite, no blocks)')
    print('Equations:', len(eqs))
    print('Variables:', len(dag)-len(eqs))
    print('Non-zeros:', dag.number_of_edges())
    print('Forbidden:', len(forbidden))
    print('Tears:', len(tears))
    print('Tear variables:', tears)
    assert bipartite.is_bipartite_node_set(dag.to_undirected(), eqs)
    assert nx.is_directed_acyclic_graph(dag)

def read_feasible_solution(problem_name, g, g_eqs, g_forbidden):
    # A serialized DAG as in total_ordering_to_dag. It is a feasible solution 
    # (a particular elimination order) of g, the undirected bipartite graph of 
    # equations and variables, without the blocks.
    filename = DATADIR + problem_name + '_TO.pkl.gz'
    dag, eqs, forbidden, tears = deserialize(filename)
    # Check if dag indeed corresponds to g
    assert eqs == g_eqs
    assert forbidden == g_forbidden
    assert sorted(dag) == sorted(g)
    info_on_total_order_dag(dag, eqs, forbidden, tears)
    return dag, tears, len(tears)
