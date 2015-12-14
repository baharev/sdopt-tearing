# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from six.moves import map as imap 
from networkx import DiGraph, Graph, relabel_nodes
from networkx.algorithms import bipartite
from expression_tree import add_recursively, defines_var_alias, \
                            gen_var_node_dicts, get_varnames, to_symbolic_form  
from flatten import DATADIR
from plot import dummy as plot # for plotting the process graph
from utils import deserialize

################################################################################

class Equation:
    UNIT       = 'unit'
    CONNECTION = 'connection'
    NESTING    = 'nesting'

    def __init__(self, names, expression_tree):
        self.names = names # sorted list of the variable names in the equation
        self.expression_tree = expression_tree       
        self.symbolic_form = None
        self.kind   = None
        self.in_out = None # If a CONNECTION equation
        self.unit   = None # If a UNIT equation
        # elims = { var name : symbolic solution }. The symbolic solutions is of
        # form v = f(v other), if possible, otherwise var name is not a key in 
        # the dictionary.
        self.elims  = { }
        self.solved_for = None # var eliminated
        # Modelica cannot give a unique ID for the equations, we generate one
        self.id = None

def gen_nonnesting_eqs(equations):
    return (eq for eq in equations if eq.kind != Equation.NESTING)

def gen_connection_eqs(equations):
    return (eq for eq in equations if eq.kind == Equation.CONNECTION)

def gen_unit_eqs(equations):
    return (eq for eq in equations if eq.kind == Equation.UNIT)

def is_connection(eq):
    return eq.kind==Equation.CONNECTION

def get_referenced_vars(equations, var_names, parameters):
    # parameters: list of (name, value) as returned by parameter_assignments
    referenced = set()
    for eq in gen_nonnesting_eqs(equations):
        referenced.update(eq.names)
    params = { t[0] for t in parameters }
    referenced -= params
    return [ name for name in var_names if name in referenced ]

################################################################################
# Not nice that the equations parse XML but couldn't find a better alternative.

def create_equations(xml_etree):
    equations = [ ]
    for xml_element_subtree in xml_etree.iter(tag='Equation'):
        expression_tree = DiGraph()
        root = next(iter(xml_element_subtree))
        add_recursively(expression_tree, root)
        #plot_tree(expression_tree) # implementation moved to in unused code
        names = get_varnames(expression_tree)
        equations.append(Equation(names, expression_tree))
    return equations

def parameter_assignments(tree):
    # Returns: list of (name, value) 
    # This ugliness must stay: binding equations (essentially parameter 
    # definitions) are structured differently than variable - number = 0.
    parameters = [ ]
    for bindeq in tree.iter(tag='BindingEquation'):
        # Assumption: binding equations look like 
        # parameter name = numeric literal
        itr_eq = iter(bindeq)
        param = next(itr_eq) # assert: subtree for a single var name
        name = get_full_name(param)
        itr_lit = iter(next(itr_eq)) # assert: subtree for a numeric literal
        assert_empty(itr_eq)
        literal = next(itr_lit)
        assert_empty(itr_lit)        
        assert hasattr(literal, 'tag'), literal
        assert literal.tag.endswith('Literal'), literal.tag
        value = literal.text # Note: int-s will be converted to floats later on
        #print(name, '=', value)    # Just store the literal as string?
        parameters.append((name, value))
    return parameters

def assert_empty(itr):
    s = next(itr, None)
    assert s is None, 'Expected an empty iterator, got: {}'.format(s)

def get_full_name(e):
    #segments = [ part.attrib['name'] for part in e ]
    # TODO This function is a temporary workaround; the JModelica 
    # seems to lose the array subscripts in a fuzzy manner.
    segments = [ ]
    for part in e:
        name = part.attrib['name']
        indices = [ i.text for i in part.iter(tag='IntegerLiteral') ]    
        if indices:
            name = name + '[' + ','.join(indices) + ']'
        segments.append(name)
    return '.'.join(segments)

################################################################################
 
def rewrite_equations(var_names, equations, parameters):
    # Returns a mapping: mapping[name] is the real unaliased name of the var, 
    # and the raw_connections containing only connector variable names.
    #
    # collect_aliasing also sets eq.kind to NESTING, CONNECTION or UNIT on each
    # equation, and eq.in_out on CONNECTIONS.
    raw_connections, in_aliases, out_aliases = collect_aliasing(equations)
    #
    mapping = resolve_aliases(equations, var_names, in_aliases, out_aliases)
    #
    # Resolving the in-in and out-out aliasing (nesting) by rewriting variable
    # names according to the bottom level (true) unit. The nesting equations
    # are superfluous after this pass, nevertheless they are kept for debugging.
    for eq in gen_nonnesting_eqs(equations): # nesting would become identity
        rewrite_names(eq, mapping)
    #
    for eq in equations:
        eq.symbolic_form = to_symbolic_form(eq.expression_tree)
    # Only connection equations are considered when collecting the unit names
    unit_names = collect_unit_names( gen_connection_eqs(equations) )
    #
    # The connections between atomic units used to be generated here, now done
    # in get_process_graph: Symbolic processing rewrites the equations again.
    #
    for eq in gen_unit_eqs(equations):
        eq.unit = get_unit_name(eq, mapping, unit_names, parameters)
    #dbg_print(equations)
    #
    # The minimum degree needs the eq.id to be set. The variables and the 
    # equations *must* have different names (IDs). By giving int-s as IDs for 
    # the equations, we are guaranteed not to have a name collision; var names 
    # are strings.
    for counter, eq in enumerate(equations):
        eq.id = counter
    return mapping, raw_connections

################################################################################

def collect_aliasing(equations):
    # Categorizing the aliasing equations: nesting, connection and unit. The 
    # actual connections will be resolved later by get_process_graph.
    in_aliases, out_aliases = DiGraph(), DiGraph()
    raw_connections = DiGraph()
    for i, eq in enumerate(equations):
        if defines_var_alias(eq.expression_tree):
            record_aliases(eq, i, in_aliases, out_aliases, raw_connections)
        else:
            eq.kind = Equation.UNIT
    dbg_print_aliasing( in_aliases, '\nInlet aliasing (nesting):')
    dbg_print_aliasing(out_aliases, '\nOutlet aliasing (nesting):')
    return raw_connections, in_aliases, out_aliases

# HACK The connector type and prefix (input or output) is encoded in the name
# of the connector instance.
INLET_TYPE  = '.inlet'
OUTLET_TYPE = '.outlet'

def record_aliases(eq, i, in_aliases, out_aliases, raw_connections):
    # HACK The corresponding variables of the connectors should also match, 
    # currently ignored
    ai = eq.names[0].partition(INLET_TYPE)
    bi = eq.names[1].partition(INLET_TYPE)
    ao = eq.names[0].partition(OUTLET_TYPE)
    bo = eq.names[1].partition(OUTLET_TYPE)
    # If the IN/OUTLET_TYPE is not found, t[1] is an empty string. 
    if   ai[1] and bi[1]: # in-in aliasing (nesting)
        add_nesting(in_aliases,  ai, bi, i)
        eq.kind = Equation.NESTING
    elif ao[1] and bo[1]: # out-out aliasing (nesting)
        add_nesting(out_aliases, ao, bo, i)
        eq.kind = Equation.NESTING        
    elif ai[1] and bo[1]: # b -> a or accidental aliasing within a unit
        add_connection_or_unit_equation(raw_connections, ai, bo, eq, i)
    elif ao[1] and bi[1]: # a -> b or accidental aliasing within a unit
        add_connection_or_unit_equation(raw_connections, bi, ao, eq, i)
    else:
        raise AssertionError(eq.names)

def add_nesting(digraph, a, b, i):
    # a and b are tuples, for example (cascade, .inlet, [1].f[2]) and 
    # (cascade.stages[1], .inlet, [1].f[2]).
    if len(a[0]) > len(b[0]):
        a, b = b, a
    # The name of the nested unit starts with the name of the parent unit       
    assert b[0].startswith(a[0]+'.'), (a, b, i)
    # The edge is directed towards the nested unit: a -> b.
    add_connection(digraph, a, b, i)

def add_connection(digraph, a, b, i):
    # t[1] is either the INLET_TYPE or the OUTLET_TYPE so it must not be empty.
    assert a[1] and b[1], (a, b, i)    
    src, dst = ''.join(a), ''.join(b)
    if not digraph.has_edge(src, dst):
        digraph.add_edge(src, dst, eq=[i])
    else:
        digraph[src][dst]['eq'].append(i)

def add_connection_or_unit_equation(raw_connections, inlet, outlet, eq, i):
    if inlet[0] != outlet[0]:
        eq.kind   = Equation.CONNECTION
        eq.in_out = ''.join(inlet), ''.join(outlet)
        add_connection(raw_connections, outlet, inlet, i)
    else: # (Accidental) alias variables within a unit
        eq.kind = Equation.UNIT

################################################################################

def resolve_aliases(equations, var_names, in_aliases, out_aliases):
    return {name: true_name(name,in_aliases,out_aliases) for name in var_names}

def true_name(name, in_aliases, out_aliases):
    if name in in_aliases:
        return walk_alias_chain(in_aliases, name)
    if name in out_aliases:
        return walk_alias_chain(out_aliases, name)
    return name

def walk_alias_chain(aliases, parent):
    children = aliases.edge[parent]
    while children:          
        assert len(children)==1, (parent, children)
        (parent,) = children
        children = aliases.edge[parent]
    return parent

################################################################################

def rewrite_names(eq, mapping):
    eq.names = sorted( mapping[name] for name in eq.names )
    assert len(eq.names)==len(set(eq.names)), eq.names
    # rewrite the names in the expression tree
    for d in gen_var_node_dicts(eq.expression_tree):
        old_name = d['value']
        new_name = mapping[old_name] 
        if new_name!=old_name:
            d['value']    = new_name
            d['origname'] = old_name
    # If a connection equation, map the inlet, outlet names as well
    if eq.in_out:
        i, o = eq.in_out
        eq.in_out = mapping[i], mapping[o]
        
################################################################################

def collect_unit_names(connection_eqs):
    unit_names = set()
    for eq in connection_eqs:
        unit_names.add(get_unit_prefix(eq.in_out[0]))
        unit_names.add(get_unit_prefix(eq.in_out[1]))
    return unit_names

# TODO Try to eliminate the need for calling this function
def get_unit_prefix(iolet):
    i, isep, _ = iolet.partition(INLET_TYPE)
    if isep:
        return i 
    o, osep, _ = iolet.partition(OUTLET_TYPE)
    if osep:
        return o
    raise AssertionError('Not an inlet or an outlet: "{}"'.format(iolet))

################################################################################

def get_unit_name(eq, mapping, unit_names, parameters):
    # All inlets and outlets must belong to the same unit. The other variables
    # (if there is any) must be within this unit. Parameters must be ignored.    
    units = set()
    var_prefixes = set()
    parameter_names = { t[0] for t in parameters } # To filter out constants        
    for name in eq.names:
        if name in parameter_names:
            continue
        if INLET_TYPE in name or OUTLET_TYPE in name:
            units.add(get_unit_prefix(mapping[name]))
        else:
            var_prefixes.add( name.rpartition('.')[0] )
    #          
    if len(units) == 1: # OK, the inlets and outlets belong to the same unit
        (unit,) = units
        assert unit in unit_names, (unit,eq.names)
        for p in var_prefixes:
            assert is_same_or_subunit_of(p, unit), (p,unit,eq.names)
        return unit
    # There are only variables in this unit equation. All vars must be embedded
    # in a unit (or a subunit of a unit) that we have seen in the connection
    # equations. TODO We would have to walk the unit's hierarchy to find the 
    # lowest level unit. Currently, we only check the shortest prefix level one.
    if len(units) == 0:
        prefs = sorted(var_prefixes)
        unit = prefs[0]
        assert unit in unit_names, (unit, prefs, 'See comment above')
        unit_with_dot = unit + '.'
        for p in prefs[1:]: 
            # p must be longer than unit=prefs[0] due to sorting var_prefixes
            assert p.startswith(unit_with_dot), (prefs, p, 'See comment above')
        return unit
    raise AssertionError((units, var_prefixes, eq.names))

def is_same_or_subunit_of(prefix, unit):
    if   len(prefix) > len(unit):
        return prefix.startswith(unit+'.')
    elif len(prefix) == len(unit):
        return prefix == unit
    else:
        raise AssertionError((prefix, unit))

################################################################################

def get_process_graph(equations, raw_connections, mapping, parameters):
    '''Returns a directed, weighted graph, describing the process graph; in 
    other words, the connections of the atomic units. The returned graph holds a
    reference to the equations. The node and edge weights are set according to 
    the number of equations and variables, such that block ordering 
    (see run_elimination) can be invoked with the returned graph.'''
    plot(raw_connections)
    process_graph = get_connections_of_atomic_units(raw_connections, mapping)
    initialize_node_and_edge_dict(process_graph)
    record_all_unit_eqs(process_graph, equations, parameters)
    record_connection_eqs(process_graph, equations)
    finalize_nodes(process_graph)
    return process_graph

def get_connections_of_atomic_units(raw_connections, mapping):
    m = { n : get_unit_prefix(mapping[n]) for n in raw_connections }
    connections = relabel_nodes(raw_connections, m, copy=True)
    plot(connections, prog='sfdp')
    return connections

def initialize_node_and_edge_dict(connections):
    for n in connections:
        d = connections.node[n]        
        d['vars'] = set()
        d['eqs']  = [ ]
        # d['weight'] will be assigned to only later, in finalize_nodes
    for _, _, d in connections.edges_iter(data=True):
        d['weight'] = 0
        d['eqs'] = [ ]

def record_all_unit_eqs(connections, equations, parameters): 
    param_names = { t[0] for t in parameters } # TODO Duplication get_unit_name
    for eq in gen_unit_eqs(equations):
        unit = eq.unit
        assert unit in connections, (unit, eq.names)
        d = connections.node[unit] 
        true_vars = set(name for name in eq.names if name not in param_names)
        d['vars'].update(true_vars)
        d['eqs'].append(eq)

def record_connection_eqs(connections, equations):
    for eq in gen_connection_eqs(equations):
        inlet, outlet = eq.in_out
        src, dst = get_unit_prefix(outlet), get_unit_prefix(inlet)
        assert connections.has_edge(src, dst), (inlet, outlet)
        d = connections[src][dst]
        d['weight'] += 1
        d['eqs'].append(eq)
        connections.node[src]['vars'].add(outlet)
        connections.node[dst]['vars'].add(inlet)

def finalize_nodes(connections):
    for n in connections:
        d = connections.node[n] 
        d['vars'] = sorted(d['vars'])
        nvars = len(d['vars'])
        neqs  = len(d['eqs'])
        d['weight'] = nvars - neqs

################################################################################

def gen_unit_conn(blocks_in_order):
    '''Returns: counter, unit_eqs, connection triples (y, x, equation id) where
    y = x and x belongs to the unit. See also dbg_ordered_blocks which uses it 
    and prints the blocks.'''
    itr = iter(blocks_in_order)
    counter = 1
    unit_eqs = next(itr)
    while unit_eqs:
        unit = unit_eqs[0].unit
        assert unit is not None, unit_eqs[0].names
        conns, next_block = get_connections_eliminated_after(unit, itr)
        yield counter, unit_eqs, conns
        counter += 1
        unit_eqs = next_block

def get_connections_eliminated_after(unit, itr):
    conns = [ ]    
    for next_block in itr:
        if next_block[0].kind==Equation.CONNECTION:
            conns += extract_iolets_in_correct_order(next_block, unit)
        else:
            assert next_block[0].kind==Equation.UNIT
            return conns, next_block
    return conns, None

def extract_iolets_in_correct_order(next_block, unit):
    # The iolets must be swapped such that the iolets of the just eliminated 
    # unit are on the right hand side.
    if get_unit_prefix(next_block[0].in_out[1]) == unit:
        func = lambda eq: (eq.in_out[0], eq.in_out[1], eq.id)
    else:
        func = lambda eq: (eq.in_out[1], eq.in_out[0], eq.id)
    return list(imap(func, next_block))

################################################################################

def dbg_ordered_blocks(blocks_in_order):
    'Illustrates how to iterate over blocks_in_order, see also gen_unit_conn.'
    for _, unit_block, conn_triples in gen_unit_conn(blocks_in_order):
        print('\nUnit')
        for eq in unit_block:
            print(eq.names)
            assert eq.kind == Equation.UNIT, eq.kind
        if conn_triples:
            print('Connections')
        for y, x, _ in conn_triples:
            # y = x, and equation id, if any, is ignored
            print('{} <- {}'.format(y, x))

################################################################################

def dbg_print_aliasing(aliases, header_str):
    print(header_str)
    src_nodes = sorted(n for n in aliases if aliases.in_degree(n)==0)
    for parent in src_nodes:
        dbg_bottom_level_alias(aliases, parent)

def dbg_bottom_level_alias(aliases, parent):
    # Somewhat duplicate of walk_alias_chain
    if parent not in aliases:          # unaliased connectors are stored only 
        print(parent, '(not aliased)') # in connections
        return parent
    children = aliases.edge[parent]
    while children:
        print(parent)            
        assert len(children)==1, (parent, children)        
        (parent,) = children
        children = aliases.edge[parent]
    print(parent)
    return parent

def dbg_print(equations):
    print('\nConnections:')
    for eq in gen_connection_eqs(equations):
        print(eq.in_out[1], '->', eq.in_out[0])
    #
    print('\nUnit equations:')
    for eq in gen_unit_eqs(equations):
        print(eq.unit, eq.names)

################################################################################

def to_bipartite_graph(equations):
    '''Returns: g, eqs, forbidden. Here, g represents the system of equations as
    an undirected bipartite graph; eqs gives the node ids of the equations; the 
    forbidden set contains the (eq,var) pairs, where variable var cannot be 
    safely or explicitly eliminated from equation eq.'''
    g, eqs, forbidden = Graph(), set(), set()
    #
    for eq in gen_unit_eqs(equations):
        eqs.add(eq.id)
        for var in eq.names:
            g.add_edge(eq.id, var)
            if var not in eq.elims:
                forbidden.add((eq.id,var))
    #
    for eq in gen_connection_eqs(equations):
        eqs.add(eq.id)
        for var in eq.names:
            g.add_edge(eq.id, var)
    #
    info_on_bipartite_graph(g, eqs, forbidden)
    #from plot import plot
    #plot(g, prog='sfdp')
    #from utils import serialize
    #serialize((g, eqs, forbidden), 'JacobsenILOSimpBounds.pkl.gz')
    return g, eqs, forbidden

def info_on_bipartite_graph(g, eqs, forbidden, log=print):
    log()
    log('Unordered equations (bipartite, no blocks)')
    log('Equations:', len(eqs))
    log('Variables:', len(g)-len(eqs))
    log('Non-zeros:', g.number_of_edges())
    log('Forbidden:', len(forbidden))
    assert bipartite.is_bipartite_node_set(g, eqs)    

def read_bipartite_graph(problem_name):
    # A serialized undirected bipartite graph equations and variables, without 
    # the blocks and unordered, as in to_bipartite_graph. 
    g, eqs, forbidden = deserialize(DATADIR+problem_name+'.pkl.gz')
    info_on_bipartite_graph(g, eqs, forbidden)
    return g, eqs, forbidden
