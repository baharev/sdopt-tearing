# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from copy import deepcopy
from codegen import dump_ampl, dump_pycode, to_unstructured_ampl
from equations import create_equations, parameter_assignments, \
                      rewrite_equations, get_referenced_vars, \
                      get_process_graph, to_bipartite_graph, to_symbolic_form, \
                      gen_nonnesting_eqs
from flatten import get_etree_without_namespace
from heap_md import matching_to_dag, min_degree as heap_mindegree
from min_degree import min_degree as lookahead_mindegree
from order_util import get_row_col_perm, deterministic_topological_sort
from plot_ordering import plot_ordering, plot_bipartite
from sympy_tree import set_symbolic_eliminations
from tearing import tearing as classic_tearing
from total_ordering import total_order as hierarchical_heuristic_tearing, \
                           to_one_block, to_spiked_form
from variable import extract_all_vars, show_missing_bounds


__all__ = ('ModelWithInletsAndOutlets', 'ModelWithNoReconstruction')


class ModelWithInletsAndOutlets:
    
    def __init__(self, xml_filename):
        # These come from the FMUX XML file:
        self.var_names  = [ ]
        self.all_var_bounds = { } # { name : (lb, ub) }
        self.parameters = [ ] # elements: (name, value) pairs
        self.equations  = [ ] # type of the elements: Equation
        # The followings come from reconstruction:
        self.mapping = None # { name : the real unaliased name of the var }
        self.raw_connections = None # nx.DiGraph, only needed temporarily
        self.bounds = None # only for referenced variables: { name : (lb, ub) }
        self.process_graph = None # nx.DiGraph, see get_process_graph
        #
        setup_model_with_iolets(xml_filename, self)
    
    def run_ilp_tearing(self):
        # Ordering the whole bipartite graph optimally with an algorithm based 
        # on integer linear programming (ILP), but without the block structure.
        # This function requires a working Gurobi installation.
        from ilp_tear import solve_problem as ilp_based_tearing
        return self.__run_bipartite_tearing( ilp_based_tearing )
    
    def run_mindegree(self):
        # FIXME A hack to match the old interface; the other ordering algorithms
        #       should be fixed instead
        def hacked_heap_md(g, eqs, forbidden):
            rowp, colp, matches, tears, sinks = heap_mindegree(g, eqs, forbidden)
            dag = matching_to_dag(g, eqs, forbidden, rowp, colp, matches, tears, sinks)
            nbunch = list(colp)
            nbunch.extend(rowp)
            order = deterministic_topological_sort(dag, nbunch)
            return dag, tears, sinks, order
        return self.__run_bipartite_tearing( hacked_heap_md )
    
    def run_mindegree_with_lookahead(self):
        return self.__run_bipartite_tearing( lookahead_mindegree )
    
    def run_classic_tearing(self):
        equations, g, eqs, forbidden = self.__create_bipartite_repr()
        torn_blocks = classic_tearing(g, eqs, forbidden)
        func = lambda: to_spiked_form(equations, torn_blocks)
        return self.__get_spiked_form_with_blocks( func ) 
    
    def run_hierarchical_heuristic_tearing(self):
        process_graph = deepcopy(self.process_graph)
        func = lambda: hierarchical_heuristic_tearing(process_graph)
        return self.__get_spiked_form_with_blocks( func )
    
    def __run_bipartite_tearing(self, tearing_algorithm):
        equations, g, eqs, forbidden = self.__create_bipartite_repr()
        dag, tears, sinks, order = tearing_algorithm(g, eqs, forbidden)
        # Get the spiked form (row and column permutation) from the ordering:
        row_perm, col_perm = get_row_col_perm(eqs, dag, tears, sinks, order)
        # A hack, turn the whole system into one big block:
        blk = to_one_block(row_perm, col_perm, equations, dag, tears, sinks)
        # blk holds references to the equations
        # Stuff for plotting and for code generation:
        return BipartiteTearingResult( (g, eqs, forbidden), row_perm, col_perm,\
                                        tears, blk, self.bounds.copy() )    
    
    def __create_bipartite_repr(self):
        # The ordering algorithms may write the equations, so do a deepcopy
        equations = deepcopy(self.equations)
        equations = list(gen_nonnesting_eqs(equations))
        g, eqs, forbidden =  to_bipartite_graph(equations)
        return equations, g, eqs, forbidden

    def __get_spiked_form_with_blocks(self, func):
        alltears, allresids, blocks = func()
        bounds = self.bounds.copy()
        return BlockTearingResult(alltears, allresids, blocks, bounds)


#-------------------------------------------------------------------------------

class ModelWithNoReconstruction:
    
    def __init__(self, xml_filename):
        self.var_names  = [ ]
        self.all_var_bounds = { }
        self.parameters = [ ]
        self.equations  = [ ]
        read_model(xml_filename, self)
        for eq in self.equations:
            eq.symbolic_form = to_symbolic_form(eq.expression_tree)
    
    def translate_to_ampl(self):
        m = self
        lines = to_unstructured_ampl(m.var_names, m.parameters, m.equations)
        for line in lines:
            print(line)

#-------------------------------------------------------------------------------

class BipartiteTearingResult:
    
    def __init__(self, g_eqs_forbidden, row_perm, col_perm, tears, blk, bounds):
        # Holds references to the equations through blk
        self.g_eqs_forbidden = g_eqs_forbidden # see to_bipartite_graph
        self.row_perm = row_perm # row permutation in the optimal ordering
        self.col_perm = col_perm # column permutation in the optimal ordering 
        self.tears = tears # list of variables names that were torn
        self.as_one_block = blk # a hack, see the comment in to_one_block
        self.bounds = bounds  # only for referenced variables: {name : (lb, ub)}
    
    def generate_ampl_code(self):
        blk = self.as_one_block
        dump_ampl(self.tears, blk.resids, [blk], self.bounds)
    
    def generate_python_code(self):
        blk = self.as_one_block
        dump_pycode(self.tears, blk.resids, [blk])
    
    def plot(self):
        g, _, forbidden = self.g_eqs_forbidden
        plot_bipartite(g, forbidden, self.row_perm, self.col_perm)

#-------------------------------------------------------------------------------

class BlockTearingResult:
    
    def __init__(self, alltears, allresids, blocks, bounds):
        self.alltears = alltears
        self.allresids = allresids
        self.blocks = blocks
        self.bounds = bounds
    
    def generate_ampl_code(self):
        dump_ampl(self.alltears, self.allresids, self.blocks, self.bounds)
    
    def generate_python_code(self):
        dump_pycode(self.alltears, self.allresids, self.blocks)
    
    def plot(self):
        plot_ordering(self.blocks)

################################################################################
# The following functions are only interesting on the implementation level.

def setup_model_with_iolets(xml_filename, m):
    read_model(xml_filename, m)
    reconstruct_connections_partially(m)
    # full reconstruction of the connection happens only later: symbolic 
    # processing rewrites the equations and eliminates the parameters
    read_variable_bounds(m)
    do_symbolic_processing(m)
    finish_connections(m)

def read_model(xml_filename, m):
    tree = get_etree_without_namespace(xml_filename)
    variables = extract_all_vars(tree) # list, element type: Variable 
    m.var_names = [ var.name for var in variables ]
    m.all_var_bounds = { v.name : (v.lb, v.ub) for v in variables } 
    m.parameters = parameter_assignments(tree)
    m.equations = create_equations(tree)

def reconstruct_connections_partially(m):
    m.mapping, m.raw_connections = rewrite_equations(m.var_names, m.equations, \
                                                     m.parameters)

def read_variable_bounds(m):
    referenced_vrs = get_referenced_vars(m.equations, m.var_names, m.parameters)
    show_missing_bounds(referenced_vrs, m.all_var_bounds)
    m.bounds = { name : m.all_var_bounds[name] for name in referenced_vrs }

def do_symbolic_processing(m):
    # Side-effect: rewrites the equations and eliminates the parameters too!
    set_symbolic_eliminations(m.equations, m.parameters, m.bounds)
    m.parameters = [ ]

def finish_connections(m):
    m.process_graph = get_process_graph( m.equations, m.raw_connections, \
                                         m.mapping, m.parameters)
    m.raw_connections = None
