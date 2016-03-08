# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>

#===============================================================================
# This is the demo application presented at http://sdopt-tearing.readthedocs.org

from __future__ import print_function
from sys import stderr
from benchmarks import create_testproblem
from mfes import run_mfes_heuristic  # MFES: Minimum Feedback Edge Set
from model_fmux import ModelWithInletsAndOutlets, ModelWithNoReconstruction
from plot import plot_tearing_result
from utils import has_gurobi, has_matplotlib, info


def main():
    
    # Read in the flattened Modelica model (FMUX). Alias variable elimination 
    # and tearing was disabled during flattening. Assumption: the connectors are
    # called inlet or outlet.
    model = ModelWithInletsAndOutlets('demo.xml.gz')
    
    # The just created model contains the process graph: a directed graph in 
    # which the devices correspond to the vertices, the edges to the material
    # flows. Furthermore, each equation has been examined and the safe and 
    # explicit eliminations have been recorded.
    
    if has_gurobi():
        # Run the exact method, based on integer linear programming (ILP). The 
        # optimal ordering is returned, with the least spike columns.
        result = model.run_ilp_tearing() # ILP: Integer Linear Programming
        show(result)
    else:
        stderr.write('The exact tearing algorithm requires Gurobi.\n')
    
    # Hierarchical tearing: first order the blocks, then the equations within 
    # each block. The blocks correspond to the devices. 
    result = model.run_hierarchical_heuristic_tearing()
    show(result)
    
    # A greedy tearing heuristic has been implemented, inspired by algorithm 
    # (2.3) of Fletcher and Hall, see http://dx.doi.org/10.1007/BF02025533
    # The heuristic resembles the minimum degree algorithm, hence the name. 
    result = model.run_mindegree()
    show(result)
    
    # The same as above but with a lookahead step when facing ties in the greedy
    # choice. The result happens to be optimal for this example.
    result = model.run_mindegree_with_lookahead()
    show(result)
    
    # Tearing in the chemical engineering sense: It is equivalent to the minimum
    # feedback edge set (MFES) problem, also known as the maximum acyclic 
    # subgraph problem. A plot will pop up if matplotlib is installed; the 3 red
    # edges form a minimum feedback edge set of the directed graph. If Gurobi is
    # installed, it will be verified that the tear set is a minimum feedback
    # edge set. 
    tearing_as_in_chemical_engineering()
    
    # Even if the inlet and outlet naming convention for the connectors is not 
    # followed, the flattened Modelica model can be still translated to an AMPL 
    # model. However, the process graph (where the devices correspond to the 
    # vertices, the edges to the material flows) cannot be reconstructed.
    model = ModelWithNoReconstruction('demo.xml.gz')
    model.translate_to_ampl()


def show(result):
    result.plot()
    result.generate_ampl_code()
    result.generate_python_code()


def tearing_as_in_chemical_engineering():
    
    # The test problems in the benchmarks module are taken from: 
    # http://dx.doi.org/10.4173/mic.1983.3.2   (available free of charge)
    directed_graph = create_testproblem('Problem 7 (opt=3)')
    info(directed_graph)
    
    # The greedy heuristic runs with lookahead:
    cost, torn_edges = run_mfes_heuristic(directed_graph)
    
    # The heuristic is generic: The graph can have edge weights, so the cost 
    # does not necessarily equal len(torn_edges) in the general case; however,
    # for Problem 7 it does, as this test problem has unit weights.
    print('It is sufficient to remove', cost, 'edges to make the graph acyclic')
    
    if has_matplotlib():
        plot_tearing_result(directed_graph, torn_edges)
    else:
        stderr.write('Plotting requires matplotlib.\n')
    
    if has_gurobi():
        minimum_cost = run_exact_mfes(directed_graph)
        print('It is necessary to remove', minimum_cost, 'edges.')


def run_exact_mfes(directed_graph):
    from grb_pcm import solve_problem
    _, cost = solve_problem(directed_graph)  # we ignore the list of torn edges
    return cost


if __name__ == '__main__':
    main()
