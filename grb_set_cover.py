# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
import six
from gurobipy import LinExpr, GRB, Model, setParam
from benchmarks import gen_benchmark_digraphs
from utils import double_check, get_all_cycles, solve_ilp

def main():
    setParam("LogFile", "/tmp/gurobi.log")
    for g in gen_benchmark_digraphs():
        success, edgelist_per_cycle = get_all_cycles(g, cutoff=10000000)
        if not success:
            continue # too many simple cycles
        solve_cm(g, edgelist_per_cycle)

# Stats only contains counters that will be mutated in-place if passed as the 
# keyword argument stats.

def solve_cm(g, edgelist_per_cycle, stats=None):
    #setParam("LogFile", "/tmp/gurobi.log")
    #setParam("OutputFlag", 0)
    m, vrs = build_ilp(g, edgelist_per_cycle)
    success = solve_ilp(m, stats)
    if not success:
        return None, None
    solution =  [ e for e, y in six.iteritems(vrs) if int(round(y.x)) == 1 ]
    objective = int(round(m.getObjective().getValue()))
    return solution, objective

def build_ilp(g, edgelist_per_cycle):
    m = Model()
    binary = GRB.BINARY
    edges = g.edges(data=True)
    # add all variables first, keep vrs= {edge:gurobi var} for the constraints 
    vrs = { (u,v) : m.addVar(vtype=binary,obj=d['weight']) for u,v,d in edges }
    m.update()
    for edgelist in edgelist_per_cycle:
        # add the set covering constraint of each simple cycle
        y = [ vrs[edge] for edge in edgelist ]
        lhs = LinExpr([1.0]*len(y), y)
        m.addConstr(lhs, GRB.GREATER_EQUAL, 1.0)
    return m, vrs

#-------------------------------------------------------------------------------

def rigorous_mfes(subgraph, cutoff):
    # It is just a convenience wrapper around solve_cm. Returns: 
    # (error msg, elims, obj, ncyc). The error message is empty iff successful.
    success, edges_per_cycle = get_all_cycles(subgraph, cutoff)
    if not success:
        return 'Too many simple cycles in relaxed SCC', None, None, None
    elims, objective = solve_cm(subgraph, edges_per_cycle)
    if elims is None:
        return 'Solver failure, giving up', None, None, None
    assert objective > 0 # in an SCC at least...
    double_check(subgraph, objective, elims)
    return '', elims, objective, len(edges_per_cycle) 

if __name__ == '__main__':
    main()
