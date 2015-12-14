# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from networkx import DiGraph

class Variable:
    ALGEBRAIC   = 'algebraic'
    INDEP_PARAM = 'indep_param'

    def __init__(self, name, lb, ub, category_str, linear):
        self.name     = name
        self.lb       = lb
        self.ub       = ub        
        self.category = _get_var_category(category_str)
        self.linear   = linear

    def __str__(self):
        return self.name

def _get_var_category(category_str):
    # Replace with a map if grows
    if category_str=='algebraic':
        return Variable.ALGEBRAIC
    if category_str=='independentParameter' or \
       category_str=='independentConstant':
        return Variable.INDEP_PARAM
    raise NotImplementedError(category_str)

#-------------------------------------------------------------------------------
# Not nice that the variables parse XML but couldn't find a better alternative.

def extract_all_vars(tree):
    return [ create_var(var) for var in tree.iter(tag='ScalarVariable') ]

def create_var(element):
    name = element.attrib['name']
    lb, ub = extract_bounds(element)
    category_str = next(element.iter(tag='VariableCategory')).text
    linear = bool( next(element.iter(tag='isLinear')).text )
    return Variable(name, lb, ub, category_str, linear)

def extract_bounds(element):
    e = next(element.iter('Real'), None)
    attrib = e.attrib if e is not None else { }
    lb = float(attrib.get('min', '-inf'))
    ub = float(attrib.get('max',  'inf'))
    return lb, ub

def show_missing_bounds(referenced_vars, all_var_bounds):
    NegInf, PosInf = float('-inf'), float('inf')
    vars_missing_bound = [ ]
    for name in referenced_vars:
        lb, ub = all_var_bounds[name]
        if lb==NegInf or ub==PosInf:
            vars_missing_bound.append(name)
    print()
    print('Missing bounds:')
    for name in vars_missing_bound:
        lb, ub = all_var_bounds[name]
        print(name, ' in [', lb, ', ', ub, ']', sep='')
    print()

#-------------------------------------------------------------------------------

def build_hierarchy(all_names):
    hierarchy = DiGraph()
    root = ' ' # root sentinel, could be '.' as well, but '' breaks pygraphviz
    for name in sorted(all_names): # Sorting seems necessary
        for prefix, node_name in walk_hierarchy(name, root):
            hierarchy.add_edge(prefix, node_name)
    return hierarchy

def walk_hierarchy(name, root):
    # a.b.c.d -> a.b.c, a.b.c.d -> a.b, a.b.c -> a, a.b -> root, a
    prefix = get_prefix(name)
    while prefix:
        yield prefix, name
        name = prefix
        prefix = get_prefix(name)
    yield root, name

def get_prefix(name):
    return name.rpartition('.')[0]
