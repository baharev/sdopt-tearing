# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import division, print_function
from copy import deepcopy
from string import Template
import six
import sympy as sp
from expression_tree import ntype, to_symbolic_form, fold_constants, \
                            get_linear_vars
from equations import gen_unit_eqs
from codegen import modelica_names_to_ampl
from utils import import_code

log = print
#def log(*args, **kwargs): pass

def set_symbolic_eliminations(equations, parameters, bounds):
    params = { name : value for name, value in parameters }    
    for eq in gen_unit_eqs(equations):
        expr_tree = replace_params(eq.expression_tree, params)
        expr_tree = fold_constants(expr_tree)
        linear_vars = get_linear_vars(expr_tree)
        varname_bnds = { name : bounds[name] for name in eq.names 
                                              if name not in params }
        # Special case aliasing and linear combinations?
        solutions = symbolic_sols(expr_tree, varname_bnds, params, linear_vars)
        eq.names = sorted(name for name in eq.names if name not in params)
        eq.expression_tree = expr_tree
        eq.symbolic_form  = to_symbolic_form(expr_tree)
        eq.elims = solutions

def symbolic_sols(expr_tree, varname_bnds, params, linear_vars):
    sympy_tree = to_sympy_expression_tree(expr_tree)
    log()
    log(sympy_tree)
    # Modelica variable name -> generated name (v0, v1, ...)
    to_ampl = { name : 'v%d' % i for i, name in enumerate(varname_bnds) 
                                  if name not in params }
    variables = [ v for v in sympy_tree.atoms(sp.Symbol) 
                     if str(v) not in params ] 
    log('Vars:', variables)
    solutions = { } 
    for v in variables:
        sol = get_solution(sympy_tree, v, to_ampl, varname_bnds, linear_vars)
        if sol:
            solutions[str(v)] = sol
    return solutions

def get_solution(sympy_tree, v, to_ampl, varname_bnds, linear_vars):
    try:
        sol = sp.solve(sympy_tree, v, rational=False)
    except NotImplementedError as nie:
        log('<<<\n', nie, '\n>>>', sep='')
        return
    if len(sol)!=1: # Either no solution or multiple solutions
        return
    # Unique and explicit solution
    solution = str(sol[0])
    log(v, '=', solution)
    py_eq = modelica_names_to_ampl(solution,to_ampl.keys(),to_ampl)
    varname = str(v)
    log(to_ampl[varname],'=',py_eq)
    safe = varname in linear_vars or check_safety(py_eq, to_ampl, varname_bnds)
    log('Is safe?', safe)
    return varname+' = '+solution if safe else None

eval_code = Template('''
try:
    from sympy.mpmath import iv
except ImportError:
    from mpmath import iv

iv.dps = 15

def is_safe():
    $varnames = $varbounds
    try:
        res = $expression
    except:
        return False # e.g. ComplexResult: logarithm of a negative number
    return res in iv.mpf([-10**15, 10**15])
''')

def check_safety(py_eq, to_ampl, varname_bnds):
    names, ivbounds = [ ], [ ]
    bound_template = Template('iv.mpf(($l, $u))')
    NegInf, PosInf = float('-inf'), float('inf')
    for name, bounds in six.iteritems(varname_bnds):
        names.append(to_ampl[name])
        lb = str(bounds[0]) if bounds[0]!=NegInf else "'-inf'"
        ub = str(bounds[1]) if bounds[1]!=PosInf else "'inf'"
        ivbounds.append(bound_template.substitute(l=lb, u=ub))
    py_eq = py_eq.replace('exp', 'iv.exp')
    py_eq = py_eq.replace('log', 'iv.log')    
    code = eval_code.substitute( varnames   = ', '.join(names),
                                 varbounds  = ', '.join(ivbounds),
                                 expression = py_eq )
    #log(code)
    m = import_code(code)
    return m.is_safe()

#-------------------------------------------------------------------------------

def replace_params(orig_tree, params):
    g = deepcopy(orig_tree)
    param_dicts = ( d for n,d in g.nodes_iter(data=True)
                       if d['kind']==ntype.VAR and d['value'] in params )
    for d in param_dicts:
        name = d['value']
        d['value'] = params[name]        
        d['kind'] = ntype.NUM
    return g

NTYPE_TO_OP = { ntype.ADD: sp.Add,
                ntype.DIV: lambda x,y: x/y,
                ntype.EXP: sp.exp,
                ntype.LOG: sp.log,
                ntype.MUL: sp.Mul,
                ntype.NEG: lambda x: -x,
                ntype.POW: sp.Pow,
                ntype.SQR: lambda x: x**2,
                ntype.SUB: lambda x,y: x-y,
              }

# networkx.DiGraph -> SymPy expression tree
# This is a no-win situation: either the expression tree will know about SymPy
# or this module will know about the networkx.DiGraph.
def to_sympy_expression_tree(expr_tree):
    root = len(expr_tree)
    return recurse(expr_tree, root)

# Quasi-duplicate of expression_tree.infix()
def recurse(dag, n):
    # Assumption: the node id order corresponds to the children order because 
    # the expression tree was constructed that way
    args = sorted(dag.edge[n])
    nargs = len(args)
    d = dag.node[n]
    kind = d['kind']
    if   nargs == 0:
        return leaf_node(kind, d['value'])
    elif nargs == 1:
        return unary_operation(kind, recurse(dag, args[0]))
    else:
        assert nargs == 2, d
        return binary_operation(kind, recurse(dag,args[0]),recurse(dag,args[1]))

def leaf_node(kind, value):
    return sp.Symbol(value, real=True) if kind==ntype.VAR else sp.Number(value)

def unary_operation(kind, arg):
    op = NTYPE_TO_OP[kind]
    return op(arg)

def binary_operation(kind, arg1, arg2):
    op = NTYPE_TO_OP[kind]
    return op(arg1, arg2)
