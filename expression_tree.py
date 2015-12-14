# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from copy import deepcopy
import math
import six
from networkx import convert_node_labels_to_integers
from py3compat import irange
from utils import pairwise

# The expression_tree is a rooted, ordered binary tree. Unfortunately, the 
# networkx.DiGraph uses hashing, so the children order has to be restored by 
# sorting the children by their node id. A post-order traversal gives the node
# ids in sorted order, from 1 to n. The root's id is n; i.e. root = len(tree).
# This is guaranteed by construction. The leaf nodes are either numbers or 
# variables, both are stored as string.

class ntype:
    ADD = 'add' # sum in SDOPT    
    DIV = 'div'
    EXP = 'exp'
    LOG = 'log'
    MUL = 'mul'
    NEG = 'neg'
    NUM = 'num'
    POW = 'pow'
    SQR = 'sqr'
    SUB = 'sub' # Not in SDOPT
    VAR = 'var'

#===============================================================================

def add_node(dag, kind):
    my_id = len(dag) + 1
    dag.add_node(my_id, kind=kind)
    return my_id

def add_binary(dag, kind, e):
    itr_children  = iter(e)
    id_left  = add_recursively(dag, next(itr_children))
    id_right = add_recursively(dag, next(itr_children))
    my_id = add_node(dag, kind)
    dag.add_edge(my_id, id_left)
    dag.add_edge(my_id, id_right)
    return my_id

def add_unary(dag, kind, e):
    child = next(iter(e))
    id_child = add_recursively(dag, child)
    my_id = add_node(dag, kind)
    dag.add_edge(my_id, id_child)
    return my_id

def add_leaf(dag, kind, e):
    my_id = add_node(dag, kind)
    if   kind == ntype.NUM:
        value = e.text 
    elif kind == ntype.VAR:
        # TODO Seems like a bug in the JModelica compiler
        # lazy import due to cyclic import
        from equations import get_full_name
        value = get_full_name(e)
    else:
        raise AssertionError(kind)
    dag.node[my_id]['value'] = value
    return my_id

NODE_BUILDER = { 'Add' : (add_binary, ntype.ADD),
                 'Sub' : (add_binary, ntype.SUB),
                 'Mul' : (add_binary, ntype.MUL),
                 'Div' : (add_binary, ntype.DIV),
                 'Pow' : (add_binary, ntype.POW),
                 'Neg' : (add_unary,  ntype.NEG),
                 'Exp' : (add_unary,  ntype.EXP),
                 'Log' : (add_unary,  ntype.LOG),
                 'Identifier' :     (add_leaf, ntype.VAR),
                 'IntegerLiteral' : (add_leaf, ntype.NUM),
                 'RealLiteral' :    (add_leaf, ntype.NUM),        
               }

def add_recursively(dag, e):
    adder, kind = NODE_BUILDER[e.tag]
    return adder(dag, kind, e)

#===============================================================================

def gen_var_node_dicts(expression_tree):
    g = expression_tree
    return (d for _, d in g.nodes_iter(data=True) if d['kind']==ntype.VAR)

def get_varnames(g):
    var_name_gen = ( d['value'] for _, d in g.nodes_iter(data=True) 
                                 if d['kind']==ntype.VAR )
    return sorted(set(var_name_gen)) 

def to_symbolic_form(dag):
    root = len(dag)
    symb = infix(dag, root)
    symb = cleanup_a_bit(symb)
    #print(symb)
    assert symb is not None
    return symb

def infix(dag, n):
    # Assumption: the node id order corresponds to the children order because 
    # the expression tree was constructed that way
    args = sorted(dag.edge[n])
    nargs = len(args)
    d = dag.node[n]
    kind = d['kind']
    if   nargs == 0:
        return leaf_printer(kind, d['value'])
    elif nargs == 1:
        return unary_printer(kind, infix(dag, args[0]))
    else:
        assert nargs == 2, d
        return binary_printer(kind, infix(dag,args[0]), infix(dag,args[1]))
    
def leaf_printer(kind, val):
    # Only negative numbers need to be parenthesized
    return '('+val+')' if kind==ntype.NUM and float(val) < 0.0 else val

NTYPE_OP_TO_AMPL = { ntype.ADD : '+',
                     ntype.DIV : '/',
                     ntype.EXP : 'exp',
                     ntype.LOG : 'log',
                     ntype.MUL : '*',
                     ntype.NEG : '-',
                     ntype.POW : 'pow',
                     ntype.SUB : '-',
                   }

def unary_printer(kind, arg):
    op = NTYPE_OP_TO_AMPL[kind]
    if kind == ntype.NEG:
        return '({}({}))'.format(op, arg)
    return '{}({})'.format(op, arg)

def binary_printer(kind, arg1, arg2):
    op = NTYPE_OP_TO_AMPL[kind]
    if kind == ntype.POW:
        return '{}({},{})'.format(op, arg1, arg2)
    return '({}{}{})'.format(arg1, op, arg2)

def cleanup_a_bit(symb):
    assert symb[0] == '(' and symb[-1] == ')', symb # Assumes last node was Sub
    symb = symb[1:-1]
    end = len(symb)
    if   symb.endswith('-0'):
        end = -2
    elif symb.endswith('-0.0'):
        end = -4
    return symb[:end]

#===============================================================================

def defines_var_alias(expression_tree):
    # v1 - v2 = 0 
    if len(expression_tree)==3:
        arg1_t = expression_tree.node[1]['kind']
        arg2_t = expression_tree.node[2]['kind']
        op_t   = expression_tree.node[3]['kind']
        if arg1_t==ntype.VAR and arg2_t==ntype.VAR and op_t==ntype.SUB:
            return True
    return False

def rearrange_if_assignment(expression_tree):
    # v42 - expression = 0  OR  expression - v42 = 0  ->  v42 = expression
    root_id   = len(expression_tree)
    root_node = expression_tree.node[root_id]
    children  = sorted(expression_tree.edge[root_id])
    if len(children) != 2 or root_node['kind'] != ntype.SUB:
        return
    left_id, right_id = children
    left_child  = expression_tree.node[left_id]
    right_child = expression_tree.node[right_id]  
    if left_child['kind'] != ntype.VAR and right_child['kind'] != ntype.VAR:
        return
    if left_child['kind'] != ntype.VAR:
        left_id, right_id = right_id, left_id
        left_child, right_child = right_child, left_child
    rhs = infix(expression_tree, right_id)
    if rhs[0] == '(' and rhs[-1] == ')':
        rhs = rhs[1:-1]
    return left_child['value'] + ' = ' + rhs

#===============================================================================

def fold_constants(expr_tree):
    # Side effect: converts all ntype.NUM to float and then back to str.
    tree  = deepcopy(expr_tree)
    nodes_in_postorder = irange(1,len(tree)+1)
    for n in nodes_in_postorder:
        fold_if_possible(tree, n)
    # If nothing happened, just return the input
    if len(tree)==len(expr_tree):
        return expr_tree
    # Undo str -> float conversion on numbers
    for n, d in tree.nodes_iter(data=True):
        if d['kind'] == ntype.NUM:
            d['value'] = str(d['value'])
    # Relabel nodes to restore invariant: node ids are 1..n in post-order
    return convert_node_labels_to_integers(tree, 1, 'sorted')

def fold_if_possible(tree, n):
    children = sorted(tree.edge[n])
    if not children:
        return
    all_number = all(tree.node[child]['kind']==ntype.NUM for child in children)
    if not all_number:
        return
    d = tree.node[n]
    op =  NTYPE_OP_TO_PY[d['kind']]
    args = [ float(tree.node[child]['value']) for child in children ]
    d['value'] = op(*args)
    d['kind']  = ntype.NUM
    tree.remove_nodes_from(children)
    
NTYPE_OP_TO_PY = { ntype.ADD : lambda x, y: x+y,
                   ntype.DIV : lambda x, y: x/y,
                   ntype.EXP : lambda x: math.exp(x),
                   ntype.LOG : lambda x: math.log(x),
                   ntype.MUL : lambda x, y: x*y,
                   ntype.NEG : lambda x: -x,
                   ntype.POW : lambda x, y: math.pow(x, y),
                   ntype.SUB : lambda x, y: x-y,
                 }

def get_linear_vars(tree):
    # Assumes that constant folding has already been performed
    linear_vars = set()
    for name, occurances in six.iteritems(collect_varnode_occurances(tree)):
        if all( path_to_root_is_linear(tree,n) for n in occurances ):
            linear_vars.add(name)
    return linear_vars

def collect_varnode_occurances(tree):
    # Returns { Modelica name : [ node IDs] }
    occurances = { }
    for n, d in tree.nodes_iter(data=True):
        if d['kind'] == ntype.VAR:
            occurances.setdefault(d['value'], [ ]).append(n)
    return occurances    

def path_to_root_is_linear(tree, n):
    # Is n an input of any nonlinear operation on the path to root?
    path = path_to_root(tree, n)
    edges = pairwise(path)
    return all(is_parent_linear(tree, child, parent) for child, parent in edges)

def path_to_root(tree, n):
    path = [ n ]
    in_edges = tree.pred[n]
    while in_edges:
        (parent,) = in_edges
        path.append(parent)
        in_edges = tree.pred[parent]
    assert path[-1] == len(tree)
    return path

def is_parent_linear(tree, child, parent):
    kind = tree.node[parent]['kind'] 
    if kind==ntype.ADD or kind==ntype.SUB:
        return True
    if kind==ntype.MUL:
        sibling = get_sibling(tree, child, parent)
        return tree.node[sibling]['kind'] == ntype.NUM
    if kind==ntype.DIV:
        left, right = sorted(tree.edge[parent])
        return child==left and tree.node[right]['kind'] == ntype.NUM
    return False

def get_sibling(tree, child, parent):
    left, right = sorted(tree.edge[parent])
    return left if left!=child else right
