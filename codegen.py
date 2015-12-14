# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from equations import gen_nonnesting_eqs
from expression_tree import rearrange_if_assignment

################################################################################
# AMPL specific code

def dump_ampl(alltears, allresids, blocks, bounds):
    var_names, pyname_map = get_all_names(blocks)
    lines = [ ]
    add_ampl_preamble(lines, var_names, alltears, bounds, pyname_map)
    add_equations(lines, blocks, pyname_map, end_of_stmt=';')
    add_ampl_postamble(lines, alltears, pyname_map)
    print()
    for line in lines:
        print(line)

def add_ampl_preamble(lines, var_names, alltears, bounds, pyname_map):
    # Write variable declarations and bounds, but only for tears
    tear_set = set(alltears)
    NegInf, PosInf = float('-inf'), float('inf')
    for name in var_names:
        if name not in tear_set:
            lines.append('var {};  # {}'.format(pyname_map[name], name))
        else:
            l, u = bounds[name]
            lb = ' >= {}'.format(l) if l != NegInf else ''
            ub = ' <= {}'.format(u) if u != PosInf else ''
            sep = ',' if lb and ub else ''
            lines.append('var {}{}{}{};  # {}'.format(pyname_map[name], 
                                                      lb, sep, ub, name))
            # E.g.: var v18;  # condenser.divider.outlet[2].f[2]
            # or:   var v19 >= 0.0, <= 1.0;  # condenser.divider.zeta
    lines.append('')

def add_ampl_postamble(lines, alltears, pyname_map):
    lines.append('solve;')
    for tear in alltears:
        name = pyname_map[tear]
        lines.append('printf "{} = %f # %s\\n",{},"{}";'.format(name,name,tear))
        # E.g. printf "v228 = %f # %s\n",v228,"condenser.noVaporFlash.L";
        # Then AMPL will print: v228 = 3.115617 # condenser.noVaporFlash.L
        # in the command line. This can then be copied into the Python code.
    lines.append('')

################################################################################
# Python specific code

def dump_pycode(alltears, allresids, blocks):
    _, pyname_map = get_all_names(blocks)
    lines = [ ]
    add_py_preamble(lines, alltears, pyname_map)
    add_equations(lines, blocks, pyname_map)
    add_py_postamble(lines, allresids, pyname_map)
    print()
    for line in lines:
        print(line)

def add_py_preamble(lines, alltears, pyname_map):
    lines.append('#!/usr/bin/env python')
    lines.append('from math import exp')
    if alltears:
        lines.append('# All tear variables')
    for tear in alltears:
        lines.append('{} = 0.1  # {}'.format(pyname_map[tear], tear))
        # Writes for example: v228 = 0.1  # condenser.noVaporFlash.L
        # The magic number 0.1 can then be replaced with the solution from AMPL.
    lines.append('')

def add_py_postamble(lines, allresids, pyname_map):
    residuals = ', '.join('r_%d' % eq.id for eq in allresids)
    lines.append('residuals = [{}]'.format(residuals))
    lines.append('print residuals')

################################################################################
# Code shared by both AMPL and Python

def get_all_names(blocks):
    var_names = list(gen_column_permutation(blocks))
    pyname_map = { name : 'v%d' % i for i, name in enumerate(var_names) }
    return var_names, pyname_map

def gen_column_permutation(blocks):
    for blk in blocks:
        for v in blk.vars:
            yield v
        for y, _, _ in blk.conn_triples:
            yield y

def add_equations(lines, blocks, pyname_map, end_of_stmt=''):
    counter = 0
    for blk in blocks:
        if blk.eqs:
            lines.append('# Block')
        if blk.tears:
            lines.append('# Tears: '+tears_to_pynames(blk.tears, pyname_map))
            # E.g. # Tears: condenser.divider.zeta (v19)
        for eq in blk.eqs:
            # A variable is eliminated (non-residual equation):
            if eq.solved_for:
                eq_str = eq.elims[eq.solved_for]
                # E.g. v33 = v31 + v32 but with Modelica names
            # Or residual equation, either in AMPL or in Python form
            elif end_of_stmt: # AMPL
                eq_str = eq.symbolic_form + ' = 0'
                # E.g. v30-(v33*v36) = 0 but with Modelica names
            else:             # Python
                eq_str = ('r_%d = ' % eq.id) + eq.symbolic_form
                # E.g. r_27 = v30-(v33*v36) but with Modelica names
            # Now replace Modelica names, and move the original eq to comment
            code = replace_names(eq_str, eq.names, pyname_map, end_of_stmt)
            # AMPL requires a name for each equation
            prefix, counter = eq_name(counter, end_of_stmt)
            lines.append(prefix + code)
            # An AMPL and a Python example line:
            # eq_42: v30-(v33*v36) = 0;  # cascade.stages[1].flash...
            # r_27 = v30-(v33*v36)  # r_27 = cascade.stages[1].flash.outlet...
        if blk.conn_triples:
            lines.append('# Connections')
        for y, x, _ in blk.conn_triples:
            # AMPL requires a name for each equation
            prefix, counter = eq_name(counter, end_of_stmt)
            lines.append(prefix + conn_to_code(y, x, pyname_map, end_of_stmt))
            # An AMPL and a Python example line:
            # eq_229: v229 = v225;  # liquidSink.inlet.H = reboiler.outlet[2].H
            # v229 = v225  # liquidSink.inlet.H = reboiler.outlet[2].H
        lines.append('')

def tears_to_pynames(names, pyname_map):
    # E.g. condenser.L (v228), condenser.outlet[1].f[1] (v224)
    return ', '.join('{} ({})'.format(name, pyname_map[name]) for name in names)

def replace_names(eq_str, names, pyname_map, end_of_stmt):
    code = modelica_names_to_ampl(eq_str, names, pyname_map)
    # E.g. v33 = v31 + v32;  # cascade.stages[1].flash.V = cascade.stages[1]...
    return '{}{}  # {}'.format(code, end_of_stmt, eq_str)

def conn_to_code(y, x, pyname_map, end_of_stmt):
    y_py, x_py = pyname_map[y], pyname_map[x]
    # E.g. v229 = v225;  # liquidSink.inlet.H = reboiler.outlet[2].H
    return '{} = {}{}  # {} = {}'.format(y_py, x_py, end_of_stmt, y, x)

def eq_name(counter, end_of_stmt):
    # AMPL requires a name for each equation
    if end_of_stmt:
        prefix = 'eq_%d: ' % counter
        counter += 1
        return prefix, counter
    return '', 0

def modelica_names_to_ampl(symbolic_form, names, to_ampl):
    # To avoid replacing only a prefix, start from the longest name
    assert symbolic_form is not None, names
    for name in sorted(names, key=len, reverse=True):
        symbolic_form = symbolic_form.replace(name, to_ampl[name])
    return symbolic_form

################################################################################

def to_unstructured_ampl(true_names, parameters, equations):
    lines = [ ]
    # Modelica variable name -> generated AMPL name (v0, v1, ...)
    to_ampl = { name : 'v%d' % i for i, name in enumerate(true_names) }
    # Write the variable declarations in AMPL
    for name in true_names:
        lines.append('var %s;  # %s' % (to_ampl[name], name))
    lines.append('\n######################################################\n')
    # Variables that are essentially parameters
    for i, t in enumerate(parameters):
        name, value = t
        lines.append('bind_%d: %s = %s;  # %s'%(i,to_ampl[name], value, name))
    lines.append('\n######################################################\n')
    # Write the equations in AMPL form
    for i, eq in enumerate(gen_nonnesting_eqs(equations)):
        expr = rearrange_if_assignment(eq.expression_tree)
        if not expr: # not an assignment
            expr = eq.symbolic_form + ' = 0'
        ampl_eq = modelica_names_to_ampl(expr, eq.names, to_ampl)
        lines.append('eq_%d: %s;  # %s' % (i, ampl_eq, eq.symbolic_form))
    return lines

