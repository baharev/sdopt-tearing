# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from networkx import Graph
from py3compat import izip_longest

__all__ = [ 'gen_testproblems' ]

def gen_testproblems():
    # yields: (undirected bipartite graph, equations, forbidden edges)
    testproblems = sorted(p for p in dir(TestProblems) if not p.startswith('_'))
    obj = TestProblems()
    print('===============================================================')    
    for problem_name in testproblems:
        print("Creating test problem '%s'" % problem_name)
        edge_list = getattr(obj, problem_name)()
        yield to_graph_eqs_forbidden(edge_list) 
        print('===============================================================')

def to_graph_eqs_forbidden(edge_list, equation_name_length=3):
    W = 'weight'
    edges = [ (e,v,{W:int(w)}) for e,v,w in grouper(edge_list.split(),3) ]
    forbidden = {(e,v) for e,v,d in edges if d[W] < 0}
    g = Graph(edges)
    eqs = [n for n in g if len(n)==equation_name_length]
    return g, eqs, forbidden

def grouper(iterable, n, fillvalue=None):
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

class TestProblems:

    def bug(self):
        return '''
                eq0  a  1    eq0  b  1    eq0  c  1
                eq1  a  1    eq1  c -1
                eq2  b  1    eq2  c -1
        '''

    def condenser_divider(self):
        return '''
                eq1  e  1    eq1  f  1
                eq2  a -1    eq2  b  1
                eq3  a -1    eq3  c  1
                eq4  a -1    eq4  d  1
                eq5  b  1    eq5  e  1
                eq6  c  1    eq6  f  1
                eq7  d  1    eq7  g  1
        '''

    def condenser_flash(self):
        return '''
                eq0  a  1    eq0  b  1    eq0  c  1
                eq1  a  1    eq1  c -1    eq1  d -1
                eq2  b  1    eq2  c -1    eq2  e -1
                eq3  d  1    eq3  f  1
                eq4  c -1    eq4  f -1    eq4  g  1
                eq5  a  1    eq5  h  1
                eq6  b  1    eq6  i  1
                eq7  f  1    eq7  j  1
        '''
    
    def reboiler(self):
        return '''
                e01  a  1
                e02  b  1
                                
                e03  a  1    e03  c  1
                e04  b  1    e04  d  1
                e05  c  1    e05  e  1
                e06  e  1    e06  f  1
                e07  e -1    e07  g  1
                e08  e -1    e08  h  1
                                
                e09  a  1    e09  i  1
                e10  b  1    e10  j  1
                                
                e11  e -1    e11  i  1    e11  k -1
                                 
                e12  g  1    e12  l  1
                e13  h  1    e13  k  1    e13  m  1
                e14  l  1    e14  m  1    e14  n  1
                                
                e15  a  1    e15  b  1
                e16  f  1    e16  j  1    e16  k  1
                                
                e17  m  1    e17  o  1
                e18  i  1    e18  p  1
                e19  j  1    e19  q  1
        '''

    def two_SCCs(self):
        return '''
                eq0  x1  1    eq0  x3  1    eq0  x4  1   
                eq1  x1  1    eq1  x2  1
                eq2  x3  1    eq2  x4  1
                eq3  x1  1    eq3  x2  1
        '''

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    from tearing import run_tests
    run_tests()
