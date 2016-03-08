#!/usr/bin/env python
#
# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from importlib import import_module
from os import listdir
from os.path import splitext

def run(): # MUST not be called main() or we will get infinite recursion!!!
    #
    # Just a hackish way of running everything in this directory that can be 
    # run.
    # For vanilla installations:
    import sys
    sys.path.append('/home/ali/gurobi/build/lib.linux-x86_64-2.7')
    ignored = { 'dm_to_pdf.py', 'fmux_creator.py', 'json_io.py' }
    modules = sorted( f for f in listdir('.') if f.endswith('.py') \
                                              and f not in ignored )
    check_copyright_notice(modules)
    for module in modules:
        print('About to import:', module)
        try:
            m = import_module(splitext(module)[0])
        except ImportError as import_error:
            msg = str(import_error)
            if msg=="No module named 'gurobipy'" or \
               msg=='No module named gurobipy': # The former Py3, the latter Py2
                continue
            else:
                raise import_error
        if hasattr(m, 'main'):
            m.main()
        if hasattr(m, 'run_tests'):
            m.run_tests()

def check_copyright_notice(modules):
    notice_missing = [m for m in modules if notice_is_missing(m)]
    if notice_missing:
        print('The following files lack the copyright notice:')
        print(notice_missing)
        quit()

def notice_is_missing(module):
    with open(module, 'r') as f:
        for line in f:
            if '# Copyright (C) ' in line:
                return False
            if 'import ' in line:
                return True
        return True

if __name__ == '__main__':
    run()
