# Only called by flatten.py as a subprocess.
#
# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from pymodelica import compile_fmux
import sys

opt_dict = { "automatic_tearing" : False, 
             "eliminate_alias_variables": False, 
             "variability_propagation" : True, 
             "generate_mof_files" : True 
           }

class_name = sys.argv[2]
file_name  = sys.argv[1]
fmux_file = compile_fmux(class_name, file_name, compiler_options=opt_dict)
