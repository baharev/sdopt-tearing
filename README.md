Exact and heuristic methods for tearing
=======================================

Documentation of the demo application is available at 
[sdopt-tearing.readthedocs.org](https://sdopt-tearing.readthedocs.org). 


Some of the algorithms are documented in the following academic papers:

  - [An exact method for the minimum feedback arc set problem](http://reliablecomputing.eu/baharev_minimum_feedback_arc_set.pdf)
  - Exact methods for optimal tearing (to be submitted soon)

**The code is pretty much a work in progress; expect breaking changes.**

Some of the code will be contributed back to 
[NetworkX](http://networkx.github.io/documentation/latest/overview.html)
when it is appropriate. The remaing part of the code will be released 
as a standalone Python package on PyPI.


Installation
------------

The code has been tested with Python 2.7 and 3.5. The `six`, `networkx`, 
and `sympy` packages are necessary; `matplotlib` and `pydot` (or 
`pydot-ng`) are recommended but not required. If you wish to run the 
exact algorithms based on integer programming, you will also need 
[Gurobi](http://www.gurobi.com/). If you do not have Gurobi installed, 
the demo application will detect its absence, and simply skips those 
steps that would require the integer programming solver.


Installing on Windows with PyGraphviz
-------------------------------------

Only Python 2.7 was tested on Dec 13, 2015. 
Installed Python 2.7 with Miniconda-3.9.1-Windows-x86.exe from the
[Miniconda installer archive](https://repo.continuum.io/miniconda/). 
(The latest Miniconda failed to install.)
Then [graphviz](http://www.graphviz.org/Download_windows.php) was installed 
with the `graphviz-2.38.msi` installer. After installation the bin 
directory of graphviz was added to the `PATH`. Next, `pygraphviz` was 
installed from 
[Unofficial Windows Binaries for Python Extension Packages](http://www.lfd.uci.edu/~gohlke/pythonlibs/#pygraphviz)
as `pip install pygraphviz-1.3.1-cp27-none-win32.whl`.
The `matplotlib` and `pydot-ng` packages were installed only *after* 
`pygraphviz`.
