
Exact and heuristic methods for tearing
=======================================

Documentation of the demo application `demo.py` is available at 
[sdopt-tearing.readthedocs.org](https://sdopt-tearing.readthedocs.org). 

Most of the implemented algorithms are described in the following 
academic papers:

  - [An exact method for the minimum feedback arc set problem](http://reliablecomputing.eu/baharev_minimum_feedback_arc_set.pdf) (submitted)
  - Exact methods for optimal tearing (to be submitted soon)

**The code is pretty much a work in progress.**

Some of the code will be contributed back to 
[NetworkX](http://networkx.github.io/documentation/latest/overview.html)
wherever it is appropriate. The remaining part of the code will be released 
as a Python package on PyPI.


Installation
------------

The demo application has been tested with Python 2.7 and 3.5. The `six`,
`networkx`, and `sympy` packages are necessary; `matplotlib` and `pydot`
(or `pydot-ng`) are recommended but not required. If you wish to run the 
exact algorithms based on integer programming, you will also need 
[Gurobi](http://www.gurobi.com/). If you do not have Gurobi installed, 
the demo application will detect its absence, and simply skips those 
steps that would require the integer programming solver.


Installing on Windows with PyGraphviz
-------------------------------------

Only Python 2.7 was tested on Dec 13, 2015. Python 2.7 was installed 
with Miniconda-3.9.1-Windows-x86.exe from the
[Miniconda installer archive](https://repo.continuum.io/miniconda/). 
(The latest Miniconda failed to install.) Then 
[graphviz](http://www.graphviz.org/Download_windows.php) was installed 
with the `graphviz-2.38.msi` installer. After installation the bin 
directory of graphviz was added to the `PATH`. Next, `pygraphviz` was 
installed from 
[Unofficial Windows Binaries for Python Extension Packages](http://www.lfd.uci.edu/~gohlke/pythonlibs/#pygraphviz)
as `pip install pygraphviz-1.3.1-cp27-none-win32.whl`.
The `matplotlib` and `pydot-ng` packages were installed only after
`pygraphviz`. The demo application works as expected.


Reproducing the results of the academic papers
----------------------------------------------

The test graphs are given in `benchmarks.py`. The tables of the paper
[An exact method for the minimum feedback arc set problem](http://reliablecomputing.eu/baharev_minimum_feedback_arc_set.pdf) 
can be reproduced as follows. 

 - The `grb_lop.py` module implements the *Integer programming 
 formulation with triangle inequalities* of Section 3.1. The `grb_` 
 prefix refers to Gurobi, the `lop` part to the linear ordering problem 
 since these constraints were developed for solving that problem. 
 - The `grb_set_cover.py` module implements the *Integer programming 
 formulation as minimum set cover* of Section 3.2. 
 - The proposed method of the paper, *An integer programming approach 
 with lazy constraint generation* of Section 4, can be found in 
 `grb_pcm.py` where `pcm` refers to the partial cycle matrix. 
 - The procedure called *Safely removing edges* in Appendix A is 
 implemented in `grb_simplifier.py`.

**Please keep in mind that the code is a work in progress.**
