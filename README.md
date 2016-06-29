
Exact and heuristic methods for tearing
=======================================

Many of the implemented algorithms are described in the following 
academic papers:

  - [An exact method for the minimum feedback arc set problem](http://reliablecomputing.eu/baharev_minimum_feedback_arc_set.pdf) (submitted)
  - [Tearing systems of nonlinear equations I. A survey.](http://reliablecomputing.eu/baharev_tearing_survey.pdf) (submitted)
  - [Tearing systems of nonlinear equations II. A practical exact algorithm](http://reliablecomputing.eu/baharev_tearing_exact_algorithm.pdf) (submitted)

See also [Reproducing the results of the academic papers](#reproducing-the-results-of-the-minimum-feedback-arc-set-paper)
below.


The code is a work in progress
------------------------------

Some of the code will be contributed back to 
[NetworkX](https://networkx.readthedocs.io)
wherever it is appropriate. The remaining part of the code will be released 
as a Python package on PyPI. In the meantime, the `rpc_api.py` is a good place 
to start looking. (`rpc` stands for remote procedure call; it can be called from
Java or C++ through the `json_io.py`). The API in `rpc_api.py` takes a sparse
matrix in coordinate format and returns the row and column permutation vectors.
Documentation of the demo application `demo.py` is available at 
[sdopt-tearing.readthedocs.io](https://sdopt-tearing.readthedocs.io).

While reading the code, **please keep in mind that the code is a work in 
progress.**


Reproducing the results of the minimum feedback arc set paper
-------------------------------------------------------------

The results of the paper
[An exact method for the minimum feedback arc set problem](http://reliablecomputing.eu/baharev_minimum_feedback_arc_set.pdf) 
can be reproduced as follows. 

**Algorithms**

 - The `grb_lop.py` module implements the *Integer programming 
 formulation with triangle inequalities* of Section 3.1. The `grb_` 
 prefix refers to Gurobi, the `lop` part to the linear ordering problem 
 since these constraints were developed for solving that problem. 
 - The `grb_set_cover.py` module implements the *Integer programming 
 formulation as minimum set cover* of Section 3.2. 
 - The `grb_lazy.py` provides the implementation of the proposed method, 
 *An integer programming approach with lazy constraint generation* of Section 4.
 - The procedure called *Safely removing edges* in Appendix A is 
 implemented in `grb_simplifier.py`.

**Input graphs**
 
The test graphs are given in plain text files. The `*.edges` file 
gives the edge list of the graph; the `*.mfes` file gives the minimum feedback
edge set. The `*.cycles` file gives those cycles that are sufficient to include 
in the cycle matrix in order to prove the optimality of the minimum feedback 
edge set, and the `*.lp` encodes the corresponding integer linear program (a 
minimum set cover problem) in CPLEX LP file format. The `*.mst` file gives the 
optimal solution vector of this integer program, and can be used as a starting 
point as well. 

The cycles in the `*.cycles` file are given one per line, and each line gives 
the edge list of one cycle. For example the line  
`1 6 6 8 8 1` encodes the 
cycle of the edges `(1, 6)`, `(6, 8)`, `(8, 1)`.

Since only the median execution time is necessary to create the median execution 
time vs. `n` plots, it was possible to give up on certain long running 
computations. These graphs are given in the `*_aborted.zip` files. They do not 
contain any `*.mfes`, `*.cycles`, `*.lp`, or `*.mst` since these problems were 
not solved to optimality. The only claim is that computing the minimum feedback 
edge set of these graphs with the proposed method takes longer than the median 
execution time for the corresponding group of graphs.

The graphs used for benchmarking are given in the following files.
 
 - Section 5.1. The easy test graphs for cross-checking correctness are given in
 the Python module `benchmarks.py`. 
 - Section 5.2. The sparse random graphs are given in 
 `benchmark_mfes/erdos_renyi.zip` and 
 `benchmark_mfes/erdos_renyi_aborted.zip`. The seeds `61` and `78` were skipped 
 since they yielded random graphs with more than one strongly connected 
 components for some of the `n`s.
 - Section 5.3. The random tournaments for testing the worst-case behavior are 
 given in `benchmark_mfes/tournament.zip` and 
 `benchmark_mfes/tournament_aborted.zip`.
 - Section 5.4. The challenging sparse graphs are given in 
 `benchmark_mfes/de_Bruijn.zip` and in `benchmark_mfes/Imase_Itoh.zip`. The 
 self-loops, if any, have been removed.
 
 
Reproducing the results of the paper on optimal tearing
-------------------------------------------------------

The algorithms are documented in the academic paper
[Tearing systems of nonlinear equations -- II. A practical exact algorithm](http://reliablecomputing.eu/baharev_tearing_exact_algorithm.pdf).

 - The `heap_md.py` module implements a greedy heuristic for ordering to lower
 Hessenberg form.
 - The `ilp_tear.py` module implements optimal tearing with integer 
 programming.
 - The `bb4_tear.py` module provides the custom branch and bound 
 algorithm for optimal tearing.
 - The names of selected test problems from the 
 [COCONUT Benchmark](http://www.mat.univie.ac.at/~neum/glopt/coconut/Benchmark/Benchmark.html) 
 are under `data/benchmark_tearing/`. The results of the 12 runs are plotted in the 
 following PDF files (but only for those problems that have at most 
 500 non-zero entries):  
 [Ordering the Jacobian of the equality constraints](http://reliablecomputing.eu/constraint_jacobian_nz_500.pdf) (33.8 MB) and  
 [Ordering the Jacobian of the first-order optimality conditions](http://reliablecomputing.eu/first_order_opt_cond_nz_500.pdf) (49.3 MB).

The tests for checking correctness are in the `test_<module name>.py` 
modules. Cross-checking the ILP-based and the custom branch and bound 
algorithm is implemented in `test_bb_ilp.py`. Running these tests 
require [Hypothesis](https://hypothesis.readthedocs.io), the 
[QuickCheck](https://en.wikipedia.org/wiki/QuickCheck) for Python.
(I am a big fan of property-based testing.)


Dependencies
------------

The demo application has been tested with Python 2.7 and 3.5. The `six`,
`networkx`, `sympy`, and `namedlist` packages are necessary; 
`matplotlib`, `pydot` (or `pydot-ng`), and `pygraphviz` are 
recommended but not required. However, if matplotlib is installed, then
`pydot` (or `pydot-ng`) and `pygraphviz` are also required to run the 
`demo.py` application.

Computing the bipartite matching can easily become the bottleneck. In 
that case, [MC21](http://www.hsl.rl.ac.uk/catalogue/mc21.html) from 
the Harwell Subroutine Library (HSL) is recommended. The wrapper code 
can be found under `data/mc21/` with compilation instructions.

If you wish to run the exact algorithms based on integer programming, 
you will also need [Gurobi](http://www.gurobi.com/). If you do not have 
Gurobi installed, the `demo.py` application will detect its absence, and 
simply skips those steps that would require the integer programming 
solver. The algorithms that use Gurobi have only been tested with Python 
2.7.


Installing on Windows with PyGraphviz
-------------------------------------

Only Python 2.7 was tested on Mar 08, 2016. Python 2.7 and the 
[dependencies](#dependencies) were installed with 
[Miniconda](http://conda.pydata.org/miniconda.html) 
(released on Dec 17, 2015). Then 
[graphviz](http://www.graphviz.org/Download_windows.php) was installed 
with the `graphviz-2.38.msi` installer. After installation the bin 
directory of graphviz was added to the `PATH`. Next, `pygraphviz` was 
installed from 
[Unofficial Windows Binaries for Python Extension Packages](http://www.lfd.uci.edu/~gohlke/pythonlibs/#pygraphviz)
as `pip install pygraphviz-1.3.1-cp27-none-win32.whl`.
The `matplotlib` and `pydot-ng` packages were installed only after
`pygraphviz`. The demo application works as expected.

