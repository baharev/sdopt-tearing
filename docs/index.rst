
=======================================
Exact and heuristic methods for tearing
=======================================

A *prototype* implementation of various tearing algorithms is presented by 
walking through a demo example. The `source code is available on GitHub 
<https://github.com/baharev/sdopt-tearing>`_ under the 3-clause BSD license.

The picture below shows a sparse matrix ordered to the so-called spiked form.
The original matrix is of size `76 x 76`; this can be reduced to a `5 x 5` 
matrix (where `5` is the number of spike columns). The blue lines correspond to 
the equipment boundaries in the technical system; the red squares are above the 
diagonal; the grey squares are "forbidden" variables (no explicit elimination 
possible).

.. image:: ./pics/SpikedForm.png
   :alt: A sparse matrix ordered to the so-called spiked form.
   :align: center
   :scale: 50%

--------------------------------------------------------------------------------

Requirements
============

The :mod:`six`, :mod:`networkx`, and :mod:`sympy` packages are necessary; 
:mod:`matplotlib` is recommended but not required. If you wish to run the exact 
algorithms based on integer programming, you will also need 
`Gurobi <http://www.gurobi.com/>`_. If you do not have Gurobi installed, the 
demo application will detect it and simply skip those steps that would require 
the integer programming solver.

--------------------------------------------------------------------------------

Demo application
================


1. Input
--------

The Modelica model :file:`data/demo.mo` has already been 
flattened with the JModelica compiler (by calling :func:`compile_fmux`; the 
relevant modules are :mod:`flatten` and :mod:`fmux_creator`). The demo 
application takes this flattened model as input.


2. Recovering the process graph
-------------------------------

A directed graph is recovered from the flattened model: The equipments 
correspond to the vertices of the process graph, the edges correspond to the 
material flows.

.. image:: ./pics/Cascade.png
   :alt: Digraph representation of a distillation column.
   :align: center
   :scale: 75%

Currently, recovering the directed edges is only possible if the input 
connectors of the equipments are called `inlet`, and their output connectors are
called `outlet`. There is an ongoing discussion with the JModelica developers on 
reconstructing this information in a generic way, without assuming any naming 
convention.


3. Symbolic manipulation of the equations
-----------------------------------------

The `expression tree <http://docs.sympy.org/latest/tutorial/manipulation.html>`_ of 
the equations are symbolically manipulated to determine which variables can be 
explicitly and safely eliminated from which equations.

.. image:: ./pics/ExprTree.png
   :alt: Expression Tree in SymPy.
   :align: center
   :scale: 75%

4. Optimal tearing
------------------

If Gurobi is installed, the system of equations is ordered 
optimally, with an exact method. 

.. image:: ./pics/OptimalTearing.png
   :alt: Optimal tearing, obtained with integer programming.
   :align: center


5. A tearing heuristic exploiting the natural block structure
-------------------------------------------------------------

Technical systems can be partitioned into smaller blocks along the equipment 
boundaries in a fairly natural way. We call this partitioning the natural block 
structure. The implemented tearing heuristic first orders the blocks, then the 
equations within each block.

.. image:: ./pics/TearingWithBlocks.png
   :alt: Tearing with the block structure.
   :align: center

The above picture shows the so-called spiked form.


6. Code generation after tearing
--------------------------------

The `AMPL <http://en.wikipedia.org/wiki/AMPL>`_
code is written out in such a way that the variables can be eliminated as 
desired. The reduced system will have as many variables and equations as the 
number of spike columns in the spiked form.
Executable Python code is also emitted: It only serves for cross-checking 
correctness. For efficient computations, templated C++ code will be emitted in
the future the Jacobian will be obtained with reverse mode automatic 
differentiation.
