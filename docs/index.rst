
=======================================
Exact and heuristic methods for tearing
=======================================

A *prototype* implementation of various tearing algorithms is presented by 
walking through a demo example. The technical details will be published in an
academic paper. The `source code is available on GitHub 
<https://github.com/baharev/sdopt-tearing>`_ under the 3-clause BSD license.

.. _spiked-form:

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

The code has been tested with Python 2.7 and 3.4. 
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

The equations are represented as binary expression trees in the flat Modelica
model. The picture below shows the expression tree for ::

    y[1] = alpha*x[1]/(1.0+(alpha-1.0)*x[1])

.. image:: ./pics/ExprTree.png
   :alt: Expression Tree in SymPy.
   :align: center
   :scale: 75%

The `expression tree <http://docs.sympy.org/latest/tutorial/manipulation.html>`_ of 
the equations are symbolically manipulated with `SymPy <http://www.sympy.org/>`_
to determine which variables can be explicitly and safely eliminated from which 
equations. An example for unsafe elimination is rearreanging `x*y=1` to `y=1/x`
if `x` may potentially take on the value `0`.


4. Optimal tearing
------------------

If Gurobi is installed, the Jacobian is ordered optimally, with an exact method.
For same system that was shown on the top (the very 
:ref:`first picture <spiked-form>`), we get an
ordering that yields a `4 x 4` reduced system (compared to `5 x 5` obtained with
the heuristic method). The integer programming approach does not need the block
structure (those were the blue lines in the first picture).

.. image:: ./pics/OptimalTearing.png
   :alt: Optimal tearing, obtained with integer programming.
   :align: center
   :scale: 75%


5. A tearing heuristic exploiting the natural block structure
-------------------------------------------------------------

Technical systems can be partitioned into smaller blocks along the equipment 
boundaries in a fairly natural way. We call this partitioning the natural block 
structure. The implemented tearing heuristic first orders the blocks, then the 
equations within each block. This is how the :ref:`first picture <spiked-form>`
with the spiked form was obtained. It is also repeated here:

.. image:: ./pics/SpikedForm.png
   :alt: Tearing with the block structure.
   :align: center
   :scale: 50%

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
