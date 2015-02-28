.. sdopt-tearing documentation master file, created by
   sphinx-quickstart on Sat Feb 28 23:04:04 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


=======================================
Exact and heuristic methods for tearing
=======================================

An example is presented here, showing the capabilities of the research prototype
of various tearing algorithms. The technical details will be published in an
academic paper. The `source code is available on GitHub 
<https://github.com/baharev/sdopt-tearing>`_ under the 3-clause BSD license.

.. _spiked-form:

The picture below shows a sparse matrix ordered to the so-called 
**spiked form**.
The matrix is of size 76x76; this can be reduced to a 5x5 
matrix, where 5 is determined by the number of spike columns -- columns with red
entries. The blue lines correspond to the equipment boundaries in the technical 
system; the entries above the diagonal are red; the gray squares are "forbidden"
variables (no explicit elimination possible).

.. image:: ./pics/SpikedForm.png
   :alt: A sparse matrix ordered to the so-called spiked form.
   :align: center
   :scale: 50%

This plot has been automatically generated with the prototype implementation.

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

Steps of the demo application
=============================


1. Input: flattened Modelica model
----------------------------------

The Modelica model :file:`data/demo.mo` has already been 
flattened with the `JModelica <http://www.jmodelica.org/>`_ compiler 
(by calling :func:`compile_fmux`; the relevant modules are :mod:`flatten` and 
:mod:`fmux_creator`). The demo application takes this flattened model as input.


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
connectors of the equipments are called ``inlet``, and their output connectors 
are called ``outlet``. There is an ongoing discussion with the JModelica 
developers on reconstructing this information in a generic way, without assuming
any naming convention.


3. Symbolic manipulation of the equations
-----------------------------------------

The equations are represented as binary expression trees in the flat Modelica
model. The picture below shows the expression tree for ::

    y[1] = alpha*x[1]/(1.0+(alpha-1.0)*x[1])

.. image:: ./pics/ExprTree2.png
   :alt: Expression Tree in SymPy.
   :align: center
   :scale: 75%

The `expression tree <http://docs.sympy.org/latest/tutorial/manipulation.html>`_ of 
the equations are symbolically manipulated with `SymPy <http://www.sympy.org/>`_
to determine which variables can be explicitly and safely eliminated from which 
equations. An example for unsafe elimination is rearranging `x*y=1` to `y=1/x`
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
   :scale: 50%


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
number of spike columns in the spiked form. An example code snippet is shown
below. ::

    # Unit
    # Tears: condenser.divider.zeta (v19)
    eq_14: v14 = v12*v19;  # condenser.divider.outlet[1].f[1] = condenser.divider.inlet[1].f[1]*condenser.divider.zeta
    eq_15: v15 = v13*v19;  # condenser.divider.outlet[1].f[2] = condenser.divider.inlet[1].f[2]*condenser.divider.zeta
    eq_16: v16 = v11*v19;  # condenser.divider.outlet[1].H = condenser.divider.inlet[1].H*condenser.divider.zeta
    eq_17: v17 = v12 - v14;  # condenser.divider.outlet[2].f[1] = condenser.divider.inlet[1].f[1] - condenser.divider.outlet[1].f[1]
    eq_18: v18 = v13 - v15;  # condenser.divider.outlet[2].f[2] = condenser.divider.inlet[1].f[2] - condenser.divider.outlet[1].f[2]
    eq_19: ((v17*32.04)+(v18*60.1))-96.0 = 0;  # ((condenser.divider.outlet[2].f[1]*32.04)+(condenser.divider.outlet[2].f[2]*60.1))-96.0 = 0
    eq_20: v20 = v11 - v16;  # condenser.divider.outlet[2].H = condenser.divider.inlet[1].H - condenser.divider.outlet[1].H
    # Connections
    eq_21: v21 = v20;  # cascade.stages[1].mixer.inlet[1].H = condenser.divider.outlet[2].H
    eq_22: v22 = v17;  # cascade.stages[1].mixer.inlet[1].f[1] = condenser.divider.outlet[2].f[1]
    eq_23: v23 = v18;  # cascade.stages[1].mixer.inlet[1].f[2] = condenser.divider.outlet[2].f[2]
    eq_24: v24 = v16;  # distillateSink.inlet.H = condenser.divider.outlet[1].H
    eq_25: v25 = v14;  # distillateSink.inlet.f[1] = condenser.divider.outlet[1].f[1]
    eq_26: v26 = v15;  # distillateSink.inlet.f[2] = condenser.divider.outlet[1].f[2]


Executable Python code is also emitted: It only serves for cross-checking 
correctness. 

--------------------------------------------------------------------------------

Future work: code generation for reverse mode automatic differentiation
=======================================================================

For efficient computations, templated C++ code will be emitted in
the future the Jacobian will be obtained with reverse mode automatic 
differentiation.

Source code is generated from the DAG representation of the expressions 
in order to compute the 
`Jacobian <http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant>`_
with reverse mode 
`automatic differentiation <http://en.wikipedia.org/wiki/Automatic_differentiation>`_. 
Currently only Python code is emitted, in the near future, templated C++ code 
will also be generated. For example, for the above example `exp(3*x+2*y)+4*z`
the following Python code is generated (hand-edited to improve readability)::

    # f = exp(3*x+2*y)+z
    # Forward sweep
    t1 = 3.0*x + 2.0*y
    t2 = exp(t1)
    f = 4.0*z + t2 - 1.0
    # Backward sweep
    u0 = 1.0
    u1 = 4.0 * u0  # df/dz = 4
    u2 = u0
    u3 = t2 * u2
    u4 = 3.0 * u3  # df/dx = 3*exp(3*x+2*y)
    u5 = 2.0 * u3  # df/dy = 2*exp(3*x+2*y)

The templated C++ version of this code will greatly benefit from code 
optimization performed by the C++ compiler; I expect the generated code to be 
as good as hand-written.

--------------------------------------------------------------------------------

Contents:

.. toctree::
   :maxdepth: 2



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

