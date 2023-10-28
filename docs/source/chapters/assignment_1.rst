.. _ch:riemann_solver:

Riemann Solver
==============

This chapter implements and tests the most basic part of our project: The *f-wave* solver for the one-dimensional shallow water equations.
The shallow water equations are a system of nonlinear hyperbolic conservations laws with an optional source term:

.. math:: \begin{bmatrix} h \\ hu \end{bmatrix}_t + \begin{bmatrix} hu \\ hu^2 + \frac{1}{2}gh^2 \end{bmatrix}_x = S(x,t).
   :label: eq:swe1d


:math:`x \in \mathbb{R}` is the location in space and :math:`t \in \mathbb{R}^+` time.
The quantities :math:`q=[h, hu]^T` are given by the height of the water column :math:`h(x,t)` and the momentum :math:`hu(x,t)`, :math:`u` is the particle velocity.
:math:`g` is the used gravity, typically we use :math:`g:=9.80665\text{m}/\text{s}^2`.
:math:`f := [hu, hu^2 + \frac{1}{2}gh^2]^T` is the flux function.

.. figure:: data_1/tsunami_quantities.svg
  :name: fig:swe_quantities

  Sketch of the quantities appearing in the one-dimensional shallow water equations.

Our source term :math:`S(x,t)` will include the effect of space-dependent bathymetry (topography of the ocean) in future parts of the project.
The introduction of additional forces, e.g., friction or the coriolis effect, as part of the term is possible.
:numref:`fig:swe_quantities` illustrates all involved variables.

**Note:** As units we use meters (m) and seconds (s) for all computations.


Literature
----------

We discuss the basics with respect to the numerical formulation, software, and development strategies in our meetings.
Nevertheless, not all details can be covered in such a short time.
For now, the following list of books, papers and guides are recommended for your personal studies:

-  `Riemann Problems and Jupyter Solutions <https://www.clawpack.org/riemann_book/>`_,
   D. I. Ketcheson, R. J. LeVeque, M. J. del Razo, 2020

-  *Finite volume methods for hyperbolic problems*, R. J. LeVeque, 2002

-  *Riemann solvers and numerical methods for fluid dynamics*, E. F.
   Toro, 2009

-  *A wave propagation method for conservation laws and balance laws
   with spatially varying flux functions*, D. S. Bale et. al., 2003

-  `Thinking in C++ <https://archive.org/details/TICPP2ndEdVolOne>`_,
   Bruce Eckel, 2000

-  `git Documentation <https://git-scm.com/doc>`_

-  `Doxygen Manual <https://www.doxygen.nl/manual/index.html>`_

-  `Catch2 Tutorial
   <https://github.com/catchorg/Catch2/blob/master/docs/tutorial.md>`_

-  `SCons User Guide
   <https://www.scons.org/doc/production/HTML/scons-user.html>`_

-  `Information on ParaView <https://www.paraview.org/>`_

-  `Information on VisIt
   <https://wci.llnl.gov/simulation/computer-codes/visit>`_


.. _ch:start:

Getting Started
---------------

Our project starts with a given code to jump-start your developments.
First, we make sure that you can use the code and work on it collaboratively.

.. admonition:: Tasks

   #. The project’s initial code is available from `Github <https://github.com/breuera/tsunami_lab>`_.
      Fork the git-repository for your work.
      You may do this in your preferred way, e.g., by using your own resources or by using a provider of your choice.
      Make sure that all team members have access to your fork.

   #. The initial code ships with an implementation of the Roe solver.
      The Roe solver is similar to the targeted f-wave solver which you will develop as part of :numref:`ch:fwave`.
      Make yourself familiar with the initial code.
      The software uses `SCons <https://www.scons.org>`_ as a build tool.
      Build the code, run the unit tests, run the solver, and visualize the results.

   #. Generate a `Doxygen <https://www.doxygen.nl/>`_ documentation from the code.
      In future, comment your code using Doxygen syntax, especially functions and function parameters.

   #. Use git’s version control features for all changes and write meaningful commit messages.
      Consider going through pull requests and code reviews for changes.
      Make yourself familiar with git.

   #. It is paramount to test new features of our software before going to production.
      For this, we'll introduce unit tests through `Catch2 <https://github.com/catchorg/Catch2>`_ whenever possible.
      Make yourself familiar with Catch2.

.. warning::
   The last three tasks are ongoing and have to be fulfilled throughout the entire project.
   This means that you have to write meaningful documentation, meaningful commit messages and unit tests for every single part of the project!

.. _ch:fwave:

F-wave Solver
-------------

The f-wave solver approximately solves the following Initial Value Problem (IVP) for the shallow water equations :eq:`eq:swe1d` over time:

.. math::
   :label: eq:rp

   q(x,0) =
         \begin{cases}
           q_{l} \quad   &\text{if } x < 0 \\
           q_{r} \quad &\text{if }   x > 0
         \end{cases} \qquad q_l, q_r \in \mathbb{R}^+ \times \mathbb{R}.


Theory shows that the solution, arising from the discontinuity at :math:`x=0`, consist of two waves.
Each wave is either a shock or a rarefaction wave.
The f-wave solver uses two shock waves to approximate the true solution.

First, we use the Roe eigenvalues :math:`\lambda^{\text{Roe}}_{1/2}` in terms of the left and right quantities :math:`q_l` and :math:`q_r` with respect to position :math:`x=0` to approximate the true wave speeds:

.. math::
   :label: eq:eigenvalues

   \begin{aligned}
         \lambda^{\text{Roe}}_{1}(q_l, q_r) &= u^{\text{Roe}}(q_l, q_r) - \sqrt{gh^{\text{Roe}}(q_l, q_r)}, \\
         \lambda^{\text{Roe}}_{2}(q_l, q_r) &= u^{\text{Roe}}(q_l, q_r) + \sqrt{gh^{\text{Roe}}(q_l, q_r)},
   \end{aligned}


where the height :math:`h^{\text{Roe}}` and particle velocity :math:`u^{\text{Roe}}` are given as:

.. math::

   \begin{aligned}
         h^{\text{Roe}}(q_l, q_r) &= \frac{1}{2} (h_l + h_r), \\
         u^{\text{Roe}}(q_l, q_r) &=  \frac{u_l \sqrt{h_l} + u_r \sqrt{h_r}}{\sqrt{h_l}+\sqrt{h_r}}.
       \end{aligned}

Using the Roe eigenvalues we can define corresponding eigenvectors :math:`r_{1/2}^{\text{Roe}}`:

.. math::

   \begin{aligned}
         r_1^{\text{Roe}} &=
           \begin{bmatrix}
             1 \\ \lambda^{\text{Roe}}_1
           \end{bmatrix}, \\
         r_2^{\text{Roe}} &=
           \begin{bmatrix}
             1 \\ \lambda^{\text{Roe}}_2
           \end{bmatrix}.
       \end{aligned}

The decomposition of the jump in the flux function :math:`f`, :math:`\Delta f := f(q_r) - f(q_l)`, into the eigenvectors gives the waves :math:`Z_{1/2}`:

.. math:: \Delta f = \sum_{p=1}^2 \alpha_p r_p \equiv  \sum_{p=1}^2 Z_p, \qquad \alpha_p \in \mathbb{R}.
   :label: eq:jumpdec

The left “cell” :math:`\mathcal{C}_{l}` is influenced by the left-going waves (:math:`\lambda_p < 0`) and the right “cell” :math:`\mathcal{C}_r` by the right-going waves (:math:`\lambda_p > 0`).
This leads to the definition of net updates which summarize the net effect of the waves to the left and right “cell”:

.. math::
   :label: eq:netupdates

   \begin{split}
         A^- \Delta Q := \sum_{p:\{ \lambda_p^\text{Roe} < 0 \}} Z_p \\
         A^+ \Delta Q := \sum_{p:\{ \lambda_p^\text{Roe} > 0 \}} Z_p
   \end{split}

**Note:** The eigencoefficients :math:`\alpha_p` in Equation
(:eq:`eq:jumpdec`) are obtained by multiplying the `inverse
<https://mathworld.wolfram.com/MatrixInverse.html>`_ of the matrix of right
eigenvectors :math:`R=[r_1^\text{Roe}, r_2^\text{Roe}]` with the jump in
fluxes:

.. math::

   \begin{bmatrix}
         \alpha_1 \\
         \alpha_2
       \end{bmatrix} =
       \begin{bmatrix}
         1 & 1 \\
         \lambda^{\text{Roe}}_1 & \lambda^{\text{Roe}}_2
       \end{bmatrix}^{-1} \Delta f.

.. admonition:: Tasks

   #. Implement the f-wave solver for the homogeneous, i.e., :math:`S(x,t)=0`, shallow water equations.
      This implementation should be independent of the main code.
      The already implemented Roe solver in the files ``Roe.h``, ``Roe.cpp`` and ``Roe.test.cpp`` might serve as a guideline.

      *  Input values are the left state :math:`q_l = [h_l, (hu)_l]^T` and the right state :math:`q_r = [h_r, (hu)_r]^T`.

      *  Output values are the left- and right-going net updates: :math:`A^- \Delta Q` and :math:`A^+ \Delta Q`.


   #. As discussed above, write meaningful unit-tests using Catch2 for all implemented features.
      Examples for the basic f-wave solver are:

      *  Verification of the eigenvalue computation: You may derive a basic set of eigenvalues for given input values :math:`q_l` and :math:`q_r` manually.

      *  Zero -- with respect to machine precision -- net updates in the case of steady states, e.g., :math:`q_l = q_r`.

      *  Correctness tests for supersonic problems :math:`\lambda_{1/2}^\text{Roe} < 0 \vee \lambda_{1/2}^\text{Roe} > 0`:
         You may derive requirements using :eq:`eq:eigenvalues`.
         Remark: This implies that one of the net updates is zero as stated in :eq:`eq:netupdates`.
