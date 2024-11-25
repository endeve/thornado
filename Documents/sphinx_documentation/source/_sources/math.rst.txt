Mathematical Overview
=======================

These pages provide a mathematical overview of the numerical methods used in
thornado.

Continuous Galerkin
-------------------

Discontinuous Galerkin
----------------------

The use of the discontinuous Galerkin method dates back to the work of ...
This overview will assume a single spatial dimension, but see the addendum
for information about the generalization to multiple dimensions.
As seen in (link to relevant sections in Physical Overview), our goal is to
solve a system of hyperbolic conservation equations of the form

.. math::

   \frac{\partial \bf{U}}{\partial t} + \frac{1}{\sqrt{\gamma}}\frac{\partial}{\partial x}\left( \sqrt{\gamma} \bf{F}(\bf{U}) \right) = \bf{S}(\bf{U}).

We start by thinking about our spatial domain :math:`\mathcal{D}`, which, in 1D,
is just a line.
We can then chop this line into a series of non-overlapping segments called
elements. On a single element, which we call :math:`\bf{K}`, we have the
lower boundary :math:`x_{L}` and the upper boundary :math:`x_{U}`.
Since we do not want to worry with the particulars of each individual element
we instead represent all of them by a reference element
:math:`\bf{I}`, which uses the coordinate :math:`\eta \in (-\frac{1}{2},
\frac{1}{2})`. Now, any well-behaved (define this!) function :math:`f(\eta)`
on :math:`\bf{I}` can be written as the infinite series

.. math::

   f(\eta) = \sum_{i = 1}^{\infty}c_{i}\phi_{i}(\eta),

where :math:`\phi_{i}(\eta)` is a basis function. Truncating this series at
some finite index i = n, we are left with the approximation space V of
functions that we can then make out of the remaining basis functions.
