Physical Overview
=================

This solver solves the Euler equations of non-relativistic hydrodynamics.
These equations represent the conservation of mass, momentum, and energy in a
non-viscous fluid. In a curvilinear coordinate system with metric tensor
components :math:`\gamma_{ij}`, they are written as

.. math::

	\frac{\partial \rho}{\partial t}
	+
	\frac{1}{\sqrt{\gamma}} \frac{\partial}{\partial x^{j}}
	\left(\sqrt{\gamma} \rho v^{j}\right) = 0,

.. math::

	\frac{\partial}{\partial t}\left(\rho v_{i}\right)
	+
	\frac{1}{\sqrt{\gamma}} \frac{\partial}{\partial x^{j}}
	\left(
		\sqrt{\gamma} \left(\rho v_{i} v^{j} + p \delta^{j}_{\, i}\right)
	\right)
    =
	\frac{1}{2}\left(v^{j} v^{k} + p \gamma^{jk}\right)
	\frac{\partial}{\partial x^{i}}\left(\gamma_{jk}\right)
	-
	\rho \frac{\partial \Phi}{\partial x^{i}},

.. math::

	\frac{\partial E}{\partial t}
	+
	\frac{1}{\sqrt{\gamma}} \frac{\partial}{\partial x^{j}}
	\left(
		\sqrt{\gamma} \left(E + p\right) v^{j}
	\right)
	=
    -\rho v^{j} \frac{\partial \Phi}{\partial x^{j}},

where :math:`\rho` is the mass density, :math:`v^{i}` is the i-th component of
the velocity, :math:`p` is the pressure, :math:`E = \rho \epsilon + \frac{1}{2}
\rho v^{i} v_{i}` is the energy density with specific internal energy
:math:`\epsilon`, :math:`\Phi` is the gravitational potential,
and :math:`\sqrt{\gamma}` is the determinant of the metric. This system
is closed by both Poisson's equation

.. math::

	\frac{1}{\sqrt{\gamma}} \frac{\partial}{\partial x^{i}}
	\left(
		\sqrt{\gamma} \gamma^{ij} \frac{\partial \Phi}{\partial x^{j}}
	\right)
	= 4 \pi G \rho,

with the gravitational constant :math:`G`,
and a general equation of state of the
form :math:`p = p(\rho, T, Y_{e})`, where :math:`T` is the temperature
and :math:`Y_{e}` is the electron fraction.
If :math:`Y_{e} \neq 0`, then we take :math:`\rho` to be the
baryon mass density, and evolve the electron mass density
:math:`\rho_{e} = \rho Y_{e}` separately with the equation

.. math::

	\frac{\partial \rho_{e}}{\partial t}
	+
	\frac{1}{\sqrt{\gamma}} \frac{\partial}{\partial x^{j}}
	\left(\sqrt{\gamma} \rho_{e} v^{j}\right) = 0
