..
  Following the convention given in
  https://lpn-doc-sphinx-primer.readthedocs.io/en/stable/concepts/
  heading.html#:~:text=Section%20headers%20are%20created%20by,
  from%20the%20succession%20of%20headings.

##################################
General Relativistic Hydrodynamics
##################################

******************
Mathematical Model
******************

The equations of general relativistic hydrodynamics under the xCFC approximation
can be expressed as a system of hyperbolic conservation laws,
representing the conservation of baryon number, momentum, and energy, all
as measured by an Eulerian observer (as defined in the 3+1 decomposition).

The conserved variables are nonlinear functions of the primitive variables:
the baryon mass density measured by a comoving observer, :math:`\rho`,
the fluid three-velocity measured by an Eulerian observer, :math:`v^{i}`,
and the internal energy density measured by a comoving observer, :math:`e`.
Collecting the primitive variables into a column vector, :math:`\boldsymbol{V}`,
we have

.. math::

  \boldsymbol{V} := \left( \rho \, , v^{i} \, , e \right)^{T} \, .

These are related to the conserved baryon mass density, :math:`D`, as

.. math::

  D := \rho \, W \, ,

where :math:`W` is the Lorentz factor of the fluid measured by an Eulerian
observer; the momentum density, :math:`S_{j}`, as

.. math::

  S_{j} := \rho \, h \, W^{2} \, v_{j} \, ,

where :math:`h` is the specific enthalpy, defined as
:math:`h := 1 + \left( e + p \right) / \rho`, with :math:`p` the thermal
pressure; and a conserved energy density, :math:`\tau`, defined as

.. math::

  \tau := \rho \, h \, W^{2} - p - \rho \, W \, .

We collect the conserved variables into a column vector, :math:`\boldsymbol{U}`,

.. math::

  \boldsymbol{U} := \left( D \, , S_{j} \, , \tau \right)^{T} \, .

With that, we can write the GRHD equations as

.. math::

  \frac{\partial \left( \psi^{6} \, \boldsymbol{U} \right)}{\partial t}
  + \frac{\psi^{6}}{\sqrt{\gamma}}
  \frac{\partial \left( \alpha \, \sqrt{\gamma} \,
  \boldsymbol{F}^{i}\left(\boldsymbol{U}\right) \right)}
  {\partial x^{i}}
  = \alpha \, \psi^{6} \, \boldsymbol{S}\left(\boldsymbol{U}\right) \, ,

where :math:`\alpha` is the lapse function,
:math:`\psi` is the conformal factor,
and :math:`\sqrt{\gamma}` is the square root of the determinant of the spatial
three-metric.
(TODO: define the fluxes and sources)
