#!/usr/bin/env python3
"""Algebraic stability audit for Thornado's limiting mixed-DG operator.

The one-dimensional radial operator below uses the same reference interval,
Gauss collocation, Lagrange face evaluation, weak volume/surface signs,
mass lumping, face averaging of scale factors, and weak scale-factor
derivatives as the OrderV streaming implementation.  It linearizes the
static isotropic moment equations, eliminates the stiff first moment, and
reports the largest eigenvalue of the mass-weighted symmetric part.  The
implemented quadratures support one through five nodes, i.e. polynomial
degrees zero through four.
"""

from __future__ import annotations

import argparse
import numpy as np


def lobatto(n: int) -> np.ndarray:
    if n == 1:
        return np.array([0.0])
    if n == 2:
        return np.array([-0.5, 0.5])
    if n == 3:
        return np.array([-0.5, 0.0, 0.5])
    if n == 4:
        return np.array([-0.5, -np.sqrt(5.0) / 10.0,
                         np.sqrt(5.0) / 10.0, 0.5])
    if n == 5:
        return np.array([-0.5, -np.sqrt(21.0) / 14.0, 0.0,
                         np.sqrt(21.0) / 14.0, 0.5])
    raise ValueError("Thornado implements one through five spatial nodes")


def lagrange_values(nodes: np.ndarray, x: float) -> np.ndarray:
    n = len(nodes)
    values = np.ones(n)
    for j in range(n):
        for k in range(n):
            if k != j:
                values[j] *= (x - nodes[k]) / (nodes[j] - nodes[k])
    return values


def derivative_matrix(nodes: np.ndarray) -> np.ndarray:
    n = len(nodes)
    matrix = np.zeros((n, n))
    for i, x in enumerate(nodes):
        for j in range(n):
            total = 0.0
            for m in range(n):
                if m == j:
                    continue
                product = 1.0 / (nodes[j] - nodes[m])
                for k in range(n):
                    if k != j and k != m:
                        product *= (x - nodes[k]) / (nodes[j] - nodes[k])
                total += product
            matrix[i, j] = total
    return matrix


def geometry_scale_factors(name: str, x: np.ndarray) -> np.ndarray:
    if name.startswith("cartesian"):
        return np.vstack((np.ones_like(x), np.ones_like(x), np.ones_like(x)))
    if name == "cylindrical-r":
        return np.vstack((np.ones_like(x), np.ones_like(x), x))
    if name in ("cylindrical-z", "cylindrical-phi"):
        radius = 1.4
        return np.vstack((np.ones_like(x), np.ones_like(x),
                          radius * np.ones_like(x)))
    if name == "spherical-r":
        return np.vstack((np.ones_like(x), x, x))
    if name == "spherical-schwarzschild-r":
        mass = 0.4
        psi2 = (1.0 + 0.5 * mass / np.abs(x)) ** 2
        return np.vstack((psi2, psi2 * x, psi2 * x))
    if name == "spherical-schwarzschild-r-strong":
        mass = 1.0
        psi2 = (1.0 + 0.5 * mass / np.abs(x)) ** 2
        return np.vstack((psi2, psi2 * x, psi2 * x))
    if name == "spherical-theta":
        radius = 1.4
        return np.vstack((np.ones_like(x), radius * np.ones_like(x),
                          radius * np.sin(x)))
    if name == "spherical-phi":
        radius = 1.4
        theta = 0.9
        return np.vstack((np.ones_like(x), radius * np.ones_like(x),
                          radius * np.sin(theta) * np.ones_like(x)))
    raise ValueError(name)


def direction_index(name: str) -> int:
    if name.endswith("-theta") or name.endswith("-z"):
        return 1
    if name.endswith("-phi"):
        return 2
    return 0


def build_operator(n: int, n_cells: int, geometry: str, eta: float,
                   r_left: float = 1.0, r_right: float = 2.0,
                   stretch: float = 1.0):
    # Thornado's reference interval is [-1/2,1/2].
    x_gl, w_gl = np.polynomial.legendre.leggauss(n)
    x_gl *= 0.5
    w_gl *= 0.5
    x_lob = lobatto(n)
    interp_lob_to_gl = np.vstack(
        [lagrange_values(x_lob, x) for x in x_gl]
    )
    deriv = derivative_matrix(x_gl)
    lower = lagrange_values(x_gl, -0.5)
    upper = lagrange_values(x_gl, 0.5)

    logical_edges = np.linspace(0.0, 1.0, n_cells + 1)
    edges = r_left + (r_right - r_left) * logical_edges ** stretch
    widths = np.diff(edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    size = n_cells * n

    # Geometry is evaluated in Lobatto points and interpolated to Gauss
    # points, exactly as GeometryComputationModule does.
    h = np.empty((n_cells, 3, n))
    for cell in range(n_cells):
        x_l = centers[cell] + widths[cell] * x_lob
        h_l = geometry_scale_factors(geometry, x_l)
        h[cell] = h_l @ interp_lob_to_gl.T
    jac = np.prod(h, axis=1)
    mass_diag = (widths[:, None] * w_gl[None, :] * jac).reshape(-1)
    mass = np.diag(mass_diag)

    # A single face geometry is constructed by averaging scale-factor
    # traces and recomputing the determinant.
    h_face = np.empty((n_cells + 1, 3))
    # Thornado fills geometry ghost cells before averaging their traces with
    # the interior trace.  For these analytic flat-coordinate metrics that
    # gives the metric at the physical boundary, including for p=0.
    h_face[0] = geometry_scale_factors(
        geometry, np.array([edges[0]])
    )[:, 0]
    h_face[-1] = geometry_scale_factors(
        geometry, np.array([edges[-1]])
    )[:, 0]
    for face in range(1, n_cells):
        h_face[face] = 0.5 * (h[face - 1] @ upper + h[face] @ lower)
    jac_face = np.prod(h_face, axis=1)
    direction = direction_index(geometry)

    def weak_scalar_derivative(values: np.ndarray,
                               common: np.ndarray) -> np.ndarray:
        """Thornado weak derivative used for the metric scale factors."""
        values = values.reshape(n_cells, n)
        result = np.empty_like(values)
        for cell in range(n_cells):
            raw = (-lower * common[cell] + upper * common[cell + 1]
                   - deriv.T @ (w_gl * values[cell]))
            result[cell] = raw / (w_gl * widths[cell])
        return result

    dh = np.empty_like(h)
    for component in range(3):
        dh[:, component] = weak_scalar_derivative(
            h[:, component], h_face[:, component]
        )

    def first_moment_rhs(density_flat: np.ndarray) -> np.ndarray:
        """Static isotropic G_1 streaming RHS including metric source."""
        density = density_flat.reshape(n_cells, n)
        traces_l = density @ lower
        traces_u = density @ upper
        d_hat = np.empty(n_cells + 1)
        d_hat[0] = traces_l[0]
        d_hat[-1] = traces_u[-1]
        d_hat[1:-1] = 0.5 * (traces_u[:-1] + traces_l[1:])
        rhs = np.empty_like(density)
        a2 = 1.0 / 3.0
        for cell in range(n_cells):
            flux_volume = a2 * jac[cell] * density[cell]
            raw = (lower * jac_face[cell] * a2 * d_hat[cell]
                   - upper * jac_face[cell + 1] * a2 * d_hat[cell + 1]
                   + deriv.T @ (w_gl * flux_volume))
            rhs[cell] = raw / (w_gl * widths[cell] * jac[cell])

        # The static isotropic part of the metric source used by the
        # observer-correction routine.  g_ii=h_i^2 and
        # S^{ii}=a^2 D/(2 g_ii), hence -dFlux is
        # a^2 D sum_i (partial h_i)/h_i.
        rhs += a2 * density * np.sum(dh / h, axis=1)
        return rhs.reshape(-1)

    # Constant scalar opacity.  With h_1=1 in all three flat coordinate
    # systems, the radial contravariant and covariant first moments agree.
    basis = np.eye(size)
    first_rhs_matrix = np.column_stack(
        [first_moment_rhs(basis[:, column]) for column in range(size)]
    )
    first_moment = first_rhs_matrix / (h[:, direction, :] ** 2).reshape(-1, 1)

    def density_rhs(first_flat: np.ndarray, density_flat: np.ndarray) -> np.ndarray:
        first = first_flat.reshape(n_cells, n)
        density = density_flat.reshape(n_cells, n)
        i_l = first @ lower
        i_u = first @ upper
        d_l = density @ lower
        d_u = density @ upper
        flux = np.zeros(n_cells + 1)  # reflecting boundaries
        for face in range(1, n_cells):
            width_face = 0.5 * (widths[face - 1] + widths[face])
            jump = d_l[face] - d_u[face - 1]
            normal_metric = 1.0 / h_face[face, direction] ** 2
            flux[face] = 0.5 * (i_u[face - 1] + i_l[face]) \
                         - eta * (1.0 / 3.0) * normal_metric * jump / width_face
        rhs = np.empty_like(density)
        for cell in range(n_cells):
            raw = (lower * jac_face[cell] * flux[cell]
                   - upper * jac_face[cell + 1] * flux[cell + 1]
                   + deriv.T @ (w_gl * jac[cell] * first[cell]))
            rhs[cell] = raw / (w_gl * widths[cell] * jac[cell])
        return rhs.reshape(-1)

    operator = np.column_stack([
        density_rhs(first_moment[:, column], basis[:, column])
        for column in range(size)
    ])
    symmetric = 0.5 * (mass @ operator + operator.T @ mass)
    eig_symmetric = np.linalg.eigvalsh(symmetric)
    eig_operator = np.linalg.eigvals(operator)
    spectral_scale = max(1.0, np.max(np.abs(eig_operator)))
    scale = max(1.0, np.linalg.norm(mass @ operator, ord=2))
    return {
        "max_symmetric": eig_symmetric[-1],
        "min_symmetric": eig_symmetric[0],
        "relative_max_symmetric": eig_symmetric[-1] / scale,
        "max_real_operator": np.max(eig_operator.real),
        "relative_max_real_operator": np.max(eig_operator.real) / spectral_scale,
        "nullity": int(np.count_nonzero(np.abs(eig_operator) < 1.0e-10 * spectral_scale)),
        "constant_residual": np.linalg.norm(operator @ np.ones(size), np.inf),
        "constant_first_moment": np.linalg.norm(
            first_moment @ np.ones(size), np.inf
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", type=int, default=8)
    parser.add_argument("--eta", type=float, default=1.0)
    parser.add_argument("--stretch", type=float, default=1.0)
    parser.add_argument(
        "--include-curved", action="store_true",
        help="also print exploratory conformal-Schwarzschild radial blocks",
    )
    args = parser.parse_args()
    print("geometry          nodes degree max_sym/scale   max_Re(L)/rho nullity  ||L1||inf        ||G1||inf")
    geometries = [
        ("cartesian-x", 1.0, 2.0),
        ("cylindrical-r", 1.0, 2.0),
        ("cylindrical-z", 0.0, 1.0),
        ("cylindrical-phi", 0.0, 1.0),
        ("spherical-r", 1.0, 2.0),
        ("spherical-theta", 0.2, np.pi - 0.2),
        ("spherical-phi", 0.0, 1.0),
    ]
    if args.include_curved:
        geometries[5:5] = [
            ("spherical-schwarzschild-r", 1.0, 2.0),
            ("spherical-schwarzschild-r-strong", 0.55, 2.0),
        ]
    for geometry, lower_bound, upper_bound in geometries:
        for n in range(1, 6):
            result = build_operator(
                n, args.cells, geometry, args.eta,
                r_left=lower_bound, r_right=upper_bound,
                stretch=args.stretch,
            )
            print(f"{geometry:17s} {n:5d} {n-1:6d} "
                  f"{result['relative_max_symmetric']: .6e} "
                  f"{result['relative_max_real_operator']: .6e} "
                  f"{result['nullity']:7d} "
                  f"{result['constant_residual']: .6e} "
                  f"{result['constant_first_moment']: .6e}")


if __name__ == "__main__":
    main()
