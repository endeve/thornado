#!/usr/bin/env python3
"""Matrix audit of Thornado's frozen-matter absorption--diffusion limit.

The matrices reproduce the one-dimensional Cartesian and flat-spherical,
static, isotropic specializations of the production weak DG streaming
operator.  Both LLF jump terms are retained.  The stationary first moment is
then eliminated exactly from the semi-discrete equations, rather than
replaced by pointwise division by the transport opacity.

This script is an algebraic audit.  It does not call Thornado and it does
not approximate the IMEX time integration.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def lagrange_values(nodes: np.ndarray, point: float) -> np.ndarray:
    values = np.ones(nodes.size)
    for j in range(nodes.size):
        for k in range(nodes.size):
            if k != j:
                values[j] *= (point - nodes[k]) / (nodes[j] - nodes[k])
    return values


def derivative_matrix(nodes: np.ndarray) -> np.ndarray:
    matrix = np.empty((nodes.size, nodes.size))
    for i, point in enumerate(nodes):
        for j in range(nodes.size):
            total = 0.0
            for m in range(nodes.size):
                if m == j:
                    continue
                term = 1.0 / (nodes[j] - nodes[m])
                for k in range(nodes.size):
                    if k != j and k != m:
                        term *= (point - nodes[k]) / (nodes[j] - nodes[k])
                total += term
            matrix[i, j] = total
    return matrix


@dataclass
class Operators:
    mass: np.ndarray
    x: np.ndarray
    opacity: np.ndarray
    density_from_flux: np.ndarray
    density_penalty: np.ndarray
    flux_from_density: np.ndarray
    flux_penalty: np.ndarray
    schur: np.ndarray
    pointwise_schur: np.ndarray


def lobatto_nodes(n_nodes: int) -> np.ndarray:
    if n_nodes == 1:
        return np.array([0.0])
    if n_nodes == 2:
        return np.array([-0.5, 0.5])
    if n_nodes == 3:
        return np.array([-0.5, 0.0, 0.5])
    if n_nodes == 4:
        return np.array([
            -0.5, -np.sqrt(5.0) / 10.0,
            +np.sqrt(5.0) / 10.0, +0.5,
        ])
    if n_nodes == 5:
        return np.array([
            -0.5, -np.sqrt(21.0) / 14.0, 0.0,
            +np.sqrt(21.0) / 14.0, +0.5,
        ])
    raise ValueError("Thornado implements one through five spatial nodes")


def build_operators(
    n_cells: int,
    n_nodes: int,
    mean_optical_depth: float,
    opacity_contrast: float = 0.0,
    source_aware: bool = True,
    alpha: float = 1.0,
    trace_pair: str = "centered",
) -> Operators:
    """Build the exact periodic Cartesian matrices used in the audit."""

    if trace_pair not in {"centered", "left-right", "right-left"}:
        raise ValueError(f"unknown trace pair: {trace_pair}")

    # Thornado uses Gauss nodes on [-1/2,1/2].
    nodes, weights = np.polynomial.legendre.leggauss(n_nodes)
    nodes *= 0.5
    weights *= 0.5
    lower = lagrange_values(nodes, -0.5)
    upper = lagrange_values(nodes, +0.5)
    derivative = derivative_matrix(nodes)

    dx = 1.0 / n_cells
    centers = (np.arange(n_cells) + 0.5) * dx
    x = (centers[:, None] + dx * nodes[None, :]).reshape(-1)
    size = n_cells * n_nodes
    mass = np.diag(np.tile(dx * weights, n_cells))

    kappa_mean = mean_optical_depth / dx
    opacity = kappa_mean * (
        1.0 + opacity_contrast * np.sin(2.0 * np.pi * x)
    )

    opacity_cells = opacity.reshape(n_cells, n_nodes)
    opacity_lower = opacity_cells @ lower
    opacity_upper = opacity_cells @ upper
    opacity_face = np.empty(n_cells)
    for face in range(n_cells):
        left = (face - 1) % n_cells
        right = face
        # This is the production equal-width specialization of the
        # width-weighted average of the two one-sided opacity traces.
        opacity_face[face] = 0.5 * (
            max(0.0, opacity_upper[left])
            + max(0.0, opacity_lower[right])
        )

    if source_aware:
        theta = 1.0 / (1.0 + 1.5 * alpha * opacity_face * dx)
    else:
        theta = np.ones(n_cells)

    def weak_rhs(volume: np.ndarray, face_flux: np.ndarray) -> np.ndarray:
        volume = volume.reshape(n_cells, n_nodes)
        result = np.empty_like(volume)
        for cell in range(n_cells):
            face_left = cell
            face_right = (cell + 1) % n_cells
            raw = (
                lower * face_flux[face_left]
                - upper * face_flux[face_right]
                + derivative.T @ (weights * volume[cell])
            )
            result[cell] = raw / (weights * dx)
        return result.reshape(-1)

    def traces(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        cells = values.reshape(n_cells, n_nodes)
        return cells @ lower, cells @ upper

    def density_rhs(first: np.ndarray, density: np.ndarray) -> np.ndarray:
        i_lower, i_upper = traces(first)
        d_lower, d_upper = traces(density)
        face_flux = np.empty(n_cells)
        for face in range(n_cells):
            left = (face - 1) % n_cells
            right = face
            jump_density = d_lower[right] - d_upper[left]
            if trace_pair == "centered":
                physical_flux = 0.5 * (i_upper[left] + i_lower[right])
            elif trace_pair == "left-right":
                physical_flux = i_upper[left]
            else:
                physical_flux = i_lower[right]
            face_flux[face] = physical_flux \
                - 0.5 * alpha * theta[face] * jump_density
        return weak_rhs(first, face_flux)

    def first_rhs(density: np.ndarray, first: np.ndarray) -> np.ndarray:
        d_lower, d_upper = traces(density)
        i_lower, i_upper = traces(first)
        face_flux = np.empty(n_cells)
        for face in range(n_cells):
            left = (face - 1) % n_cells
            right = face
            jump_first = i_lower[right] - i_upper[left]
            if trace_pair == "centered":
                physical_flux = (d_upper[left] + d_lower[right]) / 6.0
            elif trace_pair == "left-right":
                physical_flux = d_lower[right] / 3.0
            else:
                physical_flux = d_upper[left] / 3.0
            face_flux[face] = (
                physical_flux
                - 0.5 * alpha * jump_first
            )
        return weak_rhs(density / 3.0, face_flux)

    basis = np.eye(size)
    zero = np.zeros(size)
    density_from_flux = np.column_stack(
        [density_rhs(basis[:, column], zero) for column in range(size)]
    )
    density_penalty = np.column_stack(
        [density_rhs(zero, basis[:, column]) for column in range(size)]
    )
    flux_from_density = np.column_stack(
        [first_rhs(basis[:, column], zero) for column in range(size)]
    )
    flux_penalty = np.column_stack(
        [first_rhs(zero, basis[:, column]) for column in range(size)]
    )

    relaxation = np.diag(opacity) - flux_penalty
    exact_first = np.linalg.solve(relaxation, flux_from_density)
    pointwise_first = flux_from_density / opacity[:, None]
    schur = density_penalty + density_from_flux @ exact_first
    pointwise_schur = density_penalty + density_from_flux @ pointwise_first

    return Operators(
        mass=mass,
        x=x,
        opacity=opacity,
        density_from_flux=density_from_flux,
        density_penalty=density_penalty,
        flux_from_density=flux_from_density,
        flux_penalty=flux_penalty,
        schur=schur,
        pointwise_schur=pointwise_schur,
    )

def build_spherical_operators(
    n_cells: int,
    n_nodes: int,
    mean_optical_depth: float,
    opacity_contrast: float = 0.0,
    source_aware: bool = True,
    alpha: float = 1.0,
    stretch: float = 1.3,
) -> Operators:
    """Build the flat-spherical radial operator with reflecting boundaries."""

    nodes, weights = np.polynomial.legendre.leggauss(n_nodes)
    nodes *= 0.5
    weights *= 0.5
    lower = lagrange_values(nodes, -0.5)
    upper = lagrange_values(nodes, +0.5)
    derivative = derivative_matrix(nodes)

    logical_edges = np.linspace(0.0, 1.0, n_cells + 1)
    edges = 1.0 + logical_edges**stretch
    widths = np.diff(edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    x = (centers[:, None] + widths[:, None] * nodes[None, :])
    x_flat = x.reshape(-1)
    size = n_cells * n_nodes

    # GeometryComputationModule evaluates the scale factors at Lobatto
    # points and interpolates them to the Gauss points.  For flat spherical
    # coordinates h=(1,r,r), this interpolation is exact for n_nodes >= 2.
    lobatto = lobatto_nodes(n_nodes)
    interp = np.vstack([lagrange_values(lobatto, point) for point in nodes])
    radius_lobatto = centers[:, None] + widths[:, None] * lobatto[None, :]
    radius_gauss = radius_lobatto @ interp.T
    jacobian = radius_gauss**2
    jacobian_face = edges**2
    mass = np.diag((widths[:, None] * weights[None, :] * jacobian).reshape(-1))

    mean_width = float(np.mean(widths))
    kappa_mean = mean_optical_depth / mean_width
    phase = 2.0 * np.pi * (x_flat - 1.0)
    opacity = kappa_mean * (1.0 + opacity_contrast * np.cos(phase))
    opacity_cells = opacity.reshape(n_cells, n_nodes)
    opacity_lower = opacity_cells @ lower
    opacity_upper = opacity_cells @ upper
    opacity_face = np.zeros(n_cells + 1)
    opacity_face[0] = max(0.0, opacity_lower[0])
    opacity_face[-1] = max(0.0, opacity_upper[-1])
    for face in range(1, n_cells):
        face_width = 0.5 * (widths[face - 1] + widths[face])
        opacity_face[face] = (
            widths[face - 1] * max(0.0, opacity_upper[face - 1])
            + widths[face] * max(0.0, opacity_lower[face])
        ) / (2.0 * face_width)

    theta = np.ones(n_cells + 1)
    if source_aware:
        for face in range(1, n_cells):
            face_width = 0.5 * (widths[face - 1] + widths[face])
            theta[face] = 1.0 / (
                1.0 + 1.5 * alpha * opacity_face[face] * face_width
            )

    # Weak derivatives of h=(1,r,r), used by the production observer
    # correction.  The common face scale factors are the exact face radii.
    def weak_scalar_derivative(values: np.ndarray, common: np.ndarray) -> np.ndarray:
        result = np.empty_like(values)
        for cell in range(n_cells):
            raw = (
                -lower * common[cell]
                + upper * common[cell + 1]
                - derivative.T @ (weights * values[cell])
            )
            result[cell] = raw / (weights * widths[cell])
        return result

    dh_radius = weak_scalar_derivative(radius_gauss, edges)

    def weak_rhs(volume: np.ndarray, face_flux: np.ndarray) -> np.ndarray:
        volume = volume.reshape(n_cells, n_nodes)
        result = np.empty_like(volume)
        for cell in range(n_cells):
            raw = (
                lower * jacobian_face[cell] * face_flux[cell]
                - upper * jacobian_face[cell + 1] * face_flux[cell + 1]
                + derivative.T @ (weights * jacobian[cell] * volume[cell])
            )
            result[cell] = raw / (
                weights * widths[cell] * jacobian[cell]
            )
        return result.reshape(-1)

    def traces(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        cells = values.reshape(n_cells, n_nodes)
        return cells @ lower, cells @ upper

    def density_rhs(first: np.ndarray, density: np.ndarray) -> np.ndarray:
        i_lower, i_upper = traces(first)
        d_lower, d_upper = traces(density)
        face_flux = np.zeros(n_cells + 1)
        for face in range(1, n_cells):
            jump_density = d_lower[face] - d_upper[face - 1]
            face_flux[face] = 0.5 * (i_upper[face - 1] + i_lower[face]) \
                - 0.5 * alpha * theta[face] * jump_density
        return weak_rhs(first, face_flux)

    def first_rhs(density: np.ndarray, first: np.ndarray) -> np.ndarray:
        d_lower, d_upper = traces(density)
        i_lower, i_upper = traces(first)
        face_flux = np.empty(n_cells + 1)
        # Reflecting radiation boundary conditions: D is even and I is odd.
        face_flux[0] = density.reshape(n_cells, n_nodes)[0] @ lower / 3.0 \
            - alpha * i_lower[0]
        face_flux[-1] = density.reshape(n_cells, n_nodes)[-1] @ upper / 3.0 \
            + alpha * i_upper[-1]
        for face in range(1, n_cells):
            jump_first = i_lower[face] - i_upper[face - 1]
            face_flux[face] = (
                (d_upper[face - 1] + d_lower[face]) / 6.0
                - 0.5 * alpha * jump_first
            )
        result = weak_rhs(density.reshape(-1) / 3.0, face_flux)
        # Static isotropic observer correction: +2 D/(3r), constructed
        # from the same weak metric derivative as the production routine.
        metric_source = (
            (2.0 / 3.0) * density.reshape(n_cells, n_nodes)
            * dh_radius / radius_gauss
        )
        return result + metric_source.reshape(-1)

    basis = np.eye(size)
    zero = np.zeros(size)
    density_from_flux = np.column_stack(
        [density_rhs(basis[:, column], zero) for column in range(size)]
    )
    density_penalty = np.column_stack(
        [density_rhs(zero, basis[:, column]) for column in range(size)]
    )
    flux_from_density = np.column_stack(
        [first_rhs(basis[:, column], zero) for column in range(size)]
    )
    flux_penalty = np.column_stack(
        [first_rhs(zero, basis[:, column]) for column in range(size)]
    )
    relaxation = np.diag(opacity) - flux_penalty
    exact_first = np.linalg.solve(relaxation, flux_from_density)
    pointwise_first = flux_from_density / opacity[:, None]
    schur = density_penalty + density_from_flux @ exact_first
    pointwise_schur = density_penalty + density_from_flux @ pointwise_first
    return Operators(
        mass=mass,
        x=x_flat,
        opacity=opacity,
        density_from_flux=density_from_flux,
        density_penalty=density_penalty,
        flux_from_density=flux_from_density,
        flux_penalty=flux_penalty,
        schur=schur,
        pointwise_schur=pointwise_schur,
    )


def mass_norm(vector: np.ndarray, mass: np.ndarray) -> float:
    return float(np.sqrt(max(0.0, vector @ mass @ vector)))


def relative_mass_norm(
    difference: np.ndarray, reference: np.ndarray, mass: np.ndarray
) -> float:
    return mass_norm(difference, mass) / max(
        mass_norm(reference, mass), np.finfo(float).tiny
    )


def element_means(values: np.ndarray, n_cells: int, n_nodes: int) -> np.ndarray:
    _, weights = np.polynomial.legendre.leggauss(n_nodes)
    weights *= 0.5
    return values.reshape(n_cells, n_nodes) @ weights


def mass_weighted_element_means(
    values: np.ndarray, mass: np.ndarray, n_cells: int, n_nodes: int
) -> np.ndarray:
    """Element averages in the physical volume measure stored in mass."""
    cell_values = values.reshape(n_cells, n_nodes)
    cell_mass = np.diag(mass).reshape(n_cells, n_nodes)
    return np.sum(cell_mass * cell_values, axis=1) / np.sum(cell_mass, axis=1)


def checkerboard_projection(means: np.ndarray) -> float:
    alternating = (-1.0) ** np.arange(means.size)
    centered = means - np.mean(means)
    return float(abs(np.dot(centered, alternating)) / means.size)


def gridscale_rms(values: np.ndarray) -> float:
    """RMS of a high-pass filter normalized to unity on (-1)^element."""
    filtered = (
        0.5 * values
        - 0.25 * np.roll(values, 1)
        - 0.25 * np.roll(values, -1)
    )
    return float(np.sqrt(np.mean(filtered**2)))


def gridscale_rms_nonperiodic(values: np.ndarray) -> float:
    filtered = (
        0.5 * values[1:-1]
        - 0.25 * values[:-2]
        - 0.25 * values[2:]
    )
    return float(np.sqrt(np.mean(filtered**2)))


def continuum_diffusion_rhs(x: np.ndarray, kappa_mean: float, contrast: float) -> np.ndarray:
    wave = 2.0 * np.pi
    amplitude = 0.1
    density_prime = -amplitude * wave * np.sin(wave * x)
    density_second = -amplitude * wave**2 * np.cos(wave * x)
    opacity = kappa_mean * (1.0 + contrast * np.sin(wave * x))
    opacity_prime = kappa_mean * contrast * wave * np.cos(wave * x)
    diffusion = 1.0 / (3.0 * opacity)
    diffusion_prime = -opacity_prime / (3.0 * opacity**2)
    return diffusion_prime * density_prime + diffusion * density_second


def continuum_spherical_diffusion_rhs(
    radius: np.ndarray, kappa_mean: float, contrast: float
) -> np.ndarray:
    coordinate = radius - 1.0
    density_prime = -0.1 * np.pi * np.sin(np.pi * coordinate)
    density_second = -0.1 * np.pi**2 * np.cos(np.pi * coordinate)
    opacity = kappa_mean * (
        1.0 + contrast * np.cos(2.0 * np.pi * coordinate)
    )
    opacity_prime = (
        -2.0 * np.pi * kappa_mean * contrast
        * np.sin(2.0 * np.pi * coordinate)
    )
    diffusion = 1.0 / (3.0 * opacity)
    diffusion_prime = -opacity_prime / (3.0 * opacity**2)
    return (
        diffusion * density_second
        + (diffusion_prime + 2.0 * diffusion / radius) * density_prime
    )


def symbol_branches(operator: np.ndarray, n_cells: int, n_nodes: int) -> np.ndarray:
    branches = np.empty((n_cells, n_nodes), dtype=complex)
    for mode in range(n_cells):
        phase = np.exp(2.0j * np.pi * mode * np.arange(n_cells) / n_cells)
        lifting = np.zeros((n_cells * n_nodes, n_nodes), dtype=complex)
        for cell in range(n_cells):
            lifting[cell * n_nodes:(cell + 1) * n_nodes, :] \
                = phase[cell] * np.eye(n_nodes) / np.sqrt(n_cells)
        symbol = lifting.conj().T @ operator @ lifting
        branches[mode] = np.linalg.eigvals(symbol)
    return branches


def audit_p1_limiting_symbol(n_cells: int = 64) -> dict[str, float]:
    """Verify the analytic two-branch symbol of the p=1 thick-limit operator."""

    # source_aware=False exposes the unit-coefficient density jump matrix.
    operators = build_operators(
        n_cells, 2, mean_optical_depth=1.0,
        opacity_contrast=0.0, source_aware=False,
    )
    dx = 1.0 / n_cells
    limiting = (
        operators.density_from_flux @ operators.flux_from_density
        + (2.0 / (3.0 * dx)) * operators.density_penalty
    )
    numerical = symbol_branches(limiting, n_cells, 2)

    maximum_error = 0.0
    for mode in range(n_cells):
        phase = 2.0 * np.pi * mode / n_cells
        cosine = np.cos(phase)
        trace = (4.0 / (3.0 * dx**2)) * (
            cosine**2 + 2.0 * cosine - 6.0
        )
        determinant = (4.0 / (3.0 * dx**4)) * (
            1.0 - cosine
        ) * (5.0 - 3.0 * cosine)
        analytic = np.roots((1.0, -trace, determinant))
        error = np.max(np.abs(
            np.sort_complex(numerical[mode]) - np.sort_complex(analytic)
        ))
        maximum_error = max(maximum_error, float(error * dx**2))

    zero_mode = np.sort(numerical[0].real) * dx**2
    nyquist = np.sort(numerical[n_cells // 2].real) * dx**2
    eigenvalues = np.linalg.eigvals(limiting)
    scale = max(float(np.max(np.abs(eigenvalues))), 1.0)
    return {
        "formula_error_h2": maximum_error,
        "max_real_h2": float(np.max(eigenvalues.real) * dx**2),
        "nullity": float(np.count_nonzero(
            np.abs(eigenvalues) < 1.0e-11 * scale
        )),
        "zero_internal_h2": float(zero_mode[0]),
        "zero_physical_h2": float(zero_mode[1]),
        "nyquist_lower_h2": float(nyquist[0]),
        "nyquist_upper_h2": float(nyquist[1]),
    }


def audit_case(n_cells: int, n_nodes: int, tau: float, contrast: float) -> dict[str, float]:
    operators = build_operators(n_cells, n_nodes, tau, contrast, True)
    size = n_cells * n_nodes
    ones = np.ones(size)
    eig = np.linalg.eigvals(operators.schur)
    spectral_scale = max(np.max(np.abs(eig)), 1.0)
    nullity = int(np.count_nonzero(np.abs(eig) < 1.0e-11 * spectral_scale))

    density_equilibrium = 1.0 + 0.1 * np.cos(2.0 * np.pi * operators.x)
    discrete_rhs = operators.schur @ density_equilibrium
    dx = 1.0 / n_cells
    kappa_mean = tau / dx
    continuum_rhs = continuum_diffusion_rhs(
        operators.x, kappa_mean, contrast
    )
    exact_means = element_means(discrete_rhs, n_cells, n_nodes)
    continuum_means = element_means(continuum_rhs, n_cells, n_nodes)
    mean_error = exact_means - continuum_means

    # Verify the Schur complement against the complete frozen-matter steady
    # block system, including absorption.  The absorption fraction matches
    # the reduced Thornado absorption--diffusion test.
    absorption = 0.5 * operators.opacity
    absorption_matrix = np.diag(absorption)
    full_matrix = np.block([
        [operators.density_penalty - absorption_matrix,
         operators.density_from_flux],
        [operators.flux_from_density,
         operators.flux_penalty - np.diag(operators.opacity)],
    ])
    full_forcing = np.concatenate((-absorption * density_equilibrium,
                                   np.zeros(size)))
    full_solution = np.linalg.solve(full_matrix, full_forcing)
    density_schur = np.linalg.solve(
        absorption_matrix - operators.schur,
        absorption * density_equilibrium,
    )
    block_schur_error = relative_mass_norm(
        full_solution[:size] - density_schur,
        density_schur,
        operators.mass,
    )

    branches = symbol_branches(operators.schur, n_cells, n_nodes) \
        if contrast == 0.0 else None
    nyquist_max = np.nan
    if branches is not None:
        nyquist_max = float(np.max(branches[n_cells // 2].real))

    central_schur = operators.schur - operators.density_penalty
    central_eig = np.linalg.eigvals(central_schur)
    central_scale = max(np.max(np.abs(central_eig)), 1.0)
    central_nullity = int(
        np.count_nonzero(np.abs(central_eig) < 1.0e-11 * central_scale)
    )
    central_branches = symbol_branches(central_schur, n_cells, n_nodes) \
        if contrast == 0.0 else None
    central_nyquist_max = np.nan
    if central_branches is not None:
        central_nyquist_max = float(
            np.max(central_branches[n_cells // 2].real)
        )

    mass_symmetry_residual = np.linalg.norm(
        operators.mass @ operators.schur
        - operators.schur.T @ operators.mass
    ) / max(
        np.linalg.norm(operators.mass @ operators.schur),
        np.finfo(float).tiny,
    )

    adjoint_residual = np.linalg.norm(
        operators.mass @ operators.flux_from_density
        + (1.0 / 3.0)
        * operators.density_from_flux.T @ operators.mass
    ) / max(
        np.linalg.norm(operators.mass @ operators.flux_from_density),
        np.finfo(float).tiny,
    )

    return {
        "constant_residual": float(np.linalg.norm(operators.schur @ ones, np.inf)),
        "max_real_eigenvalue": float(np.max(eig.real)),
        "nullity": float(nullity),
        "central_nullity": float(central_nullity),
        "adjoint_residual": float(adjoint_residual),
        "mass_symmetry_residual": float(mass_symmetry_residual),
        "block_schur_error": block_schur_error,
        "pointwise_elimination_error": relative_mass_norm(
            (operators.schur - operators.pointwise_schur) @ density_equilibrium,
            discrete_rhs,
            operators.mass,
        ),
        "continuum_relative_error": relative_mass_norm(
            discrete_rhs - continuum_rhs, continuum_rhs, operators.mass
        ),
        "mean_continuum_relative_error": float(
            np.linalg.norm(exact_means - continuum_means)
            / max(np.linalg.norm(continuum_means), np.finfo(float).tiny)
        ),
        "checkerboard_absolute": checkerboard_projection(exact_means),
        "checkerboard_relative": checkerboard_projection(exact_means)
        / max(float(np.max(np.abs(exact_means))), np.finfo(float).tiny),
        "gridscale_error_relative": gridscale_rms(mean_error)
        / max(float(np.sqrt(np.mean(continuum_means**2))), np.finfo(float).tiny),
        "nyquist_max_real": nyquist_max,
        "central_nyquist_max_real": central_nyquist_max,
    }


def audit_spherical_case(
    n_cells: int, n_nodes: int, tau: float, contrast: float
) -> dict[str, float]:
    operators = build_spherical_operators(
        n_cells, n_nodes, tau, contrast, True
    )
    size = n_cells * n_nodes
    eig = np.linalg.eigvals(operators.schur)
    spectral_scale = max(np.max(np.abs(eig)), 1.0)
    central = operators.schur - operators.density_penalty
    central_eig = np.linalg.eigvals(central)
    central_scale = max(np.max(np.abs(central_eig)), 1.0)

    density_equilibrium = 1.0 + 0.1 * np.cos(
        np.pi * (operators.x - 1.0)
    )
    discrete_rhs = operators.schur @ density_equilibrium
    kappa_mean = tau * n_cells
    continuum_rhs = continuum_spherical_diffusion_rhs(
        operators.x, kappa_mean, contrast
    )
    discrete_means = mass_weighted_element_means(
        discrete_rhs, operators.mass, n_cells, n_nodes
    )
    continuum_means = mass_weighted_element_means(
        continuum_rhs, operators.mass, n_cells, n_nodes
    )
    interior_cells = slice(2, n_cells - 2)
    interior_nodes = np.zeros(size, dtype=bool)
    interior_nodes.reshape(n_cells, n_nodes)[interior_cells, :] = True
    interior_mass = operators.mass[np.ix_(interior_nodes, interior_nodes)]

    absorption = 0.5 * operators.opacity
    absorption_matrix = np.diag(absorption)
    full_matrix = np.block([
        [operators.density_penalty - absorption_matrix,
         operators.density_from_flux],
        [operators.flux_from_density,
         operators.flux_penalty - np.diag(operators.opacity)],
    ])
    full_forcing = np.concatenate((-absorption * density_equilibrium,
                                   np.zeros(size)))
    full_solution = np.linalg.solve(full_matrix, full_forcing)
    density_schur = np.linalg.solve(
        absorption_matrix - operators.schur,
        absorption * density_equilibrium,
    )

    adjoint_residual = np.linalg.norm(
        operators.mass @ operators.flux_from_density
        + (1.0 / 3.0)
        * operators.density_from_flux.T @ operators.mass
    ) / max(
        np.linalg.norm(operators.mass @ operators.flux_from_density),
        np.finfo(float).tiny,
    )

    mass_symmetry_residual = np.linalg.norm(
        operators.mass @ operators.schur
        - operators.schur.T @ operators.mass
    ) / max(
        np.linalg.norm(operators.mass @ operators.schur),
        np.finfo(float).tiny,
    )

    return {
        "constant_residual": float(
            np.linalg.norm(operators.schur @ np.ones(size), np.inf)
        ),
        "max_real_eigenvalue": float(np.max(eig.real)),
        "nullity": float(
            np.count_nonzero(np.abs(eig) < 1.0e-11 * spectral_scale)
        ),
        "central_nullity": float(
            np.count_nonzero(
                np.abs(central_eig) < 1.0e-11 * central_scale
            )
        ),
        "adjoint_residual": float(adjoint_residual),
        "mass_symmetry_residual": float(mass_symmetry_residual),
        "block_schur_error": relative_mass_norm(
            full_solution[:size] - density_schur,
            density_schur,
            operators.mass,
        ),
        "pointwise_elimination_error": relative_mass_norm(
            (operators.schur - operators.pointwise_schur)
            @ density_equilibrium,
            discrete_rhs,
            operators.mass,
        ),
        "continuum_relative_error": relative_mass_norm(
            discrete_rhs - continuum_rhs,
            continuum_rhs,
            operators.mass,
        ),
        "interior_continuum_relative_error": relative_mass_norm(
            (discrete_rhs - continuum_rhs)[interior_nodes],
            continuum_rhs[interior_nodes],
            interior_mass,
        ),
        "mean_continuum_relative_error": float(
            np.linalg.norm(discrete_means - continuum_means)
            / max(np.linalg.norm(continuum_means), np.finfo(float).tiny)
        ),
        "gridscale_error_relative": gridscale_rms(
            discrete_means - continuum_means
        ) / max(
            float(np.sqrt(np.mean(continuum_means**2))),
            np.finfo(float).tiny,
        ),
        "interior_mean_relative_error": float(
            np.linalg.norm(
                discrete_means[interior_cells]
                - continuum_means[interior_cells]
            )
            / max(
                np.linalg.norm(continuum_means[interior_cells]),
                np.finfo(float).tiny,
            )
        ),
        "interior_gridscale_error_relative": gridscale_rms_nonperiodic(
            (discrete_means - continuum_means)[interior_cells]
        ) / max(
            float(np.sqrt(np.mean(continuum_means[interior_cells]**2))),
            np.finfo(float).tiny,
        ),
    }


def write_figure(output: Path, tau: float, contrast: float) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), constrained_layout=True)
    for n_nodes, linestyle in ((2, "-"), (3, "--"), (4, ":")):
        operators = build_operators(32, n_nodes, tau, contrast, True)
        density = 1.0 + 0.1 * np.cos(2.0 * np.pi * operators.x)
        rhs = operators.schur @ density
        rhs_pointwise = operators.pointwise_schur @ density
        kappa_mean = tau * 32.0
        rhs_continuum = continuum_diffusion_rhs(
            operators.x, kappa_mean, contrast
        )
        axes[0, 0].plot(operators.x, rhs, linestyle, marker=".", ms=2,
                        label=f"{n_nodes} nodes")
        axes[0, 1].plot(operators.x, rhs - rhs_continuum, linestyle,
                        marker=".", ms=2, label=f"{n_nodes} nodes")
        axes[1, 0].plot(operators.x, rhs - rhs_pointwise, linestyle,
                        marker=".", ms=2, label=f"{n_nodes} nodes")
        means = element_means(rhs, 32, n_nodes)
        axes[1, 1].plot(np.arange(32), means, linestyle, marker=".", ms=3,
                        label=f"{n_nodes} nodes")

    axes[0, 0].set_title("Exact discrete Schur-complement diffusion RHS")
    axes[0, 1].set_title("Discrete minus continuum RHS")
    axes[1, 0].set_title("Exact minus pointwise first-moment elimination")
    axes[1, 1].set_title("Element-average discrete RHS")
    axes[0, 0].set_xlabel("x")
    axes[0, 1].set_xlabel("x")
    axes[1, 0].set_xlabel("x")
    axes[1, 1].set_xlabel("element")
    for axis in axes.flat:
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    fig.suptitle(
        f"Frozen-matter DG Schur audit: tau={tau:.0e}, contrast={contrast:.1f}"
    )
    fig.savefig(output)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=Path("frozen_matter_schur_audit.pdf"))
    parser.add_argument("--tau", type=float, default=1.0e5)
    parser.add_argument("--contrast", type=float, default=0.9)
    parser.add_argument("--no-figure", action="store_true")
    args = parser.parse_args()

    p1_symbol = audit_p1_limiting_symbol()
    print("p=1 opacity-scaled limiting symbol")
    print(
        "  max h^2 analytic error="
        f"{p1_symbol['formula_error_h2']:.3e}, "
        "max h^2 Re(lambda)="
        f"{p1_symbol['max_real_h2']:.3e}, "
        f"nullity={int(p1_symbol['nullity'])}"
    )
    print(
        "  h^2 lambda(phi=0)="
        f"({p1_symbol['zero_internal_h2']:.12g}, "
        f"{p1_symbol['zero_physical_h2']:.12g})"
    )
    print(
        "  h^2 lambda(phi=pi)="
        f"({p1_symbol['nyquist_lower_h2']:.12g}, "
        f"{p1_symbol['nyquist_upper_h2']:.12g})"
    )
    print()

    print(
        "cells nodes contrast const_res max_Re null/null0 adjoint symmetry block_err "
        "point_elim_err nodal_err mean_err grid_err cb_proj "
        "nyquist_Re/central"
    )
    for contrast in (0.0, args.contrast):
        for n_nodes in (2, 3, 4):
            for n_cells in (16, 32, 64):
                result = audit_case(n_cells, n_nodes, args.tau, contrast)
                print(
                    f"{n_cells:5d} {n_nodes:5d} {contrast:8.3f} "
                    f"{result['constant_residual']:9.2e} "
                    f"{result['max_real_eigenvalue']:9.2e} "
                    f"{int(result['nullity']):2d}/{int(result['central_nullity']):2d} "
                    f"{result['adjoint_residual']:9.2e} "
                    f"{result['mass_symmetry_residual']:9.2e} "
                    f"{result['block_schur_error']:9.2e} "
                    f"{result['pointwise_elimination_error']:12.4e} "
                    f"{result['continuum_relative_error']:9.2e} "
                    f"{result['mean_continuum_relative_error']:9.2e} "
                    f"{result['gridscale_error_relative']:9.2e} "
                    f"{result['checkerboard_relative']:9.2e} "
                    f"{result['nyquist_max_real']:10.3e}/"
                    f"{result['central_nyquist_max_real']:10.3e}"
                )

    print()
    print(
        "spherical: cells nodes contrast const_res max_Re null/null0 "
        "adjoint symmetry block_err point_elim_err nodal_err(int) "
        "mean_err(int) grid_err(int)"
    )
    for contrast in (0.0, args.contrast):
        for n_nodes in (2, 3, 4):
            for n_cells in (16, 32, 64):
                result = audit_spherical_case(
                    n_cells, n_nodes, args.tau, contrast
                )
                print(
                    f"{n_cells:5d} {n_nodes:5d} {contrast:8.3f} "
                    f"{result['constant_residual']:9.2e} "
                    f"{result['max_real_eigenvalue']:9.2e} "
                    f"{int(result['nullity']):2d}/{int(result['central_nullity']):2d} "
                    f"{result['adjoint_residual']:9.2e} "
                    f"{result['mass_symmetry_residual']:9.2e} "
                    f"{result['block_schur_error']:9.2e} "
                    f"{result['pointwise_elimination_error']:12.4e} "
                    f"{result['continuum_relative_error']:9.2e}"
                    f"({result['interior_continuum_relative_error']:9.2e}) "
                    f"{result['mean_continuum_relative_error']:9.2e}"
                    f"({result['interior_mean_relative_error']:9.2e}) "
                    f"{result['gridscale_error_relative']:9.2e}"
                    f"({result['interior_gridscale_error_relative']:9.2e})"
                )

    print()
    print(
        "Cartesian optical-depth sweep: cells=32 nodes=2 contrast="
        f"{args.contrast:g}"
    )
    print("tau_cell point_elim_err nodal_err mean_err grid_err")
    for tau in (1.0, 3.0, 10.0, 30.0, 1.0e2, 1.0e3, 1.0e4, 1.0e5):
        result = audit_case(32, 2, tau, args.contrast)
        print(
            f"{tau:9.1e} "
            f"{result['pointwise_elimination_error']:14.6e} "
            f"{result['continuum_relative_error']:9.3e} "
            f"{result['mean_continuum_relative_error']:9.3e} "
            f"{result['gridscale_error_relative']:9.3e}"
        )

    if not args.no_figure:
        write_figure(args.output, args.tau, args.contrast)


if __name__ == "__main__":
    main()
