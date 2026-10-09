#!/usr/bin/env python3
"""Modal and roundoff audit of a saved frozen-matter radiation operator.

The audit uses the same static, isotropic spherical DG matrices assembled by
``analyze_frozen_transient_operator.py``.  It eliminates the first moment
exactly at finite opacity, decomposes every spatial field into element-average
and within-element parts, and advances the complete affine frozen system with
a matrix exponential.  No additional Thornado evolution is required.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

from analyze_frozen_transient_operator import (
    assemble_isotropic_spherical_operator,
    column_map,
    read_table,
    select_block,
)


def relative_norm(a: np.ndarray, b: np.ndarray) -> float:
    scale = max(float(np.linalg.norm(b)), np.finfo(float).tiny)
    return float(np.linalg.norm(a - b) / scale)


def mass_norm(a: np.ndarray, mass: np.ndarray) -> float:
    value = np.dot(np.asarray(a).conj(), mass * np.asarray(a)).real
    return float(np.sqrt(max(value, 0.0)))


def projectors(
    mass: np.ndarray, n_cells: int, n_nodes: int
) -> tuple[np.ndarray, np.ndarray]:
    """Return the mass-orthogonal element-average and internal projectors."""

    mean = np.zeros((mass.size, mass.size))
    for cell in range(n_cells):
        local = slice(cell * n_nodes, (cell + 1) * n_nodes)
        weights = mass[local]
        mean[local, local] = np.outer(
            np.ones(n_nodes), weights / np.sum(weights)
        )
    internal = np.eye(mass.size) - mean
    return mean, internal


def weighted_frobenius(operator: np.ndarray, mass: np.ndarray) -> float:
    scale = np.sqrt(mass)
    transformed = scale[:, None] * operator / scale[None, :]
    return float(np.linalg.norm(transformed, ord="fro"))


def longdouble_norm(values: np.ndarray) -> np.longdouble:
    values = np.asarray(values, dtype=np.longdouble).ravel()
    return np.sqrt(np.sum(values * values, dtype=np.longdouble))


def longdouble_relative(a: np.ndarray, b: np.ndarray) -> float:
    a_ld = np.asarray(a, dtype=np.longdouble)
    b_ld = np.asarray(b, dtype=np.longdouble)
    scale = max(longdouble_norm(b_ld), np.finfo(np.longdouble).tiny)
    return float(longdouble_norm(a_ld - b_ld) / scale)


def fsum_matvec(matrix: np.ndarray, vector: np.ndarray) -> np.ndarray:
    """Matrix-vector product with correctly rounded scalar accumulation."""

    return np.fromiter(
        (
            math.fsum(float(a) * float(b) for a, b in zip(row, vector))
            for row in matrix
        ),
        dtype=float,
        count=matrix.shape[0],
    )


def kahan_matvec(matrix: np.ndarray, vector: np.ndarray) -> np.ndarray:
    result = np.empty(matrix.shape[0])
    for i, row in enumerate(matrix):
        total = 0.0
        correction = 0.0
        for a, b in zip(row, vector):
            value = float(a) * float(b) - correction
            updated = total + value
            correction = (updated - total) - value
            total = updated
        result[i] = total
    return result


def multiprecision_matvec(
    matrix: np.ndarray, vector: np.ndarray, decimal_digits: int
) -> np.ndarray:
    """Reference action retaining the exact binary64 matrix coefficients."""

    try:
        import mpmath as mp
    except ImportError as error:
        raise RuntimeError("--mp-dps requires mpmath") from error
    mp.mp.dps = decimal_digits
    result = np.empty(matrix.shape[0])
    for i, row in enumerate(matrix):
        nonzero = np.flatnonzero(row)
        total = mp.fsum(
            mp.mpf(float(row[j])) * mp.mpf(float(vector[j]))
            for j in nonzero
        )
        result[i] = float(total)
    return result


def refined_solve(
    matrix: np.ndarray, rhs: np.ndarray, iterations: int = 3
) -> tuple[np.ndarray, list[float]]:
    """Double solve with residuals accumulated in extended precision."""

    solution = np.linalg.solve(matrix, rhs)
    matrix_ld = np.asarray(matrix, dtype=np.longdouble)
    rhs_ld = np.asarray(rhs, dtype=np.longdouble)
    residual_history: list[float] = []
    rhs_scale = max(longdouble_norm(rhs_ld), np.finfo(np.longdouble).tiny)
    for _ in range(iterations + 1):
        residual_ld = rhs_ld - matrix_ld @ np.asarray(
            solution, dtype=np.longdouble
        )
        residual_history.append(float(longdouble_norm(residual_ld) / rhs_scale))
        if len(residual_history) > iterations:
            break
        correction = np.linalg.solve(matrix, np.asarray(residual_ld, dtype=float))
        solution += correction
    return solution, residual_history


def modal_metrics(
    values: np.ndarray,
    mean_projector: np.ndarray,
    internal_projector: np.ndarray,
    mass: np.ndarray,
    radius: np.ndarray,
    n_cells: int,
    n_nodes: int,
) -> dict[str, float]:
    mean = mean_projector @ values
    internal = internal_projector @ values
    cells = np.asarray(values).real.reshape(n_cells, n_nodes)
    radius_cells = radius.reshape(n_cells, n_nodes)
    mean_cells = mean.real.reshape(n_cells, n_nodes)[:, 0]
    centers = np.sum(
        mass.reshape(n_cells, n_nodes) * radius_cells, axis=1
    ) / np.sum(mass.reshape(n_cells, n_nodes), axis=1)
    if n_nodes == 1:
        return {
            "internal_fraction": 0.0,
            "backward_fraction": 0.0,
            "slope_correlation": np.nan,
            "face_jump_fraction": 0.0,
        }

    local_slope = (cells[:, -1] - cells[:, 0]) / (
        radius_cells[:, -1] - radius_cells[:, 0]
    )
    envelope_slope = np.gradient(mean_cells, centers)
    scale = max(float(np.max(np.abs(envelope_slope))), np.finfo(float).tiny)
    active = np.abs(envelope_slope) > 1.0e-10 * scale
    slope_weight = np.abs(local_slope[active])
    backward = local_slope[active] * envelope_slope[active] < 0.0
    backward_fraction = (
        float(np.sum(slope_weight[backward]) / np.sum(slope_weight))
        if np.sum(slope_weight) > 0.0
        else 0.0
    )
    if np.count_nonzero(active) > 1:
        slope_correlation = float(
            np.corrcoef(local_slope[active], envelope_slope[active])[0, 1]
        )
    else:
        slope_correlation = np.nan

    jumps = cells[1:, 0] - cells[:-1, -1]
    face_scale = max(float(np.linalg.norm(np.diff(mean_cells))), np.finfo(float).tiny)
    return {
        "internal_fraction": mass_norm(internal, mass)
        / max(mass_norm(values, mass), np.finfo(float).tiny),
        "backward_fraction": backward_fraction,
        "slope_correlation": slope_correlation,
        "face_jump_fraction": float(np.linalg.norm(jumps) / face_scale),
    }


def print_metrics(label: str, metrics: dict[str, float]) -> None:
    print(
        f"  {label:25s} internal={metrics['internal_fraction']:.6e}  "
        f"backward={metrics['backward_fraction']:.6e}  "
        f"slope corr={metrics['slope_correlation']:.6e}  "
        f"face jump={metrics['face_jump_fraction']:.6e}"
    )


def regional_slope_metrics(
    values: np.ndarray,
    radius_km: np.ndarray,
    radial_range: tuple[float, float],
    n_cells: int,
    n_nodes: int,
) -> tuple[float, float]:
    """Return slope-weighted backward fraction and slope correlation."""

    if n_nodes == 1:
        return 0.0, np.nan
    cells = np.asarray(values).real.reshape(n_cells, n_nodes)
    radii = radius_km.reshape(n_cells, n_nodes)
    centers = np.mean(radii, axis=1)
    averages = np.mean(cells, axis=1)
    local_slope = (cells[:, -1] - cells[:, 0]) / (
        radii[:, -1] - radii[:, 0]
    )
    envelope_slope = np.gradient(averages, centers)
    selected = (
        (centers >= radial_range[0])
        & (centers <= radial_range[1])
    )
    scale = max(
        float(np.max(np.abs(envelope_slope[selected]))),
        np.finfo(float).tiny,
    )
    selected &= np.abs(envelope_slope) > 1.0e-10 * scale
    weights = np.abs(local_slope[selected])
    backward = local_slope[selected] * envelope_slope[selected] < 0.0
    fraction = (
        float(np.sum(weights[backward]) / np.sum(weights))
        if np.sum(weights) > 0.0
        else 0.0
    )
    correlation_value = (
        float(np.corrcoef(
            local_slope[selected], envelope_slope[selected]
        )[0, 1])
        if np.count_nonzero(selected) > 1
        else np.nan
    )
    return fraction, correlation_value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=Path)
    parser.add_argument("--species", type=int, default=1)
    parser.add_argument("--energy", type=int, default=25)
    parser.add_argument("--plot", type=Path)
    parser.add_argument("--radius-min-km", type=float, default=10.0)
    parser.add_argument("--radius-max-km", type=float, default=45.0)
    parser.add_argument(
        "--full-spectrum",
        action="store_true",
        help="also diagonalize the complete 2N by 2N transient generator",
    )
    parser.add_argument(
        "--mp-dps",
        type=int,
        default=0,
        help="also recompute the spatial action with this many decimal digits",
    )
    args = parser.parse_args()

    metadata, columns, table = read_table(args.file)
    data = column_map(columns, table)
    block = select_block(data, args.species, args.energy)
    order = np.lexsort((block["spatial_node"], block["element"]))
    for name in block:
        block[name] = block[name][order]

    if max(np.max(np.abs(block[name])) for name in ("v1", "v2", "v3")) > 1.0e-12:
        raise ValueError("this audit currently requires a static fluid")

    n_cells = int(np.max(block["element"]))
    n_nodes = int(np.max(block["spatial_node"]))
    n = n_cells * n_nodes
    source_aware = bool(round(metadata.get("source_aware_llf", 0.0)))
    a_dd, a_di, a_id, a_ii = assemble_isotropic_spherical_operator(
        block, source_aware
    )

    chi = block["chi_emab"]
    kappa = block["kappa_transport"]
    relaxation = np.diag(kappa) - a_ii
    fick_map, refinement = refined_solve(relaxation, a_id)
    schur = a_dd - np.diag(chi) + a_di @ fick_map
    forcing_density = chi * block["D0"]

    radius = block["r_code"]
    mass = block["wx"] * block["dr_code"] * block["sqrtgm"]
    mean, internal = projectors(mass, n_cells, n_nodes)

    density = block["N_before"]
    first = block["G1_before"]
    first_equilibrated = fick_map @ density
    first_equilibrated_refined, vector_refinement = refined_solve(
        relaxation, a_id @ density
    )
    first_from_means = fick_map @ (mean @ density)
    first_from_internal = fick_map @ (internal @ density)

    print(
        f"Frozen real-profile Schur audit: species={args.species}, "
        f"energy={args.energy}, cells={n_cells}, nodes={n_nodes}"
    )
    print(f"  extended precision: {np.finfo(np.longdouble).nmant + 1} significand bits")
    print(f"  cond_2(K-J_I)={np.linalg.cond(relaxation):.6e}")
    print(
        "  iterative-refinement residuals: "
        + " ".join(f"{value:.3e}" for value in refinement)
    )
    print(
        "  physical-vector refinement residuals: "
        + " ".join(f"{value:.3e}" for value in vector_refinement)
    )
    print(
        "  refined change in R^-1 B D="
        f"{relative_norm(first_equilibrated_refined, first_equilibrated):.6e}"
    )

    stream = np.block([[a_dd, a_di], [a_id, a_ii]])
    initial = np.concatenate((density, first))
    action_64 = stream @ initial
    action_ld = np.asarray(
        np.asarray(stream, dtype=np.longdouble)
        @ np.asarray(initial, dtype=np.longdouble),
        dtype=float,
    )
    action_fsum = fsum_matvec(stream, initial)
    action_kahan = kahan_matvec(stream, initial)
    print("Precision sensitivity of the frozen spatial action:")
    print(f"  float64 vs extended accumulation={relative_norm(action_64, action_ld):.6e}")
    print(f"  float64 vs fsum accumulation   ={relative_norm(action_64, action_fsum):.6e}")
    print(f"  float64 vs Kahan accumulation  ={relative_norm(action_64, action_kahan):.6e}")
    if args.mp_dps > 0:
        action_mp = multiprecision_matvec(stream, initial, args.mp_dps)
        print(
            f"  float64 vs {args.mp_dps}-digit accumulation="
            f"{relative_norm(action_64, action_mp):.6e}"
        )

    print("First-moment spatial representation:")
    for label, values in (
        ("initial I", first),
        ("R^-1 B D", first_equilibrated),
        ("from element means", first_from_means),
        ("from internal D modes", first_from_internal),
    ):
        print_metrics(
            label,
            modal_metrics(values, mean, internal, mass, radius, n_cells, n_nodes),
        )
    print(
        "  contribution norms: "
        f"||I(mean D)||/||I||={mass_norm(first_from_means, mass) / max(mass_norm(first_equilibrated, mass), np.finfo(float).tiny):.6e}, "
        f"||I(internal D)||/||I||={mass_norm(first_from_internal, mass) / max(mass_norm(first_equilibrated, mass), np.finfo(float).tiny):.6e}"
    )
    radius_km = radius / metadata["centimeter_code"] / 1.0e5
    radial_range = (args.radius_min_km, args.radius_max_km)
    print(
        f"Slope audit on {radial_range[0]:g}--{radial_range[1]:g} km "
        "(weighted by the absolute within-element slope):"
    )
    for label, values in (
        ("initial I", first),
        ("R^-1 B D", first_equilibrated),
        ("from element means", first_from_means),
        ("from internal D modes", first_from_internal),
    ):
        backward, slope_correlation = regional_slope_metrics(
            values, radius_km, radial_range, n_cells, n_nodes
        )
        print(
            f"  {label:25s} backward={backward:.6e}  "
            f"slope corr={slope_correlation:.6e}"
        )

    # A global interpolant is not an independent reference solution, but it
    # is a useful diagnostic of whether the backwards element slope is already
    # implied by the smooth nodal envelope.  Cubic and monotone interpolants
    # are both reported so this conclusion does not depend on one interpolant.
    try:
        from scipy.interpolate import CubicSpline, PchipInterpolator
    except ImportError as error:
        raise RuntimeError("this audit requires scipy") from error
    selected_nodes = (
        (radius_km >= radial_range[0])
        & (radius_km <= radial_range[1])
    )
    print("Smooth-envelope Fick comparison on the same radial interval:")
    for label, interpolant in (
        ("cubic", CubicSpline(radius, density)),
        ("monotone", PchipInterpolator(radius, density)),
    ):
        reference = -interpolant.derivative()(radius) / (3.0 * kappa)
        if np.dot(reference[selected_nodes], first_equilibrated[selected_nodes]) < 0.0:
            reference = -reference
        backward, slope_correlation = regional_slope_metrics(
            reference, radius_km, radial_range, n_cells, n_nodes
        )
        print(
            f"  {label:8s} rel={relative_norm(first_equilibrated[selected_nodes], reference[selected_nodes]):.6e}  "
            f"corr={np.corrcoef(first_equilibrated[selected_nodes], reference[selected_nodes])[0, 1]:.6e}  "
            f"backward={backward:.6e}  slope corr={slope_correlation:.6e}"
        )

    # Block norms quantify transfer between element averages and internal
    # modes.  These are properties of the frozen operator, independent of the
    # particular initial data.
    schur_blocks = {
        "mean <- mean": mean @ schur @ mean,
        "mean <- internal": mean @ schur @ internal,
        "internal <- mean": internal @ schur @ mean,
        "internal <- internal": internal @ schur @ internal,
    }
    schur_norm = max(weighted_frobenius(schur, mass), np.finfo(float).tiny)
    print("Mass-scaled Schur block Frobenius norms / ||S||:")
    for label, operator in schur_blocks.items():
        print(f"  {label:20s} {weighted_frobenius(operator, mass) / schur_norm:.6e}")

    eigenvalues = np.linalg.eigvals(schur)
    maximum_real = float(np.max(eigenvalues.real))
    spectral_scale = max(float(np.max(np.abs(eigenvalues))), 1.0)
    unstable = int(np.count_nonzero(eigenvalues.real > 1.0e-11 * spectral_scale))
    mass_sqrt = np.sqrt(mass)
    schur_mass = mass_sqrt[:, None] * schur / mass_sqrt[None, :]
    numerical_abscissa = float(
        np.max(np.linalg.eigvalsh(0.5 * (schur_mass + schur_mass.T)))
    )
    nonnormality = float(
        np.linalg.norm(schur_mass.T @ schur_mass - schur_mass @ schur_mass.T)
        / max(np.linalg.norm(schur_mass) ** 2, np.finfo(float).tiny)
    )
    print("Frozen Schur spectrum:")
    print(f"  max Re(lambda)={maximum_real:.6e}  unstable modes={unstable}")
    print(f"  numerical abscissa={numerical_abscissa:.6e}")
    print(f"  normalized nonnormality={nonnormality:.6e}")

    # The complete affine system supplies both the exact one-step state and
    # its stationary limit.  Comparing a single exponential with two half
    # exponentials checks the numerical matrix-exponential evaluation.
    try:
        from scipy.sparse import bmat, csr_matrix
        from scipy.sparse.linalg import expm_multiply
    except ImportError as error:
        raise RuntimeError("this audit requires scipy") from error

    generator = np.block(
        [
            [a_dd - np.diag(chi), a_di],
            [a_id, a_ii - np.diag(kappa)],
        ]
    )
    forcing = np.concatenate((forcing_density, np.zeros(n)))
    augmented = bmat(
        [
            [csr_matrix(generator), csr_matrix(forcing[:, None])],
            [None, csr_matrix((1, 1))],
        ],
        format="csr",
    )
    augmented_initial = np.concatenate((initial, [1.0]))
    dt = metadata["timestep_code"]
    exact_one = expm_multiply(augmented * dt, augmented_initial)
    exact_half = expm_multiply(augmented * (0.5 * dt), augmented_initial)
    exact_two_half = expm_multiply(augmented * (0.5 * dt), exact_half)
    print(
        "  exp(dt A) vs exp(dt A/2)^2="
        f"{relative_norm(exact_one[:-1], exact_two_half[:-1]):.6e}"
    )

    stationary, stationary_refinement = refined_solve(generator, -forcing)
    stationary_residual = np.asarray(generator, dtype=np.longdouble) @ np.asarray(
        stationary, dtype=np.longdouble
    ) + np.asarray(forcing, dtype=np.longdouble)
    stationary_scale = max(
        longdouble_norm(
            np.asarray(generator, dtype=np.longdouble)
            @ np.asarray(stationary, dtype=np.longdouble)
        ),
        longdouble_norm(np.asarray(forcing, dtype=np.longdouble)),
        np.finfo(np.longdouble).tiny,
    )
    print(
        "  stationary affine residual (extended accumulation)="
        f"{float(longdouble_norm(stationary_residual) / stationary_scale):.6e}"
    )
    print(
        "  stationary refinement residuals: "
        + " ".join(f"{value:.3e}" for value in stationary_refinement)
    )
    print_metrics(
        "exact I after one dt",
        modal_metrics(
            exact_one[n:2 * n], mean, internal, mass, radius, n_cells, n_nodes
        ),
    )
    print_metrics(
        "stationary I",
        modal_metrics(
            stationary[n:], mean, internal, mass, radius, n_cells, n_nodes
        ),
    )
    for label, values in (
        ("exact I after one dt", exact_one[n:2 * n]),
        ("stationary I", stationary[n:]),
    ):
        backward, slope_correlation = regional_slope_metrics(
            values, radius_km, radial_range, n_cells, n_nodes
        )
        print(
            f"  {label:25s} on {radial_range[0]:g}--{radial_range[1]:g} km: "
            f"backward={backward:.6e}  slope corr={slope_correlation:.6e}"
        )

    if args.full_spectrum:
        from scipy.linalg import eig

        transient_eigenvalues, transient_eigenvectors = eig(
            generator, check_finite=False, overwrite_a=False
        )
        i_mode = int(np.argmax(transient_eigenvalues.real))
        mode = transient_eigenvectors[:, i_mode]
        mode_power = (
            np.abs(mode[:n].reshape(n_cells, n_nodes)) ** 2
            + np.abs(mode[n:].reshape(n_cells, n_nodes)) ** 2
        ).sum(axis=1)
        centers_km = np.mean(radius_km.reshape(n_cells, n_nodes), axis=1)
        total_power = max(float(np.sum(mode_power)), np.finfo(float).tiny)
        print("Complete transient spectrum:")
        print(
            f"  max Re(lambda)={transient_eigenvalues[i_mode].real:.6e}  "
            f"Im(lambda)={transient_eigenvalues[i_mode].imag:.6e}"
        )
        print(
            "  leading-mode power fractions: "
            f"r < {radial_range[1]:g} km="
            f"{np.sum(mode_power[centers_km < radial_range[1]]) / total_power:.6e}, "
            f"r > 100 km={np.sum(mode_power[centers_km > 100.0]) / total_power:.6e}"
        )

    if args.plot:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
        axes[0].plot(radius_km, density, ".-", label=r"$D$")
        axes[0].plot(radius_km, mean @ density, "-", label="element mean")
        axes[0].set_ylabel("radiation density")
        axes[0].legend(loc="best")

        axes[1].plot(radius_km, first_equilibrated, ".-", label=r"$R^{-1}BD$")
        axes[1].plot(radius_km, first_from_means, "-", label="from element means")
        axes[1].plot(
            radius_km,
            first_from_internal,
            "-",
            label="from internal density modes",
        )
        axes[1].set_ylabel("stationary first moment")
        axes[1].legend(loc="best")

        axes[2].plot(
            radius_km,
            exact_one[n:2 * n],
            ".-",
            label="exact after one timestep",
        )
        axes[2].plot(radius_km, stationary[n:], "-", label="stationary")
        axes[2].set_ylabel("first moment")
        axes[2].set_xlabel("radius [km]")
        axes[2].legend(loc="best")
        for axis in axes:
            axis.set_xlim(radial_range)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=180)


if __name__ == "__main__":
    main()
