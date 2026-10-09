#!/usr/bin/env python3
"""Compare the frozen DG first moment with a projected smooth Fick flux.

For each spectral node the saved nodal density and transport opacity are
continued smoothly in radius.  The continuum static diffusion relation

    I(r,E) = -[3 kappa(r,E)]^{-1} dD(r,E)/dr

is evaluated at high-order quadrature points and projected into the same
spatial DG space using the spherical volume measure.  The audit separates
the weak pressure-gradient operator, pointwise opacity division, and exact
finite-opacity elimination of the first-moment LLF operator.  It also forms
the energy-weighted grey first moment used by the existing diagnostics.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.linalg import solve

from analyze_frozen_transient_operator import (
    assemble_isotropic_spherical_operator,
    column_map,
    lagrange_values,
    read_table,
    select_block,
)


@dataclass
class SmoothField:
    value: object
    derivative: object


def smooth_positive_field(
    radius: np.ndarray, values: np.ndarray, method: str
) -> SmoothField:
    if np.any(values <= 0.0):
        raise ValueError(f"{method} positive-field interpolation received nonpositive data")
    if method == "pchip-linear":
        interpolation = PchipInterpolator(radius, values)
        return SmoothField(interpolation, interpolation.derivative())
    if method == "pchip-log":
        interpolation = PchipInterpolator(radius, np.log(values))
    elif method == "cubic-log":
        interpolation = CubicSpline(radius, np.log(values))
    else:
        raise ValueError(f"unknown interpolation method {method}")

    def value(points: np.ndarray) -> np.ndarray:
        return np.exp(interpolation(points))

    def derivative(points: np.ndarray) -> np.ndarray:
        return value(points) * interpolation.derivative()(points)

    return SmoothField(value, derivative)


def project_function(
    function,
    radius: np.ndarray,
    widths: np.ndarray,
    nodes: np.ndarray,
    quadrature_order: int,
) -> np.ndarray:
    """Full spherical-volume L2 projection into the nodal DG basis."""

    n_cells = widths.size
    n_nodes = nodes.size
    centers = np.mean(radius.reshape(n_cells, n_nodes), axis=1)
    quadrature_nodes, quadrature_weights = np.polynomial.legendre.leggauss(
        quadrature_order
    )
    quadrature_nodes *= 0.5
    quadrature_weights *= 0.5
    basis = np.column_stack(
        [lagrange_values(nodes, point) for point in quadrature_nodes]
    ).T
    result = np.empty((n_cells, n_nodes))
    for cell in range(n_cells):
        points = centers[cell] + widths[cell] * quadrature_nodes
        measure = quadrature_weights * points**2
        mass = widths[cell] * basis.T @ (measure[:, None] * basis)
        rhs = widths[cell] * basis.T @ (measure * function(points))
        result[cell] = np.linalg.solve(mass, rhs)
    return result.reshape(-1)


def element_decomposition(
    values: np.ndarray, mass: np.ndarray, n_cells: int, n_nodes: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cells = values.reshape(n_cells, n_nodes)
    weights = mass.reshape(n_cells, n_nodes)
    averages = np.sum(weights * cells, axis=1) / np.sum(weights, axis=1)
    mean = np.repeat(averages[:, None], n_nodes, axis=1)
    return averages, mean.reshape(-1), (cells - mean).reshape(-1)


def weighted_norm(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sqrt(max(np.sum(weights * np.abs(values) ** 2), 0.0)))


def relative_error(
    numerical: np.ndarray,
    reference: np.ndarray,
    weights: np.ndarray,
) -> float:
    return weighted_norm(numerical - reference, weights) / max(
        weighted_norm(reference, weights), np.finfo(float).tiny
    )


def radial_mask(
    radius_km: np.ndarray, lower: float, upper: float
) -> np.ndarray:
    return (radius_km >= lower) & (radius_km <= upper)


def slope_diagnostics(
    values: np.ndarray,
    radius_km: np.ndarray,
    selected_cells: np.ndarray,
    n_cells: int,
    n_nodes: int,
) -> tuple[np.ndarray, float, float]:
    cells = values.reshape(n_cells, n_nodes)
    radii = radius_km.reshape(n_cells, n_nodes)
    centers = np.mean(radii, axis=1)
    averages = np.mean(cells, axis=1)
    if n_nodes == 1:
        return np.zeros(n_cells), 0.0, np.nan
    slopes = (cells[:, -1] - cells[:, 0]) / (
        radii[:, -1] - radii[:, 0]
    )
    envelope = np.gradient(averages, centers)
    scale = max(
        float(np.max(np.abs(envelope[selected_cells]))),
        np.finfo(float).tiny,
    )
    active = selected_cells & (np.abs(envelope) > 1.0e-10 * scale)
    slope_weight = np.abs(slopes[active])
    backward = slopes[active] * envelope[active] < 0.0
    backward_fraction = (
        float(np.sum(slope_weight[backward]) / np.sum(slope_weight))
        if np.sum(slope_weight) > 0.0
        else 0.0
    )
    correlation = (
        float(np.corrcoef(slopes[active], envelope[active])[0, 1])
        if np.count_nonzero(active) > 1
        else np.nan
    )
    return slopes, backward_fraction, correlation


def mev_code() -> float:
    speed_of_light_mks = 2.99792458e8
    gravitational_constant_mks = 6.673e-11
    electron_volt_mks = 1.602176462e-19
    kilogram = gravitational_constant_mks / speed_of_light_mks**2
    joule = kilogram / speed_of_light_mks**2
    return 1.0e6 * electron_volt_mks * joule


def audit_energy(
    block: dict[str, np.ndarray],
    source_aware: bool,
    method: str,
    quadrature_order: int,
    spatial_operators: tuple[np.ndarray, np.ndarray] | None = None,
) -> dict[str, np.ndarray]:
    order = np.lexsort((block["spatial_node"], block["element"]))
    block = {name: values[order] for name, values in block.items()}
    n_cells = int(np.max(block["element"]))
    n_nodes = int(np.max(block["spatial_node"]))
    radius = block["r_code"]
    widths = block["dr_code"].reshape(n_cells, n_nodes)[:, 0]
    centers = np.mean(radius.reshape(n_cells, n_nodes), axis=1)
    nodes = (radius[:n_nodes] - centers[0]) / widths[0]
    mass = block["wx"] * block["dr_code"] * block["sqrtgm"]

    if spatial_operators is None:
        _, _, pressure_gradient, first_llf = assemble_isotropic_spherical_operator(
            block, source_aware
        )
    else:
        pressure_gradient, first_llf = spatial_operators
    density = block["N_before"]
    kappa = block["kappa_transport"]
    discrete_gradient = pressure_gradient @ density
    pointwise_first = discrete_gradient / kappa
    discrete_first = solve(
        np.diag(kappa) - first_llf,
        discrete_gradient,
        check_finite=False,
        assume_a="gen",
    )

    density_smooth = smooth_positive_field(radius, density, method)
    kappa_smooth = smooth_positive_field(radius, kappa, method)

    def continuum_gradient(points: np.ndarray) -> np.ndarray:
        return -density_smooth.derivative(points) / 3.0

    def continuum_first(points: np.ndarray) -> np.ndarray:
        return continuum_gradient(points) / kappa_smooth.value(points)

    projected_gradient = project_function(
        continuum_gradient, radius, widths, nodes, quadrature_order
    )
    projected_first = project_function(
        continuum_first, radius, widths, nodes, quadrature_order
    )
    sampled_first = continuum_first(radius)

    return {
        "radius": radius,
        "mass": mass,
        "density": density,
        "kappa": kappa,
        "discrete_gradient": discrete_gradient,
        "projected_gradient": projected_gradient,
        "pointwise_first": pointwise_first,
        "discrete_first": discrete_first,
        "projected_first": projected_first,
        "sampled_first": sampled_first,
        "initial_first": block["G1_before"],
        "energy": np.array([block["E_code"][0]]),
        "energy_weight": np.array([
            block["wE"][0] * block["dE_code"][0] * block["E_code"][0] ** 3
        ]),
    }


def print_comparison(
    name: str,
    numerical: np.ndarray,
    reference: np.ndarray,
    mass: np.ndarray,
    selected_nodes: np.ndarray,
    radius_km: np.ndarray,
    selected_cells: np.ndarray,
    n_cells: int,
    n_nodes: int,
) -> None:
    numerical_average, _, numerical_internal = element_decomposition(
        numerical, mass, n_cells, n_nodes
    )
    reference_average, _, reference_internal = element_decomposition(
        reference, mass, n_cells, n_nodes
    )
    cell_weights = mass.reshape(n_cells, n_nodes).sum(axis=1)
    _, backward, slope_correlation = slope_diagnostics(
        numerical, radius_km, selected_cells, n_cells, n_nodes
    )
    numerical_slopes, _, _ = slope_diagnostics(
        numerical, radius_km, selected_cells, n_cells, n_nodes
    )
    reference_slopes, reference_backward, _ = slope_diagnostics(
        reference, radius_km, selected_cells, n_cells, n_nodes
    )
    print(
        f"  {name:27s} nodal={relative_error(numerical[selected_nodes], reference[selected_nodes], mass[selected_nodes]):.6e}  "
        f"mean={relative_error(numerical_average[selected_cells], reference_average[selected_cells], cell_weights[selected_cells]):.6e}  "
        f"internal={relative_error(numerical_internal[selected_nodes], reference_internal[selected_nodes], mass[selected_nodes]):.6e}  "
        f"slope={relative_error(numerical_slopes[selected_cells], reference_slopes[selected_cells], cell_weights[selected_cells]):.6e}  "
        f"backward={backward:.3f}/{reference_backward:.3f}  corr={slope_correlation:.3f}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=Path)
    parser.add_argument("--species", type=int, default=1)
    parser.add_argument("--energy", type=int, default=25)
    parser.add_argument(
        "--method",
        choices=("pchip-linear", "pchip-log", "cubic-log"),
        default="pchip-log",
    )
    parser.add_argument("--quadrature-order", type=int, default=16)
    parser.add_argument("--radius-min-km", type=float, default=10.0)
    parser.add_argument("--radius-max-km", type=float, default=45.0)
    parser.add_argument("--table", type=Path)
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args()

    metadata, columns, values = read_table(args.file)
    data = column_map(columns, values)
    source_aware = bool(round(metadata.get("source_aware_llf", 0.0)))
    energies = np.unique(data["flat_energy"].astype(int))

    selected_block = select_block(data, args.species, args.energy)
    selected_order = np.lexsort((
        selected_block["spatial_node"], selected_block["element"]
    ))
    selected_block_ordered = {
        name: field[selected_order] for name, field in selected_block.items()
    }
    _, _, pressure_gradient, first_llf = assemble_isotropic_spherical_operator(
        selected_block_ordered, source_aware
    )
    spatial_operators = (pressure_gradient, first_llf)
    selected = audit_energy(
        selected_block,
        source_aware,
        args.method,
        args.quadrature_order,
        spatial_operators,
    )
    n_nodes = int(metadata["n_spatial_nodes"])
    n_cells = selected["radius"].size // n_nodes
    radius_km = (
        selected["radius"] / metadata["centimeter_code"] / 1.0e5
    )
    centers_km = np.mean(radius_km.reshape(n_cells, n_nodes), axis=1)
    selected_nodes = radial_mask(
        radius_km, args.radius_min_km, args.radius_max_km
    )
    selected_cells = radial_mask(
        centers_km, args.radius_min_km, args.radius_max_km
    )

    print(
        f"Projection-consistency audit: species={args.species}, "
        f"energy-node={args.energy}, method={args.method}, "
        f"quadrature={args.quadrature_order}"
    )
    print(
        f"  interval={args.radius_min_km:g}--{args.radius_max_km:g} km, "
        f"E={selected['energy'][0] / mev_code():.8e} MeV"
    )
    print_comparison(
        "B_h D vs projected -D'/3",
        selected["discrete_gradient"],
        selected["projected_gradient"],
        selected["mass"],
        selected_nodes,
        radius_km,
        selected_cells,
        n_cells,
        n_nodes,
    )
    print_comparison(
        "(B_h D)/kappa vs Fick",
        selected["pointwise_first"],
        selected["projected_first"],
        selected["mass"],
        selected_nodes,
        radius_km,
        selected_cells,
        n_cells,
        n_nodes,
    )
    print_comparison(
        "saved first moment vs Fick",
        selected["initial_first"],
        selected["projected_first"],
        selected["mass"],
        selected_nodes,
        radius_km,
        selected_cells,
        n_cells,
        n_nodes,
    )
    print_comparison(
        "(K-J_I)^-1 B_h D vs Fick",
        selected["discrete_first"],
        selected["projected_first"],
        selected["mass"],
        selected_nodes,
        radius_km,
        selected_cells,
        n_cells,
        n_nodes,
    )
    print_comparison(
        "sampled smooth Fick vs projection",
        selected["sampled_first"],
        selected["projected_first"],
        selected["mass"],
        selected_nodes,
        radius_km,
        selected_cells,
        n_cells,
        n_nodes,
    )

    energy_rows: list[list[float]] = []
    grey_initial = np.zeros_like(selected["discrete_first"])
    grey_discrete = np.zeros_like(grey_initial)
    grey_reference = np.zeros_like(grey_initial)
    energy_scale = max(data["E_code"])
    width_scale = max(data["dE_code"])
    for energy in energies:
        audit = audit_energy(
            select_block(data, args.species, int(energy)),
            source_aware,
            args.method,
            args.quadrature_order,
            spatial_operators,
        )
        weight = (
            audit["energy_weight"][0] / (energy_scale**3 * width_scale)
        )
        grey_initial += weight * audit["initial_first"]
        grey_discrete += weight * audit["discrete_first"]
        grey_reference += weight * audit["projected_first"]
        _, backward, slope_correlation = slope_diagnostics(
            audit["discrete_first"],
            radius_km,
            selected_cells,
            n_cells,
            n_nodes,
        )
        _, reference_backward, _ = slope_diagnostics(
            audit["projected_first"],
            radius_km,
            selected_cells,
            n_cells,
            n_nodes,
        )
        energy_rows.append([
            float(energy),
            audit["energy"][0] / mev_code(),
            weight,
            backward,
            reference_backward,
            slope_correlation,
            relative_error(
                audit["discrete_first"][selected_nodes],
                audit["projected_first"][selected_nodes],
                audit["mass"][selected_nodes],
            ),
            weighted_norm(
                (weight * audit["discrete_first"])[selected_nodes],
                audit["mass"][selected_nodes],
            ),
        ])

    print("Energy-weighted grey first moment:")
    print_comparison(
        "initial grey vs projected Fick",
        grey_initial,
        grey_reference,
        selected["mass"],
        selected_nodes,
        radius_km,
        selected_cells,
        n_cells,
        n_nodes,
    )
    print_comparison(
        "discrete grey vs projected Fick",
        grey_discrete,
        grey_reference,
        selected["mass"],
        selected_nodes,
        radius_km,
        selected_cells,
        n_cells,
        n_nodes,
    )

    energy_table = np.asarray(energy_rows)
    contribution = energy_table[:, 7]
    important = contribution > 1.0e-3 * np.max(contribution)
    print(
        "  spectrally important nodes with predominantly backward slopes: "
        f"{np.count_nonzero(important & (energy_table[:, 3] > 0.5))}/"
        f"{np.count_nonzero(important)}"
    )
    if args.table:
        np.savetxt(
            args.table,
            energy_table,
            header=(
                "flat_energy E_MeV normalized_grey_weight "
                "discrete_backward_fraction reference_backward_fraction "
                "discrete_slope_correlation relative_nodal_error "
                "weighted_grey_contribution_norm"
            ),
        )

    if args.plot:
        import matplotlib.pyplot as plt

        def plot_elements(axis, values, **kwargs):
            cells = values.reshape(n_cells, n_nodes)
            radii = radius_km.reshape(n_cells, n_nodes)
            label = kwargs.pop("label", None)
            for cell in np.flatnonzero(selected_cells):
                axis.plot(
                    radii[cell], cells[cell],
                    label=label if label is not None else None,
                    **kwargs,
                )
                label = None

        fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
        plot_elements(
            axes[0, 0], selected["discrete_first"],
            color="tab:blue", marker="o", label=r"$(K-J_I)^{-1}B_hD$",
        )
        plot_elements(
            axes[0, 0], selected["initial_first"],
            color="tab:orange", marker=".", label="saved first moment",
        )
        plot_elements(
            axes[0, 0], selected["projected_first"],
            color="black", label="projected smooth Fick flux",
        )
        axes[0, 0].set_title(f"spectral node {args.energy}")
        axes[0, 0].set_ylabel("first moment [code]")
        axes[0, 0].legend(loc="best")

        discrete_slopes, _, _ = slope_diagnostics(
            selected["discrete_first"], radius_km, selected_cells,
            n_cells, n_nodes,
        )
        reference_slopes, _, _ = slope_diagnostics(
            selected["projected_first"], radius_km, selected_cells,
            n_cells, n_nodes,
        )
        axes[1, 0].plot(
            centers_km[selected_cells], discrete_slopes[selected_cells],
            label="discrete",
        )
        initial_slopes, _, _ = slope_diagnostics(
            selected["initial_first"], radius_km, selected_cells,
            n_cells, n_nodes,
        )
        axes[1, 0].plot(
            centers_km[selected_cells], initial_slopes[selected_cells],
            label="saved",
        )
        axes[1, 0].plot(
            centers_km[selected_cells], reference_slopes[selected_cells],
            label="projected smooth",
        )
        axes[1, 0].set_ylabel("within-element radial slope")
        axes[1, 0].set_xlabel("radius [km]")
        axes[1, 0].legend(loc="best")

        plot_elements(
            axes[0, 1], grey_discrete,
            color="tab:blue", marker="o", label="discrete grey",
        )
        plot_elements(
            axes[0, 1], grey_initial,
            color="tab:orange", marker=".", label="saved grey",
        )
        plot_elements(
            axes[0, 1], grey_reference,
            color="black", label="projected smooth grey",
        )
        axes[0, 1].set_title("energy-weighted grey first moment")
        axes[0, 1].legend(loc="best")

        grey_slopes, _, _ = slope_diagnostics(
            grey_discrete, radius_km, selected_cells, n_cells, n_nodes
        )
        grey_reference_slopes, _, _ = slope_diagnostics(
            grey_reference, radius_km, selected_cells, n_cells, n_nodes
        )
        axes[1, 1].plot(
            centers_km[selected_cells], grey_slopes[selected_cells],
            label="discrete grey",
        )
        grey_initial_slopes, _, _ = slope_diagnostics(
            grey_initial, radius_km, selected_cells, n_cells, n_nodes
        )
        axes[1, 1].plot(
            centers_km[selected_cells], grey_initial_slopes[selected_cells],
            label="saved grey",
        )
        axes[1, 1].plot(
            centers_km[selected_cells], grey_reference_slopes[selected_cells],
            label="projected smooth grey",
        )
        axes[1, 1].set_xlabel("radius [km]")
        axes[1, 1].legend(loc="best")
        for axis in axes.flat:
            axis.set_xscale("log")
            axis.set_xlim(args.radius_min_km, args.radius_max_km)
            axis.grid(alpha=0.2)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=180)


if __name__ == "__main__":
    main()
