#!/usr/bin/env python3
"""Exact p=1 modal audit of Thornado's centered mixed-DG derivative.

The Cartesian periodic calculation derives and verifies the two Fourier
branches of the centered weak derivative.  An optional frozen-profile audit
then measures how the production density divergence acts on the
within-element first-moment component found by
``projection_consistency_audit.py``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from analyze_frozen_transient_operator import (
    assemble_isotropic_spherical_operator,
    column_map,
    read_table,
    select_block,
)
from frozen_matter_schur_audit import build_operators
from projection_consistency_audit import (
    audit_energy,
    element_decomposition,
    radial_mask,
    relative_error,
    weighted_norm,
)


def analytic_centered_symbol(theta: float, width: float) -> np.ndarray:
    """Symbol of the negative centered DG derivative in [mean,slope]."""

    sine = np.sin(theta)
    cosine = np.cos(theta)
    return np.array(
        [
            [-1.0j * sine, cosine - 1.0],
            [3.0 * (1.0 - cosine), 3.0j * sine],
        ],
        dtype=complex,
    ) / width


def modal_symbol(
    operator: np.ndarray, theta: float, n_cells: int
) -> np.ndarray:
    """Extract a two-by-two Fourier symbol in the p=1 modal basis."""

    gauss = np.array([-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)])
    modal_to_nodal = np.column_stack((np.ones(2), gauss))
    nodal_to_modal = np.linalg.inv(modal_to_nodal)
    phase = np.exp(1.0j * theta * np.arange(n_cells))
    symbol = np.empty((2, 2), dtype=complex)
    for column in range(2):
        modal = np.zeros((n_cells, 2), dtype=complex)
        modal[:, column] = phase
        nodal = (modal @ modal_to_nodal.T).reshape(-1)
        result_nodal = (operator @ nodal).reshape(n_cells, 2)
        result_modal = result_nodal @ nodal_to_modal.T
        symbol[:, column] = np.mean(
            result_modal * phase.conj()[:, None], axis=0
        )
    return symbol


def mass_internal_fraction(vector: np.ndarray) -> float:
    """Internal modal fraction for mass diag(1,1/3)."""

    norm2 = abs(vector[0]) ** 2 + abs(vector[1]) ** 2 / 3.0
    return float(
        np.sqrt((abs(vector[1]) ** 2 / 3.0) / max(norm2, np.finfo(float).tiny))
    )


def branch_data(theta_values: np.ndarray, width: float) -> dict[str, np.ndarray]:
    physical = np.empty(theta_values.size, dtype=complex)
    internal = np.empty(theta_values.size, dtype=complex)
    physical_fraction = np.empty(theta_values.size)
    internal_fraction = np.empty(theta_values.size)
    for index, theta in enumerate(theta_values):
        symbol = analytic_centered_symbol(theta, width)
        eigenvalues, eigenvectors = np.linalg.eig(symbol)
        target = -1.0j * theta / width
        i_physical = int(np.argmin(np.abs(eigenvalues - target)))
        i_internal = 1 - i_physical
        physical[index] = eigenvalues[i_physical]
        internal[index] = eigenvalues[i_internal]
        physical_fraction[index] = mass_internal_fraction(
            eigenvectors[:, i_physical]
        )
        internal_fraction[index] = mass_internal_fraction(
            eigenvectors[:, i_internal]
        )
    return {
        "physical": physical,
        "internal": internal,
        "physical_fraction": physical_fraction,
        "internal_fraction": internal_fraction,
    }


def p1_symbol_audit(n_cells: int) -> dict[str, object]:
    operators_centered = build_operators(
        n_cells,
        2,
        mean_optical_depth=1.0,
        opacity_contrast=0.0,
        source_aware=False,
    )
    width = 1.0 / n_cells
    maximum_a_error = 0.0
    maximum_b_error = 0.0
    maximum_ab_error = 0.0
    for mode in range(n_cells):
        theta = 2.0 * np.pi * mode / n_cells
        centered = analytic_centered_symbol(theta, width)
        numerical_a = modal_symbol(
            operators_centered.density_from_flux, theta, n_cells
        )
        numerical_b = modal_symbol(
            operators_centered.flux_from_density, theta, n_cells
        )
        maximum_a_error = max(
            maximum_a_error,
            float(np.linalg.norm(width * (numerical_a - centered))),
        )
        maximum_b_error = max(
            maximum_b_error,
            float(np.linalg.norm(width * (numerical_b - centered / 3.0))),
        )
        maximum_ab_error = max(
            maximum_ab_error,
            float(np.linalg.norm(
                width**2 * (numerical_a @ numerical_b - centered @ centered / 3.0)
            )),
        )

    gauss = np.array([-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)])
    sawtooth = np.tile(gauss, n_cells)
    constant = np.ones(2 * n_cells)
    norm = lambda vector: float(np.linalg.norm(vector, np.inf))

    operators_thick = build_operators(
        n_cells,
        2,
        mean_optical_depth=1.0e5,
        opacity_contrast=0.0,
        source_aware=True,
    )
    return {
        "width": width,
        "operators_centered": operators_centered,
        "operators_thick": operators_thick,
        "maximum_a_error": maximum_a_error,
        "maximum_b_error": maximum_b_error,
        "maximum_ab_error": maximum_ab_error,
        "a_constant": norm(operators_centered.density_from_flux @ constant),
        "b_constant": norm(operators_centered.flux_from_density @ constant),
        "a_sawtooth": norm(operators_centered.density_from_flux @ sawtooth),
        "b_sawtooth": norm(operators_centered.flux_from_density @ sawtooth),
        "centered_schur_sawtooth": norm(
            (
                operators_centered.density_from_flux
                @ operators_centered.flux_from_density
            )
            @ sawtooth
        ),
        "penalty_sawtooth": norm(
            operators_thick.density_penalty @ sawtooth
        ),
        "thick_schur_sawtooth": norm(operators_thick.schur @ sawtooth),
    }


def projector(
    mass: np.ndarray, n_cells: int, n_nodes: int
) -> tuple[np.ndarray, np.ndarray]:
    mean = np.zeros((mass.size, mass.size))
    for cell in range(n_cells):
        local = slice(cell * n_nodes, (cell + 1) * n_nodes)
        weights = mass[local]
        mean[local, local] = np.outer(
            np.ones(n_nodes), weights / np.sum(weights)
        )
    return mean, np.eye(mass.size) - mean


def print_response_comparison(
    label: str,
    numerical: np.ndarray,
    reference: np.ndarray,
    mass: np.ndarray,
    selected_nodes: np.ndarray,
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
    cell_mass = mass.reshape(n_cells, n_nodes).sum(axis=1)
    print(
        f"  {label:25s} nodal={relative_error(numerical[selected_nodes], reference[selected_nodes], mass[selected_nodes]):.6e}  "
        f"mean={relative_error(numerical_average[selected_cells], reference_average[selected_cells], cell_mass[selected_cells]):.6e}  "
        f"internal={relative_error(numerical_internal[selected_nodes], reference_internal[selected_nodes], mass[selected_nodes]):.6e}"
    )


def frozen_profile_audit(
    path: Path,
    species: int,
    energy: int,
    radius_min_km: float,
    radius_max_km: float,
) -> dict[str, np.ndarray]:
    metadata, columns, table = read_table(path)
    data = column_map(columns, table)
    block = select_block(data, species, energy)
    order = np.lexsort((block["spatial_node"], block["element"]))
    block = {name: field[order] for name, field in block.items()}
    source_aware = bool(round(metadata.get("source_aware_llf", 0.0)))
    a_dd, a_di, _, _ = assemble_isotropic_spherical_operator(
        block, source_aware
    )
    audit = audit_energy(
        block, source_aware, "pchip-log", 16
    )
    n_nodes = int(metadata["n_spatial_nodes"])
    n_cells = audit["radius"].size // n_nodes
    radius_km = audit["radius"] / metadata["centimeter_code"] / 1.0e5
    centers_km = np.mean(radius_km.reshape(n_cells, n_nodes), axis=1)
    selected_nodes = radial_mask(radius_km, radius_min_km, radius_max_km)
    selected_cells = radial_mask(centers_km, radius_min_km, radius_max_km)
    mean, internal = projector(audit["mass"], n_cells, n_nodes)

    discrete_first = audit["discrete_first"]
    reference_first = audit["projected_first"]
    first_error = discrete_first - reference_first
    divergence_discrete = a_di @ discrete_first
    divergence_reference = a_di @ reference_first
    divergence_error = a_di @ first_error
    divergence_from_mean = a_di @ (mean @ discrete_first)
    divergence_from_internal = a_di @ (internal @ discrete_first)
    density_penalty = a_dd @ audit["density"]
    spatial_discrete = density_penalty + divergence_discrete
    spatial_reference = density_penalty + divergence_reference

    print(
        f"Frozen-profile action of A_h: species={species}, energy={energy}, "
        f"interval={radius_min_km:g}--{radius_max_km:g} km"
    )
    print_response_comparison(
        "A_h I_h vs A_h I_proj",
        divergence_discrete,
        divergence_reference,
        audit["mass"],
        selected_nodes,
        selected_cells,
        n_cells,
        n_nodes,
    )
    print_response_comparison(
        "J_D D+A_h I_h",
        spatial_discrete,
        spatial_reference,
        audit["mass"],
        selected_nodes,
        selected_cells,
        n_cells,
        n_nodes,
    )
    response_scale = max(
        weighted_norm(divergence_discrete[selected_nodes], audit["mass"][selected_nodes]),
        np.finfo(float).tiny,
    )
    print(
        "  response norms / ||A_h I_h||: "
        f"A_h(mean I)={weighted_norm(divergence_from_mean[selected_nodes], audit['mass'][selected_nodes]) / response_scale:.6e}, "
        f"A_h(internal I)={weighted_norm(divergence_from_internal[selected_nodes], audit['mass'][selected_nodes]) / response_scale:.6e}, "
        f"A_h(delta I)={weighted_norm(divergence_error[selected_nodes], audit['mass'][selected_nodes]) / response_scale:.6e}, "
        f"J_D D={weighted_norm(density_penalty[selected_nodes], audit['mass'][selected_nodes]) / response_scale:.6e}"
    )
    return {
        "radius_km": radius_km,
        "selected_nodes": selected_nodes,
        "discrete_first": discrete_first,
        "reference_first": reference_first,
        "divergence_discrete": divergence_discrete,
        "divergence_reference": divergence_reference,
        "divergence_from_mean": divergence_from_mean,
        "divergence_from_internal": divergence_from_internal,
        "density_penalty": density_penalty,
    }


def frozen_profile_grey_audit(
    path: Path,
    species: int,
    radius_min_km: float,
    radius_max_km: float,
) -> None:
    """Apply the same divergence audit after energy integration."""

    metadata, columns, table = read_table(path)
    data = column_map(columns, table)
    source_aware = bool(round(metadata.get("source_aware_llf", 0.0)))
    energies = np.unique(data["flat_energy"].astype(int))
    first_block = select_block(data, species, int(energies[0]))
    order = np.lexsort((first_block["spatial_node"], first_block["element"]))
    first_block = {name: field[order] for name, field in first_block.items()}
    _, divergence, pressure_gradient, first_llf = (
        assemble_isotropic_spherical_operator(first_block, source_aware)
    )
    first_audit = audit_energy(
        first_block,
        source_aware,
        "pchip-log",
        16,
        (pressure_gradient, first_llf),
    )
    grey_discrete = np.zeros_like(first_audit["discrete_first"])
    grey_reference = np.zeros_like(grey_discrete)
    energy_scale = max(data["E_code"])
    width_scale = max(data["dE_code"])
    for energy in energies:
        audit = audit_energy(
            select_block(data, species, int(energy)),
            source_aware,
            "pchip-log",
            16,
            (pressure_gradient, first_llf),
        )
        weight = audit["energy_weight"][0] / (energy_scale**3 * width_scale)
        grey_discrete += weight * audit["discrete_first"]
        grey_reference += weight * audit["projected_first"]

    n_nodes = int(metadata["n_spatial_nodes"])
    n_cells = grey_discrete.size // n_nodes
    radius_km = (
        first_audit["radius"] / metadata["centimeter_code"] / 1.0e5
    )
    centers_km = np.mean(radius_km.reshape(n_cells, n_nodes), axis=1)
    selected_nodes = radial_mask(radius_km, radius_min_km, radius_max_km)
    selected_cells = radial_mask(centers_km, radius_min_km, radius_max_km)
    mass = first_audit["mass"]
    print(f"Energy-integrated action of A_h: species={species}")
    print_response_comparison(
        "A_h I_grey vs reference",
        divergence @ grey_discrete,
        divergence @ grey_reference,
        mass,
        selected_nodes,
        selected_cells,
        n_cells,
        n_nodes,
    )
    mean, internal = projector(mass, n_cells, n_nodes)
    response_scale = max(
        weighted_norm((divergence @ grey_discrete)[selected_nodes], mass[selected_nodes]),
        np.finfo(float).tiny,
    )
    print(
        "  response norms / ||A_h I_grey||: "
        f"A_h(mean I)={weighted_norm((divergence @ (mean @ grey_discrete))[selected_nodes], mass[selected_nodes]) / response_scale:.6e}, "
        f"A_h(internal I)={weighted_norm((divergence @ (internal @ grey_discrete))[selected_nodes], mass[selected_nodes]) / response_scale:.6e}, "
        f"A_h(delta I)={weighted_norm((divergence @ (grey_discrete - grey_reference))[selected_nodes], mass[selected_nodes]) / response_scale:.6e}"
    )


def write_figure(
    output: Path,
    n_cells: int,
    frozen: dict[str, np.ndarray] | None,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    width = 1.0 / n_cells
    theta = np.linspace(1.0e-6, np.pi, 512)
    branches = branch_data(theta, width)
    n_columns = 2 if frozen is not None else 1
    fig, axes = plt.subplots(2, n_columns, figsize=(7 * n_columns, 8))
    if n_columns == 1:
        axes = np.asarray(axes)[:, None]

    axes[0, 0].plot(
        theta, width * branches["physical"].imag, label="physical branch"
    )
    axes[0, 0].plot(
        theta, width * branches["internal"].imag, label="internal branch"
    )
    axes[0, 0].plot(
        theta, -theta, "k--", label=r"continuous derivative $-\varphi$"
    )
    axes[0, 0].set_ylabel(
        r"$h\,\mathrm{Im}(\lambda)$ of $\widehat C_h$"
    )
    axes[0, 0].set_xlabel(r"element wavenumber $\varphi$")
    axes[0, 0].set_title(
        r"(a) Eigenvalues of centered derivative $C_h$"
    )
    axes[0, 0].legend(loc="best")

    axes[1, 0].plot(
        theta, branches["physical_fraction"], label="physical branch"
    )
    axes[1, 0].plot(
        theta, branches["internal_fraction"], label="internal branch"
    )
    axes[1, 0].set_ylabel(
        r"fraction of eigenvector norm in slope $b_j$"
    )
    axes[1, 0].set_xlabel(r"element wavenumber $\varphi$")
    axes[1, 0].set_title(
        r"(b) Element-internal content of each branch"
    )
    axes[1, 0].legend(loc="best")

    if frozen is not None:
        selected = frozen["selected_nodes"]
        radius = frozen["radius_km"][selected]
        axes[0, 1].plot(
            radius,
            frozen["discrete_first"][selected],
            ".-",
            label=r"$(K-J_I)^{-1}B_hD_h$",
        )
        axes[0, 1].plot(
            radius,
            frozen["reference_first"][selected],
            "-",
            label=r"$\Pi_h[-(\partial_rD)/(3\kappa)]$",
        )
        axes[0, 1].set_title(
            r"(c) Frozen $\nu_e$, $E=119.19$ MeV first moment"
        )
        axes[0, 1].set_ylabel(r"$I_h^r$ [code units]")
        axes[0, 1].set_xlabel("radius [km]")
        axes[0, 1].legend(loc="best")

        axes[1, 1].plot(
            radius,
            frozen["divergence_from_mean"][selected],
            label=r"$A_h(\Pi_0 I_h)$",
        )
        axes[1, 1].plot(
            radius,
            frozen["divergence_from_internal"][selected],
            label=r"$A_h[(1-\Pi_0)I_h]$",
        )
        axes[1, 1].plot(
            radius,
            frozen["divergence_discrete"][selected],
            "k-",
            label=r"$A_hI_h$",
        )
        axes[1, 1].set_title(
            r"(d) Decomposition of discrete divergence $A_hI_h$"
        )
        axes[1, 1].set_ylabel(r"contribution to $\partial_tD_h$ [code units]")
        axes[1, 1].set_xlabel("radius [km]")
        axes[1, 1].legend(loc="best")

    for axis in axes.flat:
        axis.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output, dpi=180)


def write_trace_null_mode_figure(output: Path) -> None:
    """Show explicitly why the repeated p=1 slope is invisible to C_h."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_cells = 4
    width = 1.0
    gauss = np.array([-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)])
    repeated_slope = np.tile(gauss, n_cells)
    operators = build_operators(
        n_cells,
        2,
        mean_optical_depth=1.0,
        opacity_contrast=0.0,
        source_aware=False,
    )
    discrete_derivative = (
        operators.density_from_flux @ repeated_slope
    ).reshape(n_cells, 2)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.3))
    local_xi = np.linspace(-1.0, 1.0, 64)
    for cell in range(n_cells):
        physical_x = cell + 0.5 * (local_xi + 1.0)
        axes[0].plot(physical_x, local_xi, color="C0", linewidth=2)
    for face in range(1, n_cells):
        axes[0].plot(face, 1.0, "o", color="C1", label=None)
        axes[0].plot(face, -1.0, "o", color="C2", label=None)
        axes[0].plot(face, 0.0, "ks", markersize=5, label=None)
        axes[0].vlines(face, -1.0, 1.0, colors="0.75", linestyles=":")
    axes[0].plot([], [], "o", color="C1", label=r"left trace $q^-=+1$")
    axes[0].plot([], [], "o", color="C2", label=r"right trace $q^+=-1$")
    axes[0].plot([], [], "ks", markersize=5,
                 label=r"centered trace $\widehat q=0$")
    axes[0].set_xlim(0.0, float(n_cells))
    axes[0].set_ylim(-1.2, 1.2)
    axes[0].set_xlabel(r"$x/h$")
    axes[0].set_ylabel(r"$q_h$")
    axes[0].set_title(r"(a) Repeated $p=1$ slope and its face traces")
    axes[0].legend(loc="upper right")

    node_x = np.concatenate([
        cell + 0.5 * (gauss + 1.0) for cell in range(n_cells)
    ])
    axes[1].axhline(
        -2.0 / width,
        color="C1",
        linestyle="--",
        label=r"elementwise $-\partial_x q_h=-2/h$",
    )
    axes[1].plot(
        node_x,
        discrete_derivative.reshape(-1),
        "o-",
        color="C0",
        label=r"centered weak derivative $C_hq_h$",
    )
    for face in range(1, n_cells):
        axes[1].axvline(face, color="0.75", linestyle=":")
    axes[1].set_xlim(0.0, float(n_cells))
    axes[1].set_xlabel(r"$x/h$")
    axes[1].set_ylabel("derivative")
    axes[1].set_title(r"(b) The centered weak derivative is blind to the mode")
    axes[1].legend(loc="best")

    for axis in axes:
        axis.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output, dpi=180)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", type=int, default=64)
    parser.add_argument("--frozen-file", type=Path)
    parser.add_argument("--species", type=int, default=1)
    parser.add_argument("--energy", type=int, default=25)
    parser.add_argument("--radius-min-km", type=float, default=10.0)
    parser.add_argument("--radius-max-km", type=float, default=45.0)
    parser.add_argument("--plot", type=Path)
    parser.add_argument("--trace-plot", type=Path)
    args = parser.parse_args()

    audit = p1_symbol_audit(args.cells)
    print("Exact centered p=1 modal-symbol audit")
    print(
        f"  max dimensionless symbol errors: A={audit['maximum_a_error']:.6e}, "
        f"B={audit['maximum_b_error']:.6e}, A*B={audit['maximum_ab_error']:.6e}"
    )
    print(
        "  constant mode residuals: "
        f"||A 1||_inf={audit['a_constant']:.6e}, "
        f"||B 1||_inf={audit['b_constant']:.6e}"
    )
    print(
        "  repeated within-element slope residuals: "
        f"||A s||_inf={audit['a_sawtooth']:.6e}, "
        f"||B s||_inf={audit['b_sawtooth']:.6e}, "
        f"||A B s||_inf={audit['centered_schur_sawtooth']:.6e}"
    )
    print(
        "  source-aware thick operator on repeated slope: "
        f"||J_D s||_inf={audit['penalty_sawtooth']:.6e}, "
        f"||L_kappa s||_inf={audit['thick_schur_sawtooth']:.6e}"
    )
    theta = 2.0 * np.pi / args.cells
    branch = branch_data(np.array([theta]), audit["width"])
    print(
        "  lowest nonzero element wavenumber: "
        f"h Im(lambda_phys)={audit['width'] * branch['physical'][0].imag:.6e}, "
        f"h Im(lambda_int)={audit['width'] * branch['internal'][0].imag:.6e}, "
        f"internal fractions={branch['physical_fraction'][0]:.6e}/"
        f"{branch['internal_fraction'][0]:.6e}"
    )

    frozen = None
    if args.frozen_file:
        frozen = frozen_profile_audit(
            args.frozen_file,
            args.species,
            args.energy,
            args.radius_min_km,
            args.radius_max_km,
        )
        frozen_profile_grey_audit(
            args.frozen_file,
            args.species,
            args.radius_min_km,
            args.radius_max_km,
        )
    if args.plot:
        write_figure(args.plot, args.cells, frozen)
    if args.trace_plot:
        write_trace_null_mode_figure(args.trace_plot)


if __name__ == "__main__":
    main()
