#!/usr/bin/env python3
"""Causal modal audit of one frozen Thornado radiation step.

The audit uses the static, isotropic specialization of the production radial
DG operator and the IMEX tableau recorded with ``*.frozen_transient``.  It
does not propose a modified evolution scheme.  Instead, it repeats the saved
frozen-coefficient step after removing selected blocks of the spatial
operator.  This distinguishes the density LLF term from the pathway

    density -> nonconstant first moment -> density -> EmAb exchange.

For two spatial nodes, ``P1`` below denotes the complete nonconstant spatial
mode in each element.  For higher order it denotes all modes orthogonal to
the element-constant mode.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_frozen_transient_operator import (
    assemble_isotropic_spherical_operator,
    column_map,
    lagrange_values,
    parse_tableau,
    read_table,
    select_block,
)


@dataclass
class FrozenStep:
    density: np.ndarray
    first: np.ndarray
    stages_density: list[np.ndarray]
    stages_first: list[np.ndarray]
    matter_increment: np.ndarray


def apply_matrix(matrix: np.ndarray, values: np.ndarray) -> np.ndarray:
    if matrix.ndim == 2:
        return matrix @ values
    if matrix.ndim == 3 and matrix.shape[0] == values.shape[1]:
        return np.einsum("bij,jb->ib", matrix, values, optimize=True)
    raise ValueError("operator has an incompatible shape")


def assemble_density_llf_matrix(
    block: dict[str, np.ndarray], source_aware: bool
) -> np.ndarray:
    """Assemble only the density-jump contribution to the density equation."""

    element = block["element"].astype(int)
    node = block["spatial_node"].astype(int)
    order = np.lexsort((node, element))
    n_cells = int(element.max())
    n_nodes = int(node.max())

    def shaped(name: str) -> np.ndarray:
        return block[name][order].reshape(n_cells, n_nodes)

    radius = shaped("r_code")
    widths = shaped("dr_code")[:, 0]
    weights = shaped("wx")[0]
    centers = np.mean(radius, axis=1)
    nodes = (radius[0] - centers[0]) / widths[0]
    lower = lagrange_values(nodes, -0.5)
    upper = lagrange_values(nodes, +0.5)
    sqrtgm = shaped("sqrtgm")
    h1 = shaped("h1")
    h2 = shaped("h2")
    h3 = shaped("h3")

    h1_face = np.empty(n_cells + 1)
    h2_face = np.empty(n_cells + 1)
    h3_face = np.empty(n_cells + 1)
    h1_face[0], h2_face[0], h3_face[0] = (
        h1[0] @ lower,
        h2[0] @ lower,
        h3[0] @ lower,
    )
    h1_face[-1], h2_face[-1], h3_face[-1] = (
        h1[-1] @ upper,
        h2[-1] @ upper,
        h3[-1] @ upper,
    )
    for face in range(1, n_cells):
        h1_face[face] = 0.5 * (h1[face - 1] @ upper + h1[face] @ lower)
        h2_face[face] = 0.5 * (h2[face - 1] @ upper + h2[face] @ lower)
        h3_face[face] = 0.5 * (h3[face - 1] @ upper + h3[face] @ lower)
    sqrtgm_face = h1_face * h2_face * h3_face

    theta = np.ones(n_cells + 1)
    if source_aware:
        opacity = shaped("kappa_transport")
        opacity_lower = opacity @ lower
        opacity_upper = opacity @ upper
        boundary_opacity = 0.5 * (
            max(0.0, opacity_upper[0]) + max(0.0, opacity_lower[0])
        )
        theta[0] = 1.0 / (
            1.0
            + 1.5 * abs(h1_face[0]) ** 2 * boundary_opacity * widths[0]
        )
        for face in range(1, n_cells):
            face_width = 0.5 * (widths[face - 1] + widths[face])
            opacity_face = (
                widths[face - 1] * max(0.0, opacity_upper[face - 1])
                + widths[face] * max(0.0, opacity_lower[face])
            ) / (2.0 * face_width)
            normal = abs(h1_face[face])
            theta[face] = 1.0 / (
                1.0 + 1.5 * normal**2 * opacity_face * face_width
            )
        boundary_opacity = 0.5 * (
            max(0.0, opacity_upper[-1]) + max(0.0, opacity_lower[-1])
        )
        theta[-1] = 1.0 / (
            1.0
            + 1.5 * abs(h1_face[-1]) ** 2 * boundary_opacity * widths[-1]
        )

    size = n_cells * n_nodes
    matrix = np.zeros((size, size))

    def add_face(
        face: int,
        trace_blocks: tuple[tuple[int, np.ndarray], ...],
        test_blocks: tuple[tuple[int, np.ndarray], ...],
    ) -> None:
        for test_cell, test in test_blocks:
            rows = slice(test_cell * n_nodes, (test_cell + 1) * n_nodes)
            test_factor = (
                test
                * sqrtgm_face[face]
                / (weights * widths[test_cell] * sqrtgm[test_cell])
            )
            for trace_cell, trace in trace_blocks:
                columns = slice(
                    trace_cell * n_nodes, (trace_cell + 1) * n_nodes
                )
                matrix[rows, columns] += np.outer(
                    test_factor, -0.5 * theta[face] * trace
                )

    add_face(
        0,
        ((0, lower - upper),),
        ((0, lower),),
    )
    for face in range(1, n_cells):
        add_face(
            face,
            ((face - 1, -upper), (face, lower)),
            ((face - 1, -upper), (face, lower)),
        )
    add_face(
        n_cells,
        ((n_cells - 1, lower - upper),),
        ((n_cells - 1, -upper),),
    )
    return matrix


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path)
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("modal_causality"),
    )
    parser.add_argument("--radius-min-km", type=float, default=7.0)
    parser.add_argument("--radius-max-km", type=float, default=100.0)
    parser.add_argument("--min-optical-depth", type=float, default=10.0)
    return parser.parse_args()


def relative_norm(candidate: np.ndarray, reference: np.ndarray) -> float:
    return float(
        np.linalg.norm(candidate - reference)
        / max(np.linalg.norm(reference), np.finfo(float).tiny)
    )


def weighted_norm(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sqrt(max(np.sum(weights * values**2), 0.0)))


def weighted_block_norm(
    values: np.ndarray,
    spatial_weights: np.ndarray,
    block_weights: np.ndarray,
) -> float:
    weighted = values * block_weights[None, :]
    return float(
        np.sqrt(
            max(np.sum(spatial_weights[:, None] * weighted**2), 0.0)
        )
    )


def weighted_block_cosine(
    left: np.ndarray,
    right: np.ndarray,
    spatial_weights: np.ndarray,
    block_weights: np.ndarray,
) -> float:
    scale = np.sqrt(spatial_weights)[:, None] * block_weights[None, :]
    left_weighted = (left * scale).ravel()
    right_weighted = (right * scale).ravel()
    denominator = np.linalg.norm(left_weighted) * np.linalg.norm(right_weighted)
    if denominator == 0.0:
        return float("nan")
    return float(np.dot(left_weighted, right_weighted) / denominator)


def correlation(left: np.ndarray, right: np.ndarray) -> float:
    left = left - np.mean(left)
    right = right - np.mean(right)
    denominator = np.linalg.norm(left) * np.linalg.norm(right)
    if denominator == 0.0:
        return float("nan")
    return float(np.dot(left, right) / denominator)


def best_fit_factor(candidate: np.ndarray, reference: np.ndarray) -> float:
    denominator = float(np.dot(candidate, candidate))
    if denominator == 0.0:
        return 0.0
    return float(np.dot(candidate, reference)) / denominator


def modal_projectors(block: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    element = block["element"].astype(int)
    node = block["spatial_node"].astype(int)
    order = np.lexsort((node, element))
    n_cells = int(element.max())
    n_nodes = int(node.max())
    weights = (
        block["wx"][order]
        * block["dr_code"][order]
        * block["sqrtgm"][order]
    ).reshape(n_cells, n_nodes)
    p0 = np.zeros((n_cells * n_nodes, n_cells * n_nodes))
    for cell in range(n_cells):
        indices = slice(cell * n_nodes, (cell + 1) * n_nodes)
        normalized = weights[cell] / np.sum(weights[cell])
        p0[indices, indices] = normalized[None, :]
    return p0, np.eye(n_cells * n_nodes) - p0


def advance_frozen_batch(
    initial_density: np.ndarray,
    initial_first: np.ndarray,
    equilibrium_density: np.ndarray,
    chi: np.ndarray,
    kappa: np.ndarray,
    a_dd: np.ndarray,
    a_di: np.ndarray,
    a_id: np.ndarray,
    a_ii: np.ndarray,
    dt: float,
    tableau: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
) -> FrozenStep:
    """Advance independent energy/species columns with frozen coefficients."""

    a_ex, a_im, w_ex, w_im = tableau
    stages_density: list[np.ndarray] = []
    stages_first: list[np.ndarray] = []
    for stage in range(a_ex.shape[0]):
        rhs_density = initial_density.copy()
        rhs_first = initial_first.copy()
        for previous in range(stage):
            density_previous = stages_density[previous]
            first_previous = stages_first[previous]
            stream_density = (
                apply_matrix(a_dd, density_previous) + a_di @ first_previous
            )
            stream_first = (
                a_id @ density_previous + a_ii @ first_previous
            )
            collision_density = chi * (
                equilibrium_density - density_previous
            )
            collision_first = -kappa * first_previous
            rhs_density += dt * (
                a_ex[stage, previous] * stream_density
                + a_im[stage, previous] * collision_density
            )
            rhs_first += dt * (
                a_ex[stage, previous] * stream_first
                + a_im[stage, previous] * collision_first
            )

        diagonal = a_im[stage, stage]
        stages_density.append(
            (rhs_density + dt * diagonal * chi * equilibrium_density)
            / (1.0 + dt * diagonal * chi)
        )
        stages_first.append(
            rhs_first / (1.0 + dt * diagonal * kappa)
        )

    density = initial_density.copy()
    first = initial_first.copy()
    matter_increment = np.zeros_like(initial_density)
    for stage, (density_stage, first_stage) in enumerate(
        zip(stages_density, stages_first, strict=True)
    ):
        stream_density = apply_matrix(a_dd, density_stage) + a_di @ first_stage
        stream_first = a_id @ density_stage + a_ii @ first_stage
        collision_density = chi * (equilibrium_density - density_stage)
        collision_first = -kappa * first_stage
        density += dt * (
            w_ex[stage] * stream_density
            + w_im[stage] * collision_density
        )
        first += dt * (
            w_ex[stage] * stream_first
            + w_im[stage] * collision_first
        )
        matter_increment -= dt * w_im[stage] * collision_density

    return FrozenStep(
        density=density,
        first=first,
        stages_density=stages_density,
        stages_first=stages_first,
        matter_increment=matter_increment,
    )


def element_modes(
    values: np.ndarray,
    p0: np.ndarray,
    p1: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    return p0 @ values, p1 @ values


def slope_fraction(
    values: np.ndarray,
    radius_km: np.ndarray,
    selected_cells: np.ndarray,
    n_cells: int,
    n_nodes: int,
) -> tuple[float, float]:
    cells = values.reshape(n_cells, n_nodes)
    radius = radius_km.reshape(n_cells, n_nodes)
    centers = np.mean(radius, axis=1)
    averages = np.mean(cells, axis=1)
    if n_nodes == 1:
        return 0.0, float("nan")
    local_slope = (cells[:, -1] - cells[:, 0]) / (
        radius[:, -1] - radius[:, 0]
    )
    envelope_slope = np.gradient(averages, centers)
    scale = max(
        float(np.max(np.abs(envelope_slope[selected_cells]))),
        np.finfo(float).tiny,
    )
    active = selected_cells & (np.abs(envelope_slope) > 1.0e-10 * scale)
    slope_weight = np.abs(local_slope[active])
    backward = local_slope[active] * envelope_slope[active] < 0.0
    fraction = (
        float(np.sum(slope_weight[backward]) / np.sum(slope_weight))
        if np.sum(slope_weight) > 0.0
        else 0.0
    )
    return fraction, correlation(local_slope[active], envelope_slope[active])


def element_average(
    values: np.ndarray,
    spatial_weights: np.ndarray,
    n_cells: int,
    n_nodes: int,
) -> np.ndarray:
    cells = values.reshape(n_cells, n_nodes)
    weights = spatial_weights.reshape(n_cells, n_nodes)
    return np.sum(weights * cells, axis=1) / np.sum(weights, axis=1)


_ATTRIBUTION_CACHE: dict[
    tuple[Path, float], dict[tuple[str, str], np.ndarray]
] = {}


def read_attribution_increment(
    path: Path,
    time_s: float,
    quantity: str,
    component: str,
) -> np.ndarray | None:
    """Read one nodal attribution vector without loading its full history."""

    if not path.exists():
        return None
    cache_key = (path.resolve(), time_s)
    if cache_key not in _ATTRIBUTION_CACHE:
        columns: list[str] | None = None
        indices: dict[str, int] = {}
        rows_by_key: dict[tuple[str, str], list[tuple[int, int, float]]] = {}
        tolerance = 32.0 * np.finfo(float).eps * max(abs(time_s), 1.0)
        with path.open() as stream:
            for line in stream:
                if line.startswith("# cycle"):
                    columns = line[1:].split()
                    indices = {name: index for index, name in enumerate(columns)}
                    continue
                if line.startswith("#") or not line.strip():
                    continue
                if columns is None:
                    raise ValueError(f"{path}: data precede the column header")
                words = line.split()
                if abs(float(words[indices["time_s"]]) - time_s) > tolerance:
                    continue
                key = (
                    words[indices["quantity"]],
                    words[indices["component"]],
                )
                rows_by_key.setdefault(key, []).append(
                    (
                        int(words[indices["element"]]),
                        int(words[indices["node"]]),
                        float(words[indices["step_increment"]]),
                    )
                )
        cached: dict[tuple[str, str], np.ndarray] = {}
        for key, rows in rows_by_key.items():
            rows.sort(key=lambda item: (item[0], item[1]))
            cached[key] = np.array([item[2] for item in rows])
        _ATTRIBUTION_CACHE[cache_key] = cached
    return _ATTRIBUTION_CACHE[cache_key].get((quantity, component))


def read_recorded_stage_states(
    path: Path,
    block_keys: list[tuple[int, int]],
    n_energy_nodes: int,
) -> dict[int, dict[int, tuple[np.ndarray, np.ndarray]]]:
    """Read every saved state at every IMEX stage."""

    columns: list[str] | None = None
    with path.open() as stream:
        for line in stream:
            if line.startswith("# stage state element"):
                columns = line[1:].split()
                break
    if columns is None:
        raise ValueError(f"{path}: missing stage-state column header")
    raw = np.loadtxt(path, comments="#", ndmin=2)
    if raw.shape[1] != len(columns):
        raise ValueError(
            f"{path}: found {raw.shape[1]} columns, expected {len(columns)}"
        )
    table = {name: raw[:, index] for index, name in enumerate(columns)}
    stage_numbers = table["stage"].astype(int)
    state_numbers = table["state"].astype(int)
    stages: dict[int, dict[int, tuple[np.ndarray, np.ndarray]]] = {}
    for stage in np.unique(stage_numbers):
        stages[stage] = {}
        for state in np.unique(state_numbers[stage_numbers == stage]):
            density_columns: list[np.ndarray] = []
            first_columns: list[np.ndarray] = []
            for species, flat_energy in block_keys:
                energy_element = (flat_energy - 1) // n_energy_nodes + 1
                energy_node = (flat_energy - 1) % n_energy_nodes + 1
                selected = (
                    (stage_numbers == stage)
                    & (state_numbers == state)
                    & (table["species"].astype(int) == species)
                    & (table["energy_element"].astype(int) == energy_element)
                    & (table["energy_node"].astype(int) == energy_node)
                )
                order = np.lexsort(
                    (
                        table["spatial_node"][selected],
                        table["element"][selected],
                    )
                )
                density_columns.append(table["N"][selected][order])
                first_columns.append(table["G1"][selected][order])
            stages[stage][state] = (
                np.column_stack(density_columns),
                np.column_stack(first_columns),
            )
    return stages


def backward_p1_component(
    first: np.ndarray,
    p0: np.ndarray,
    p1: np.ndarray,
    radius_km: np.ndarray,
    n_cells: int,
    n_nodes: int,
) -> np.ndarray:
    """Retain p=1 values only where their slope opposes the mean envelope."""

    mean = p0 @ first
    high = p1 @ first
    radius = radius_km.reshape(n_cells, n_nodes)
    centers = np.mean(radius, axis=1)
    result = np.zeros_like(high)
    for column in range(first.shape[1]):
        mean_cells = mean[:, column].reshape(n_cells, n_nodes)[:, 0]
        high_cells = high[:, column].reshape(n_cells, n_nodes)
        envelope = np.gradient(mean_cells, centers)
        local = (high_cells[:, -1] - high_cells[:, 0]) / (
            radius[:, -1] - radius[:, 0]
        )
        backward = local * envelope < 0.0
        result[:, column] = (
            high_cells * backward[:, None]
        ).reshape(-1)
    return result


def main() -> None:
    args = parse_args()
    metadata, columns, values = read_table(args.file)
    data = column_map(columns, values)
    dt = metadata["timestep_code"]
    tableau_path = Path(
        str(args.file).replace(".frozen_transient", ".frozen_transient_stages")
    )
    tableau = parse_tableau(tableau_path)

    species_values = np.unique(data["species"].astype(int))
    energy_values = np.unique(data["flat_energy"].astype(int))
    block_keys = [
        (species, energy)
        for species in species_values
        for energy in energy_values
    ]
    blocks: list[dict[str, np.ndarray]] = []
    for species, energy in block_keys:
        block = select_block(data, species, energy)
        order = np.lexsort((block["spatial_node"], block["element"]))
        blocks.append({name: column[order] for name, column in block.items()})

    reference = blocks[0]
    n_cells = int(reference["element"].max())
    n_nodes = int(reference["spatial_node"].max())
    n_spatial = n_cells * n_nodes
    if n_nodes != 2:
        raise ValueError(
            "the causal interpretation as a p=1 audit requires two spatial nodes"
        )
    p0, p1 = modal_projectors(reference)

    source_aware = bool(round(metadata.get("source_aware_llf", 0.0)))
    a_dd, a_di, a_id, a_ii = assemble_isotropic_spherical_operator(
        reference, False
    )
    a_dd_direct = assemble_density_llf_matrix(reference, False)
    density_llf_assembly_error = relative_norm(a_dd_direct, a_dd)
    if density_llf_assembly_error > 1.0e-12:
        raise RuntimeError(
            "direct density-LLF assembly does not reproduce the production "
            f"matrix (relative error {density_llf_assembly_error:.6e})"
        )
    a_dd_source_aware = np.stack(
        [assemble_density_llf_matrix(block, True) for block in blocks],
        axis=0,
    )
    a_dd_production = a_dd_source_aware if source_aware else a_dd

    initial_density = np.column_stack([block["N_before"] for block in blocks])
    initial_first = np.column_stack([block["G1_before"] for block in blocks])
    accepted_density = np.column_stack([block["N_after"] for block in blocks])
    accepted_first = np.column_stack([block["G1_after"] for block in blocks])
    equilibrium_density = np.column_stack([block["D0"] for block in blocks])
    chi = np.column_stack([block["chi_emab"] for block in blocks])
    kappa = np.column_stack([block["kappa_transport"] for block in blocks])

    variants = {
        "baseline": (a_dd_production, a_di, a_id, a_ii),
        "without density LLF": (np.zeros_like(a_dd), a_di, a_id, a_ii),
        "without first-moment LLF": (
            a_dd,
            a_di,
            a_id,
            np.zeros_like(a_ii),
        ),
        "source-aware density LLF": (
            a_dd_source_aware,
            a_di,
            a_id,
            a_ii,
        ),
        "source-aware without p=1 feedback": (
            a_dd_source_aware,
            a_di @ p0,
            a_id,
            a_ii,
        ),
        "source-aware without p=1 pathway": (
            a_dd_source_aware,
            a_di @ p0,
            p0 @ a_id,
            a_ii,
        ),
        "without p=1 first-moment generation": (
            a_dd,
            a_di,
            p0 @ a_id,
            a_ii,
        ),
        "without p=1 first-moment feedback": (
            a_dd,
            a_di @ p0,
            a_id,
            a_ii,
        ),
        "without p=1 first-moment pathway": (
            a_dd,
            a_di @ p0,
            p0 @ a_id,
            a_ii,
        ),
        "without density LLF or p=1 feedback": (
            np.zeros_like(a_dd),
            a_di @ p0,
            a_id,
            a_ii,
        ),
        "without density LLF or p=1 pathway": (
            np.zeros_like(a_dd),
            a_di @ p0,
            p0 @ a_id,
            a_ii,
        ),
    }
    results = {
        name: advance_frozen_batch(
            initial_density,
            initial_first,
            equilibrium_density,
            chi,
            kappa,
            *operators,
            dt,
            tableau,
        )
        for name, operators in variants.items()
    }
    results["without initial p=1 first moment"] = advance_frozen_batch(
        initial_density,
        p0 @ initial_first,
        equilibrium_density,
        chi,
        kappa,
        a_dd_production,
        a_di,
        a_id,
        a_ii,
        dt,
        tableau,
    )

    baseline = results["baseline"]
    accepted = np.vstack((accepted_density, accepted_first))
    predicted = np.vstack((baseline.density, baseline.first))
    baseline_error = relative_norm(predicted, accepted)
    density_increment_error = relative_norm(
        baseline.density - initial_density,
        accepted_density - initial_density,
    )
    first_increment_error = relative_norm(
        baseline.first - initial_first,
        accepted_first - initial_first,
    )
    exported_stream_density = np.column_stack(
        [block["stream_total_N"] for block in blocks]
    )
    exported_stream_first = np.column_stack(
        [block["stream_total_G1"] for block in blocks]
    )
    matrix_stream_density = (
        apply_matrix(a_dd_production, initial_density) + a_di @ initial_first
    )
    matrix_stream_first = a_id @ initial_density + a_ii @ initial_first

    energy_reference = max(float(block["E_code"][0]) for block in blocks)
    qh_weights = np.array([
        block["wE"][0]
        * (block["dE_code"][0] / energy_reference)
        * (block["E_code"][0] / energy_reference) ** 3
        for block in blocks
    ])
    ql_weights = np.zeros(len(blocks))
    for index, ((species, _), block) in enumerate(zip(block_keys, blocks, strict=True)):
        lepton_number = 1.0 if species == 1 else -1.0 if species == 2 else 0.0
        ql_weights[index] = (
            lepton_number
            * block["wE"][0]
            * (block["dE_code"][0] / energy_reference)
            * (block["E_code"][0] / energy_reference) ** 2
        )

    exchanges: dict[str, dict[str, np.ndarray]] = {}
    grey_first: dict[str, np.ndarray] = {}
    for name, result in results.items():
        exchanges[name] = {
            "QH": result.matter_increment @ qh_weights,
            "QL": result.matter_increment @ ql_weights,
        }
        grey_first[name] = result.first @ qh_weights

    attribution_path = Path(
        str(args.file).replace(
            ".frozen_transient", ".exchange_attribution.nodes"
        )
    )
    attribution_time = (
        metadata.get("state_time_s", np.nan)
        + metadata.get("timestep_s", np.nan)
    )

    radius_km = (
        reference["r_code"] / metadata["centimeter_code"] / 1.0e5
    )
    spatial_weights = (
        reference["wx"] * reference["dr_code"] * reference["sqrtgm"]
    )
    selected_nodes = (
        (radius_km >= args.radius_min_km)
        & (radius_km <= args.radius_max_km)
    )
    radius_cells = element_average(
        radius_km, spatial_weights, n_cells, n_nodes
    )
    selected_cells = (
        (radius_cells >= args.radius_min_km)
        & (radius_cells <= args.radius_max_km)
    )

    optical_depth = np.stack([
        block["kappa_transport"]
        * block["dr_code"]
        * np.sqrt(block["gm11"])
        for block in blocks
    ], axis=1)
    thick_weight = np.sum(
        qh_weights[None, :] * (optical_depth >= args.min_optical_depth),
        axis=1,
    )
    total_weight = np.sum(qh_weights)
    thick_nodes = thick_weight >= 0.5 * total_weight
    complete_audit_cells = (
        selected_nodes.reshape(n_cells, n_nodes).all(axis=1)
        & thick_nodes.reshape(n_cells, n_nodes).all(axis=1)
    )
    audit_nodes = np.repeat(complete_audit_cells, n_nodes)
    if not np.any(audit_nodes):
        raise ValueError("no nodes satisfy the selected radius and optical-depth cuts")
    thick_block_entries = selected_nodes[:, None] & (
        optical_depth >= args.min_optical_depth
    )

    # The explicit density RHS at zero fluid velocity is exactly the sum of
    # the density LLF term and the centered divergence of the first moment.
    # Use the recorded post-implicit-limiter stage states, which are the
    # states supplied to the production explicit RHS, to split that update
    # without rerunning or linearizing the coupled matter solve.
    recorded_stage_states = read_recorded_stage_states(
        tableau_path,
        block_keys,
        int(round(metadata["n_energy_nodes"])),
    )
    _, a_im, w_ex, w_im = tableau
    stage_density_components = {
        name: np.zeros_like(initial_density)
        for name in (
            "standard density LLF",
            "source-aware density LLF",
            "centered p=0 first moment",
            "centered p=1 first moment",
            "centered backward p=1 subset",
            "centered forward p=1 subset",
        )
    }
    stage_first_components = {
        name: np.zeros_like(initial_first)
        for name in (
            "centered pressure from p=0 density",
            "centered pressure from p=1 density",
            "first-moment LLF",
        )
    }
    for stage, weight in enumerate(w_ex, start=1):
        if weight == 0.0:
            continue
        density_stage, first_stage = recorded_stage_states[stage][4]
        first_mean = p0 @ first_stage
        first_high = p1 @ first_stage
        first_backward = backward_p1_component(
            first_stage,
            p0,
            p1,
            radius_km,
            n_cells,
            n_nodes,
        )
        factor = dt * weight
        stage_density_components["standard density LLF"] += (
            factor * (a_dd @ density_stage)
        )
        stage_density_components["source-aware density LLF"] += (
            factor * apply_matrix(a_dd_source_aware, density_stage)
        )
        stage_density_components["centered p=0 first moment"] += (
            factor * (a_di @ first_mean)
        )
        stage_density_components["centered p=1 first moment"] += (
            factor * (a_di @ first_high)
        )
        stage_density_components["centered backward p=1 subset"] += (
            factor * (a_di @ first_backward)
        )
        stage_density_components["centered forward p=1 subset"] += (
            factor * (a_di @ (first_high - first_backward))
        )
        stage_first_components["centered pressure from p=0 density"] += (
            factor * (a_id @ (p0 @ density_stage))
        )
        stage_first_components["centered pressure from p=1 density"] += (
            factor * (a_id @ (p1 @ density_stage))
        )
        stage_first_components["first-moment LLF"] += (
            factor * (a_ii @ first_stage)
        )

    stage_density_components["standard total"] = (
        stage_density_components["standard density LLF"]
        + stage_density_components["centered p=0 first moment"]
        + stage_density_components["centered p=1 first moment"]
    )
    stage_density_components["source-aware total"] = (
        stage_density_components["source-aware density LLF"]
        + stage_density_components["centered p=0 first moment"]
        + stage_density_components["centered p=1 first moment"]
    )
    stage_first_components["total"] = sum(stage_first_components.values())
    stage_first_collision = np.zeros_like(initial_first)
    for stage, weight in enumerate(w_im, start=1):
        diagonal = a_im[stage - 1, stage - 1]
        if weight == 0.0 or diagonal == 0.0:
            continue
        first_before = recorded_stage_states[stage][2][1]
        first_after = recorded_stage_states[stage][3][1]
        stage_first_collision += (weight / diagonal) * (
            first_after - first_before
        )
    stage_first_balance = (
        stage_first_components["total"] + stage_first_collision
    )
    stage_exchange_components = {
        "QH": {
            name: values @ qh_weights
            for name, values in stage_density_components.items()
        },
        "QL": {
            name: values @ ql_weights
            for name, values in stage_density_components.items()
        },
    }
    stage_first_grey_components = {
        "QH": {
            name: values @ qh_weights
            for name, values in stage_first_components.items()
        },
        "QL": {
            name: values @ ql_weights
            for name, values in stage_first_components.items()
        },
    }

    actual_transport_response: dict[str, dict[str, np.ndarray]] = {}
    for quantity in ("QH", "QL"):
        no_stream = read_attribution_increment(
            attribution_path,
            attribution_time,
            quantity,
            "no_stream_shadow",
        )
        standard_total = read_attribution_increment(
            attribution_path,
            attribution_time,
            quantity,
            "standard_total",
        )
        ap_total = read_attribution_increment(
            attribution_path,
            attribution_time,
            quantity,
            "ap_total",
        )
        if no_stream is None or standard_total is None or ap_total is None:
            continue
        actual_transport_response[quantity] = {
            "standard": standard_total - no_stream,
            "source-aware": ap_total - no_stream,
            "standard minus source-aware": standard_total - ap_total,
        }

    report_lines = [
        "Frozen-step modal causality audit",
        f"  input: {args.file}",
        f"  saved time: {metadata.get('state_time_s', np.nan):.8e} s",
        f"  timestep: {metadata.get('timestep_s', np.nan):.8e} s",
        f"  baseline frozen prediction vs Thornado: {baseline_error:.6e}",
        f"  direct density-LLF matrix assembly: {density_llf_assembly_error:.6e}",
        (
            "  density-increment prediction vs Thornado: "
            f"{density_increment_error:.6e}"
        ),
        (
            "  first-moment-increment prediction vs Thornado: "
            f"{first_increment_error:.6e}"
        ),
        (
            "  isotropic matrix density action vs production: "
            f"{relative_norm(matrix_stream_density, exported_stream_density):.6e}"
        ),
        (
            "  isotropic matrix first-moment action vs production: "
            f"{relative_norm(matrix_stream_first, exported_stream_first):.6e}"
        ),
        (
            "  optically thick density action vs production: "
            f"{relative_norm(matrix_stream_density[thick_block_entries], exported_stream_density[thick_block_entries]):.6e}"
        ),
        (
            "  optically thick first-moment action vs production: "
            f"{relative_norm(matrix_stream_first[thick_block_entries], exported_stream_first[thick_block_entries]):.6e}"
        ),
        (
            "  optically thick density-increment prediction: "
            f"{relative_norm((baseline.density - initial_density)[thick_block_entries], (accepted_density - initial_density)[thick_block_entries]):.6e}"
        ),
        (
            "  optically thick first-moment-increment prediction: "
            f"{relative_norm((baseline.first - initial_first)[thick_block_entries], (accepted_first - initial_first)[thick_block_entries]):.6e}"
        ),
        (
            "  trusted audit region: "
            f"{args.radius_min_km:g} <= r <= {args.radius_max_km:g} km, "
            f"energy-weighted tau >= {args.min_optical_depth:g} at >= 50%"
        ),
        "",
    ]

    production_stage_name = (
        "source-aware total" if source_aware else "standard total"
    )
    report_lines.extend(
        [
            "Recorded-stage density-update composition:",
            f"  production density flux: {'source-aware LLF' if source_aware else 'standard LLF'}",
        ]
    )
    for quantity in ("QH", "QL"):
        components = stage_exchange_components[quantity]
        production_norm = weighted_norm(
            components[production_stage_name][audit_nodes],
            spatial_weights[audit_nodes],
        )
        p1_norm = weighted_norm(
            components["centered p=1 first moment"][audit_nodes],
            spatial_weights[audit_nodes],
        )
        density_name = (
            "source-aware density LLF"
            if source_aware
            else "standard density LLF"
        )
        report_lines.extend(
            [
                f"  {quantity}:",
                (
                    f"    density LLF norm / total: "
                    f"{weighted_norm(components[density_name][audit_nodes], spatial_weights[audit_nodes]) / max(production_norm, np.finfo(float).tiny):.6e}"
                ),
                (
                    "    centered p=0 first moment norm / total: "
                    f"{weighted_norm(components['centered p=0 first moment'][audit_nodes], spatial_weights[audit_nodes]) / max(production_norm, np.finfo(float).tiny):.6e}"
                ),
                (
                    "    centered p=1 first moment norm / total: "
                    f"{p1_norm / max(production_norm, np.finfo(float).tiny):.6e}"
                ),
                (
                    "    backward p=1 subset norm / p=1: "
                    f"{weighted_norm(components['centered backward p=1 subset'][audit_nodes], spatial_weights[audit_nodes]) / max(p1_norm, np.finfo(float).tiny):.6e}"
                ),
            ]
        )
    report_lines.append("")
    report_lines.append("Recorded-stage first-moment update composition:")
    for quantity in ("QH", "QL"):
        components = stage_first_grey_components[quantity]
        total_norm = weighted_norm(
            components["total"][audit_nodes], spatial_weights[audit_nodes]
        )
        report_lines.append(f"  {quantity} weighting:")
        for name in (
            "centered pressure from p=0 density",
            "centered pressure from p=1 density",
            "first-moment LLF",
        ):
            report_lines.append(
                f"    {name} norm / total: "
                f"{weighted_norm(components[name][audit_nodes], spatial_weights[audit_nodes]) / max(total_norm, np.finfo(float).tiny):.6e}"
            )
        collision_grey = stage_first_collision @ (
            qh_weights if quantity == "QH" else ql_weights
        )
        balance_grey = stage_first_balance @ (
            qh_weights if quantity == "QH" else ql_weights
        )
        collision_norm = weighted_norm(
            collision_grey[audit_nodes], spatial_weights[audit_nodes]
        )
        balance_norm = weighted_norm(
            balance_grey[audit_nodes], spatial_weights[audit_nodes]
        )
        report_lines.extend(
            [
                (
                    "    implicit first-moment collision / spatial: "
                    f"{collision_norm / max(total_norm, np.finfo(float).tiny):.6e}"
                ),
                (
                    "    completed first-moment increment / spatial: "
                    f"{balance_norm / max(total_norm, np.finfo(float).tiny):.6e}"
                ),
                (
                    "    collision--spatial correlation: "
                    f"{correlation(collision_grey[audit_nodes], -components['total'][audit_nodes]):.6e}"
                ),
                (
                    "    spectral collision / spatial: "
                    f"{weighted_block_norm(stage_first_collision[audit_nodes], spatial_weights[audit_nodes], qh_weights if quantity == 'QH' else ql_weights) / max(weighted_block_norm(stage_first_components['total'][audit_nodes], spatial_weights[audit_nodes], qh_weights if quantity == 'QH' else ql_weights), np.finfo(float).tiny):.6e}"
                ),
                (
                    "    spectral completed increment / spatial: "
                    f"{weighted_block_norm(stage_first_balance[audit_nodes], spatial_weights[audit_nodes], qh_weights if quantity == 'QH' else ql_weights) / max(weighted_block_norm(stage_first_components['total'][audit_nodes], spatial_weights[audit_nodes], qh_weights if quantity == 'QH' else ql_weights), np.finfo(float).tiny):.6e}"
                ),
                (
                    "    spectral collision--spatial cosine: "
                    f"{weighted_block_cosine(stage_first_collision[audit_nodes], -stage_first_components['total'][audit_nodes], spatial_weights[audit_nodes], qh_weights if quantity == 'QH' else ql_weights):.6e}"
                ),
            ]
        )
    report_lines.append("")

    for quantity in ("QH", "QL"):
        actual = read_attribution_increment(
            attribution_path,
            attribution_time,
            quantity,
            "emab_full",
        )
        if actual is None:
            continue
        frozen = exchanges["baseline"][quantity]
        fit = float(np.dot(frozen[audit_nodes], actual[audit_nodes])) / max(
            float(np.dot(frozen[audit_nodes], frozen[audit_nodes])),
            np.finfo(float).tiny,
        )
        report_lines.extend(
            [
                f"{quantity} frozen EmAb increment validation:",
                f"  best-fit unit factor: {fit:.8e}",
                (
                    "  shape error after unit fit: "
                    f"{relative_norm(fit * frozen[audit_nodes], actual[audit_nodes]):.6e}"
                ),
                (
                    "  spatial correlation: "
                    f"{correlation(frozen[audit_nodes], actual[audit_nodes]):.6e}"
                ),
                "",
            ]
        )

    for quantity in ("QH", "QL"):
        baseline_exchange = exchanges["baseline"][quantity]
        baseline_mean, baseline_high = element_modes(
            baseline_exchange, p0, p1
        )
        baseline_norm = weighted_norm(
            baseline_exchange[audit_nodes], spatial_weights[audit_nodes]
        )
        baseline_high_norm = weighted_norm(
            baseline_high[audit_nodes], spatial_weights[audit_nodes]
        )
        baseline_mean_norm = weighted_norm(
            baseline_mean[audit_nodes], spatial_weights[audit_nodes]
        )
        report_lines.extend(
            [
                f"{quantity}:",
                f"  baseline norm: {baseline_norm:.6e}",
                f"  baseline p=1 fraction: {baseline_high_norm / max(baseline_norm, np.finfo(float).tiny):.6e}",
                f"  baseline p=0 fraction: {baseline_mean_norm / max(baseline_norm, np.finfo(float).tiny):.6e}",
            ]
        )
        for name in variants:
            if name == "baseline":
                continue
            candidate = exchanges[name][quantity]
            candidate_mean, candidate_high = element_modes(candidate, p0, p1)
            difference = baseline_exchange - candidate
            difference_mean, difference_high = element_modes(difference, p0, p1)
            report_lines.extend(
                [
                    f"  {name}:",
                    (
                        "    full response retained: "
                        f"{weighted_norm(candidate[audit_nodes], spatial_weights[audit_nodes]) / max(baseline_norm, np.finfo(float).tiny):.6e}"
                    ),
                    (
                        "    p=1 response retained: "
                        f"{weighted_norm(candidate_high[audit_nodes], spatial_weights[audit_nodes]) / max(baseline_high_norm, np.finfo(float).tiny):.6e}"
                    ),
                    (
                        "    removed response / baseline: "
                        f"{weighted_norm(difference[audit_nodes], spatial_weights[audit_nodes]) / max(baseline_norm, np.finfo(float).tiny):.6e}"
                    ),
                    (
                        "    removed p=0 / baseline p=0: "
                        f"{weighted_norm(difference_mean[audit_nodes], spatial_weights[audit_nodes]) / max(baseline_mean_norm, np.finfo(float).tiny):.6e}"
                    ),
                    (
                        "    removed p=1 / baseline p=1: "
                        f"{weighted_norm(difference_high[audit_nodes], spatial_weights[audit_nodes]) / max(baseline_high_norm, np.finfo(float).tiny):.6e}"
                    ),
                    (
                        "    corr(removed, baseline): "
                        f"{correlation(difference[audit_nodes], baseline_exchange[audit_nodes]):.6e}"
                    ),
                ]
            )
        report_lines.append("")

        no_density_llf = exchanges["without density LLF"][quantity]
        no_density_llf_norm = weighted_norm(
            no_density_llf[audit_nodes], spatial_weights[audit_nodes]
        )
        for name in (
            "without density LLF or p=1 feedback",
            "without density LLF or p=1 pathway",
        ):
            candidate = exchanges[name][quantity]
            removed = no_density_llf - candidate
            report_lines.extend(
                [
                    f"  {name}, relative to the no-density-LLF response:",
                    (
                        "    full response retained: "
                        f"{weighted_norm(candidate[audit_nodes], spatial_weights[audit_nodes]) / max(no_density_llf_norm, np.finfo(float).tiny):.6e}"
                    ),
                    (
                        "    removed response / no-density-LLF response: "
                        f"{weighted_norm(removed[audit_nodes], spatial_weights[audit_nodes]) / max(no_density_llf_norm, np.finfo(float).tiny):.6e}"
                    ),
                    (
                        "    corr(removed, no-density-LLF response): "
                        f"{correlation(removed[audit_nodes], no_density_llf[audit_nodes]):.6e}"
                    ),
                ]
            )
        initial_p1_candidate = exchanges["without initial p=1 first moment"][quantity]
        initial_p1_removed = baseline_exchange - initial_p1_candidate
        report_lines.extend(
            [
                "  without initial p=1 first moment:",
                (
                    "    full response retained: "
                    f"{weighted_norm(initial_p1_candidate[audit_nodes], spatial_weights[audit_nodes]) / max(baseline_norm, np.finfo(float).tiny):.6e}"
                ),
                (
                    "    removed response / baseline: "
                    f"{weighted_norm(initial_p1_removed[audit_nodes], spatial_weights[audit_nodes]) / max(baseline_norm, np.finfo(float).tiny):.6e}"
                ),
                (
                    "    corr(removed, baseline): "
                    f"{correlation(initial_p1_removed[audit_nodes], baseline_exchange[audit_nodes]):.6e}"
                ),
                "",
            ]
        )
        source_aware = exchanges["source-aware density LLF"][quantity]
        source_aware_norm = weighted_norm(
            source_aware[audit_nodes], spatial_weights[audit_nodes]
        )
        for name in (
            "source-aware without p=1 feedback",
            "source-aware without p=1 pathway",
        ):
            candidate = exchanges[name][quantity]
            removed = source_aware - candidate
            report_lines.extend(
                [
                    f"  {name}, relative to the source-aware response:",
                    (
                        "    full response retained: "
                        f"{weighted_norm(candidate[audit_nodes], spatial_weights[audit_nodes]) / max(source_aware_norm, np.finfo(float).tiny):.6e}"
                    ),
                    (
                        "    removed response / source-aware response: "
                        f"{weighted_norm(removed[audit_nodes], spatial_weights[audit_nodes]) / max(source_aware_norm, np.finfo(float).tiny):.6e}"
                    ),
                    (
                        "    corr(removed, source-aware response): "
                        f"{correlation(removed[audit_nodes], source_aware[audit_nodes]):.6e}"
                    ),
                ]
            )

    report_lines.extend(
        [
            "",
            "Actual-stage density-equation decomposition:",
            (
                "  The standard split uses the exact post-implicit-limiter "
                "states supplied to each explicit stage."
            ),
        ]
    )
    for quantity in ("QH", "QL"):
        if quantity not in actual_transport_response:
            continue
        components = stage_exchange_components[quantity]
        actual_standard = actual_transport_response[quantity]["standard"]
        actual_ap = actual_transport_response[quantity]["source-aware"]
        actual_difference = actual_transport_response[quantity][
            "standard minus source-aware"
        ]
        fitted_standard = best_fit_factor(
            components["standard total"][audit_nodes],
            actual_standard[audit_nodes],
        )
        fitted_ap = best_fit_factor(
            components["source-aware total"][audit_nodes],
            actual_ap[audit_nodes],
        )
        fitted_difference = best_fit_factor(
            (
                components["standard total"]
                - components["source-aware total"]
            )[audit_nodes],
            actual_difference[audit_nodes],
        )
        actual_standard_norm = weighted_norm(
            actual_standard[audit_nodes], spatial_weights[audit_nodes]
        )
        actual_ap_norm = weighted_norm(
            actual_ap[audit_nodes], spatial_weights[audit_nodes]
        )
        p1_norm = weighted_norm(
            components["centered p=1 first moment"][audit_nodes],
            spatial_weights[audit_nodes],
        )
        stage_standard_norm = weighted_norm(
            components["standard total"][audit_nodes],
            spatial_weights[audit_nodes],
        )
        report_lines.extend(
            [
                f"  {quantity}:",
                (
                    "    standard stage split vs standard shadow: "
                    f"shape error={relative_norm(fitted_standard * components['standard total'][audit_nodes], actual_standard[audit_nodes]):.6e}, "
                    f"corr={correlation(components['standard total'][audit_nodes], actual_standard[audit_nodes]):.6e}"
                ),
                (
                    "    source-aware stage split vs source-aware shadow: "
                    f"shape error={relative_norm(fitted_ap * components['source-aware total'][audit_nodes], actual_ap[audit_nodes]):.6e}, "
                    f"corr={correlation(components['source-aware total'][audit_nodes], actual_ap[audit_nodes]):.6e}"
                ),
                (
                    "    LLF change vs shadow standard-source-aware change: "
                    f"shape error={relative_norm(fitted_difference * (components['standard total'] - components['source-aware total'])[audit_nodes], actual_difference[audit_nodes]):.6e}, "
                    f"corr={correlation((components['standard total'] - components['source-aware total'])[audit_nodes], actual_difference[audit_nodes]):.6e}"
                ),
            ]
        )
        for name in (
            "standard density LLF",
            "source-aware density LLF",
            "centered p=0 first moment",
            "centered p=1 first moment",
        ):
            reference_norm = (
                actual_ap_norm if name == "source-aware density LLF" else actual_standard_norm
            )
            fit = fitted_ap if name == "source-aware density LLF" else fitted_standard
            actual = actual_ap if name == "source-aware density LLF" else actual_standard
            report_lines.append(
                f"    {name}: fitted norm/actual={weighted_norm(fit * components[name][audit_nodes], spatial_weights[audit_nodes]) / max(reference_norm, np.finfo(float).tiny):.6e}, "
                f"norm/standard-stage-total={weighted_norm(components[name][audit_nodes], spatial_weights[audit_nodes]) / max(stage_standard_norm, np.finfo(float).tiny):.6e}, "
                f"corr(actual)={correlation(components[name][audit_nodes], actual[audit_nodes]):.6e}"
            )
        for name in (
            "centered backward p=1 subset",
            "centered forward p=1 subset",
        ):
            report_lines.append(
                f"    {name}: norm/p1={weighted_norm(components[name][audit_nodes], spatial_weights[audit_nodes]) / max(p1_norm, np.finfo(float).tiny):.6e}, "
                f"corr(standard shadow)={correlation(components[name][audit_nodes], actual_standard[audit_nodes]):.6e}"
            )
        report_lines.append("")

    baseline_first = grey_first["baseline"]
    baseline_backward, baseline_slope_correlation = slope_fraction(
        baseline_first,
        radius_km,
        selected_cells,
        n_cells,
        n_nodes,
    )
    no_generation_first = grey_first["without p=1 first-moment generation"]
    no_generation_backward, no_generation_slope_correlation = slope_fraction(
        no_generation_first,
        radius_km,
        selected_cells,
        n_cells,
        n_nodes,
    )
    report_lines.extend(
        [
            "Energy-weighted first moment:",
            f"  baseline backward-slope fraction: {baseline_backward:.6e}",
            f"  baseline local/envelope slope correlation: {baseline_slope_correlation:.6e}",
            (
                "  without p=1 generation backward-slope fraction: "
                f"{no_generation_backward:.6e}"
            ),
            (
                "  without p=1 generation local/envelope correlation: "
                f"{no_generation_slope_correlation:.6e}"
            ),
            "",
        ]
    )

    # Direct block actions at the beginning of the captured step.  These do
    # not alter the evolution; they identify which modal blocks generate and
    # return the nonconstant first moment.
    p0_density = p0 @ initial_density
    p1_density = p1 @ initial_density
    p1_first = p1 @ initial_first
    b10 = p1 @ (a_id @ p0_density)
    b11 = p1 @ (a_id @ p1_density)
    a01 = p0 @ (a_di @ p1_first)
    a11 = p1 @ (a_di @ p1_first)
    b10_grey = b10 @ qh_weights
    b11_grey = b11 @ qh_weights
    a01_grey = a01 @ qh_weights
    a11_grey = a11 @ qh_weights
    report_lines.extend(
        [
            "Direct modal actions at the saved state:",
            (
                "  ||P1 B P0 D|| / ||P1 B D||: "
                f"{weighted_norm(b10_grey[audit_nodes], spatial_weights[audit_nodes]) / max(weighted_norm((b10_grey + b11_grey)[audit_nodes], spatial_weights[audit_nodes]), np.finfo(float).tiny):.6e}"
            ),
            (
                "  ||P1 B P1 D|| / ||P1 B D||: "
                f"{weighted_norm(b11_grey[audit_nodes], spatial_weights[audit_nodes]) / max(weighted_norm((b10_grey + b11_grey)[audit_nodes], spatial_weights[audit_nodes]), np.finfo(float).tiny):.6e}"
            ),
            (
                "  corr(P1 B P0 D, P1 B P1 D): "
                f"{correlation(b10_grey[audit_nodes], b11_grey[audit_nodes]):.6e}"
            ),
            (
                "  ||P0 A P1 I|| / ||A P1 I||: "
                f"{weighted_norm(a01_grey[audit_nodes], spatial_weights[audit_nodes]) / max(weighted_norm((a01_grey + a11_grey)[audit_nodes], spatial_weights[audit_nodes]), np.finfo(float).tiny):.6e}"
            ),
            (
                "  ||P1 A P1 I|| / ||A P1 I||: "
                f"{weighted_norm(a11_grey[audit_nodes], spatial_weights[audit_nodes]) / max(weighted_norm((a01_grey + a11_grey)[audit_nodes], spatial_weights[audit_nodes]), np.finfo(float).tiny):.6e}"
            ),
            (
                "  corr(P0 A P1 I, P1 A P1 I): "
                f"{correlation(a01_grey[audit_nodes], a11_grey[audit_nodes]):.6e}"
            ),
        ]
    )

    summary_path = args.output_prefix.with_suffix(".txt")
    summary_path.write_text("\n".join(report_lines) + "\n")
    print("\n".join(report_lines))

    plot_names = (
        "baseline",
        "without density LLF",
        "without p=1 first-moment pathway",
    )
    colors = {
        "baseline": "black",
        "without density LLF": "#0072B2",
        "without p=1 first-moment pathway": "#D55E00",
    }
    fig, axes = plt.subplots(3, 1, figsize=(11, 11), sharex=True)
    for axis, quantity in zip(axes[:2], ("QH", "QL"), strict=True):
        scale = max(
            np.max(np.abs(exchanges["baseline"][quantity][selected_nodes])),
            np.finfo(float).tiny,
        )
        for name in plot_names:
            axis.plot(
                radius_km[selected_nodes],
                exchanges[name][quantity][selected_nodes] / scale,
                ".-",
                color=colors[name],
                linewidth=0.8,
                markersize=2.0,
                label=name,
            )
        axis.set_ylabel(f"{quantity} / max |baseline|")
        axis.grid(alpha=0.2)
        axis.legend(loc="best")

    first_scale = max(
        np.max(np.abs(baseline_first[selected_nodes])),
        np.finfo(float).tiny,
    )
    axes[2].plot(
        radius_km[selected_nodes],
        baseline_first[selected_nodes] / first_scale,
        ".-",
        color=colors["baseline"],
        linewidth=0.8,
        markersize=2.0,
        label="baseline first moment",
    )
    axes[2].plot(
        radius_km[selected_nodes],
        no_generation_first[selected_nodes] / first_scale,
        ".-",
        color="#CC79A7",
        linewidth=0.8,
        markersize=2.0,
        label="without p=1 first-moment generation",
    )
    axes[2].set_ylabel("energy-weighted first moment / max |baseline|")
    axes[2].set_xlabel("radius [km]")
    axes[2].grid(alpha=0.2)
    axes[2].legend(loc="best")
    axes[2].set_xscale("log")
    axes[2].set_xlim(args.radius_min_km, args.radius_max_km)
    fig.suptitle("Frozen-step modal causality audit")
    fig.tight_layout()
    fig.savefig(args.output_prefix.with_suffix(".png"), dpi=180)

    if actual_transport_response:
        fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
        component_colors = {
            "standard density LLF": "#D55E00",
            "source-aware density LLF": "#E69F00",
            "centered p=0 first moment": "#0072B2",
            "centered p=1 first moment": "#009E73",
            "centered backward p=1 subset": "#CC79A7",
        }
        for row, quantity in enumerate(("QH", "QL")):
            if quantity not in actual_transport_response:
                continue
            components = stage_exchange_components[quantity]
            actual_standard = actual_transport_response[quantity]["standard"]
            fit_standard = best_fit_factor(
                components["standard total"][audit_nodes],
                actual_standard[audit_nodes],
            )
            scale = max(
                float(np.max(np.abs(actual_standard[selected_nodes]))),
                np.finfo(float).tiny,
            )
            axis = axes[row, 0]
            axis.plot(
                radius_km[selected_nodes],
                actual_standard[selected_nodes] / scale,
                color="black",
                linewidth=1.2,
                label="standard shadow response",
            )
            axis.plot(
                radius_km[selected_nodes],
                fit_standard * components["standard total"][selected_nodes] / scale,
                color="#56B4E9",
                linewidth=0.9,
                label="recorded-stage standard split",
            )
            axis.set_ylabel(f"{quantity} / max |standard shadow|")
            axis.grid(alpha=0.2)

            axis = axes[row, 1]
            floor = 1.0e-14
            for name in component_colors:
                normalized = np.maximum(
                    np.abs(fit_standard * components[name][selected_nodes]) / scale,
                    floor,
                )
                axis.plot(
                    radius_km[selected_nodes],
                    normalized,
                    color=component_colors[name],
                    linewidth=1.0,
                    label=(
                        "source-aware density LLF (same stages)"
                        if name == "source-aware density LLF"
                        else name
                    ),
                )
            axis.set_yscale("log")
            axis.set_ylabel(f"|fitted {quantity} component| / max |shadow|")
            axis.grid(alpha=0.2)
        for axis in axes[-1, :]:
            axis.set_xlabel("radius [km]")
        for axis in axes.flat:
            axis.set_xscale("log")
            axis.set_xlim(args.radius_min_km, args.radius_max_km)
        axes[0, 0].legend(loc="upper left", fontsize=8)
        axes[0, 1].legend(loc="upper left", fontsize=8)
        fig.suptitle("Recorded-stage density-update attribution")
        fig.tight_layout()
        for suffix in ("png", "pdf"):
            fig.savefig(
                Path(str(args.output_prefix) + f".stage_attribution.{suffix}"),
                dpi=180,
            )

    fig, axes = plt.subplots(2, 2, figsize=(13, 7), sharex=True)
    for row, quantity in enumerate(("QH", "QL")):
        weights = qh_weights if quantity == "QH" else ql_weights
        spatial = stage_first_components["total"] @ weights
        collision = stage_first_collision @ weights
        completed = stage_first_balance @ weights
        scale = max(
            float(np.max(np.abs(spatial[selected_nodes]))),
            np.finfo(float).tiny,
        )
        axis = axes[row, 0]
        axis.plot(
            radius_km[selected_nodes],
            spatial[selected_nodes] / scale,
            color="#0072B2",
            linewidth=1.0,
            label="explicit spatial increment",
        )
        axis.plot(
            radius_km[selected_nodes],
            collision[selected_nodes] / scale,
            color="#D55E00",
            linewidth=1.0,
            linestyle="--",
            label="implicit collision increment",
        )
        axis.plot(
            radius_km[selected_nodes],
            completed[selected_nodes] / scale,
            color="black",
            linewidth=1.0,
            label="sum",
        )
        axis.set_ylabel(f"{quantity}-weighted increment\n/ max |spatial|")
        axis.grid(alpha=0.2)
        axis = axes[row, 1]
        floor = 1.0e-16
        for values, color, linestyle, label in (
            (spatial, "#0072B2", "-", "explicit spatial increment"),
            (collision, "#D55E00", "--", "implicit collision increment"),
            (completed, "black", "-", "sum"),
        ):
            axis.plot(
                radius_km[selected_nodes],
                np.maximum(np.abs(values[selected_nodes]) / scale, floor),
                color=color,
                linewidth=1.0,
                linestyle=linestyle,
                label=label,
            )
        axis.set_yscale("log")
        axis.set_ylabel(f"|{quantity}-weighted increment|\n/ max |spatial|")
        axis.grid(alpha=0.2)
    axes[0, 0].legend(loc="upper right", fontsize=8)
    axes[0, 1].legend(loc="upper right", fontsize=8)
    for axis in axes.flat:
        axis.set_xscale("log")
        axis.set_xlim(args.radius_min_km, args.radius_max_km)
    for axis in axes[-1, :]:
        axis.set_xlabel("radius [km]")
    fig.suptitle("Recorded-stage first-moment balance")
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(
            Path(str(args.output_prefix) + f".first_moment_balance.{suffix}"),
            dpi=180,
        )


if __name__ == "__main__":
    main()
