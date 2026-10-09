#!/usr/bin/env python3
"""Audit one frozen-matter Thornado radiation step at nodal resolution.

The input is ``matter_diagnostics.dat.frozen_transient`` written by
``StandaloneMatterDiagnosticsModule``.  The file contains the state at the
start of the captured step, the collision coefficients evaluated at that
state, separated applications of the production streaming operator, and the
state accepted by Thornado at the end of the step.

The default audit uses those production operator applications directly.  The
optional matrix audit assembles the static, radial, isotropic two-moment DG
operator on the saved mesh for one energy/species block and advances it with
the saved IMEX tableau.  Its streaming-action mismatch is reported before the
one-step comparison; a large mismatch means that the isotropic matrix is not a
faithful representation of that particular saved state.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np


def read_table(path: Path) -> tuple[dict[str, float], list[str], np.ndarray]:
    metadata: dict[str, float] = {}
    columns: list[str] | None = None
    with path.open() as stream:
        for line in stream:
            if not line.startswith("#"):
                continue
            words = line[1:].split()
            if not words:
                continue
            if words[0] == "element":
                columns = words
                break
            if len(words) == 2:
                try:
                    metadata[words[0]] = float(words[1])
                except ValueError:
                    pass
    if columns is None:
        raise ValueError(f"no column header in {path}")
    values = np.loadtxt(path, comments="#", ndmin=2)
    if values.shape[1] != len(columns):
        raise ValueError(
            f"{path}: {values.shape[1]} values but {len(columns)} columns"
        )
    return metadata, columns, values


def column_map(columns: list[str], values: np.ndarray) -> dict[str, np.ndarray]:
    return {name: values[:, index] for index, name in enumerate(columns)}


def relative_norm(a: np.ndarray, b: np.ndarray) -> float:
    scale = max(np.linalg.norm(b), np.finfo(float).tiny)
    return float(np.linalg.norm(a - b) / scale)


def lagrange_values(nodes: np.ndarray, point: float) -> np.ndarray:
    result = np.ones(nodes.size)
    for j in range(nodes.size):
        for k in range(nodes.size):
            if j != k:
                result[j] *= (point - nodes[k]) / (nodes[j] - nodes[k])
    return result


def derivative_matrix(nodes: np.ndarray) -> np.ndarray:
    matrix = np.empty((nodes.size, nodes.size))
    for i, point in enumerate(nodes):
        for j in range(nodes.size):
            value = 0.0
            for m in range(nodes.size):
                if m == j:
                    continue
                term = 1.0 / (nodes[j] - nodes[m])
                for k in range(nodes.size):
                    if k != j and k != m:
                        term *= (point - nodes[k]) / (nodes[j] - nodes[k])
                value += term
            matrix[i, j] = value
    return matrix


def assemble_isotropic_spherical_operator(
    block: dict[str, np.ndarray], source_aware: bool
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Assemble the static radial p_N DG streaming matrices on a saved mesh."""

    element = block["element"].astype(int)
    node = block["spatial_node"].astype(int)
    order = np.lexsort((node, element))
    n_cells = int(element.max())
    n_nodes = int(node.max())
    size = n_cells * n_nodes
    if order.size != size:
        raise ValueError("selected energy/species block is not a complete spatial grid")

    def shaped(name: str) -> np.ndarray:
        return block[name][order].reshape(n_cells, n_nodes)

    radius = shaped("r_code")
    widths = shaped("dr_code")[:, 0]
    weights = shaped("wx")[0]
    centers = np.mean(radius, axis=1)
    nodes = (radius[0] - centers[0]) / widths[0]
    lower = lagrange_values(nodes, -0.5)
    upper = lagrange_values(nodes, +0.5)
    derivative = derivative_matrix(nodes)
    sqrtgm = shaped("sqrtgm")
    h1 = shaped("h1")

    # Thornado interpolates h_i to a face and then reconstructs sqrt(gamma).
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

    opacity = shaped("kappa_transport")
    opacity_lower = opacity @ lower
    opacity_upper = opacity @ upper
    theta = np.ones(n_cells + 1)
    if source_aware:
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
            tau = opacity_face * normal * face_width
            theta[face] = 1.0 / (1.0 + 1.5 * normal * tau)
        boundary_opacity = 0.5 * (
            max(0.0, opacity_upper[-1]) + max(0.0, opacity_lower[-1])
        )
        theta[-1] = 1.0 / (
            1.0
            + 1.5 * abs(h1_face[-1]) ** 2 * boundary_opacity * widths[-1]
        )

    def traces(vector: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        cells = vector.reshape(n_cells, n_nodes)
        return cells @ lower, cells @ upper

    def weak_rhs(volume: np.ndarray, face_flux: np.ndarray) -> np.ndarray:
        cells = volume.reshape(n_cells, n_nodes)
        result = np.empty_like(cells)
        for cell in range(n_cells):
            raw = (
                lower * sqrtgm_face[cell] * face_flux[cell]
                - upper * sqrtgm_face[cell + 1] * face_flux[cell + 1]
                + derivative.T @ (weights * sqrtgm[cell] * cells[cell])
            )
            result[cell] = raw / (weights * widths[cell] * sqrtgm[cell])
        return result.reshape(-1)

    def density_rhs(first: np.ndarray, density: np.ndarray) -> np.ndarray:
        i_lower, i_upper = traces(first)
        d_lower, d_upper = traces(density)
        face_flux = np.zeros(n_cells + 1)
        # bcX=31 copies the nodal polynomial into each ghost element.  At the
        # inner boundary G1 is reflected; at the outer boundary this linear
        # representation applies while the copied nodal G1 values are >= 0.
        face_flux[0] = -0.5 * theta[0] * (d_lower[0] - d_upper[0])
        face_flux[-1] = 0.5 * (i_upper[-1] + i_lower[-1]) \
            - 0.5 * theta[-1] * (d_lower[-1] - d_upper[-1])
        for face in range(1, n_cells):
            face_flux[face] = 0.5 * (i_upper[face - 1] + i_lower[face]) \
                - 0.5 * theta[face] * (d_lower[face] - d_upper[face - 1])
        return weak_rhs(first, face_flux)

    def first_rhs(density: np.ndarray, first: np.ndarray) -> np.ndarray:
        d_lower, d_upper = traces(density)
        i_lower, i_upper = traces(first)
        face_flux = np.empty(n_cells + 1)
        face_flux[0] = (d_upper[0] + d_lower[0]) / 6.0 - i_lower[0]
        face_flux[-1] = (d_upper[-1] + d_lower[-1]) / 6.0 \
            - 0.5 * (i_lower[-1] - i_upper[-1])
        for face in range(1, n_cells):
            face_flux[face] = (d_upper[face - 1] + d_lower[face]) / 6.0 \
                - 0.5 * (i_lower[face] - i_upper[face - 1])
        result = weak_rhs(density / 3.0, face_flux)

        # Static isotropic spherical-coordinate source.  The metric
        # derivative is evaluated with the same weak derivative used above.
        metric_derivative = np.empty_like(h2)
        for cell in range(n_cells):
            raw = (
                -lower * h2_face[cell]
                + upper * h2_face[cell + 1]
                - derivative.T @ (weights * h2[cell])
            )
            metric_derivative[cell] = raw / (weights * widths[cell])
        result += (
            (2.0 / 3.0)
            * density.reshape(n_cells, n_nodes)
            * metric_derivative
            / h2
        ).reshape(-1)
        return result

    basis = np.eye(size)
    zero = np.zeros(size)
    a_di = np.column_stack([density_rhs(basis[:, j], zero) for j in range(size)])
    a_dd = np.column_stack([density_rhs(zero, basis[:, j]) for j in range(size)])
    a_id = np.column_stack([first_rhs(basis[:, j], zero) for j in range(size)])
    a_ii = np.column_stack([first_rhs(zero, basis[:, j]) for j in range(size)])
    return a_dd, a_di, a_id, a_ii


def parse_tableau(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    a_ex_rows: dict[int, list[float]] = {}
    a_im_rows: dict[int, list[float]] = {}
    weights: dict[int, tuple[float, float]] = {}
    pattern = re.compile(r"^# stage\s+(\d+)\s+(\S+)\s+(.*)$")
    with path.open() as stream:
        for line in stream:
            match = pattern.match(line.strip())
            if not match:
                continue
            stage = int(match.group(1))
            kind = match.group(2)
            numbers = [float(value) for value in match.group(3).split()]
            if kind == "aEX":
                a_ex_rows[stage] = numbers
            elif kind == "aIM":
                a_im_rows[stage] = numbers
            elif kind == "aii_wEX_wIM":
                weights[stage] = (numbers[1], numbers[2])
    n_stages = max(a_ex_rows)
    a_ex = np.array([a_ex_rows[i] for i in range(1, n_stages + 1)])
    a_im = np.array([a_im_rows[i] for i in range(1, n_stages + 1)])
    w_ex = np.array([weights[i][0] for i in range(1, n_stages + 1)])
    w_im = np.array([weights[i][1] for i in range(1, n_stages + 1)])
    return a_ex, a_im, w_ex, w_im


def imex_step(
    initial: np.ndarray,
    streaming: np.ndarray,
    collision: np.ndarray,
    forcing: np.ndarray,
    dt: float,
    tableau: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
) -> np.ndarray:
    a_ex, a_im, w_ex, w_im = tableau
    identity = np.eye(initial.size)
    stages: list[np.ndarray] = []
    for i in range(a_ex.shape[0]):
        rhs = initial.copy()
        for j in range(i):
            rhs += dt * (
                a_ex[i, j] * (streaming @ stages[j])
                + a_im[i, j] * (collision @ stages[j] + forcing)
            )
        rhs += dt * a_im[i, i] * forcing
        stages.append(np.linalg.solve(identity - dt * a_im[i, i] * collision, rhs))
    result = initial.copy()
    for i, stage in enumerate(stages):
        result += dt * (
            w_ex[i] * (streaming @ stage)
            + w_im[i] * (collision @ stage + forcing)
        )
    return result


def select_block(data: dict[str, np.ndarray], species: int, energy: int) -> dict[str, np.ndarray]:
    mask = (data["species"].astype(int) == species) & (
        data["flat_energy"].astype(int) == energy
    )
    if not np.any(mask):
        raise ValueError(f"no rows for species={species}, flat_energy={energy}")
    return {name: value[mask] for name, value in data.items()}


def apply_centered_density_operator(
    block: dict[str, np.ndarray], first: np.ndarray
) -> np.ndarray:
    """Apply the centered spatial density operator to a radial first moment."""

    element = block["element"].astype(int)
    node = block["spatial_node"].astype(int)
    order = np.lexsort((node, element))
    n_cells = int(element.max())
    n_nodes = int(node.max())
    if order.size != n_cells * n_nodes:
        raise ValueError("selected block is not a complete spatial grid")

    def shaped(name: str) -> np.ndarray:
        return block[name][order].reshape(n_cells, n_nodes)

    radius = shaped("r_code")
    widths = shaped("dr_code")[:, 0]
    weights = shaped("wx")[0]
    centers = np.mean(radius, axis=1)
    nodes = (radius[0] - centers[0]) / widths[0]
    lower = lagrange_values(nodes, -0.5)
    upper = lagrange_values(nodes, +0.5)
    derivative = derivative_matrix(nodes)
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

    first_cells = first[order].reshape(n_cells, n_nodes)
    first_lower = first_cells @ lower
    first_upper = first_cells @ upper
    face_flux = np.zeros(n_cells + 1)
    face_flux[-1] = 0.5 * (first_upper[-1] + first_lower[-1])
    for face in range(1, n_cells):
        face_flux[face] = 0.5 * (first_upper[face - 1] + first_lower[face])

    result = np.empty_like(first_cells)
    for cell in range(n_cells):
        raw = (
            lower * sqrtgm_face[cell] * face_flux[cell]
            - upper * sqrtgm_face[cell + 1] * face_flux[cell + 1]
            + derivative.T @ (weights * sqrtgm[cell] * first_cells[cell])
        )
        result[cell] = raw / (weights * widths[cell] * sqrtgm[cell])

    inverse_order = np.empty_like(order)
    inverse_order[order] = np.arange(order.size)
    return result.reshape(-1)[inverse_order]


def apply_centered_first_moment_operator(
    block: dict[str, np.ndarray], density: np.ndarray
) -> np.ndarray:
    """Apply the centered isotropic pressure-gradient operator."""

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
    derivative = derivative_matrix(nodes)
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

    density_cells = density[order].reshape(n_cells, n_nodes)
    density_lower = density_cells @ lower
    density_upper = density_cells @ upper
    face_flux = np.empty(n_cells + 1)
    face_flux[0] = (density_upper[0] + density_lower[0]) / 6.0
    face_flux[-1] = (density_upper[-1] + density_lower[-1]) / 6.0
    for face in range(1, n_cells):
        face_flux[face] = (
            density_upper[face - 1] + density_lower[face]
        ) / 6.0

    result = np.empty_like(density_cells)
    for cell in range(n_cells):
        raw = (
            lower * sqrtgm_face[cell] * face_flux[cell]
            - upper * sqrtgm_face[cell + 1] * face_flux[cell + 1]
            + derivative.T
            @ (weights * sqrtgm[cell] * density_cells[cell] / 3.0)
        )
        result[cell] = raw / (weights * widths[cell] * sqrtgm[cell])

    metric_derivative = np.empty_like(h2)
    for cell in range(n_cells):
        raw = (
            -lower * h2_face[cell]
            + upper * h2_face[cell + 1]
            - derivative.T @ (weights * h2[cell])
        )
        metric_derivative[cell] = raw / (weights * widths[cell])
    result += (2.0 / 3.0) * density_cells * metric_derivative / h2

    inverse_order = np.empty_like(order)
    inverse_order[order] = np.arange(order.size)
    return result.reshape(-1)[inverse_order]


def apply_first_moment_llf_operator(
    block: dict[str, np.ndarray], first: np.ndarray
) -> np.ndarray:
    """Apply only the unscaled LLF first-moment jump term."""

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

    first_cells = first[order].reshape(n_cells, n_nodes)
    first_lower = first_cells @ lower
    first_upper = first_cells @ upper
    face_flux = np.empty(n_cells + 1)
    face_flux[0] = -first_lower[0]
    face_flux[-1] = -0.5 * (first_lower[-1] - first_upper[-1])
    for face in range(1, n_cells):
        face_flux[face] = -0.5 * (
            first_lower[face] - first_upper[face - 1]
        )

    result = np.empty_like(first_cells)
    for cell in range(n_cells):
        raw = (
            lower * sqrtgm_face[cell] * face_flux[cell]
            - upper * sqrtgm_face[cell + 1] * face_flux[cell + 1]
        )
        result[cell] = raw / (weights * widths[cell] * sqrtgm[cell])

    inverse_order = np.empty_like(order)
    inverse_order[order] = np.arange(order.size)
    return result.reshape(-1)[inverse_order]


def remove_element_means(
    values: np.ndarray, element: np.ndarray, node: np.ndarray, weights: np.ndarray
) -> np.ndarray:
    result = values.copy()
    for i_element in np.unique(element):
        mask = element == i_element
        order = np.argsort(node[mask])
        local = values[mask][order]
        local_weights = weights[mask][order]
        mean = np.sum(local_weights * local) / np.sum(local_weights)
        indices = np.flatnonzero(mask)[order]
        result[indices] -= mean
    return result


def correlation(a: np.ndarray, b: np.ndarray) -> float:
    if np.linalg.norm(a) == 0.0 or np.linalg.norm(b) == 0.0:
        return np.nan
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def first_moment_audit(
    data: dict[str, np.ndarray],
    metadata: dict[str, float],
    min_tau: float,
    state: str,
) -> dict[str, np.ndarray]:
    """Form the Fick residual and its density/heating response for all blocks."""

    flat_spatial = data["flat_spatial"].astype(int)
    n_spatial = int(flat_spatial.max())
    representative = np.zeros(n_spatial, dtype=int)
    for i_spatial in range(1, n_spatial + 1):
        representative[i_spatial - 1] = np.flatnonzero(flat_spatial == i_spatial)[0]

    result = {
        "element": data["element"][representative].astype(int),
        "spatial_node": data["spatial_node"][representative].astype(int),
        "r_km": data["r_code"][representative]
        / metadata["centimeter_code"]
        / 1.0e5,
        "wx": data["wx"][representative],
    }
    for name in (
        "matter_heating",
        "A_delta_I",
        "centered_density",
        "density_LLF",
        "delta_I",
        "I_Fick",
        "R_I",
        "J_I_I",
        "thick_weight",
        "total_weight",
    ):
        result[name] = np.zeros(n_spatial)

    species_values = np.unique(data["species"].astype(int))
    energy_values = np.unique(data["flat_energy"].astype(int))
    energy_reference = np.max(data["E_code"])

    for species in species_values:
        for energy in energy_values:
            block = select_block(data, species, energy)
            order = np.argsort(block["flat_spatial"].astype(int))
            for name in block:
                block[name] = block[name][order]

            kappa = block["kappa_transport"]
            optical_depth = (
                kappa * block["dr_code"] * np.sqrt(block["gm11"])
            )
            if state == "before":
                density = block["N_before"]
                first = block["G1_before"]
                centered_first = block["stream_centered_G1"]
                centered_density = block["stream_centered_N"]
                collision_density = block["coll_N"]
                collision_first = block["coll_G1"]
                j_i_i = (
                    block["stream_total_G1"] - block["stream_noLLF_G1"]
                )
            else:
                density = block["N_after"]
                first = block["G1_after"]
                centered_first = apply_centered_first_moment_operator(
                    block, density
                )
                centered_density = apply_centered_density_operator(block, first)
                collision_density = block["chi_emab"] * (
                    block["D0"] - density
                )
                collision_first = -kappa * first
                j_i_i = apply_first_moment_llf_operator(block, first)

            active = kappa > max(np.max(kappa) * 1.0e-14, np.finfo(float).tiny)
            i_fick = first.copy()
            i_fick[active] = (
                centered_first[active] / kappa[active]
            )
            delta_i = first - i_fick
            a_delta_i = apply_centered_density_operator(block, delta_i)
            r_i = centered_first + j_i_i + collision_first

            weight = (
                block["wE"]
                * block["dE_code"]
                * (block["E_code"] / energy_reference) ** 3
            )
            indices = block["flat_spatial"].astype(int) - 1
            result["matter_heating"][indices] += weight * (-collision_density)
            result["A_delta_I"][indices] += weight * a_delta_i
            result["centered_density"][indices] += (
                weight * centered_density
            )
            if state == "before":
                result["density_LLF"][indices] += weight * (
                    block["stream_total_N"] - block["stream_noLLF_N"]
                )
            result["delta_I"][indices] += weight * delta_i
            result["I_Fick"][indices] += weight * i_fick
            result["R_I"][indices] += weight * r_i
            result["J_I_I"][indices] += weight * j_i_i
            result["thick_weight"][indices] += weight * (optical_depth >= min_tau)
            result["total_weight"][indices] += weight

    result["thick_fraction"] = result["thick_weight"] / np.maximum(
        result["total_weight"], np.finfo(float).tiny
    )
    for name in ("matter_heating", "A_delta_I", "R_I", "J_I_I"):
        result[f"{name}_highpass"] = remove_element_means(
            result[name],
            result["element"],
            result["spatial_node"],
            result["wx"],
        )
    return result


def write_first_moment_table(path: Path, audit: dict[str, np.ndarray]) -> None:
    names = list(audit)
    values = np.column_stack([audit[name] for name in names])
    np.savetxt(path, values, header=" ".join(names))


def write_exact_first_moment_plot(
    path: Path,
    block: dict[str, np.ndarray],
    metadata: dict[str, float],
    exact_first: np.ndarray,
    accepted_first: np.ndarray,
    species: int,
    energy: int,
    radius_min_km: float,
    radius_max_km: float,
) -> None:
    """Plot the nodal and element-average first moment for one saved block."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_nodes = int(block["spatial_node"].max())
    n_cells = int(block["element"].max())
    radius_km = block["r_code"] / metadata["centimeter_code"] / 1.0e5
    radius = radius_km.reshape(n_cells, n_nodes)
    weights = block["wx"].reshape(n_cells, n_nodes).copy()
    weights /= np.sum(weights, axis=1)[:, None]

    exact = exact_first.reshape(n_cells, n_nodes)
    accepted = accepted_first.reshape(n_cells, n_nodes)
    radius_average = np.sum(weights * radius, axis=1)
    exact_average = np.sum(weights * exact, axis=1)
    accepted_average = np.sum(weights * accepted, axis=1)

    selected_nodes = (
        (radius_km >= radius_min_km) & (radius_km <= radius_max_km)
    )
    selected_averages = (
        (radius_average >= radius_min_km) & (radius_average <= radius_max_km)
    )
    if not np.any(selected_nodes) or not np.any(selected_averages):
        raise ValueError(
            "first-moment plot contains no data in "
            f"{radius_min_km} <= r <= {radius_max_km} km"
        )

    fig, axis = plt.subplots(figsize=(11, 6))
    axis.plot(
        radius_km[selected_nodes],
        exact_first[selected_nodes],
        "o-",
        markersize=2.0,
        linewidth=0.8,
        label="exact frozen semidiscrete solution: nodal",
    )
    axis.plot(
        radius_km[selected_nodes],
        accepted_first[selected_nodes],
        "o-",
        markersize=2.0,
        linewidth=0.8,
        label="Thornado PDARS: nodal",
    )
    axis.plot(
        radius_average[selected_averages],
        exact_average[selected_averages],
        "o-",
        color="black",
        markersize=3.0,
        linewidth=1.2,
        label="exact element average",
    )
    axis.plot(
        radius_average[selected_averages],
        accepted_average[selected_averages],
        "o-",
        color="0.5",
        markersize=3.0,
        linewidth=1.2,
        label="PDARS element average",
    )
    plotted_values = np.concatenate(
        (
            exact_first[selected_nodes],
            accepted_first[selected_nodes],
            exact_average[selected_averages],
            accepted_average[selected_averages],
        )
    )
    positive = plotted_values[plotted_values > 0.0]
    if positive.size == 0:
        raise ValueError("first-moment plot contains no positive values")
    positive_radius = np.concatenate(
        (radius_km[selected_nodes], radius_average[selected_averages])
    )
    positive_radius = positive_radius[positive_radius > 0.0]
    if positive_radius.size == 0:
        raise ValueError("first-moment plot contains no positive radii")
    axis.set_xscale("log")
    axis.set_xlim(
        radius_min_km,
        radius_max_km,
    )
    axis.set_yscale("log")
    axis.set_ylim(0.7 * np.min(positive), 1.4 * np.max(positive))
    axis.set_xlabel("radius [km]")
    axis.set_ylabel(r"$G_1$ [code units]")
    axis.set_title(
        f"species {species}, flattened energy node {energy}, first timestep"
    )
    axis.grid(alpha=0.2)
    axis.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=180)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=Path)
    parser.add_argument("--species", type=int, default=1)
    parser.add_argument("--energy", type=int, default=1, help="flattened energy-node index")
    parser.add_argument("--matrix-audit", action="store_true")
    parser.add_argument(
        "--exact-time-audit",
        action="store_true",
        help="compare Thornado with the exact-in-time frozen semidiscrete system",
    )
    parser.add_argument(
        "--exact-first-moment-plot",
        type=Path,
        help="plot the exact and accepted first moment; requires --exact-time-audit",
    )
    parser.add_argument(
        "--exact-first-moment-radius-min-km",
        type=float,
        default=7.0,
        help="minimum radius shown in the exact first-moment plot",
    )
    parser.add_argument(
        "--exact-first-moment-radius-max-km",
        type=float,
        default=100.0,
        help="maximum radius shown in the exact first-moment plot",
    )
    parser.add_argument("--plot", type=Path)
    parser.add_argument(
        "--first-moment-audit",
        action="store_true",
        help="test the Fick residual, its density response, and first-moment LLF term",
    )
    parser.add_argument("--first-moment-table", type=Path)
    parser.add_argument("--first-moment-plot", type=Path)
    parser.add_argument("--min-optical-depth", type=float, default=10.0)
    parser.add_argument(
        "--audit-state",
        choices=("before", "after"),
        default="after",
        help="state at which to evaluate the frozen first-moment balance",
    )
    args = parser.parse_args()
    if args.exact_first_moment_plot and not args.exact_time_audit:
        parser.error("--exact-first-moment-plot requires --exact-time-audit")

    metadata, columns, values = read_table(args.file)
    data = column_map(columns, values)
    dt = metadata["timestep_code"]

    collision_formula_n = data["chi_emab"] * (data["D0"] - data["D_before"])
    collision_formula_g1 = -data["kappa_transport"] * data["gm11"] * data["I1_before"]
    collision_n_error = relative_norm(data["coll_N"], collision_formula_n)
    collision_g1_error = relative_norm(data["coll_G1"], collision_formula_g1)
    print(f"rows={values.shape[0]}  saved_state_time_s={metadata.get('state_time_s', np.nan):.8e}")
    print(f"collision N vs EmAb formula: rel={collision_n_error:.6e}")
    print(f"collision G1 vs -kappa I_1: rel={collision_g1_error:.6e}")
    print(
        "streaming split closure N: "
        f"rel={relative_norm(data['stream_centered_N'], data['stream_surface_N'] + data['stream_volume_N']):.6e}"
    )
    print(
        "streaming split closure G1: "
        f"rel={relative_norm(data['stream_centered_G1'], data['stream_surface_G1'] + data['stream_volume_G1']):.6e}"
    )
    print(
        "stored semidiscrete closure G1: "
        f"rel={relative_norm(data['semidiscrete_G1'], data['stream_total_G1'] + data['coll_G1']):.6e}"
    )
    print(
        "forward Euler vs accepted full-step state: "
        f"N rel={relative_norm(data['forward_euler_N'], data['N_after']):.6e}, "
        f"G1 rel={relative_norm(data['forward_euler_G1'], data['G1_after']):.6e}"
    )

    block = select_block(data, args.species, args.energy)
    order = np.lexsort((block["spatial_node"], block["element"]))
    for name in block:
        block[name] = block[name][order]
    first_residual = block["stream_total_G1"] + block["coll_G1"]
    first_scale = np.abs(block["stream_total_G1"]) + np.abs(block["coll_G1"])
    print(
        f"selected block s={args.species}, e={args.energy}: "
        f"max |R_I|/(|stream_I|+|collision_I|)="
        f"{np.max(np.abs(first_residual) / np.maximum(first_scale, np.finfo(float).tiny)):.6e}"
    )

    matrix_prediction = None
    exact_prediction = None
    if args.matrix_audit or args.exact_time_audit:
        max_velocity = max(
            np.max(np.abs(data["v1"])),
            np.max(np.abs(data["v2"])),
            np.max(np.abs(data["v3"])),
        )
        if max_velocity > 1.0e-12:
            raise ValueError(
                "--matrix-audit currently represents only a static fluid; "
                f"max |v|={max_velocity:.6e}"
            )
        if collision_n_error > 1.0e-8 or collision_g1_error > 1.0e-8:
            raise ValueError(
                "--matrix-audit currently represents only EmAb+Iso collision "
                "terms; the exported collision operator does not have that form "
                f"(density rel={collision_n_error:.6e}, "
                f"first-moment rel={collision_g1_error:.6e})"
            )
        outer_element = int(block["element"].max())
        outer_flux = block["G1_before"][
            block["element"].astype(int) == outer_element
        ]
        if np.any(outer_flux < 0.0):
            raise ValueError(
                "--matrix-audit assumes the bcX=31 outer-boundary MAX is on "
                "its linear, nonnegative-G1 branch"
            )
        source_aware = bool(round(metadata.get("source_aware_llf", 0.0)))
        a_dd, a_di, a_id, a_ii = assemble_isotropic_spherical_operator(
            block, source_aware
        )
        n = block["D_before"].size
        streaming = np.block([[a_dd, a_di], [a_id, a_ii]])
        initial = np.concatenate((block["N_before"], block["G1_before"]))
        exported_stream = np.concatenate((block["stream_total_N"], block["stream_total_G1"]))
        action = streaming @ initial
        interior = np.ones(2 * n, dtype=bool)
        n_nodes = int(block["spatial_node"].max())
        interior[:n_nodes] = False
        interior[n - n_nodes:n] = False
        interior[n:n + n_nodes] = False
        interior[-n_nodes:] = False
        print(
            "isotropic matrix vs production streaming action: "
            f"all rel={relative_norm(action, exported_stream):.6e}, "
            f"interior rel={relative_norm(action[interior], exported_stream[interior]):.6e}"
        )

        chi = block["chi_emab"]
        kappa = block["kappa_transport"]
        collision = np.diag(np.concatenate((-chi, -kappa)))
        forcing = np.concatenate((chi * block["D0"], np.zeros(n)))
        tableau_path = Path(str(args.file).replace(".frozen_transient", ".frozen_transient_stages"))
        tableau = parse_tableau(tableau_path)
        matrix_prediction = imex_step(initial, streaming, collision, forcing, dt, tableau)
        accepted = np.concatenate((block["N_after"], block["G1_after"]))
        print(
            "frozen isotropic IMEX prediction vs Thornado accepted state: "
            f"all rel={relative_norm(matrix_prediction, accepted):.6e}, "
            f"interior rel={relative_norm(matrix_prediction[interior], accepted[interior]):.6e}"
        )

        if args.exact_time_audit:
            try:
                from scipy.sparse import bmat, csr_matrix
                from scipy.sparse.linalg import expm_multiply
            except ImportError as error:
                raise RuntimeError(
                    "--exact-time-audit requires scipy"
                ) from error

            generator = streaming + collision
            augmented = bmat(
                [
                    [csr_matrix(generator), csr_matrix(forcing[:, None])],
                    [None, csr_matrix((1, 1))],
                ],
                format="csr",
            )
            exact_prediction = expm_multiply(
                augmented * dt, np.concatenate((initial, [1.0]))
            )[:-1]
            exact_first = exact_prediction[n:]
            accepted_first = accepted[n:]
            exact_highpass = remove_element_means(
                exact_first,
                block["element"].astype(int),
                block["spatial_node"].astype(int),
                block["wx"],
            )
            accepted_highpass = remove_element_means(
                accepted_first,
                block["element"].astype(int),
                block["spatial_node"].astype(int),
                block["wx"],
            )
            print(
                "exact frozen semidiscrete solution vs Thornado: "
                f"state rel={relative_norm(exact_prediction, accepted):.6e}, "
                f"first moment rel={relative_norm(exact_first, accepted_first):.6e}"
            )
            print(
                "selected-block within-element first-moment fraction: "
                f"exact={np.linalg.norm(exact_highpass) / max(np.linalg.norm(exact_first), np.finfo(float).tiny):.6e}, "
                f"Thornado={np.linalg.norm(accepted_highpass) / max(np.linalg.norm(accepted_first), np.finfo(float).tiny):.6e}; "
                f"component difference rel={relative_norm(exact_highpass, accepted_highpass):.6e}"
            )
            if args.exact_first_moment_plot:
                write_exact_first_moment_plot(
                    args.exact_first_moment_plot,
                    block,
                    metadata,
                    exact_first,
                    accepted_first,
                    args.species,
                    args.energy,
                    args.exact_first_moment_radius_min_km,
                    args.exact_first_moment_radius_max_km,
                )

    first_audit = None
    if args.first_moment_audit or args.first_moment_table or args.first_moment_plot:
        first_audit = first_moment_audit(
            data, metadata, args.min_optical_depth, args.audit_state
        )
        thick = first_audit["thick_fraction"] >= 0.5
        if not np.any(thick):
            raise ValueError(
                "no spatial points satisfy the requested optically thick criterion"
            )
        heating_hp = first_audit["matter_heating_highpass"][thick]
        response_hp = first_audit["A_delta_I_highpass"][thick]
        residual_hp = first_audit["R_I_highpass"][thick]
        llf_hp = first_audit["J_I_I_highpass"][thick]
        print(
            "optically thick energy-integrated element-scale audit "
            f"at the {args.audit_state}-step state "
            f"(weighted thick fraction >= 0.5, tau >= {args.min_optical_depth:g}):"
        )
        print(
            "  A delta-I vs matter heating: "
            f"rel={relative_norm(response_hp, heating_hp):.6e}, "
            f"correlation={correlation(response_hp, heating_hp):.6e}, "
            f"norm ratio={np.linalg.norm(response_hp) / max(np.linalg.norm(heating_hp), np.finfo(float).tiny):.6e}"
        )
        print(
            "  first-moment LLF vs full first-moment residual: "
            f"norm ratio={np.linalg.norm(llf_hp) / max(np.linalg.norm(residual_hp), np.finfo(float).tiny):.6e}, "
            f"correlation={correlation(llf_hp, residual_hp):.6e}"
        )
        print(
            "  delta-I vs I_Fick: "
            f"norm ratio={np.linalg.norm(first_audit['delta_I'][thick]) / max(np.linalg.norm(first_audit['I_Fick'][thick]), np.finfo(float).tiny):.6e}"
        )
        if args.first_moment_table:
            write_first_moment_table(args.first_moment_table, first_audit)

        if args.first_moment_plot:
            import matplotlib.pyplot as plt

            radius = first_audit["r_km"]
            thick_radius = radius[thick]
            radial_limit = 1.1 * np.max(thick_radius)

            def normalized(values: np.ndarray) -> np.ndarray:
                return values / max(np.max(np.abs(values[thick])), np.finfo(float).tiny)

            fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
            axes[0].plot(
                radius,
                normalized(first_audit["matter_heating_highpass"]),
                label="matter heating",
            )
            axes[0].plot(
                radius,
                normalized(first_audit["A_delta_I_highpass"]),
                label=r"$A\,\delta I$",
            )
            axes[0].set_ylabel("normalized element-scale part")
            axes[0].legend(loc="best")

            axes[1].plot(radius, first_audit["delta_I"], label=r"$I-I_{\rm Fick}$")
            axes[1].plot(radius, first_audit["I_Fick"], label=r"$I_{\rm Fick}$")
            axes[1].set_ylabel("energy-weighted first moment")
            axes[1].legend(loc="best")

            axes[2].plot(
                radius,
                normalized(first_audit["R_I_highpass"]),
                label=r"$R_I$",
            )
            axes[2].plot(
                radius,
                normalized(first_audit["J_I_I_highpass"]),
                label=r"$J_I I$",
            )
            axes[2].set_ylabel("normalized element-scale part")
            axes[2].set_xlabel("radius [km]")
            axes[2].legend(loc="best")
            axes[2].set_xlim(0.0, radial_limit)
            fig.tight_layout()
            fig.savefig(args.first_moment_plot, dpi=180)

    if args.plot:
        import matplotlib.pyplot as plt

        radius_km = block["r_code"] / metadata["centimeter_code"] / 1.0e5
        fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
        axes[0].plot(radius_km, block["stream_total_G1"], label="streaming")
        axes[0].plot(radius_km, block["coll_G1"], label="collision")
        axes[0].plot(radius_km, first_residual, label="sum")
        axes[0].set_ylabel("first-moment rate [code]")
        axes[0].legend(loc="best")

        axes[1].plot(radius_km, block["actual_rate_N"], label="actual step")
        axes[1].plot(radius_km, block["semidiscrete_N"], label="initial RHS")
        axes[1].set_ylabel("density rate [code]")
        axes[1].legend(loc="best")

        axes[2].plot(radius_km, block["N_after"], label="Thornado")
        axes[2].plot(radius_km, block["forward_euler_N"], label="forward Euler")
        if matrix_prediction is not None:
            axes[2].plot(radius_km, matrix_prediction[: radius_km.size], label="frozen IMEX")
        axes[2].set_ylabel("N")
        axes[2].set_xlabel("radius [km]")
        axes[2].legend(loc="best")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=180)


if __name__ == "__main__":
    main()
