#!/usr/bin/env python3
"""Plot the standalone transport-versus-collision diagnostic ledger."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import numpy as np


def read_named_table(path: Path) -> dict[str, np.ndarray]:
    with path.open() as stream:
        for line in stream:
            if line.startswith("# r_cm") or line.startswith("# element"):
                names = line[1:].split()
                break
        else:
            raise ValueError(f"{path}: missing commented column header")
    values = np.atleast_2d(np.loadtxt(path, comments="#"))
    if values.shape[1] != len(names):
        raise ValueError(
            f"{path}: found {values.shape[1]} columns, expected {len(names)}"
        )
    return dict(zip(names, values.T, strict=True))


def symlog_axis(axis, values: list[np.ndarray]) -> None:
    maximum = max(float(np.max(np.abs(value))) for value in values)
    if maximum > 0.0:
        axis.set_yscale("symlog", linthresh=maximum * 1.0e-8)


def center_edges(centers: np.ndarray, logarithmic: bool) -> np.ndarray:
    if centers.size == 1:
        width = 0.1 * centers[0] if centers[0] > 0.0 else 0.5
        return np.asarray([centers[0] - width, centers[0] + width])
    if logarithmic:
        if np.any(centers <= 0.0):
            raise ValueError("logarithmic centers must be positive")
        interior = np.sqrt(centers[:-1] * centers[1:])
        return np.concatenate(
            ([centers[0] ** 2 / interior[0]], interior, [centers[-1] ** 2 / interior[-1]])
        )
    interior = 0.5 * (centers[:-1] + centers[1:])
    return np.concatenate(
        ([2.0 * centers[0] - interior[0]], interior, [2.0 * centers[-1] - interior[-1]])
    )


def spectral_element_average(
    spectrum: dict[str, np.ndarray],
    ledger: dict[str, np.ndarray],
    field: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    elements = np.unique(spectrum["element"].astype(int))
    energies = np.unique(spectrum["E_MeV"])
    result = np.zeros((energies.size, elements.size))

    weight_by_node = {
        (int(element), int(node)): weight
        for element, node, weight, species in zip(
            ledger["element"],
            ledger["spatial_node"],
            ledger["volume_weight_cm3"],
            ledger["species"],
            strict=True,
        )
        if int(species) == 0
    }

    for i_element, element in enumerate(elements):
        element_rows = spectrum["element"].astype(int) == element
        nodes = np.unique(spectrum["spatial_node"][element_rows].astype(int))
        total_weight = sum(weight_by_node[(element, node)] for node in nodes)
        for i_energy, energy in enumerate(energies):
            value = 0.0
            for node in nodes:
                selected = (
                    element_rows
                    & (spectrum["spatial_node"].astype(int) == node)
                    & np.isclose(spectrum["E_MeV"], energy, rtol=1.0e-13, atol=0.0)
                )
                value += weight_by_node[(element, node)] * np.sum(
                    spectrum[field][selected]
                )
            result[i_energy, i_element] = value / total_weight

    radius = np.asarray(
        [
            np.average(
                ledger["r_cm"][
                    (ledger["species"] == 0)
                    & (ledger["element"].astype(int) == element)
                ],
                weights=ledger["volume_weight_cm3"][
                    (ledger["species"] == 0)
                    & (ledger["element"].astype(int) == element)
                ],
            )
            for element in elements
        ]
    )
    return radius * 1.0e-5, energies, result


def add_spectral_panel(axis, radius, energy, values, title) -> None:
    maximum = float(np.max(np.abs(values)))
    if maximum == 0.0:
        maximum = 1.0
    norm = SymLogNorm(linthresh=maximum * 1.0e-8, vmin=-maximum, vmax=maximum)
    image = axis.pcolormesh(
        center_edges(radius, logarithmic=True),
        center_edges(energy, logarithmic=True),
        values,
        shading="flat",
        cmap="coolwarm",
        norm=norm,
    )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("radius [km]")
    axis.set_ylabel("energy [MeV]")
    axis.set_title(title)
    plt.colorbar(image, ax=axis, label=r"$Q_H(E)$ [erg cm$^{-3}$ s$^{-1}$]")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "diagnostic",
        type=Path,
        help="main matter_diagnostics*.dat file",
    )
    parser.add_argument("--output", type=Path, help="output PNG path")
    args = parser.parse_args()

    main_data = read_named_table(args.diagnostic)
    ledger_path = Path(str(args.diagnostic) + ".transport_collision")
    spectrum_path = Path(str(ledger_path) + "_spectrum")
    if not ledger_path.exists() or not spectrum_path.exists():
        raise FileNotFoundError(
            f"expected {ledger_path} and {spectrum_path}; rerun with the new diagnostic code"
        )
    ledger = read_named_table(ledger_path)
    spectrum = read_named_table(spectrum_path)

    radius = main_data["r_cm"] * 1.0e-5
    collision = main_data["QH_total_erg_cm3_s"]
    transport = main_data["QH_transport_equivalent_erg_cm3_s"]
    difference = main_data["QH_transport_minus_collision_erg_cm3_s"]

    total = ledger["species"] == 0
    nodal_radius = ledger["r_cm"][total] * 1.0e-5
    nodal_element = ledger["element"][total].astype(int)
    nodal_node = ledger["spatial_node"][total].astype(int)
    nodal_weight = ledger["volume_weight_cm3"][total]
    nodal_collision = ledger["QH_collision_erg_cm3_s"][total]
    nodal_transport = ledger["QH_transport_erg_cm3_s"][total]
    nodal_difference = ledger["QH_difference_erg_cm3_s"][total]

    elements = np.unique(nodal_element)
    element_radius = np.zeros(elements.size)
    collision_variation = np.zeros(elements.size)
    transport_variation = np.zeros(elements.size)
    difference_variation = np.zeros(elements.size)
    for i_element, element in enumerate(elements):
        selected = nodal_element == element
        weight = nodal_weight[selected]
        element_radius[i_element] = np.average(
            nodal_radius[selected], weights=weight
        )
        for output, values in (
            (collision_variation, nodal_collision),
            (transport_variation, nodal_transport),
            (difference_variation, nodal_difference),
        ):
            average = np.average(values[selected], weights=weight)
            output[i_element] = np.sum(
                weight * np.abs(values[selected] - average)
            ) / max(np.sum(weight * np.abs(values[selected])), np.finfo(float).tiny)

    fig, axes = plt.subplots(2, 3, figsize=(20, 11), constrained_layout=True)

    axis = axes[0, 0]
    axis.plot(radius, collision, label="collision", color="C0")
    axis.plot(radius, transport, label="streaming equivalent", color="C1")
    axis.plot(radius, difference, label="streaming - collision", color="C3")
    axis.axhline(0.0, color="k", linewidth=0.7)
    axis.set_xscale("log")
    symlog_axis(axis, [collision, transport, difference])
    axis.set_xlabel("radius [km]")
    axis.set_ylabel(r"$Q_H$ [erg cm$^{-3}$ s$^{-1}$]")
    axis.set_title("spherical-volume element averages")
    axis.legend(fontsize="small")

    axis = axes[0, 1]
    for node in np.unique(nodal_node):
        selected = nodal_node == node
        axis.plot(
            nodal_radius[selected], nodal_collision[selected],
            marker="o", markersize=2.5, linewidth=0.8,
            label=f"collision, node {node}",
        )
        axis.plot(
            nodal_radius[selected], nodal_transport[selected],
            linestyle="--", linewidth=0.9,
            label=f"transport, node {node}",
        )
    axis.axhline(0.0, color="k", linewidth=0.7)
    axis.set_xscale("log")
    symlog_axis(axis, [nodal_collision, nodal_transport])
    axis.set_xlabel("radius [km]")
    axis.set_ylabel(r"$Q_H$ [erg cm$^{-3}$ s$^{-1}$]")
    axis.set_title("DG spatial nodes")
    axis.legend(fontsize="x-small", ncol=2)

    axis = axes[0, 2]
    axis.plot(element_radius, collision_variation, label="collision")
    axis.plot(element_radius, transport_variation, label="transport")
    axis.plot(element_radius, difference_variation, label="difference")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("radius [km]")
    axis.set_ylabel(r"$\sum w|q-\bar q|/\sum w|q|$")
    axis.set_title("within-element variation")
    axis.legend(fontsize="small")

    spectral_radius, energy, spectral_collision = spectral_element_average(
        spectrum, ledger, "QH_collision_erg_cm3_s"
    )
    _, _, spectral_transport = spectral_element_average(
        spectrum, ledger, "QH_transport_erg_cm3_s"
    )
    _, _, spectral_difference = spectral_element_average(
        spectrum, ledger, "QH_difference_erg_cm3_s"
    )
    add_spectral_panel(
        axes[1, 0], spectral_radius, energy, spectral_collision,
        "collision contribution by energy",
    )
    add_spectral_panel(
        axes[1, 1], spectral_radius, energy, spectral_transport,
        "streaming-equivalent contribution by energy",
    )
    add_spectral_panel(
        axes[1, 2], spectral_radius, energy, spectral_difference,
        "streaming - collision by energy",
    )

    output = args.output
    if output is None:
        output = args.diagnostic.with_name(
            args.diagnostic.stem + "_transport_collision.png"
        )
    fig.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
