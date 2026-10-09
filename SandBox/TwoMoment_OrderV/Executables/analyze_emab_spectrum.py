#!/usr/bin/env python3
"""Inspect the energy-resolved EmAb heating residual in the standalone test."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class ElementSpectrum:
    radius_km: np.ndarray
    energy_mev: np.ndarray
    dnu_pre: np.ndarray
    dnu_post: np.ndarray
    dnu_eq: np.ndarray
    qh_pre: np.ndarray
    qh_post: np.ndarray


def read_columns(path: Path) -> dict[str, np.ndarray]:
    with path.open() as stream:
        header = next((line for line in stream if line.startswith("# iX ")), None)
    if header is None:
        raise ValueError(f"{path}: missing spectral diagnostic header")
    names = header[1:].split()
    values = np.atleast_2d(np.loadtxt(path, comments="#"))
    if values.shape[1] != len(names):
        raise ValueError(
            f"{path}: found {values.shape[1]} columns, expected {len(names)}"
        )
    return dict(zip(names, values.T, strict=True))


def element_average(columns: dict[str, np.ndarray], species: int) -> ElementSpectrum:
    selected = columns["species"].astype(int) == species
    i_x = columns["iX"][selected].astype(int)
    i_e = columns["iE"][selected].astype(int)
    i_node_e = columns["iNodeE"][selected].astype(int)
    i_node_x = columns["iNodeX"][selected].astype(int)
    radius = columns["r_cm"][selected] * 1.0e-5
    energy = columns["E_MeV"][selected]
    weight = columns["volume_weight_cm3"][selected]

    x_indices = np.unique(i_x)
    energy_pairs = np.unique(np.column_stack((i_e, i_node_e)), axis=0)
    n_x = x_indices.size
    n_energy = energy_pairs.shape[0]

    radius_element = np.empty(n_x)
    energy_nodes = np.empty(n_energy)
    fields = {
        name: np.empty((n_x, n_energy))
        for name in (
            "Dnu_pre", "Dnu_post", "Dnu_eq_post",
            "QH_pre_erg_cm3_s", "QH_post_erg_cm3_s",
        )
    }

    for ix_out, ix_value in enumerate(x_indices):
        spatial = i_x == ix_value
        first_energy = (
            (i_e == energy_pairs[0, 0])
            & (i_node_e == energy_pairs[0, 1])
        )
        radial = spatial & first_energy
        radius_element[ix_out] = np.mean(radius[radial])

        for ie_out, (ie_value, inode_value) in enumerate(energy_pairs):
            points = spatial & (i_e == ie_value) & (i_node_e == inode_value)
            if ix_out == 0:
                energy_nodes[ie_out] = energy[points][0]
            denominator = np.sum(weight[points])
            for name in fields:
                fields[name][ix_out, ie_out] = np.sum(
                    columns[name][selected][points] * weight[points]
                ) / denominator

    order = np.argsort(energy_nodes)
    return ElementSpectrum(
        radius_km=radius_element,
        energy_mev=energy_nodes[order],
        dnu_pre=fields["Dnu_pre"][:, order],
        dnu_post=fields["Dnu_post"][:, order],
        dnu_eq=fields["Dnu_eq_post"][:, order],
        qh_pre=fields["QH_pre_erg_cm3_s"][:, order],
        qh_post=fields["QH_post_erg_cm3_s"][:, order],
    )


def logarithmic_edges(centers: np.ndarray) -> np.ndarray:
    if np.any(centers <= 0.0) or np.any(np.diff(centers) <= 0.0):
        raise ValueError("plot coordinates must be positive and increasing")
    edges = np.empty(centers.size + 1)
    edges[1:-1] = np.sqrt(centers[:-1] * centers[1:])
    edges[0] = centers[0] ** 2 / edges[1]
    edges[-1] = centers[-1] ** 2 / edges[-2]
    return edges


def report(
    spectra: dict[int, ElementSpectrum], gray_file: Path | None = None
) -> None:
    for species, data in spectra.items():
        total = np.sum(data.qh_post, axis=1)
        gross = np.sum(np.abs(data.qh_post), axis=1)
        dominant = np.argmax(np.abs(data.qh_post), axis=1)
        print(f"\nspecies {species}")
        print("  r_km       QH_post          cancellation  dominant_E_MeV  fraction")
        for index in range(min(16, data.radius_km.size)):
            cancellation = abs(total[index]) / max(
                gross[index], np.finfo(float).tiny
            )
            dominant_fraction = abs(data.qh_post[index, dominant[index]]) / max(
                gross[index], np.finfo(float).tiny
            )
            print(
                f"  {data.radius_km[index]:7.2f}  {total[index]: .6e}  "
                f"{cancellation: .3e}       "
                f"{data.energy_mev[dominant[index]]:10.4f}  "
                f"{dominant_fraction: .3e}"
            )

    if gray_file is not None and gray_file.exists():
        from analyze_standalone_matter_diagnostics import read_diagnostic

        gray = read_diagnostic(gray_file)
        print("\nenergy-resolved sum vs gray EmAb ledger:")
        for species, data in spectra.items():
            spectral_sum = np.sum(data.qh_post, axis=1)
            gray_qh = gray.columns[f"EmAb_QH_s{species}_erg_cm3_s"]
            scale = max(
                float(np.max(np.abs(spectral_sum))),
                float(np.max(np.abs(gray_qh))),
                np.finfo(float).tiny,
            )
            mismatch = float(np.max(np.abs(spectral_sum - gray_qh))) / scale
            print(f"  species {species}: max relative mismatch={mismatch:.3e}")


def plot(spectra: dict[int, ElementSpectrum], output: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm

    figure, axes = plt.subplots(
        len(spectra), 3, figsize=(19, 5.2 * len(spectra)), constrained_layout=True
    )
    axes = np.atleast_2d(axes)

    for row, (species, data) in enumerate(spectra.items()):
        radius_edges = logarithmic_edges(data.radius_km)
        energy_edges = logarithmic_edges(data.energy_mev)
        denominator = np.abs(data.dnu_post) + np.abs(data.dnu_eq)
        relative_departure = np.divide(
            data.dnu_post - data.dnu_eq,
            denominator,
            out=np.zeros_like(denominator),
            where=denominator > 0.0,
        )

        axis = axes[row, 0]
        image = axis.pcolormesh(
            radius_edges,
            energy_edges,
            relative_departure.T,
            shading="flat",
            cmap="coolwarm",
            norm=SymLogNorm(linthresh=1.0e-8, vmin=-1.0, vmax=1.0),
        )
        figure.colorbar(image, ax=axis, label=r"$(D-D_{\rm eq})/(|D|+|D_{\rm eq}|)$")
        axis.set_title(f"species {species}: spectral disequilibrium")

        axis = axes[row, 1]
        maximum = float(np.max(np.abs(data.qh_post)))
        if maximum > 0.0:
            norm = SymLogNorm(
                linthresh=maximum * 1.0e-10,
                vmin=-maximum,
                vmax=maximum,
            )
            image = axis.pcolormesh(
                radius_edges,
                energy_edges,
                data.qh_post.T,
                shading="flat",
                cmap="coolwarm",
                norm=norm,
            )
            figure.colorbar(
                image, ax=axis,
                label=r"energy-node contribution to $Q_H$ [erg cm$^{-3}$ s$^{-1}$]",
            )
        axis.set_title(f"species {species}: post-limiter EmAb heating")

        axis = axes[row, 2]
        for index in range(min(8, data.radius_km.size)):
            scale = max(
                float(np.sum(np.abs(data.qh_post[index]))),
                np.finfo(float).tiny,
            )
            line = axis.plot(
                data.energy_mev,
                data.qh_post[index] / scale,
                label=f"{data.radius_km[index]:.1f} km",
            )[0]
            axis.plot(
                data.energy_mev,
                data.qh_pre[index] / scale,
                color=line.get_color(),
                linestyle="--",
                alpha=0.65,
            )
        axis.axhline(0.0, color="k", linewidth=0.8)
        axis.set_xscale("log")
        axis.set_yscale("symlog", linthresh=1.0e-8)
        axis.set_title("adjacent elements (solid: post; dashed: pre-limiter)")
        axis.set_ylabel(r"energy-node $Q_H/\sum_q|Q_{H,q}|$")
        axis.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize="small")

        for axis in axes[row]:
            axis.set_xlabel("radius [km]" if axis is not axes[row, 2] else "energy [MeV]")
        axes[row, 0].set_xscale("log")
        axes[row, 0].set_yscale("log")
        axes[row, 1].set_xscale("log")
        axes[row, 1].set_yscale("log")
        axes[row, 0].set_ylabel("energy [MeV]")
        axes[row, 1].set_ylabel("energy [MeV]")

    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=Path)
    parser.add_argument("--plot", type=Path)
    parser.add_argument("--gray", type=Path)
    args = parser.parse_args()

    columns = read_columns(args.file)
    species_values = sorted(set(columns["species"].astype(int)))
    spectra = {
        species: element_average(columns, species) for species in species_values
    }
    gray_file = args.gray
    if gray_file is None and str(args.file).endswith(".emab_spectrum"):
        gray_file = Path(str(args.file)[: -len(".emab_spectrum")])
    report(spectra, gray_file)
    if args.plot is not None:
        plot(spectra, args.plot)


if __name__ == "__main__":
    main()
