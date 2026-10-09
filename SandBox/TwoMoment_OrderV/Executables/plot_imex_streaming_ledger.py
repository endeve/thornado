#!/usr/bin/env python3

"""Plot the stage-resolved decomposition of the explicit radiation RHS."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


QUANTITIES = ("QE", "QH", "QM1")
COMPONENTS = (
    ("total", "total", "black", "-"),
    ("spatial_centered", "centered spatial DG", "tab:blue", "-"),
    ("LLF_penalty", "LLF penalty", "tab:red", "--"),
    ("observer", "observer", "tab:green", ":"),
)
DG_COMPONENTS = (
    ("spatial_centered", "centered spatial DG", "black", "-"),
    ("surface_centered", "centered surface", "tab:orange", "--"),
    ("volume", "volume", "tab:purple", ":"),
    ("spatial_closure", "spatial closure", "tab:red", "-."),
)


def read_ledger(path: Path) -> np.ndarray:
    header = None
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("# stage explicit_rhs_used"):
                header = line[2:].split()
                break
    if header is None:
        raise RuntimeError(f"Could not find column header in {path}")
    return np.genfromtxt(path, names=header, comments="#")


def symlog_axis(axis: plt.Axes, values: list[np.ndarray]) -> None:
    maximum = max((np.nanmax(np.abs(value)) for value in values), default=0.0)
    if np.isfinite(maximum) and maximum > 0.0:
        axis.set_yscale("symlog", linthresh=max(maximum * 1.0e-8, 1.0e-300))
    axis.axhline(0.0, color="0.7", linewidth=0.7)
    axis.grid(alpha=0.2)


def plot_family(
    data: np.ndarray,
    stages: np.ndarray,
    species: int,
    components: tuple[tuple[str, str, str, str], ...],
    output: Path,
) -> None:
    figure, axes = plt.subplots(
        len(QUANTITIES), len(stages),
        figsize=(5.0 * len(stages), 3.4 * len(QUANTITIES)),
        sharex="col", squeeze=False,
    )

    for column, stage in enumerate(stages):
        selected = data[(data["stage"] == stage) & (data["species"] == species)]
        selected = np.sort(selected, order="r_cm")
        radius = selected["r_cm"] / 1.0e5
        used = bool(np.max(selected["explicit_rhs_used"]))

        for row, quantity in enumerate(QUANTITIES):
            axis = axes[row, column]
            plotted = []
            for suffix, label, color, linestyle in components:
                field = f"{quantity}_{suffix}_erg_cm3_s"
                values = selected[field]
                plotted.append(values)
                axis.plot(
                    radius, values, label=label, color=color,
                    linestyle=linestyle, linewidth=1.2,
                )
            symlog_axis(axis, plotted)
            if row == 0:
                axis.set_title(f"stage {int(stage)}; explicit RHS used={used}")
            if column == 0:
                axis.set_ylabel(rf"${quantity}$ [erg cm$^{{-3}}$ s$^{{-1}}$]")
            if row == len(QUANTITIES) - 1:
                axis.set_xlabel("radius [km]")
            if row == 0 and column == len(stages) - 1:
                axis.legend(loc="best", fontsize="small")

    figure.suptitle(f"explicit radiation RHS decomposition; species={species}")
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def report_closure(data: np.ndarray, stages: np.ndarray, species: int) -> None:
    for stage in stages:
        selected = data[(data["stage"] == stage) & (data["species"] == species)]
        used = bool(np.max(selected["explicit_rhs_used"]))
        for quantity in QUANTITIES:
            spatial = selected[f"{quantity}_spatial_centered_erg_cm3_s"]
            spatial_closure = selected[f"{quantity}_spatial_closure_erg_cm3_s"]
            total = selected[f"{quantity}_total_erg_cm3_s"]
            total_closure = selected[f"{quantity}_total_closure_erg_cm3_s"]
            scale = max(np.nanmax(np.abs(spatial)), np.finfo(float).tiny)
            print(
                f"stage={int(stage)} used={used} species={species} {quantity}: "
                f"max|spatial closure|={np.nanmax(np.abs(spatial_closure)):.6e} "
                f"relative={np.nanmax(np.abs(spatial_closure)) / scale:.6e}; "
                f"max|total closure|={np.nanmax(np.abs(total_closure)):.6e} "
                f"relative={np.nanmax(np.abs(total_closure)) / max(np.nanmax(np.abs(total)), np.finfo(float).tiny):.6e}"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "diagnostic", type=Path,
        help="main matter_diagnostics*.dat file or its .imex_streaming ledger",
    )
    parser.add_argument(
        "--species", type=int, default=0,
        help="0 for the sum over species; otherwise the one-based species index",
    )
    parser.add_argument("--output-prefix", type=Path)
    args = parser.parse_args()

    ledger = args.diagnostic
    if not str(ledger).endswith(".imex_streaming"):
        ledger = Path(str(ledger) + ".imex_streaming")
    data = read_ledger(ledger)

    species_values = np.unique(data["species"]).astype(int)
    if args.species not in species_values:
        raise ValueError(
            f"species={args.species} is unavailable; choices are {species_values.tolist()}"
        )

    stages = np.unique(data["stage"]).astype(int)
    prefix = args.output_prefix
    if prefix is None:
        prefix = ledger.with_suffix("")

    plot_family(
        data, stages, args.species, COMPONENTS,
        Path(f"{prefix}_components_s{args.species}.png"),
    )
    plot_family(
        data, stages, args.species, DG_COMPONENTS,
        Path(f"{prefix}_dg_balance_s{args.species}.png"),
    )
    report_closure(data, stages, args.species)


if __name__ == "__main__":
    main()
