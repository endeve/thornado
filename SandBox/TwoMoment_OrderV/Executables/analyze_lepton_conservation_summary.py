#!/usr/bin/env python3

"""Plot the complete time-dependent electron-lepton conservation budget."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_history(path: Path) -> np.ndarray:
    data = np.atleast_1d(np.genfromtxt(path, names=True))
    if data.size < 2:
        raise RuntimeError(f"{path}: conservation record has no completed step")
    required = (
        "time_s",
        "matter_lepton",
        "radiation_lepton",
        "total_lepton",
        "cumulative_boundary_lepton",
        "step_boundary_lepton",
    )
    missing = [name for name in required if name not in data.dtype.names]
    if missing:
        raise RuntimeError(f"{path}: missing columns: {', '.join(missing)}")
    return data


def budget(data: np.ndarray) -> dict[str, np.ndarray | float]:
    initial = data[0]
    completed = data[1:]
    total_initial = float(initial["total_lepton"])
    matter_initial = float(initial["matter_lepton"])
    radiation_initial = float(initial["radiation_lepton"])

    cumulative_matter = completed["matter_lepton"] - matter_initial
    cumulative_radiation_on_grid = (
        completed["radiation_lepton"] - radiation_initial
    )
    cumulative_off_grid = completed["cumulative_boundary_lepton"]
    cumulative_radiation = cumulative_radiation_on_grid - cumulative_off_grid
    cumulative_residual = cumulative_matter + cumulative_radiation_on_grid \
        - cumulative_off_grid

    matter_all = data["matter_lepton"]
    radiation_all = data["radiation_lepton"]
    step_matter = np.diff(matter_all)
    step_radiation = (
        np.diff(radiation_all) - completed["step_boundary_lepton"]
    )
    step_residual = step_matter + step_radiation
    step_exchange = np.abs(step_matter) + np.abs(step_radiation)
    step_unclosed_fraction = np.abs(step_residual) / np.maximum(
        step_exchange, np.finfo(float).tiny
    )

    return {
        "time_ms": completed["time_s"] * 1.0e3,
        "initial": total_initial,
        "cumulative_matter": cumulative_matter,
        "cumulative_radiation_on_grid": cumulative_radiation_on_grid,
        "cumulative_off_grid": cumulative_off_grid,
        "cumulative_radiation": cumulative_radiation,
        "cumulative_residual": cumulative_residual,
        "step_matter": step_matter,
        "step_radiation": step_radiation,
        "step_residual": step_residual,
        "step_unclosed_fraction": step_unclosed_fraction,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("standard", type=Path)
    parser.add_argument("source_aware", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    results = {
        "standard LLF": budget(read_history(args.standard)),
        "source-aware AP LLF": budget(read_history(args.source_aware)),
    }
    method_colors = {
        "standard LLF": "tab:blue",
        "source-aware AP LLF": "tab:orange",
    }

    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "cm",
            "axes.grid": True,
            "grid.alpha": 0.25,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), sharex="col")

    for axis, (method, values) in zip(axes[0], results.items()):
        initial = float(values["initial"])
        time_ms = np.asarray(values["time_ms"])
        axis.plot(
            time_ms,
            np.asarray(values["cumulative_matter"]) / initial,
            color="black",
            label=r"fluid: $\Delta L_{\rm matter}$",
        )
        axis.plot(
            time_ms,
            np.asarray(values["cumulative_radiation_on_grid"]) / initial,
            color="tab:blue",
            label=r"neutrinos: $\Delta L_\nu$",
        )
        axis.plot(
            time_ms,
            np.asarray(values["cumulative_off_grid"]) / initial,
            color="tab:red",
            label=r"off-grid tally: $L_{\rm off\,grid}$",
        )
        axis.plot(
            time_ms,
            np.asarray(values["cumulative_residual"]) / initial,
            color="black",
            linestyle=":",
            linewidth=1.8,
            label=r"total: $R_L=\Delta L_{\rm matter}+\Delta L_\nu-L_{\rm off\,grid}$",
        )
        axis.axhline(0.0, color="0.4", linewidth=0.7)
        axis.set_title(method)
        axis.set_ylabel(r"cumulative change$/L_{\rm tot}(0)$")
        axis.legend(frameon=False, fontsize=8)

    for method, values in results.items():
        time_ms = np.asarray(values["time_ms"])
        initial = float(values["initial"])
        axes[1, 0].semilogy(
            time_ms,
            np.maximum(
                np.abs(np.asarray(values["cumulative_residual"]) / initial),
                np.finfo(float).tiny,
            ),
            color=method_colors[method],
            label=method,
        )
        axes[1, 1].semilogy(
            time_ms,
            np.maximum(
                np.asarray(values["step_unclosed_fraction"]),
                np.finfo(float).tiny,
            ),
            color=method_colors[method],
            label=method,
        )

    axes[1, 0].set_title("Accumulated global conservation error")
    axes[1, 0].set_ylabel(r"$|R_L(t)|/|L_{\rm tot}(0)|$")
    axes[1, 0].set_xlabel(r"$t\,[\mathrm{ms}]$")
    axes[1, 0].legend(frameon=False)
    axes[1, 1].set_title("Unclosed fraction of each timestep's exchange")
    axes[1, 1].set_ylabel(
        r"$|R_L^n|/(|\Delta L_{\rm matter}^n|"
        r"+|\Delta L_\nu^n-L_{\rm off\,grid}^n|)$"
    )
    axes[1, 1].set_xlabel(r"$t\,[\mathrm{ms}]$")
    axes[1, 1].legend(frameon=False)

    standard = results["standard LLF"]
    source_aware = results["source-aware AP LLF"]
    reduction = (
        abs(np.asarray(standard["cumulative_residual"])[-1])
        / abs(np.asarray(source_aware["cumulative_residual"])[-1])
    )
    fig.suptitle(
        r"Electron-lepton budget: "
        r"$R_L=\Delta L_{\rm matter}+\Delta L_\nu-L_{\rm off\,grid}$"
    )
    fig.text(
        0.5,
        0.012,
        r"The off-grid sign follows Thornado's tally definition: "
        r"$L_{\nu,\rm change}=\Delta L_\nu-L_{\rm off\,grid}$.  "
        f"At the final time, the source-aware AP LLF reduces "
        f"the accumulated residual by {reduction:.0f}x.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0.0, 0.045, 1.0, 0.96))
    fig.savefig(args.output)

    for method, values in results.items():
        initial = float(values["initial"])
        matter = np.asarray(values["cumulative_matter"])[-1] / initial
        radiation = np.asarray(values["cumulative_radiation"])[-1] / initial
        residual = np.asarray(values["cumulative_residual"])[-1] / initial
        unclosed = np.asarray(values["step_unclosed_fraction"])
        print(method)
        print(f"  final matter change / L0             = {matter:+.6e}")
        print(f"  final radiation-boundary change / L0 = {radiation:+.6e}")
        print(f"  final residual / L0                  = {residual:+.6e}")
        print(f"  median per-step unclosed fraction    = {np.median(unclosed):.6e}")
        print(f"  final per-step unclosed fraction     = {unclosed[-1]:.6e}")
    print(f"figure: {args.output}")


if __name__ == "__main__":
    main()
