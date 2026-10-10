#!/usr/bin/env python3

"""Check and plot the per-step matter--radiation conservation record."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


def read_history(path: Path) -> np.ndarray:
    data = np.genfromtxt(path, names=True)
    data = np.atleast_1d(data)
    if data.size < 2:
        raise RuntimeError(f"{path}: conservation record contains no completed step")

    lepton_prefix = (
        "electron_lepton"
        if "cumulative_electron_lepton_relative" in data.dtype.names
        else "lepton"
    )
    required = (
        "cycle",
        "time_s",
        "cumulative_energy_relative",
        "cumulative_energy_exchange_relative",
        "step_energy_relative",
        "step_energy_exchange_relative",
        f"cumulative_{lepton_prefix}_relative",
        f"cumulative_{lepton_prefix}_exchange_relative",
        f"step_{lepton_prefix}_relative",
        f"step_{lepton_prefix}_exchange_relative",
        "cumulative_energy_residual_erg",
        f"cumulative_{lepton_prefix}_residual",
    )
    missing = [name for name in required if name not in data.dtype.names]
    if missing:
        raise RuntimeError(f"{path}: missing columns: {', '.join(missing)}")

    for name in data.dtype.names:
        if not np.all(np.isfinite(data[name])):
            raise RuntimeError(f"{path}: non-finite values in {name}")
    return data[1:]


def lepton_prefix(data: np.ndarray) -> str:
    if "cumulative_electron_lepton_relative" in data.dtype.names:
        return "electron_lepton"
    return "lepton"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=Path, nargs="+")
    parser.add_argument(
        "--labels",
        nargs="+",
        help="labels for input histories; defaults to the file stems",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    histories = [read_history(path) for path in args.file]
    if args.labels is None:
        labels = [path.stem for path in args.file]
    else:
        if len(args.labels) != len(args.file):
            raise RuntimeError("--labels must provide one label per input file")
        labels = args.labels

    for label, completed in zip(labels, histories):
        print(f"{label}: steps={completed.size}")
        print(
            "  energy: max step/inventory="
            f"{np.max(completed['step_energy_relative']):.6e}, "
            "max cumulative/inventory="
            f"{np.max(completed['cumulative_energy_relative']):.6e}, "
            "max step/exchange="
            f"{np.max(completed['step_energy_exchange_relative']):.6e}, "
            "final cumulative/exchange="
            f"{completed['cumulative_energy_exchange_relative'][-1]:.6e}, "
            "final residual="
            f"{completed['cumulative_energy_residual_erg'][-1]:.6e} erg"
        )
        prefix = lepton_prefix(completed)
        print(
            "  lepton: max step/inventory="
            f"{np.max(completed[f'step_{prefix}_relative']):.6e}, "
            "max cumulative/inventory="
            f"{np.max(completed[f'cumulative_{prefix}_relative']):.6e}, "
            "max step/exchange="
            f"{np.max(completed[f'step_{prefix}_exchange_relative']):.6e}, "
            "final cumulative/exchange="
            f"{completed[f'cumulative_{prefix}_exchange_relative'][-1]:.6e}, "
            "final residual="
            f"{completed[f'cumulative_{prefix}_residual'][-1]:.6e}"
        )

    output = args.output
    if output is None:
        output = args.file[0].with_suffix(args.file[0].suffix + ".pdf")

    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "cm",
            "axes.grid": True,
            "grid.alpha": 0.25,
        }
    )
    fig, axes = plt.subplots(2, 1, figsize=(10.0, 7.0), sharex=True)
    line_styles = ("-", "--", ":", "-.")
    metric_specs = (
        (
            "cumulative_{quantity}_relative",
            r"$|R(t)|/$ initial total",
            "tab:blue",
        ),
        (
            "cumulative_{quantity}_exchange_relative",
            r"$|R(t)|/$ cumulative exchanged amount",
            "tab:orange",
        ),
        (
            "step_{quantity}_relative",
            r"$|R^n|/$ total at start of step",
            "tab:green",
        ),
        (
            "step_{quantity}_exchange_relative",
            r"$|R^n|/$ amount exchanged in step",
            "tab:red",
        ),
    )
    method_styles = [line_styles[index % len(line_styles)] for index in range(len(histories))]
    for axis_index, (axis, quantity) in enumerate(zip(axes, ("energy", "lepton"))):
        for completed, method_style in zip(histories, method_styles):
            time_ms = completed["time_s"] * 1.0e3
            field_quantity = (
                lepton_prefix(completed) if quantity == "lepton" else quantity
            )
            for template, _, color in metric_specs:
                field = template.format(quantity=field_quantity)
                axis.semilogy(
                    time_ms,
                    np.maximum(completed[field], np.finfo(float).tiny),
                    color=color,
                    linestyle=method_style,
                )
    axes[0].set_ylabel("relative conservation residual")
    axes[0].set_title(
        r"Energy: $R_E=\Delta(E_{\rm matter}+E_\nu)-E_{\rm off\,grid}$"
    )
    metric_handles = [
        Line2D([0], [0], color=color, label=label)
        for _, label, color in metric_specs
    ]
    method_handles = [
        Line2D([0], [0], color="black", linestyle=style, label=label)
        for label, style in zip(labels, method_styles)
    ]
    fig.legend(
        handles=metric_handles,
        title="Normalization",
        frameon=False,
        fontsize=8,
        loc="upper right",
        bbox_to_anchor=(0.995, 0.87),
    )
    fig.legend(
        handles=method_handles,
        title="Numerical flux",
        frameon=False,
        fontsize=8,
        loc="upper right",
        bbox_to_anchor=(0.995, 0.48),
    )
    axes[1].set_ylabel("relative conservation residual")
    axes[1].set_title(
        r"Electron lepton number: "
        r"$R_L=\Delta(L_{\rm matter}+L_\nu)-L_{\rm off\,grid}$"
    )
    axes[1].set_xlabel(r"$t\,[\mathrm{ms}]$")
    fig.tight_layout(rect=(0.0, 0.0, 0.79, 1.0))
    fig.savefig(output)
    print(f"figure: {output}")


if __name__ == "__main__":
    main()
