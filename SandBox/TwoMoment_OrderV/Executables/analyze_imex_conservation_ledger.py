#!/usr/bin/env python3

"""Validate and plot the per-timestep IMEX conservation decomposition."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


COMPONENTS = (
    "collision_pair_residual",
    "transport_boundary_residual",
    "matter_remainder",
    "radiation_remainder",
    "conservation_residual",
)

COLORS = {
    "collision_pair_residual": "tab:red",
    "transport_boundary_residual": "tab:blue",
    "matter_remainder": "tab:orange",
    "radiation_remainder": "tab:green",
    "conservation_residual": "black",
}

LABELS = {
    "collision_pair_residual": r"implicit matter+radiation",
    "transport_boundary_residual": r"explicit update$-$boundary tally",
    "matter_remainder": r"matter assembly/limiter remainder",
    "radiation_remainder": r"radiation assembly/limiter remainder",
    "conservation_residual": r"measured full-step residual",
}


def read_ledger(path: Path) -> np.ndarray:
    data = np.atleast_1d(
        np.genfromtxt(
            path,
            names=True,
            dtype=None,
            encoding="utf-8",
        )
    )
    required = {
        "cycle",
        "time_s",
        "dt_s",
        "quantity",
        "component",
        "step_value",
        "cumulative_value",
    }
    missing = required.difference(data.dtype.names or ())
    if missing:
        raise RuntimeError(f"{path}: missing columns: {', '.join(sorted(missing))}")
    if data.size == 0:
        raise RuntimeError(f"{path}: empty IMEX conservation ledger")
    return data


def series(data: np.ndarray, quantity: str, component: str) -> np.ndarray:
    selected = data[
        (data["quantity"] == quantity) & (data["component"] == component)
    ]
    if selected.size == 0:
        raise RuntimeError(f"missing {quantity}/{component} ledger rows")
    return selected


def validate_ledger(path: Path, data: np.ndarray) -> None:
    for quantity in ("energy", "lepton"):
        closure = series(data, quantity, "algebraic_closure")["step_value"]
        scales = np.zeros_like(closure)
        for component in COMPONENTS:
            scales += np.abs(series(data, quantity, component)["step_value"])
        relative = np.max(np.abs(closure) / np.maximum(scales, np.finfo(float).tiny))
        if relative > 5.0e-12:
            raise RuntimeError(
                f"{path}: {quantity} algebraic closure failed: {relative:.3e}"
            )


def validate_inventory(path: Path, data: np.ndarray, history_path: Path) -> None:
    history = np.atleast_1d(np.genfromtxt(history_path, names=True))[1:]
    cycles = series(data, "energy", "conservation_residual")["cycle"]
    if history.size != cycles.size or not np.array_equal(history["cycle"], cycles):
        raise RuntimeError(f"{path}: ledger and inventory histories have different cycles")

    checks = (
        ("energy", "step_energy_residual_erg", "cumulative_energy_residual_erg"),
        ("lepton", "step_lepton_residual", "cumulative_lepton_residual"),
    )
    for quantity, step_field, cumulative_field in checks:
        residual = series(data, quantity, "conservation_residual")
        step_scale = np.maximum(
            np.maximum(np.abs(residual["step_value"]), np.abs(history[step_field])),
            np.finfo(float).tiny,
        )
        cumulative_scale = np.maximum(
            np.maximum(
                np.abs(residual["cumulative_value"]),
                np.abs(history[cumulative_field]),
            ),
            np.finfo(float).tiny,
        )
        step_error = np.max(
            np.abs(residual["step_value"] - history[step_field]) / step_scale
        )
        cumulative_error = np.max(
            np.abs(residual["cumulative_value"] - history[cumulative_field])
            / cumulative_scale
        )
        if step_error > 5.0e-11 or cumulative_error > 5.0e-11:
            raise RuntimeError(
                f"{path}: {quantity} ledger disagrees with inventory history: "
                f"step={step_error:.3e}, cumulative={cumulative_error:.3e}"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("ledger", type=Path, nargs="+")
    parser.add_argument("--labels", nargs="+")
    parser.add_argument("--inventory-history", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.labels is not None and len(args.labels) != len(args.ledger):
        raise RuntimeError("--labels must provide one label per ledger")
    if args.inventory_history is not None and len(args.inventory_history) != len(
        args.ledger
    ):
        raise RuntimeError(
            "--inventory-history must provide one history per ledger"
        )

    labels = args.labels or [path.stem for path in args.ledger]
    ledgers = [read_ledger(path) for path in args.ledger]
    for index, (path, data) in enumerate(zip(args.ledger, ledgers)):
        validate_ledger(path, data)
        if args.inventory_history is not None:
            validate_inventory(path, data, args.inventory_history[index])

        print(labels[index])
        for quantity in ("energy", "lepton"):
            print(f"  {quantity} final cumulative decomposition:")
            for component in COMPONENTS:
                value = series(data, quantity, component)["cumulative_value"][-1]
                unit = "erg" if quantity == "energy" else ""
                print(f"    {component:29s} {value:+.8e} {unit}")

    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "cm",
            "axes.grid": True,
            "grid.alpha": 0.25,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.0), sharex="col")
    method_styles = ("-", "--", ":", "-.")

    for row, quantity in enumerate(("energy", "lepton")):
        for data, method_style in zip(ledgers, method_styles):
            for component in COMPONENTS:
                values = series(data, quantity, component)
                time_ms = values["time_s"] * 1.0e3
                rate = values["step_value"] / values["dt_s"]
                axes[row, 0].plot(
                    time_ms,
                    rate,
                    color=COLORS[component],
                    linestyle=method_style,
                    linewidth=1.1,
                )
                axes[row, 1].plot(
                    time_ms,
                    values["cumulative_value"],
                    color=COLORS[component],
                    linestyle=method_style,
                    linewidth=1.1,
                )

        rate_values = np.concatenate(
            [
                np.abs(
                    series(data, quantity, component)["step_value"]
                    / series(data, quantity, component)["dt_s"]
                )
                for data in ledgers
                for component in COMPONENTS
            ]
        )
        cumulative_values = np.concatenate(
            [
                np.abs(series(data, quantity, component)["cumulative_value"])
                for data in ledgers
                for component in COMPONENTS
            ]
        )
        positive_rate = rate_values[rate_values > 0.0]
        positive_cumulative = cumulative_values[cumulative_values > 0.0]
        if positive_rate.size:
            axes[row, 0].set_yscale(
                "symlog", linthresh=max(np.percentile(positive_rate, 5), 1.0e-30)
            )
        if positive_cumulative.size:
            axes[row, 1].set_yscale(
                "symlog",
                linthresh=max(np.percentile(positive_cumulative, 5), 1.0e-30),
            )

    axes[0, 0].set_title("Energy: contribution to each timestep")
    axes[0, 1].set_title("Energy: accumulated contribution")
    axes[1, 0].set_title("Electron lepton number: contribution to each timestep")
    axes[1, 1].set_title("Electron lepton number: accumulated contribution")
    axes[0, 0].set_ylabel(r"rate $[\mathrm{erg\,s^{-1}}]$")
    axes[0, 1].set_ylabel(r"change $[\mathrm{erg}]$")
    axes[1, 0].set_ylabel(r"rate $[\mathrm{s^{-1}}]$")
    axes[1, 1].set_ylabel(r"change")
    axes[1, 0].set_xlabel(r"$t\,[\mathrm{ms}]$")
    axes[1, 1].set_xlabel(r"$t\,[\mathrm{ms}]$")

    component_handles = [
        Line2D([0], [0], color=COLORS[name], label=LABELS[name])
        for name in COMPONENTS
    ]
    method_handles = [
        Line2D([0], [0], color="0.25", linestyle=style, label=label)
        for style, label in zip(method_styles, labels)
    ]
    fig.legend(
        handles=component_handles,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.38, 0.995),
        ncol=2,
        fontsize=8,
    )
    fig.legend(
        handles=method_handles,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.82, 0.995),
        ncol=1,
        fontsize=8,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.89))
    fig.savefig(args.output)
    print(f"figure: {args.output}")


if __name__ == "__main__":
    main()
