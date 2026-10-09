#!/usr/bin/env python3

"""Compare native-grid heating diagnostics at two spatial resolutions."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_discrete_step_budget import (
    field,
    high_pass,
    read_ledger,
    rows,
    validate,
    weighted_stage_sum,
)


def diagnostics(path: Path) -> dict[str, object]:
    ledger = read_ledger(path)
    stages = validate(ledger, 5.0e-12)
    _, radius_cm, volume_cm3 = rows(ledger, "FULL_STEP_MATTER_ACTUAL")
    return {
        "ledger": ledger,
        "stages": stages,
        "radius_km": radius_cm / 1.0e5,
        "volume_cm3": volume_cm3,
        "heating": field(
            ledger, "FULL_STEP_MATTER_ACTUAL", "QH_exact"
        ),
        "surface": weighted_stage_sum(
            ledger, stages, "STREAMING_SURFACE", "QH", "explicit"
        ),
        "llf": weighted_stage_sum(
            ledger, stages, "STREAMING_LLF", "QH", "explicit"
        ),
        "radiation": field(
            ledger, "FULL_STEP_RADIATION_ACTUAL", "QH"
        ),
        "streaming": field(
            ledger, "FULL_WEIGHTED_STREAMING", "QH"
        ),
    }


def plot_global_balance(
    datasets: tuple[dict[str, object], ...],
    labels: tuple[str, ...],
    output: Path,
    rmax_km: float,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(15, 9), sharex=True)
    for data, label in zip(datasets, labels):
        radius = data["radius_km"]
        volume = data["volume_cm3"]
        matter_shell = data["heating"] * volume
        radiation_shell = data["radiation"] * volume
        streaming_shell = data["streaming"] * volume

        axes[0, 0].plot(radius, np.cumsum(matter_shell), label=label)
        axes[0, 1].plot(
            radius,
            np.cumsum(np.maximum(matter_shell, 0.0)),
            label=f"heating, {label}",
        )
        axes[0, 1].plot(
            radius,
            np.cumsum(np.minimum(matter_shell, 0.0)),
            linestyle="--",
            label=f"cooling, {label}",
        )
        axes[1, 0].plot(
            radius,
            np.cumsum(matter_shell + radiation_shell),
            label=f"matter + radiation, {label}",
        )
        axes[1, 0].plot(
            radius,
            np.cumsum(streaming_shell),
            linestyle="--",
            label=f"streaming, {label}",
        )
        axes[1, 1].plot(
            radius,
            np.cumsum(np.abs(matter_shell)),
            label=label,
        )

    titles = (
        "A. Cumulative signed fluid heating",
        "B. Cumulative heating and cooling",
        "C. Total-system change and streaming",
        r"D. Cumulative $L^1$ fluid heating",
    )
    for axis, title in zip(axes.flat, titles):
        axis.axhline(0.0, color="0.5", linewidth=0.8)
        axis.set_xlim(0.0, rmax_km)
        axis.set_title(title)
        axis.set_ylabel("integrated rate [erg s$^{-1}$]")
        axis.grid(alpha=0.2)
        axis.legend(loc="best", fontsize=8)
    for axis in axes[-1, :]:
        axis.set_xlabel("radius [km]")
    fig.suptitle("Native-grid cumulative matter-radiation energy budget", fontsize=14)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_stage_coupling(
    datasets: tuple[dict[str, object], ...],
    labels: tuple[str, ...],
    output: Path,
    rmax_km: float,
    explicit_stage: int,
    implicit_stage: int,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(15, 9), sharex=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for index, (data, label) in enumerate(zip(datasets, labels)):
        ledger = data["ledger"]
        radius = data["radius_km"]
        volume = data["volume_cm3"]
        total = field(
            ledger, "STREAMING_TOTAL", "QH", stage=explicit_stage
        )
        surface = field(
            ledger, "STREAMING_SURFACE", "QH", stage=explicit_stage
        )
        llf = field(
            ledger, "STREAMING_LLF", "QH", stage=explicit_stage
        )
        collision = field(
            ledger,
            "COLLISION_MATTER_ACTUAL",
            "QH_exact",
            stage=implicit_stage,
        )
        color = colors[index % len(colors)]

        axes[0, 0].plot(
            radius, total, color=color, label=f"total, {label}"
        )
        axes[0, 0].plot(
            radius,
            surface,
            color=color,
            linestyle="--",
            label=f"centered surface, {label}",
        )
        axes[0, 1].plot(
            radius,
            np.cumsum(total * volume),
            color=color,
            label=f"total, {label}",
        )
        axes[0, 1].plot(
            radius,
            np.cumsum(surface * volume),
            color=color,
            linestyle="--",
            label=f"centered surface, {label}",
        )
        axes[0, 1].plot(
            radius,
            np.cumsum(llf * volume),
            color=color,
            linestyle=":",
            label=f"LLF, {label}",
        )
        axes[1, 0].plot(radius, collision, color=color, label=label)
        axes[1, 1].plot(
            radius,
            np.cumsum(collision * volume),
            color=color,
            label=label,
        )

    titles = (
        f"A. Explicit stage {explicit_stage} streaming",
        f"B. Cumulative explicit stage {explicit_stage} streaming",
        f"C. Actual implicit stage {implicit_stage} matter exchange",
        f"D. Cumulative implicit stage {implicit_stage} matter exchange",
    )
    for axis, title in zip(axes.flat, titles):
        axis.axhline(0.0, color="0.5", linewidth=0.8)
        axis.set_xlim(0.0, rmax_km)
        axis.set_title(title)
        axis.grid(alpha=0.2)
        axis.legend(loc="best", fontsize=8)
    axes[0, 0].set_ylabel(r"rate [erg cm$^{-3}$ s$^{-1}$]")
    axes[1, 0].set_ylabel(r"rate [erg cm$^{-3}$ s$^{-1}$]")
    axes[0, 1].set_ylabel("integrated rate [erg s$^{-1}$]")
    axes[1, 1].set_ylabel("integrated rate [erg s$^{-1}$]")
    for axis in axes[-1, :]:
        axis.set_xlabel("radius [km]")
    fig.suptitle(
        "Resolution dependence of selected IMEX transport and collision stages",
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("coarse", type=Path)
    parser.add_argument("fine", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--global-output",
        type=Path,
        help="write cumulative volume-integrated matter/radiation balances",
    )
    parser.add_argument(
        "--stage-coupling-output",
        type=Path,
        help="compare selected explicit-streaming and implicit-collision stages",
    )
    parser.add_argument("--explicit-stage", type=int, default=1)
    parser.add_argument("--implicit-stage", type=int, default=2)
    parser.add_argument("--rmax-km", type=float, default=50.0)
    args = parser.parse_args()

    datasets = (diagnostics(args.coarse), diagnostics(args.fine))
    labels = tuple(f"nX={len(data['radius_km'])}" for data in datasets)

    fig, axes = plt.subplots(2, 2, figsize=(15, 9), sharex=True)
    quantities = (
        ("heating", "A. Actual applied fluid heating"),
        ("surface", "B. Centered DG surface contribution"),
        ("heating_high_pass", "C. Grid-scale part of applied heating"),
        ("llf", "D. LLF contribution after source-aware scaling"),
    )

    for axis, (quantity, title) in zip(axes.flat, quantities):
        for data, label in zip(datasets, labels):
            values = (
                high_pass(data["heating"])
                if quantity == "heating_high_pass"
                else data[quantity]
            )
            axis.plot(
                data["radius_km"],
                values,
                marker="o",
                markersize=2.0,
                linewidth=1.0,
                label=label,
            )
        axis.axhline(0.0, color="0.5", linewidth=0.8)
        axis.set_xlim(0.0, args.rmax_km)
        axis.set_title(title)
        axis.set_ylabel(r"rate [erg cm$^{-3}$ s$^{-1}$]")
        axis.grid(alpha=0.2)
        axis.legend(loc="best")

    for axis in axes[-1, :]:
        axis.set_xlabel("radius [km]")
    fig.suptitle(
        "Native element-average heating profiles; no restriction or averaging",
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(args.output, dpi=180)
    plt.close(fig)
    if args.global_output is not None:
        plot_global_balance(
            datasets, labels, args.global_output, args.rmax_km
        )
    if args.stage_coupling_output is not None:
        plot_stage_coupling(
            datasets,
            labels,
            args.stage_coupling_output,
            args.rmax_km,
            args.explicit_stage,
            args.implicit_stage,
        )


if __name__ == "__main__":
    main()
