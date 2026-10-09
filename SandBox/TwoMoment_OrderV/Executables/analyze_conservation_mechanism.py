#!/usr/bin/env python3

"""Compare global conservation and its final-step collision origin."""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_budget_module(script: Path):
    spec = importlib.util.spec_from_file_location("discrete_budget_analysis", script)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {script}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def read_conservation(path: Path) -> np.ndarray:
    data = np.atleast_1d(np.genfromtxt(path, names=True))
    if data.size < 2:
        raise RuntimeError(f"{path}: no completed timesteps")
    return data[1:]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("standard_conservation", type=Path)
    parser.add_argument("ap_conservation", type=Path)
    parser.add_argument("standard_budget", type=Path)
    parser.add_argument("ap_budget", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    budget_script = Path(__file__).with_name("analyze_discrete_step_budget.py")
    budget = load_budget_module(budget_script)

    histories = {
        "standard LLF": read_conservation(args.standard_conservation),
        "source-aware AP LLF": read_conservation(args.ap_conservation),
    }
    ledgers = {
        "standard LLF": budget.read_ledger(args.standard_budget),
        "source-aware AP LLF": budget.read_ledger(args.ap_budget),
    }
    colors = {
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
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))

    component_names = (
        "FULL_WEIGHTED_STREAMING",
        "FULL_WEIGHTED_COLLISION_RADIATION",
        "FULL_WEIGHTED_COLLISION_MATTER",
    )
    component_labels = (
        "streaming\nradiation",
        "collision\nradiation",
        "collision\nmatter",
    )
    x = np.arange(len(component_names), dtype=float)
    width = 0.36
    signed_scale = 1.0e56
    collision_sums = {}

    for index, (label, ledger) in enumerate(ledgers.items()):
        _, _, volume = budget.rows(ledger, "FULL_STEP_MATTER_ACTUAL")
        signed = np.asarray(
            [
                np.sum(budget.field(ledger, name, "QL") * volume)
                for name in component_names
            ]
        )
        collision_sums[label] = signed[1] + signed[2]
        bars = axes[0].bar(
            x + (index - 0.5) * width,
            signed / signed_scale,
            width,
            color=colors[label],
            label=label,
        )
        axes[0].bar_label(bars, fmt="%.2g", fontsize=8, padding=2)

    axes[0].axhline(0.0, color="0.3", linewidth=0.8)
    axes[0].set_xticks(x, component_labels)
    axes[0].set_ylabel(r"signed rate $[10^{56}\,\mathrm{s}^{-1}]$")
    axes[0].set_title(
        "Terms in the global lepton budget\n"
        r"collision conservation requires "
        r"$\dot L_{\nu,\mathrm{coll}}+\dot L_{\mathrm{matter}}=0$"
    )
    axes[0].legend(frameon=False)

    ratio_labels = (
        "grid-scale LLF response /\n"
        "grid-scale centered response",
        "unclosed collision exchange /\n"
        "applied matter exchange",
    )
    ratio_values = {}
    for index, (label, ledger) in enumerate(ledgers.items()):
        stages = sorted(set(ledger.stage[ledger.stage > 0]))
        _, _, volume = budget.rows(ledger, "FULL_STEP_MATTER_ACTUAL")
        streaming = budget.field(ledger, "FULL_WEIGHTED_STREAMING", "QL")
        llf = budget.weighted_stage_sum(
            ledger, stages, "STREAMING_LLF", "QL", "explicit"
        )
        centered = budget.weighted_stage_sum(
            ledger, stages, "STREAMING_CENTERED", "QL", "explicit"
        )
        llf_norm = np.sum(np.abs(budget.high_pass(llf)) * volume)
        centered_norm = np.sum(np.abs(budget.high_pass(centered)) * volume)
        collision_residual = budget.field(
            ledger, "FULL_ALGEBRAIC_CLOSURE", "QL"
        )
        applied_collision = budget.field(
            ledger, "FULL_STEP_MATTER_ACTUAL", "QL"
        )
        ratios = np.asarray(
            [
                llf_norm / max(centered_norm, np.finfo(float).tiny),
                np.sum(np.abs(collision_residual) * volume)
                / max(
                    np.sum(np.abs(applied_collision) * volume),
                    np.finfo(float).tiny,
                ),
            ]
        )
        ratio_values[label] = ratios
        bars = axes[1].bar(
            np.arange(2) + (index - 0.5) * width,
            ratios,
            width,
            color=colors[label],
            label=label,
        )
        axes[1].bar_label(
            bars,
            labels=[f"{value:.3g}" for value in ratios],
            fontsize=8,
            padding=2,
        )

    axes[1].set_xticks(np.arange(2), ratio_labels)
    axes[1].set_yscale("log")
    axes[1].set_ylabel(r"volume-weighted $L^1$ ratio")
    axes[1].set_title(
        "Grid-scale streaming response and resulting\n"
        "collision-stage conservation error"
    )
    axes[1].legend(frameon=False)

    final_time_ms = histories["standard LLF"]["time_s"][-1] * 1.0e3
    fig.suptitle(
        f"Final-timestep operator diagnostic at $t={final_time_ms:.3f}$ ms"
    )
    standard = "standard LLF"
    ap = "source-aware AP LLF"
    cumulative_reduction = (
        histories[standard]["cumulative_lepton_relative"][-1]
        / histories[ap]["cumulative_lepton_relative"][-1]
    )
    llf_reduction = ratio_values[standard][0] / ratio_values[ap][0]
    collision_reduction = ratio_values[standard][1] / ratio_values[ap][1]
    fig.text(
        0.5,
        0.012,
        "Detailed operator decomposition: final timestep only.\n"
        f"AP/standard reductions: grid-scale LLF response {llf_reduction:.1f}x; "
        f"collision imbalance {collision_reduction:.0f}x.  "
        f"The inventory history gives a {cumulative_reduction:.0f}x reduction "
        "in accumulated lepton error.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0.0, 0.10, 1.0, 0.94))
    fig.savefig(args.output)
    print(f"figure: {args.output}")


if __name__ == "__main__":
    main()
