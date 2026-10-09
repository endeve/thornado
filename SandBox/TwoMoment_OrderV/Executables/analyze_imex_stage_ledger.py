#!/usr/bin/env python3

"""Summarize and plot the one-step IMEX stage ledger."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


STATE_NAMES = {
    1: "assembled",
    2: "pre-implicit",
    3: "post-implicit",
    4: "post-limiter",
    5: "completed-step",
}


def read_ledger(path: Path) -> np.ndarray:
    header = None
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("# stage state element "):
                header = line[2:].split()
                break
    if header is None:
        raise ValueError(f"No column header found in {path}")

    values = np.loadtxt(path, comments="#", ndmin=2)
    if values.shape[1] != len(header):
        raise ValueError(
            f"Header has {len(header)} columns, data have {values.shape[1]}"
        )
    return np.rec.fromarrays(values.T, names=header)


def select(data: np.ndarray, stage: int, state: int) -> np.ndarray:
    mask = (data.stage.astype(int) == stage) & (data.state.astype(int) == state)
    rows = data[mask]
    return rows[np.argsort(rows.element)]


def weighted_l1(values: np.ndarray, volumes: np.ndarray) -> float:
    return float(np.sum(np.abs(values) * volumes))


def signed_integral(values: np.ndarray, volumes: np.ndarray) -> float:
    return float(np.sum(values * volumes))


def report(data: np.ndarray) -> None:
    stages = np.unique(data.stage.astype(int))
    erad_names = [name for name in data.dtype.names if name.startswith("Erad_s")]
    previous_limited = None

    print("IMEX one-step stage ledger")
    for stage in stages:
        snapshots = {
            state: select(data, stage, state)
            for state in STATE_NAMES
            if np.any(
                (data.stage.astype(int) == stage)
                & (data.state.astype(int) == state)
            )
        }
        pre = snapshots[2]
        post = snapshots[3]
        limited = snapshots[4]
        assembled = snapshots[1]
        volume = pre.volume_cm3
        implicit_dt = float(pre.implicit_dt_s[0])

        print(f"\n  stage {stage}: implicit dt = {implicit_dt:.8e} s")
        for state, rows in snapshots.items():
            qh = rows.QH_total_erg_cm3_s
            print(
                f"    {STATE_NAMES[state]:13s} QH: "
                f"integral={signed_integral(qh, volume): .6e} erg/s  "
                f"L1={weighted_l1(qh, volume):.6e} erg/s"
            )

        if not all(
            np.all(np.isfinite(rows.eint_erg_cm3))
            for rows in snapshots.values()
        ):
            print(
                "    WARNING: non-finite stage internal energy; "
                "RHS snapshots remain usable, but state-difference diagnostics do not"
            )
            previous_limited = None
            continue

        if previous_limited is not None:
            fluid_delta = assembled.eint_erg_cm3 - previous_limited.eint_erg_cm3
            radiation_delta = np.zeros(assembled.shape)
            for name in erad_names:
                radiation_delta += assembled[name] - previous_limited[name]
            print(
                "    stage assembly/tableau: "
                f"fluid L1={weighted_l1(fluid_delta, volume):.6e} erg  "
                f"radiation L1={weighted_l1(radiation_delta, volume):.6e} erg"
            )

        if implicit_dt != 0.0:
            applied = (post.eint_erg_cm3 - pre.eint_erg_cm3) / implicit_dt
            rhs_scale = max(
                weighted_l1(applied, volume),
                weighted_l1(pre.QH_total_erg_cm3_s, volume),
                np.finfo(float).tiny,
            )
            mismatch_pre = weighted_l1(
                applied - pre.QH_total_erg_cm3_s, volume
            ) / rhs_scale
            mismatch_post = weighted_l1(
                applied - post.QH_total_erg_cm3_s, volume
            ) / rhs_scale
            print(
                "    implicit fluid increment: "
                f"integral={signed_integral(applied, volume): .6e} erg/s  "
                f"L1={weighted_l1(applied, volume):.6e} erg/s"
            )
            print(
                "    relative L1 mismatch to instantaneous QH: "
                f"pre={mismatch_pre:.6e} post={mismatch_post:.6e}"
            )

            radiation_delta = np.zeros(post.shape)
            for name in erad_names:
                radiation_delta += (post[name] - pre[name]) / implicit_dt
            closure = applied + radiation_delta
            closure_scale = max(
                weighted_l1(applied, volume),
                weighted_l1(radiation_delta, volume),
                np.finfo(float).tiny,
            )
            print(
                "    comoving fluid+radiation collision closure: "
                f"integral={signed_integral(closure, volume): .6e} erg/s  "
                f"relative L1={weighted_l1(closure, volume) / closure_scale:.6e}"
            )

        for label, before, after in (
            ("pre-stage limiter", assembled, pre),
            ("post-implicit limiter", post, limited),
        ):
            fluid_delta = after.eint_erg_cm3 - before.eint_erg_cm3
            radiation_delta = np.zeros(after.shape)
            for name in erad_names:
                radiation_delta += after[name] - before[name]
            scale = max(
                weighted_l1(fluid_delta, volume),
                weighted_l1(radiation_delta, volume),
                np.finfo(float).tiny,
            )
            print(
                f"    {label:22s}: "
                f"fluid L1={weighted_l1(fluid_delta, volume):.6e} erg  "
                f"radiation L1={weighted_l1(radiation_delta, volume):.6e} erg  "
                f"combined rel-L1={weighted_l1(fluid_delta + radiation_delta, volume) / scale:.6e}"
            )

        if 5 in snapshots:
            completed = snapshots[5]
            fluid_delta = completed.eint_erg_cm3 - limited.eint_erg_cm3
            radiation_delta = np.zeros(completed.shape)
            for name in erad_names:
                radiation_delta += completed[name] - limited[name]
            print(
                "    final assembly/limiter: "
                f"fluid L1={weighted_l1(fluid_delta, volume):.6e} erg  "
                f"radiation L1={weighted_l1(radiation_delta, volume):.6e} erg"
            )

        previous_limited = limited


def plot(data: np.ndarray, output: Path) -> None:
    stages = np.unique(data.stage.astype(int))
    fig, axes = plt.subplots(
        len(stages), 3, figsize=(16, 4.3 * len(stages)), squeeze=False
    )

    previous_limited = None
    for row, stage in enumerate(stages):
        snapshots = {
            state: select(data, stage, state)
            for state in STATE_NAMES
            if np.any(
                (data.stage.astype(int) == stage)
                & (data.state.astype(int) == state)
            )
        }
        for state, rows in snapshots.items():
            radius_km = rows.r_cm / 1.0e5
            axes[row, 0].plot(
                radius_km,
                rows.QH_total_erg_cm3_s,
                label=STATE_NAMES[state],
            )

        pre = snapshots[2]
        post = snapshots[3]
        limited = snapshots[4]
        radius_km = pre.r_cm / 1.0e5
        implicit_dt = float(pre.implicit_dt_s[0])

        if implicit_dt != 0.0:
            applied = (post.eint_erg_cm3 - pre.eint_erg_cm3) / implicit_dt
            axes[row, 1].plot(radius_km, applied, label="implicit fluid increment")
        axes[row, 1].plot(
            radius_km, pre.QH_total_erg_cm3_s, label="RHS pre-implicit"
        )
        axes[row, 1].plot(
            radius_km, post.QH_total_erg_cm3_s, label="RHS post-implicit"
        )

        axes[row, 2].plot(
            radius_km,
            pre.eint_erg_cm3 - snapshots[1].eint_erg_cm3,
            label="pre-stage limiter",
        )
        axes[row, 2].plot(
            radius_km,
            limited.eint_erg_cm3 - post.eint_erg_cm3,
            label="post-implicit limiter",
        )
        if previous_limited is not None:
            axes[row, 2].plot(
                radius_km,
                snapshots[1].eint_erg_cm3 - previous_limited.eint_erg_cm3,
                label="stage assembly/tableau",
            )

        axes[row, 0].set_title(f"stage {stage}: instantaneous collision heating")
        axes[row, 1].set_title(f"stage {stage}: collision update")
        axes[row, 2].set_title(f"stage {stage}: non-collision state changes")
        axes[row, 0].set_ylabel(r"erg cm$^{-3}$ s$^{-1}$")
        axes[row, 1].set_ylabel(r"erg cm$^{-3}$ s$^{-1}$")
        axes[row, 2].set_ylabel(r"erg cm$^{-3}$")
        for axis in axes[row]:
            axis.set_xlabel("radius [km]")
            axis.set_xscale("log")
            axis.set_yscale("symlog", linthresh=1.0)
            axis.grid(alpha=0.25)
            axis.legend(fontsize="small", loc="best")

        previous_limited = limited

    fig.tight_layout()
    fig.savefig(output, dpi=180)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("ledger", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    data = read_ledger(args.ledger)
    report(data)
    output = args.output or args.ledger.with_suffix(args.ledger.suffix + ".png")
    plot(data, output)
    print(f"\nwrote {output}")


if __name__ == "__main__":
    main()
