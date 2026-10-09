#!/usr/bin/env python3

"""Compare one-step matter exchange with an independently converged run."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_discrete_step_budget import (
    field,
    high_pass,
    infer_reference_state_path,
    norm_l1,
    origin_rates,
    read_ledger,
    read_reference_state,
    reference_reservoirs,
    rows,
    validate,
)


@dataclass(frozen=True)
class Case:
    budget_path: Path
    radius_cm: np.ndarray
    width_cm: np.ndarray
    volume_cm3: np.ndarray
    edges_cm: np.ndarray
    rho: np.ndarray
    eint: np.ndarray
    ne: np.ndarray
    dt_s: float
    origins: dict[str, dict[str, np.ndarray]]


def timestep(path: Path) -> float:
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("# timestep_s "):
                return float(line.split()[2])
    raise ValueError(f"{path}: missing timestep_s")


def cell_edges(radius: np.ndarray, width: np.ndarray) -> np.ndarray:
    left = radius - 0.5 * width
    right = radius + 0.5 * width
    scale = max(float(np.max(np.abs(right))), 1.0)
    if np.max(np.abs(right[:-1] - left[1:])) > 1.0e-10 * scale:
        raise ValueError("radial elements are not contiguous")
    edges = np.concatenate(([left[0]], right))
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError("radial element edges are not strictly increasing")
    return edges


def load_case(budget_path: Path, state_path: Path | None) -> Case:
    ledger = read_ledger(budget_path)
    stages = validate(ledger, 5.0e-12)
    if state_path is None:
        state_path = infer_reference_state_path(budget_path)
    state = read_reference_state(state_path)
    _, radius, volume = rows(ledger, "FULL_STEP_MATTER_ACTUAL")
    if len(state["r_cm"]) != len(radius) or not np.allclose(
        state["r_cm"], radius, rtol=1.0e-12, atol=0.0
    ):
        raise ValueError(f"{state_path}: state and budget grids differ")
    rho, eint, ne = reference_reservoirs(state)
    return Case(
        budget_path=budget_path,
        radius_cm=radius,
        width_cm=state["dr_cm"],
        volume_cm3=volume,
        edges_cm=cell_edges(radius, state["dr_cm"]),
        rho=rho,
        eint=eint,
        ne=ne,
        dt_s=timestep(budget_path),
        origins={
            quantity: origin_rates(ledger, stages, quantity)
            for quantity in ("QH", "QL")
        },
    )


def shell_volume(left: float, right: float) -> float:
    return 4.0 * np.pi / 3.0 * (right**3 - left**3)


def restrict_spherical_average(
    values: np.ndarray, source_edges: np.ndarray, target_edges: np.ndarray
) -> np.ndarray:
    """Conservatively restrict piecewise-constant spherical-shell averages."""

    tolerance = 1.0e-10 * max(abs(target_edges[-1]), 1.0)
    if (
        abs(source_edges[0] - target_edges[0]) > tolerance
        or abs(source_edges[-1] - target_edges[-1]) > tolerance
    ):
        raise ValueError("candidate and reference radial domains differ")
    result = np.zeros(len(target_edges) - 1)
    source = 0
    for target in range(len(result)):
        left = target_edges[target]
        right = target_edges[target + 1]
        while source + 1 < len(source_edges) and source_edges[source + 1] <= left:
            source += 1
        current = source
        integral = 0.0
        while current < len(values) and source_edges[current] < right:
            overlap_left = max(left, source_edges[current])
            overlap_right = min(right, source_edges[current + 1])
            if overlap_right > overlap_left:
                integral += values[current] * shell_volume(
                    overlap_left, overlap_right
                )
            if source_edges[current + 1] >= right:
                break
            current += 1
        result[target] = integral / shell_volume(left, right)
    return result


def mapped_reference(
    candidate: Case, reference: Case, quantity: str
) -> dict[str, np.ndarray]:
    return {
        name: restrict_spherical_average(
            values, reference.edges_cm, candidate.edges_cm
        )
        for name, values in reference.origins[quantity].items()
    }


def exchange_metrics(
    candidate: Case,
    expected: np.ndarray,
    numerical: np.ndarray,
    reservoir: np.ndarray,
    selected: np.ndarray,
) -> dict[str, float]:
    volume = candidate.volume_cm3[selected]
    expected = expected[selected]
    numerical = numerical[selected]
    reservoir = reservoir[selected]
    error = numerical - expected
    expected_l1 = norm_l1(expected, volume)
    error_l1 = norm_l1(error, volume)
    reservoir_integral = float(np.sum(reservoir * volume))
    tau_expected = (
        np.inf if expected_l1 == 0.0 else reservoir_integral / expected_l1
    )
    tau_error = np.inf if error_l1 == 0.0 else reservoir_integral / error_l1
    fractional_error = error / reservoir
    return {
        "relative_l1": error_l1 / max(expected_l1, np.finfo(float).tiny),
        "signed_error": float(np.sum(error * volume)),
        "expected_l1": expected_l1,
        "error_l1": error_l1,
        "tau_expected": tau_expected,
        "tau_error": tau_error,
        "tau_ratio": (
            np.inf if error_l1 == 0.0 else expected_l1 / error_l1
        ),
        "max_step_fraction": candidate.dt_s
        * float(np.max(np.abs(fractional_error))),
        "grid_scale_error_l1": norm_l1(high_pass(error), volume),
    }


def report(candidate: Case, reference: Case, rmax_km: float | None) -> None:
    selected = np.ones(len(candidate.radius_cm), dtype=bool)
    if rmax_km is not None:
        selected &= candidate.radius_cm <= rmax_km * 1.0e5
    if not np.any(selected):
        raise ValueError("radial selection is empty")

    print(f"candidate: {candidate.budget_path}")
    print(f"reference: {reference.budget_path}")
    print(
        "Reference rates are conservatively restricted as piecewise-constant "
        "spherical-shell averages."
    )
    for quantity, reservoir in (("QH", candidate.eint), ("QL", candidate.ne)):
        expected_origins = mapped_reference(candidate, reference, quantity)
        print(f"\n{quantity} comparison:")
        for name, numerical in candidate.origins[quantity].items():
            metrics = exchange_metrics(
                candidate,
                expected_origins[name],
                numerical,
                reservoir,
                selected,
            )
            print(
                f"  {name:25s} rel-L1={metrics['relative_l1']:.6e}  "
                f"signed-error={metrics['signed_error']:.6e}  "
                f"tau_ref={metrics['tau_expected']:.6e} s  "
                f"tau_err={metrics['tau_error']:.6e} s  "
                f"tau_err/tau_ref={metrics['tau_ratio']:.6e}  "
                f"max|dt*error/reservoir|={metrics['max_step_fraction']:.6e}"
            )


def plot_comparison(
    candidate: Case,
    reference: Case,
    output: Path,
    rmax_km: float | None,
) -> None:
    radius_km = candidate.radius_cm / 1.0e5
    selected = np.ones(len(radius_km), dtype=bool)
    if rmax_km is not None:
        selected &= radius_km <= rmax_km

    fig, axes = plt.subplots(2, 4, figsize=(23, 10), sharex=True)
    for row, (quantity, reservoir, source_unit) in enumerate(
        (
            ("QH", candidate.eint, r"erg cm$^{-3}$ s$^{-1}$"),
            ("QL", candidate.ne, r"cm$^{-3}$ s$^{-1}$"),
        )
    ):
        expected_origins = mapped_reference(candidate, reference, quantity)
        numerical = candidate.origins[quantity]["applied matter change"]
        expected = expected_origins["applied matter change"]
        error = numerical - expected
        tau_expected = np.full_like(expected, np.nan)
        tau_error = np.full_like(error, np.nan)
        tau_expected[expected != 0.0] = (
            reservoir[expected != 0.0] / np.abs(expected[expected != 0.0])
        )
        tau_error[error != 0.0] = reservoir[error != 0.0] / np.abs(error[error != 0.0])

        axes[row, 0].plot(radius_km, numerical, label="candidate")
        axes[row, 0].plot(radius_km, expected, label="restricted reference")
        axes[row, 1].plot(radius_km, error, label="candidate - reference")
        axes[row, 2].plot(radius_km, tau_expected, label="physical/reference")
        axes[row, 2].plot(radius_km, tau_error, label="numerical error")
        axes[row, 2].set_yscale("log")

        names = tuple(candidate.origins[quantity])
        relative_errors = []
        for name in names:
            metrics = exchange_metrics(
                candidate,
                expected_origins[name],
                candidate.origins[quantity][name],
                reservoir,
                selected,
            )
            relative_errors.append(metrics["relative_l1"])
        positions = np.arange(len(names))
        axes[row, 3].bar(positions, np.maximum(relative_errors, 1.0e-16))
        axes[row, 3].set_xticks(
            positions,
            [name.replace(" ", "\n") for name in names],
            rotation=35,
            ha="right",
            fontsize=7,
        )
        axes[row, 3].set_yscale("log")
        axes[row, 3].set_ylabel(r"$L^1$ error / reference $L^1$")
        axes[row, 0].set_ylabel(f"{quantity} [{source_unit}]")
        axes[row, 1].set_ylabel(f"{quantity} error [{source_unit}]")
        axes[row, 2].set_ylabel("timescale [s]")

    titles = (
        "Applied matter exchange",
        "Error relative to reference",
        "Physical and numerical-error timescales",
        "Error by numerical origin",
    )
    for row in range(2):
        for column, title in enumerate(titles):
            axes[row, column].set_title(f"{'Heating' if row == 0 else 'Lepton'}: {title}")
            axes[row, column].grid(alpha=0.2)
            if column < 3:
                axes[row, column].set_xlim(radius_km[selected][0], radius_km[selected][-1])
                axes[row, column].set_xlabel("radius [km]")
                axes[row, column].legend(fontsize=8, loc="best")
        axes[row, 0].axhline(0.0, color="0.5", lw=0.8)
        axes[row, 1].axhline(0.0, color="0.5", lw=0.8)
    fig.suptitle(
        "Matter exchange relative to an independently converged reference",
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("--candidate-state", type=Path)
    parser.add_argument("--reference-state", type=Path)
    parser.add_argument("--plot", type=Path)
    parser.add_argument("--rmax-km", type=float)
    args = parser.parse_args()

    candidate = load_case(args.candidate, args.candidate_state)
    reference = load_case(args.reference, args.reference_state)
    report(candidate, reference, args.rmax_km)
    if args.plot is not None:
        plot_comparison(candidate, reference, args.plot, args.rmax_km)


if __name__ == "__main__":
    main()
