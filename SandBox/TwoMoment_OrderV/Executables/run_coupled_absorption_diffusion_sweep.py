#!/usr/bin/env python3

"""Run the finite-capacity coupled absorption--diffusion test.

The executable must be built with
``ProgramName = 'CoupledAbsorptionDiffusion'``.  The test evolves an exact
slow Fourier eigenmode of a linear radiation--matter system.  Species one is
the transported radiation field; the zeroth moment of species two is a local,
non-streaming matter reservoir.  The implicit source conserves their sum
exactly, while the standard Thornado DG streaming operator transports species
one.

The sweep compares the first-order coupled IMEX scheme with PD-ARS while
varying spatial resolution and optical depth.  Plain ``BackwardEuler`` is
deliberately excluded: in Thornado that tableau has no explicit weights and
therefore disables spatial streaming.

The ``ap-well-balanced`` mode uses the same stage-local equilibrium state in
the implicit exchange and in the numerical trace.  The trace applies the
full LLF dissipation to the nonequilibrium departure and the source-aware
factor to the equilibrium contribution.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import re
import shutil
import subprocess


PROFILE_NAME = "coupled_absorption_diffusion.dat"

TRIPLET_PATTERNS = {
    "total": re.compile(
        r"total rel\. L1/L2/Linf\s*=\s*"
        r"([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)"
    ),
    "partition": re.compile(
        r"partition rel\. L1/L2/Linf\s*=\s*"
        r"([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)"
    ),
    "flux": re.compile(
        r"flux rel\. L1/L2/Linf\s*=\s*"
        r"([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)"
    ),
    "exchange": re.compile(
        r"exchange rel\. L1/L2/Linf\s*=\s*"
        r"([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)"
    ),
}

SCALAR_PATTERNS = {
    "efoldings": re.compile(
        r"elapsed slow-mode e-foldings\s*=\s*([+\-0-9.Ee]+)"
    ),
    "mean_total_error": re.compile(
        r"mean total conservation err\.\s*=\s*([+\-0-9.Ee]+)"
    ),
    "total_checkerboard": re.compile(
        r"total checkerboard amplitude\s*=\s*([+\-0-9.Ee]+)"
    ),
    "exchange_checkerboard": re.compile(
        r"exch\. checkerboard amplitude\s*=\s*([+\-0-9.Ee]+)"
    ),
}


def depth_label(depth: float) -> str:
    return f"tau_{depth:.6e}".replace("+", "").replace("-", "m").replace(".", "p")


def observed_order(coarse: float, fine: float, mesh_ratio: float) -> float:
    if coarse <= 0.0 or fine <= 0.0:
        return math.nan
    return math.log(coarse / fine) / math.log(mesh_ratio)


def extract_metrics(log_text: str) -> dict[str, float]:
    metrics: dict[str, float] = {}

    for name, pattern in TRIPLET_PATTERNS.items():
        match = pattern.search(log_text)
        if match is None:
            raise ValueError(f"missing {name} error summary")
        metrics[f"{name}_l1"] = float(match.group(1))
        metrics[f"{name}_l2"] = float(match.group(2))
        metrics[f"{name}_linf"] = float(match.group(3))

    for name, pattern in SCALAR_PATTERNS.items():
        match = pattern.search(log_text)
        if match is None:
            raise ValueError(f"missing {name} summary")
        metrics[name] = float(match.group(1))

    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument(
        "depths",
        nargs="*",
        type=float,
        default=(1.0e2, 3.0e2, 1.0e3),
        help="mean domain optical depths (default: 1e2 3e2 1e3)",
    )
    parser.add_argument(
        "--resolutions",
        nargs="+",
        type=int,
        default=(16, 32),
        help="spatial element counts (default: 16 32)",
    )
    parser.add_argument(
        "--schemes",
        nargs="+",
        choices=("IMEX_ARS_111", "IMEX_PDARS"),
        default=("IMEX_ARS_111", "IMEX_PDARS"),
    )
    parser.add_argument(
        "--llf-modes",
        nargs="+",
        choices=("standard", "source-aware", "ap-well-balanced"),
        default=("standard",),
        help="flux-source treatments to run (default: standard)",
    )
    parser.add_argument("--workdir", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("coupled_absorption_diffusion_sweep"),
    )
    args = parser.parse_args()

    executable = args.executable.expanduser().resolve()
    workdir = args.workdir.expanduser().resolve()
    output_dir = args.output_dir.expanduser()
    if not output_dir.is_absolute():
        output_dir = workdir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    if not executable.is_file():
        raise FileNotFoundError(executable)
    if any(depth <= 0.0 for depth in args.depths):
        raise ValueError("all optical depths must be positive")
    if any(resolution < 2 for resolution in args.resolutions):
        raise ValueError("all spatial resolutions must contain at least two elements")

    rows: list[dict[str, float | int | str]] = []

    for llf_mode in args.llf_modes:
        for scheme in args.schemes:
            for resolution in args.resolutions:
                for depth in args.depths:
                    label = (
                        f"{llf_mode}_{scheme}_nx_{resolution}_"
                        f"{depth_label(depth)}"
                    )
                    environment = os.environ.copy()
                    environment["THORNADO_COUPLED_DIFFUSION_SCHEME"] = scheme
                    environment["THORNADO_COUPLED_DIFFUSION_NX"] = str(resolution)
                    environment["THORNADO_COUPLED_DIFFUSION_OPTICAL_DEPTH"] = repr(depth)
                    environment["THORNADO_SOURCE_AWARE_LLF"] = str(
                        llf_mode in ("source-aware", "ap-well-balanced")
                    ).lower()
                    environment["THORNADO_AP_WELL_BALANCED_FLUX_SOURCE"] = str(
                        llf_mode == "ap-well-balanced"
                    ).lower()
                    log_path = output_dir / f"{label}.log"

                    print(
                        f"running {llf_mode} LLF, {scheme}, "
                        f"{resolution} elements, optical depth {depth:.6e}"
                    )
                    with log_path.open("w", encoding="utf-8") as log:
                        subprocess.run(
                            [str(executable)],
                            cwd=workdir,
                            env=environment,
                            stdout=log,
                            stderr=subprocess.STDOUT,
                            check=True,
                        )

                    source = workdir / PROFILE_NAME
                    if not source.is_file():
                        raise FileNotFoundError(
                            f"{source} was not produced; see {log_path}"
                        )
                    shutil.copy2(
                        source,
                        output_dir / f"{source.stem}_{label}{source.suffix}",
                    )

                    metrics = extract_metrics(
                        log_path.read_text(encoding="utf-8")
                    )
                    rows.append(
                        {
                            "llf_mode": llf_mode,
                            "scheme": scheme,
                            "spatial_elements": resolution,
                            "optical_depth": depth,
                            **metrics,
                        }
                    )

    fieldnames = [
        "llf_mode",
        "scheme",
        "spatial_elements",
        "optical_depth",
        "total_l1",
        "total_l2",
        "total_linf",
        "partition_l1",
        "partition_l2",
        "partition_linf",
        "flux_l1",
        "flux_l2",
        "flux_linf",
        "exchange_l1",
        "exchange_l2",
        "exchange_linf",
        "efoldings",
        "mean_total_error",
        "total_checkerboard",
        "exchange_checkerboard",
    ]
    summary_path = output_dir / "summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print()
    print(
        " LLF mode      scheme          nX  optical depth       total L2   partition L2"
        "          flux L2      exchange L2  exchange checker"
    )
    for row in rows:
        print(
            f" {str(row['llf_mode']):12s}"
            f" {str(row['scheme']):14s}"
            f" {int(row['spatial_elements']):4d}"
            f" {float(row['optical_depth']):13.6e}"
            f" {float(row['total_l2']):14.6e}"
            f" {float(row['partition_l2']):14.6e}"
            f" {float(row['flux_l2']):14.6e}"
            f" {float(row['exchange_l2']):14.6e}"
            f" {float(row['exchange_checkerboard']):17.6e}"
        )

    print("\nspatial orders at each optical depth:")
    for llf_mode in args.llf_modes:
        for scheme in args.schemes:
            for depth in args.depths:
                selected = sorted(
                    (
                        row
                        for row in rows
                        if row["llf_mode"] == llf_mode
                        and row["scheme"] == scheme
                        and float(row["optical_depth"]) == depth
                    ),
                    key=lambda row: int(row["spatial_elements"]),
                )
                for coarse, fine in zip(selected[:-1], selected[1:]):
                    mesh_ratio = int(fine["spatial_elements"]) / int(
                        coarse["spatial_elements"]
                    )
                    print(
                        f"  {llf_mode:12s} {scheme:14s} tau={depth:.6e} "
                        f"nX={int(coarse['spatial_elements'])}->{int(fine['spatial_elements'])}: "
                        f"total={observed_order(float(coarse['total_l2']), float(fine['total_l2']), mesh_ratio):.4f}, "
                        f"partition={observed_order(float(coarse['partition_l2']), float(fine['partition_l2']), mesh_ratio):.4f}, "
                        f"exchange={observed_order(float(coarse['exchange_l2']), float(fine['exchange_l2']), mesh_ratio):.4f}"
                    )

    print(f"\nsummary: {summary_path}")


if __name__ == "__main__":
    main()
