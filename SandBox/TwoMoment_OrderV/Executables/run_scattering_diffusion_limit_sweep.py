#!/usr/bin/env python3

"""Exercise Thornado's variable-opacity scattering diffusion limit.

The executable must be built with ``ProgramName = 'ScatteringDiffusionLimit'``.
At each optical depth the test evolves the same continuum steady solution for
six nominal diffusion e-folding times.  This gives the upwind DG solution time
to enter its continuous asymptotic subspace before it is compared with the
variable-coefficient diffusion solution.

The relevant asymptotic signatures are:

* the density error approaches an optical-depth-independent spatial plateau;
* the density trace jump tends to zero as optical depth increases;
* increasing optical depth does not produce the linear error growth seen when
  the operator is merely evaluated on the initial discontinuous interpolant.
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


PROFILE_NAMES = (
    "scattering_diffusion_limit.dat",
    "scattering_diffusion_limit_faces.dat",
)

TRIPLET_PATTERNS = {
    "density": re.compile(
        r"density rel\. L1/L2/Linf\s*=\s*"
        r"([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)"
    ),
    "flux": re.compile(
        r"flux rel\. L1/L2/Linf\s*=\s*"
        r"([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)"
    ),
}

PAIR_PATTERNS = {
    "jump": re.compile(
        r"density face-jump rel\. L2/Linf\s*=\s*"
        r"([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)"
    ),
    "scaled_jump": re.compile(
        r"optical-depth-scaled jumps\s*=\s*"
        r"([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)"
    ),
}

SCALAR_PATTERNS = {
    "checkerboard": re.compile(
        r"checkerboard error amplitude\s*=\s*([+\-0-9.Ee]+)"
    ),
    "efoldings": re.compile(
        r"nominal diffusion e-foldings\s*=\s*([+\-0-9.Ee]+)"
    ),
}


def depth_label(depth: float) -> str:
    return f"tau_{depth:.6e}".replace("+", "").replace("-", "m").replace(".", "p")


def safe_ratio(numerator: float, denominator: float) -> float:
    if denominator == 0.0:
        return math.nan
    return numerator / denominator


def observed_order(coarse: float, fine: float, mesh_ratio: float) -> float:
    ratio = safe_ratio(coarse, fine)
    if ratio <= 0.0 or not math.isfinite(ratio):
        return math.nan
    return math.log(ratio) / math.log(mesh_ratio)


def extract_metrics(log_text: str) -> dict[str, float]:
    metrics: dict[str, float] = {}

    for name, pattern in TRIPLET_PATTERNS.items():
        match = pattern.search(log_text)
        if match is None:
            raise ValueError(f"missing {name} error summary")
        metrics[f"{name}_l1"] = float(match.group(1))
        metrics[f"{name}_l2"] = float(match.group(2))
        metrics[f"{name}_linf"] = float(match.group(3))

    for name, pattern in PAIR_PATTERNS.items():
        match = pattern.search(log_text)
        if match is None:
            raise ValueError(f"missing {name} summary")
        metrics[f"{name}_l2"] = float(match.group(1))
        metrics[f"{name}_linf"] = float(match.group(2))

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
        "--llf-modes",
        nargs="+",
        choices=("standard", "source-aware"),
        default=("standard",),
        help="density-penalty treatments to run (default: standard)",
    )
    parser.add_argument("--workdir", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("scattering_diffusion_limit_sweep"),
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
        for resolution in args.resolutions:
            for depth in args.depths:
                label = f"{llf_mode}_nx_{resolution}_{depth_label(depth)}"
                environment = os.environ.copy()
                environment["THORNADO_SCATTERING_DIFFUSION_OPTICAL_DEPTH"] = repr(depth)
                environment["THORNADO_SCATTERING_DIFFUSION_NX"] = str(resolution)
                environment["THORNADO_SOURCE_AWARE_LLF"] = str(
                    llf_mode == "source-aware"
                ).lower()
                log_path = output_dir / f"{label}.log"

                print(
                    f"running {llf_mode} LLF, {resolution} elements "
                    f"at mean optical depth {depth:.6e}"
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

                for profile_name in PROFILE_NAMES:
                    source = workdir / profile_name
                    if not source.is_file():
                        raise FileNotFoundError(
                            f"{source} was not produced; see {log_path}"
                        )
                    destination = (
                        output_dir / f"{source.stem}_{label}{source.suffix}"
                    )
                    shutil.copy2(source, destination)

                metrics = extract_metrics(log_path.read_text(encoding="utf-8"))
                rows.append(
                    {
                        "llf_mode": llf_mode,
                        "spatial_elements": resolution,
                        "optical_depth": depth,
                        **metrics,
                    }
                )

    fieldnames = [
        "llf_mode",
        "spatial_elements",
        "optical_depth",
        "density_l1",
        "density_l2",
        "density_linf",
        "flux_l1",
        "flux_l2",
        "flux_linf",
        "jump_l2",
        "jump_linf",
        "scaled_jump_l2",
        "scaled_jump_linf",
        "checkerboard",
        "efoldings",
    ]
    summary_path = output_dir / "summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print()
    print(
        " LLF mode       nX  optical depth       density L2"
        "          flux L2     face-jump L2"
    )
    for row in rows:
        print(
            f" {str(row['llf_mode']):12s}"
            f" {int(row['spatial_elements']):4d}"
            f" {row['optical_depth']:13.6e}"
            f" {row['density_l2']:16.8e}"
            f" {row['flux_l2']:16.8e}"
            f" {row['jump_l2']:16.8e}"
        )

    print("\nfixed-grid thick-limit ratios (largest tau / preceding tau):")
    for llf_mode in args.llf_modes:
        for resolution in args.resolutions:
            resolution_rows = sorted(
                (
                    row
                    for row in rows
                    if row["llf_mode"] == llf_mode
                    and int(row["spatial_elements"]) == resolution
                ),
                key=lambda row: row["optical_depth"],
            )
            if len(resolution_rows) < 2:
                continue
            previous, largest = resolution_rows[-2:]
            print(
                f"  {llf_mode:12s} nX={resolution:4d}:"
                f" density={safe_ratio(largest['density_l2'], previous['density_l2']):.6e},"
                f" flux={safe_ratio(largest['flux_l2'], previous['flux_l2']):.6e},"
                f" jump={safe_ratio(largest['jump_l2'], previous['jump_l2']):.6e}"
            )

    largest_depth = max(args.depths)
    for llf_mode in args.llf_modes:
        thick_rows = sorted(
            (
                row
                for row in rows
                if row["llf_mode"] == llf_mode
                and row["optical_depth"] == largest_depth
            ),
            key=lambda row: row["spatial_elements"],
        )
        if len(thick_rows) >= 2:
            print(
                f"\n{llf_mode} spatial orders at tau={largest_depth:.6e}:"
            )
            for coarse, fine in zip(thick_rows[:-1], thick_rows[1:]):
                mesh_ratio = (
                    fine["spatial_elements"] / coarse["spatial_elements"]
                )
                density_order = observed_order(
                    coarse["density_l2"], fine["density_l2"], mesh_ratio
                )
                flux_order = observed_order(
                    coarse["flux_l2"], fine["flux_l2"], mesh_ratio
                )
                print(
                    f"  {int(coarse['spatial_elements'])} ->"
                    f" {int(fine['spatial_elements'])}:"
                    f" density={density_order:.6f}, flux={flux_order:.6f}"
                )
    print(f"\nsummary: {summary_path}")


if __name__ == "__main__":
    main()
