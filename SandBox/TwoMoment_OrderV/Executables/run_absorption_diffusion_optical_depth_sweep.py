#!/usr/bin/env python3

"""Run the absorption--diffusion balance test at fixed spatial resolution.

Only the mean transport optical depth is varied.  The executable must have
``ProgramName = 'AbsorptionDiffusionBalance'``.  Each run is archived before
the next run overwrites Thornado's two fixed profile filenames.  By default,
the standard and experimentally diffusion-scaled LLF fluxes are both run.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess


PROFILE_NAMES = (
    "absorption_diffusion_operator_balance.dat",
    "absorption_diffusion_balance.dat",
)


def depth_label(depth: float) -> str:
    return f"tau_{depth:.6e}".replace("+", "").replace("-", "m").replace(".", "p")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument(
        "depths",
        nargs="*",
        type=float,
        default=(1.0e2, 1.0e3, 1.0e4, 1.0e5),
        help="mean domain optical depths (default: 1e2 1e3 1e4 1e5)",
    )
    parser.add_argument("--workdir", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("absorption_diffusion_optical_depth_sweep"),
    )
    parser.add_argument(
        "--llf-mode",
        choices=("standard", "scaled", "source-aware", "both", "all"),
        default="both",
        help=(
            "LLF dissipation mode to run; 'both' retains the historical "
            "standard/scaled pair, while 'all' also runs source-aware"
        ),
    )
    parser.add_argument(
        "--nodes",
        type=int,
        help="override the number of DG nodes with THORNADO_NNODES",
    )
    parser.add_argument(
        "--elements",
        type=int,
        help="override the number of X1 elements with THORNADO_NX1",
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
    if args.nodes is not None and args.nodes < 1:
        raise ValueError("--nodes must be positive")
    if args.elements is not None and args.elements < 1:
        raise ValueError("--elements must be positive")

    if args.nodes is not None and args.output_dir == Path(
        "absorption_diffusion_optical_depth_sweep"
    ):
        output_dir = workdir / f"absorption_diffusion_optical_depth_sweep_nNodes{args.nodes}"
        output_dir.mkdir(parents=True, exist_ok=True)

    if args.llf_mode == "both":
        modes = ("standard", "scaled")
    elif args.llf_mode == "all":
        modes = ("standard", "scaled", "source-aware")
    else:
        modes = (args.llf_mode,)

    for mode in modes:
        for depth in args.depths:
            label = f"{mode}_{depth_label(depth)}"
            environment = os.environ.copy()
            environment["THORNADO_PROGRAM_NAME"] = "AbsorptionDiffusionBalance"
            environment["THORNADO_ABSORPTION_DIFFUSION_OPTICAL_DEPTH"] = repr(depth)
            environment["THORNADO_SOURCE_AWARE_LLF"] = str(
                mode == "source-aware"
            ).lower()
            environment["THORNADO_ABSORPTION_DIFFUSION_LLF_SCALING"] = (
                "standard" if mode == "source-aware" else mode
            )
            if args.nodes is not None:
                environment["THORNADO_NNODES"] = str(args.nodes)
            if args.elements is not None:
                environment["THORNADO_NX1"] = str(args.elements)
            log_path = output_dir / f"{label}.log"

            print(f"running {mode} LLF at mean optical depth {depth:.6e}")
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
                destination = output_dir / f"{source.stem}_{label}{source.suffix}"
                shutil.copy2(source, destination)


if __name__ == "__main__":
    main()
