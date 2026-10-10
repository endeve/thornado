#!/usr/bin/env python3

"""Run the DeleptonizationWave collision-solver/flux comparison matrix.

The runs differ only in the source-aware AP LLF, terminal conservative
projection, and radiation-only nonlinear-solver switches.  A single
executable is reused; all run selection is performed through environment
variables read by the standalone driver.  Grouped no-projection selections
record only each configured trajectory by default and therefore skip the
counterfactual no-stream and opposite-flux shadow integrations.
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import subprocess
import sys


RUNS = (
    ("standard_no_projection", False, False, False),
    ("standard_projection", False, True, False),
    ("source_aware_no_projection", True, False, False),
    ("source_aware_projection", True, True, False),
    ("standard_radiation_only", False, False, True),
    ("source_aware_radiation_only", True, False, True),
)


def boolean_string(value: bool) -> str:
    return "true" if value else "false"


def run_and_tee(
    executable: Path,
    workdir: Path,
    environment: dict[str, str],
    log_path: Path,
) -> None:
    """Run Thornado while copying its combined output to the terminal and log."""

    with log_path.open("w", encoding="utf-8", buffering=1) as log:
        process = subprocess.Popen(
            [str(executable)],
            cwd=workdir,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            log.write(line)
            sys.stdout.write(line)
            sys.stdout.flush()
        return_code = process.wait()

    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, [str(executable)])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("--workdir", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("deleptonization_source_aware_matrix_nX256_nNodes2"),
    )
    parser.add_argument("--elements", type=int, default=256)
    parser.add_argument("--nodes", type=int, default=2)
    parser.add_argument(
        "--t-end-ms",
        type=float,
        help="override the driver end time in milliseconds",
    )
    parser.add_argument(
        "--max-cycles",
        type=int,
        help="override the driver cycle limit; use 1 to capture the t=0 operator",
    )
    parser.add_argument(
        "--case",
        choices=(
            "all",
            "both_no_projection",
            "both_radiation_only",
            "solver_comparison",
            "standard_solver_comparison",
            "source_aware_solver_comparison",
            *(run[0] for run in RUNS),
        ),
        default="all",
        help=(
            "run one named experiment, a grouped solver/flux comparison, "
            "or the complete matrix"
        ),
    )
    parser.add_argument(
        "--attribution-mode",
        choices=("auto", "full", "actual-only"),
        default="auto",
        help=(
            "full computes the no-stream and opposite-flux shadow steps; "
            "actual-only records only the configured production trajectory; "
            "auto uses actual-only for grouped comparisons and full otherwise"
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="replace an existing log or diagnostic from this matrix",
    )
    args = parser.parse_args()

    executable = args.executable.expanduser().resolve()
    workdir = args.workdir.expanduser().resolve()
    output_dir = args.output_dir.expanduser()
    if not output_dir.is_absolute():
        output_dir = workdir / output_dir

    if not executable.is_file():
        raise FileNotFoundError(executable)
    if not workdir.is_dir():
        raise NotADirectoryError(workdir)
    if args.elements < 1:
        raise ValueError("--elements must be positive")
    if args.nodes < 1:
        raise ValueError("--nodes must be positive")
    if args.t_end_ms is not None and args.t_end_ms <= 0.0:
        raise ValueError("--t-end-ms must be positive")
    if args.max_cycles is not None and args.max_cycles < 1:
        raise ValueError("--max-cycles must be positive")

    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str | int | float | bool]] = []
    if args.case == "all":
        selected_runs = RUNS
    elif args.case == "both_no_projection":
        selected_runs = tuple(
            run for run in RUNS if not run[2] and not run[3]
        )
    elif args.case == "both_radiation_only":
        selected_runs = tuple(run for run in RUNS if run[3])
    elif args.case == "solver_comparison":
        selected_runs = tuple(run for run in RUNS if not run[2])
    elif args.case == "standard_solver_comparison":
        selected_runs = tuple(
            run for run in RUNS if not run[1] and not run[2]
        )
    elif args.case == "source_aware_solver_comparison":
        selected_runs = tuple(
            run for run in RUNS if run[1] and not run[2]
        )
    else:
        selected_runs = tuple(run for run in RUNS if run[0] == args.case)

    if args.attribution_mode == "auto":
        compute_shadows = args.case not in (
            "both_no_projection",
            "both_radiation_only",
            "solver_comparison",
            "standard_solver_comparison",
            "source_aware_solver_comparison",
        )
    else:
        compute_shadows = args.attribution_mode == "full"

    for (
        label,
        source_aware,
        conservative_projection,
        radiation_only,
    ) in selected_runs:
        log_path = output_dir / f"{label}.log"
        diagnostic_path = output_dir / f"matter_diagnostics_{label}.dat"
        attribution_prefix = Path(str(diagnostic_path) + ".exchange_attribution")
        attribution_nodes = Path(str(attribution_prefix) + ".nodes")
        attribution_elements = Path(str(attribution_prefix) + ".elements")
        cumulative_conservation = Path(
            str(diagnostic_path) + ".cumulative_conservation"
        )
        imex_conservation_ledger = Path(
            str(cumulative_conservation) + ".imex"
        )
        existing = [log_path, *output_dir.glob(f"{diagnostic_path.name}*")]
        existing = [path for path in existing if path.exists()]
        if existing and not args.overwrite:
            names = ", ".join(path.name for path in existing)
            raise FileExistsError(
                f"outputs already exist for {label}: {names}; use --overwrite"
            )

        environment = os.environ.copy()
        # Prevent GNU Fortran from block-buffering preconnected output after
        # stdout is redirected through the tee pipe below.  Other runtimes
        # ignore this variable.
        environment.setdefault("GFORTRAN_UNBUFFERED_PRECONNECTED", "y")
        # The legacy variable takes precedence in the driver.  Remove an
        # inherited value so this matrix is controlled solely by the setting
        # below.
        environment.pop("THORNADO_WELL_BALANCED_TRANSPORT", None)
        environment["THORNADO_PROGRAM_NAME"] = "DeleptonizationWave1D"
        environment["THORNADO_NX1"] = str(args.elements)
        environment["THORNADO_NNODES"] = str(args.nodes)
        environment["THORNADO_SOURCE_AWARE_LLF"] = boolean_string(source_aware)
        environment["THORNADO_CONSERVATIVE_PROJECTION"] = boolean_string(
            conservative_projection
        )
        environment["THORNADO_RADIATION_ONLY_SOLVE"] = boolean_string(
            radiation_only
        )
        environment["THORNADO_AP_WELL_BALANCED_FLUX_SOURCE"] = "false"
        environment["THORNADO_MATTER_EXCHANGE_SHADOWS"] = boolean_string(
            compute_shadows
        )
        environment["THORNADO_MATTER_DIAGNOSTICS_FILE"] = str(diagnostic_path)
        if args.t_end_ms is not None:
            environment["THORNADO_T_END_MS"] = str(args.t_end_ms)
        if args.max_cycles is not None:
            environment["THORNADO_MAX_CYCLES"] = str(args.max_cycles)

        print(
            f"running {label}: source-aware={source_aware}, "
            f"conservative-projection={conservative_projection}, "
            f"radiation-only={radiation_only}, "
            f"attribution={'full' if compute_shadows else 'actual-only'}"
        )
        run_and_tee(executable, workdir, environment, log_path)

        if radiation_only and "[Invariant-constrained solve] Failed after" in \
                log_path.read_text(encoding="utf-8"):
            raise RuntimeError(
                f"{label} used the nested-solver fallback; rebuild with the "
                "no-fallback radiation-only solver"
            )

        if not diagnostic_path.is_file():
            raise FileNotFoundError(
                f"{diagnostic_path} was not produced; see {log_path}"
            )
        missing_attribution = [
            path
            for path in (attribution_nodes, attribution_elements)
            if not path.is_file()
        ]
        if missing_attribution:
            names = ", ".join(str(path) for path in missing_attribution)
            raise FileNotFoundError(
                f"matter-exchange attribution output was not produced: {names}; "
                f"verify that the executable contains the diagnostic changes; "
                f"see {log_path}"
            )
        if not cumulative_conservation.is_file():
            raise FileNotFoundError(
                "per-step cumulative conservation output was not produced: "
                f"{cumulative_conservation}; see {log_path}"
            )
        if not imex_conservation_ledger.is_file():
            raise FileNotFoundError(
                "per-step IMEX conservation ledger was not produced: "
                f"{imex_conservation_ledger}; see {log_path}"
            )

        rows.append(
            {
                "label": label,
                "source_aware_llf": source_aware,
                "conservative_projection": conservative_projection,
                "radiation_only_solve": radiation_only,
                "matter_exchange_shadows": compute_shadows,
                "spatial_elements": args.elements,
                "dg_nodes": args.nodes,
                "t_end_ms": "driver default"
                if args.t_end_ms is None
                else args.t_end_ms,
                "log": log_path.name,
                "diagnostic": diagnostic_path.name,
                "attribution_nodes": attribution_nodes.name,
                "attribution_elements": attribution_elements.name,
                "cumulative_conservation": cumulative_conservation.name,
                "imex_conservation_ledger": imex_conservation_ledger.name,
            }
        )

    manifest_path = output_dir / "matrix.csv"
    with manifest_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"matrix complete: {output_dir}")


if __name__ == "__main__":
    main()
