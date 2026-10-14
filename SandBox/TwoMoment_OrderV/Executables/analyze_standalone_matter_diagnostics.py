#!/usr/bin/env python3
"""Analyse gray collision and inverse-mean-free-path sandbox diagnostics.

Each input is the ``matter_diagnostics.dat`` file written by
``StandaloneMatterDiagnosticsModule``.  With more than one input, --compare
reports successive-resolution differences; it deliberately does not designate
any input as a reference solution.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable

import numpy as np


PROCESS_NAMES = (
    "EmAb", "Iso", "NNS", "NES", "Pair", "NuPair", "Brem",
    "NMS", "LFE", "LFC", "IMD",
)
HEATING_PROCESS_NAMES = (
    "EmAb", "NNS", "NES", "Pair", "NuPair", "Brem",
    "NMS", "LFE", "LFC", "IMD",
)
NONISO_PROCESS_NAMES = HEATING_PROCESS_NAMES
GENERAL_LEPTON_PROCESS_NAMES = ("NMS", "LFE", "LFC", "IMD")
TAU_SURFACE = 2.0 / 3.0
SPECIES_PLOT_INFO = (
    (1, r"$\nu_e$", "C0"),
    (2, r"$\bar\nu_e$", "C1"),
    (3, r"$\nu_\mu$", "C2"),
    (4, r"$\bar\nu_\mu$", "C3"),
    (5, r"$\nu_\tau$", "C4"),
    (6, r"$\bar\nu_\tau$", "C5"),
)


@dataclass
class DiagnosticFile:
    path: Path
    time_s: float
    applied_dt_s: float
    n_energy_elements: int
    n_energy_nodes: int
    n_species: int
    columns: dict[str, np.ndarray]


def read_diagnostic(path: Path) -> DiagnosticFile:
    """Read the self-describing ASCII diagnostic written by the sandbox."""

    names: list[str] | None = None
    time_s = float("nan")
    applied_dt_s = float("nan")
    n_energy_elements = -1
    n_energy_nodes = -1

    with path.open() as stream:
        for line in stream:
            if line.startswith("# time_s"):
                time_s = float(line.split()[2])
            elif line.startswith("# applied_dt_s"):
                applied_dt_s = float(line.split()[2])
            elif line.startswith("# nE"):
                match = re.fullmatch(r"# nE\s+(\d+)\s*nDOFE\s+(\d+)\s*", line)
                if match is None:
                    raise ValueError(f"{path}: malformed '# nE ... nDOFE ...' header")
                n_energy_elements = int(match.group(1))
                n_energy_nodes = int(match.group(2))
            elif line.startswith("# r_cm"):
                names = line[1:].split()
                break

    if names is None:
        raise ValueError(f"{path}: missing '# r_cm ...' column header")

    species = sorted(
        int(match.group(1))
        for name in names
        if (match := re.fullmatch(r"T_spec_s(\d+)", name)) is not None
    )
    if not species or species != list(range(1, species[-1] + 1)):
        raise ValueError(
            f"{path}: could not infer a contiguous species range from T_spec columns"
        )

    values = np.loadtxt(path, comments="#")
    values = np.atleast_2d(values)
    if values.shape[1] != len(names):
        raise ValueError(
            f"{path}: found {values.shape[1]} numeric columns, expected {len(names)}"
        )

    columns = dict(zip(names, values.T, strict=True))

    # The decompositions are assembled from the same per-species/per-process
    # source array as QH_total.  A failure here therefore indicates a corrupt
    # or mislabeled diagnostic row, rather than a physics closure error.
    total_key = "QH_total_erg_cm3_s"
    species_keys = [
        f"QH_species_s{index}" for index in range(1, species[-1] + 1)
    ]
    process_keys = [f"QH_process_{process}" for process in PROCESS_NAMES]
    if (
        total_key in columns
        and all(key in columns for key in species_keys)
        and all(key in columns for key in process_keys)
    ):
        total = columns[total_key]
        by_species = np.sum([columns[key] for key in species_keys], axis=0)
        by_process = np.sum([columns[key] for key in process_keys], axis=0)
        scale = max(
            float(np.max(np.abs(total))),
            float(np.max(np.sum([np.abs(columns[key]) for key in species_keys], axis=0))),
            float(np.max(np.sum([np.abs(columns[key]) for key in process_keys], axis=0))),
            np.finfo(float).tiny,
        )
        mismatch = max(
            float(np.max(np.abs(by_species - total))),
            float(np.max(np.abs(by_process - total))),
        ) / scale
        if not np.isfinite(mismatch) or mismatch > 1.0e-10:
            raise ValueError(
                f"{path}: inconsistent QH decomposition (relative mismatch "
                f"{mismatch:.3e}); the header and numeric column order likely differ"
            )

    return DiagnosticFile(
        path=path,
        time_s=time_s,
        applied_dt_s=applied_dt_s,
        n_energy_elements=n_energy_elements,
        n_energy_nodes=n_energy_nodes,
        n_species=species[-1],
        columns=columns,
    )


def read_named_table(path: Path) -> dict[str, np.ndarray]:
    """Read a whitespace table whose first line is a commented header."""

    with path.open() as stream:
        header = stream.readline()
    if not header.startswith("#"):
        raise ValueError(f"{path}: missing commented column header")
    names = header[1:].split()
    values = np.atleast_2d(np.loadtxt(path, comments="#"))
    if values.shape[1] != len(names):
        raise ValueError(
            f"{path}: found {values.shape[1]} numeric columns, expected {len(names)}"
        )
    return dict(zip(names, values.T, strict=True))


def plotted_species(
    data: DiagnosticFile, maximum_species: int | None = None
) -> tuple[tuple[int, str, str], ...]:
    """Return plotting metadata only for species present in the file."""

    upper = data.n_species
    if maximum_species is not None:
        upper = min(upper, maximum_species)
    return tuple(info for info in SPECIES_PLOT_INFO if info[0] <= upper)


def weighted_integral(values: np.ndarray, volume: np.ndarray) -> float:
    return float(np.sum(values * volume))


def weighted_l1(values: np.ndarray, volume: np.ndarray) -> float:
    return float(np.sum(np.abs(values) * volume))


def cumulative_volume_integral(values: np.ndarray, volume: np.ndarray) -> np.ndarray:
    """Return the radial cumulative integral of element-averaged values."""

    return np.cumsum(values * volume)


def process_inverse_mfp(data: DiagnosticFile, species: int, process: str) -> np.ndarray:
    return data.columns[f"inv_mfp_{process}_s{species}_cm_inv"]


def collision_ledger_name(
    process: str, quantity: str, species: int, unit: str
) -> str:
    return f"{process}_{quantity}_s{species}_{unit}"


def nns_ledger_name(quantity: str, species: int, unit: str) -> str:
    return collision_ledger_name("NNS", quantity, species, unit)


def emab_ledger_name(quantity: str, species: int, unit: str) -> str:
    return collision_ledger_name("EmAb", quantity, species, unit)


def radius_at_tau(
    radius: np.ndarray,
    width: np.ndarray,
    inverse_mfp: np.ndarray,
    target: float = TAU_SURFACE,
) -> float:
    """Locate tau=target using radial element averages of inverse MFP."""

    if np.any(np.diff(radius) <= 0.0) or np.any(width <= 0.0):
        raise ValueError("radial element centers and widths must be strictly increasing/positive")
    if np.any(inverse_mfp < 0.0):
        raise ValueError("optical depth is undefined for a negative combined inverse MFP")

    tau_outer = 0.0
    for i in range(radius.size - 1, -1, -1):
        tau_inner = tau_outer + inverse_mfp[i] * width[i]
        if tau_inner >= target and tau_outer < target and inverse_mfp[i] > 0.0:
            return float(
                radius[i] + 0.5 * width[i] - (target - tau_outer) / inverse_mfp[i]
            )
        tau_outer = tau_inner

    return float("nan")


def optical_depth_at_centers(width: np.ndarray, inverse_mfp: np.ndarray) -> np.ndarray:
    """Return optical depth at each radial element center."""

    if np.any(width <= 0.0) or np.any(inverse_mfp < 0.0):
        raise ValueError("optical depth requires positive widths and nonnegative inverse MFP")

    tau = np.zeros_like(inverse_mfp)
    tau_outer = 0.0
    for i in range(inverse_mfp.size - 1, -1, -1):
        tau[i] = tau_outer + 0.5 * inverse_mfp[i] * width[i]
        tau_outer += inverse_mfp[i] * width[i]
    return tau


def nearest(values: np.ndarray, radius: np.ndarray, target_radius: float) -> float:
    if not np.isfinite(target_radius):
        return float("nan")
    return float(values[np.argmin(np.abs(radius - target_radius))])


def report_file(data: DiagnosticFile) -> None:
    c = data.columns
    radius = c["r_cm"]
    width = c["dr_cm"]
    volume = c["dV_cm3"]

    print(f"\n{data.path}")
    print(
        f"  species: {data.n_species}; "
        f"energy elements/nodes: {data.n_energy_elements}/{data.n_energy_nodes}; "
        f"radial elements: {radius.size}"
    )
    print(
        f"  time={data.time_s:.8e} s; "
        f"applied-step dt={data.applied_dt_s:.8e} s"
    )
    roundoff_dt = (
        100.0 * np.finfo(float).eps
        * max(abs(data.time_s), np.finfo(float).tiny)
    )
    if data.applied_dt_s <= roundoff_dt:
        print(
            "  WARNING: applied-step dt is at the final-time roundoff scale; "
            "the applied heating is not meaningful"
        )

    principal_heating_fields = (
        "QH_total_erg_cm3_s",
        "QH_IMEX_source_total_erg_cm3_s",
        "QH_IMEX_increment_erg_cm3_s",
        "QH_applied_last_step_erg_cm3_s",
        "QH_conservative_projection_final_stage_erg_cm3_s",
        "QH_outer_equation_defect_erg_cm3_s",
    )
    print("  principal heating integrals [erg s^-1]:")
    for name in (key for key in principal_heating_fields if key in c):
        signed = weighted_integral(c[name], volume)
        magnitude = weighted_l1(c[name], volume)
        print(f"    {name:28s} signed={signed: .8e}  L1={magnitude: .8e}")

    transport_key = "QH_transport_equivalent_erg_cm3_s"
    balance_key = "QH_transport_minus_collision_erg_cm3_s"
    if transport_key in c and balance_key in c:
        collision = c["QH_total_erg_cm3_s"]
        transport = c[transport_key]
        balance = c[balance_key]
        scale = max(
            weighted_l1(collision, volume),
            weighted_l1(transport, volume),
            np.finfo(float).tiny,
        )
        print("  final-state transport/collision balance (element averages):")
        print(
            "    L1(transport-collision) / max[L1]="
            f"{weighted_l1(balance, volume) / scale:.8e}  "
            "signed(transport-collision) / max[L1]="
            f"{weighted_integral(balance, volume) / scale:.8e}"
        )

        ledger_path = Path(str(data.path) + ".transport_collision")
        if ledger_path.exists():
            ledger = read_named_table(ledger_path)
            total = ledger["species"] == 0
            weight = ledger["volume_weight_cm3"][total]
            nodal_collision = ledger["QH_collision_erg_cm3_s"][total]
            nodal_transport = ledger["QH_transport_erg_cm3_s"][total]
            nodal_balance = ledger["QH_difference_erg_cm3_s"][total]
            nodal_scale = max(
                weighted_l1(nodal_collision, weight),
                weighted_l1(nodal_transport, weight),
                np.finfo(float).tiny,
            )
            print("  final-state transport/collision balance (spatial nodes):")
            print(
                "    L1(transport-collision) / max[L1]="
                f"{weighted_l1(nodal_balance, weight) / nodal_scale:.8e}"
            )

            element = ledger["element"][total].astype(int)
            for label, values in (
                ("collision", nodal_collision),
                ("transport", nodal_transport),
                ("difference", nodal_balance),
            ):
                variation = 0.0
                magnitude = 0.0
                for index in np.unique(element):
                    selected = element == index
                    element_weight = weight[selected]
                    average = np.sum(element_weight * values[selected]) / np.sum(
                        element_weight
                    )
                    variation += np.sum(
                        element_weight * np.abs(values[selected] - average)
                    )
                    magnitude += np.sum(element_weight * np.abs(values[selected]))
                print(
                    f"    within-element variation / nodal L1 ({label:10s})="
                    f"{variation / max(magnitude, np.finfo(float).tiny):.8e}"
                )

            element_collision = []
            element_transport = []
            for index in np.unique(element):
                selected = element == index
                element_weight = weight[selected]
                element_collision.append(
                    np.sum(element_weight * nodal_collision[selected])
                    / np.sum(element_weight)
                )
                element_transport.append(
                    np.sum(element_weight * nodal_transport[selected])
                    / np.sum(element_weight)
                )
            element_collision = np.asarray(element_collision)
            element_transport = np.asarray(element_transport)
            projection_scale = max(
                float(np.max(np.abs(collision))),
                float(np.max(np.abs(transport))),
                np.finfo(float).tiny,
            )
            projection_mismatch = max(
                float(np.max(np.abs(element_collision - collision))),
                float(np.max(np.abs(element_transport - transport))),
            ) / projection_scale
            print(
                "    nodal-to-element volume-average closure="
                f"{projection_mismatch:.8e}"
            )

            spectrum_path = Path(str(ledger_path) + "_spectrum")
            if spectrum_path.exists():
                spectrum = read_named_table(spectrum_path)
                spectral_keys = np.column_stack(
                    (
                        spectrum["element"].astype(int),
                        spectrum["spatial_node"].astype(int),
                        spectrum["species"].astype(int),
                    )
                )
                unique_keys, inverse = np.unique(
                    spectral_keys, axis=0, return_inverse=True
                )
                spectral_collision = np.zeros(unique_keys.shape[0])
                spectral_transport = np.zeros(unique_keys.shape[0])
                np.add.at(
                    spectral_collision,
                    inverse,
                    spectrum["QH_collision_erg_cm3_s"],
                )
                np.add.at(
                    spectral_transport,
                    inverse,
                    spectrum["QH_transport_erg_cm3_s"],
                )

                species_rows = ledger["species"] > 0
                nodal_keys = np.column_stack(
                    (
                        ledger["element"][species_rows].astype(int),
                        ledger["spatial_node"][species_rows].astype(int),
                        ledger["species"][species_rows].astype(int),
                    )
                )
                nodal_order = np.lexsort(nodal_keys[:, ::-1].T)
                spectral_order = np.lexsort(unique_keys[:, ::-1].T)
                if not np.array_equal(
                    nodal_keys[nodal_order], unique_keys[spectral_order]
                ):
                    raise ValueError(
                        f"{spectrum_path}: node/species keys do not match {ledger_path}"
                    )
                nodal_species_collision = ledger[
                    "QH_collision_erg_cm3_s"
                ][species_rows][nodal_order]
                nodal_species_transport = ledger[
                    "QH_transport_erg_cm3_s"
                ][species_rows][nodal_order]
                spectral_collision = spectral_collision[spectral_order]
                spectral_transport = spectral_transport[spectral_order]
                spectral_scale = max(
                    float(np.max(np.abs(nodal_species_collision))),
                    float(np.max(np.abs(nodal_species_transport))),
                    np.finfo(float).tiny,
                )
                spectral_mismatch = max(
                    float(
                        np.max(
                            np.abs(spectral_collision - nodal_species_collision)
                        )
                    ),
                    float(
                        np.max(
                            np.abs(spectral_transport - nodal_species_transport)
                        )
                    ),
                ) / spectral_scale
                print(
                    "    energy-resolved-to-gray closure="
                    f"{spectral_mismatch:.8e}"
                )

    applied_key = "QH_applied_last_step_erg_cm3_s"
    if applied_key in c:
        rhs = c["QH_total_erg_cm3_s"]
        applied = c[applied_key]
        difference = rhs - applied
        applied_l1 = weighted_l1(applied, volume)
        applied_scale = max(applied_l1, np.finfo(float).tiny)
        cumulative_difference = np.max(np.abs(
            cumulative_volume_integral(difference, volume)
        ))
        print("  instantaneous total RHS vs applied last-step heating:")
        print(
            f"    L1 difference / L1(applied)="
            f"{weighted_l1(difference, volume) / applied_scale:.8e}  "
            f"max|cumulative difference| / L1(applied)="
            f"{cumulative_difference / applied_scale:.8e}"
        )

    imex_source_key = "QH_IMEX_source_total_erg_cm3_s"
    imex_increment_key = "QH_IMEX_increment_erg_cm3_s"
    if imex_source_key in c and imex_increment_key in c:
        imex_source = c[imex_source_key]
        imex_increment = c[imex_increment_key]
        applied = c[applied_key]
        scale = max(weighted_l1(applied, volume), np.finfo(float).tiny)

        print("  IMEX-weighted collision diagnostics:")
        for label, values in (
            ("stage source vs solver increment", imex_source - imex_increment),
            ("solver increment vs applied fluid", imex_increment - applied),
            ("stage source vs applied fluid", imex_source - applied),
        ):
            cumulative = np.max(np.abs(cumulative_volume_integral(values, volume)))
            print(
                f"    {label:34s}: "
                f"L1/applied={weighted_l1(values, volume) / scale:.8e}  "
                f"cumulative/applied={cumulative / scale:.8e}"
            )

        species_sum = sum(
            c[f"QH_IMEX_source_species_s{species}_erg_cm3_s"]
            for species in range(1, data.n_species + 1)
        )
        process_sum = sum(
            c[f"QH_IMEX_source_process_{process}_erg_cm3_s"]
            for process in PROCESS_NAMES
        )
        source_scale = max(weighted_l1(imex_source, volume), np.finfo(float).tiny)
        print(
            "    closure by species/process       : "
            f"{weighted_l1(species_sum - imex_source, volume) / source_scale:.8e} / "
            f"{weighted_l1(process_sum - imex_source, volume) / source_scale:.8e}"
        )

    projected_heating_key = "QH_conservative_projection_final_stage_erg_cm3_s"
    terminal_residual_key = "QH_outer_equation_defect_erg_cm3_s"
    cached_heating_key = "QH_cached_final_stage_erg_cm3_s"
    solver_heating_key = "QH_IMEX_increment_erg_cm3_s"
    if all(
        key in c
        for key in (
            projected_heating_key,
            terminal_residual_key,
            cached_heating_key,
            solver_heating_key,
        )
    ):
        projected = c[projected_heating_key]
        terminal_residual = c[terminal_residual_key]
        cached = c[cached_heating_key]
        solver = c[solver_heating_key]
        radiation_defect = cached - projected
        matter_defect = projected - solver
        total_defect = cached - solver
        closure = total_defect - radiation_defect - matter_defect
        scale = max(weighted_l1(solver, volume), np.finfo(float).tiny)
        print("  terminal conservative solve (final-stage BE comparison):")
        print(
            "    L1 defects / L1(solver increment): "
            f"radiation={weighted_l1(radiation_defect, volume) / scale:.8e}  "
            f"matter={weighted_l1(matter_defect, volume) / scale:.8e}  "
            f"total={weighted_l1(total_defect, volume) / scale:.8e}"
        )
        print(
            "    terminal matter residual / L1(solver increment)="
            f"{weighted_l1(terminal_residual, volume) / scale:.8e}  "
            "decomposition closure="
            f"{weighted_l1(closure, volume) / scale:.8e}"
        )

    nodal_path = Path(f"{data.path}.solver_convergence")
    if nodal_path.exists():
        nodal = read_named_table(nodal_path)
        projection_key = "conservative_projection_radiation_relative_residual"
        accepted_key = "accepted_radiation_relative_residual"
        status_key = "conservative_projection_status"
        nodal_solver_key = "QH_solver_increment_erg_cm3_s"
        required = (
            projection_key,
            accepted_key,
            status_key,
            projected_heating_key,
            terminal_residual_key,
            cached_heating_key,
            nodal_solver_key,
        )
        if all(key in nodal for key in required):
            status = nodal[status_key].astype(int)
            active = status != -2
            valid = status >= 0
            passed = status == 1
            failed = status == 0
            invalid = status == -1
            conservation_failed = status == -3
            accepted_residual = nodal[accepted_key]
            projected_residual = nodal[projection_key]
            print("  terminal conservative projection:")
            print(
                "    active/valid/pass/fail/EOS-invalid/conservation-failed="
                f"{np.count_nonzero(active)}/{np.count_nonzero(valid)}/"
                f"{np.count_nonzero(passed)}/{np.count_nonzero(failed)}/"
                f"{np.count_nonzero(invalid)}/"
                f"{np.count_nonzero(conservation_failed)}"
            )
            if np.any(valid):
                print(
                    "    pre-/post-projection radiation residual max="
                    f"{np.max(accepted_residual[valid]):.8e}/"
                    f"{np.max(projected_residual[valid]):.8e}"
                )
                projected = nodal[projected_heating_key]
                cached = nodal[cached_heating_key]
                solver = nodal[nodal_solver_key]
                terminal_residual = nodal[terminal_residual_key]
                scale = max(
                    float(np.max(np.abs(solver[valid]))), np.finfo(float).tiny
                )
                print(
                    "    nodal max defects / max|solver increment|: "
                    f"radiation={np.max(np.abs(cached[valid] - projected[valid])) / scale:.8e}  "
                    f"matter={np.max(np.abs(projected[valid] - solver[valid])) / scale:.8e}  "
                    f"terminal residual={np.max(np.abs(terminal_residual[valid])) / scale:.8e}"
                )
            if np.any(failed):
                failed_radius = nodal["r_cm"][failed] * 1.0e-5
                print(
                    "    failed-projection radial range="
                    f"{np.min(failed_radius):.8e}--"
                    f"{np.max(failed_radius):.8e} km"
                )

    fresh_pre_key = "QH_fresh_pre_limiter_erg_cm3_s"
    cached_final_key = "QH_cached_final_stage_erg_cm3_s"
    if fresh_pre_key in c and cached_final_key in c:
        fresh_post = c["QH_total_erg_cm3_s"]
        fresh_pre = c[fresh_pre_key]
        cached_final = c[cached_final_key]
        applied = c[applied_key]
        scale = max(weighted_l1(applied, volume), np.finfo(float).tiny)

        print("  final-stage limiter separation:")
        for label, values in (
            ("cached source vs fresh pre-limiter", cached_final - fresh_pre),
            ("fresh pre- vs post-limiter RHS", fresh_pre - fresh_post),
            ("cached source vs fresh post-limiter", cached_final - fresh_post),
        ):
            cumulative = np.max(np.abs(cumulative_volume_integral(values, volume)))
            print(
                f"    {label:36s}: "
                f"L1/applied={weighted_l1(values, volume) / scale:.8e}  "
                f"cumulative/applied={cumulative / scale:.8e}"
            )

        if "PL_theta1" in c:
            limited = c["PL_theta1"] < 1.0 - 100.0 * np.finfo(float).eps
            limiter_difference = np.abs(fresh_post - fresh_pre)
            print(
                "    pre/post difference in theta1-limited elements: "
                f"{weighted_l1(limiter_difference[limited], volume[limited]) / scale:.8e}"
            )
            print(
                "    pre/post difference outside theta1-limited elements: "
                f"{weighted_l1(limiter_difference[~limited], volume[~limited]) / scale:.8e}"
            )

    fractional_key = "fractional_eint_change_last_step"
    if fractional_key in c:
        fractional = c[fractional_key]
        index = int(np.argmax(np.abs(fractional)))
        print(
            "  last-step internal-energy change: "
            f"max|Delta eint/eint|={abs(fractional[index]):.8e} "
            f"at r={radius[index] * 1.0e-5:.8e} km"
        )

    inner_max_key = "solver_final_inner_iterations_max"
    if inner_max_key in c:
        inner_min = c["solver_final_inner_iterations_min"]
        inner_max = c[inner_max_key]
        outer_max = c["solver_outer_iterations_max"]
        first_residual = c["solver_first_inner_relative_residual_max"]
        update_fraction = c["solver_final_inner_update_fraction"]
        active = inner_max > 0.0
        accepted_without_update = active & (inner_max == 1.0)
        mixed = active & (inner_min != inner_max)
        print("  final accepted nonlinear iteration:")
        print(
            f"    active elements={np.count_nonzero(active)}  "
            f"all nodes accepted at k_inner=1="
            f"{np.count_nonzero(accepted_without_update)}  "
            f"mixed nodal iteration counts={np.count_nonzero(mixed)}"
        )
        if np.any(active):
            print(
                f"    final k_inner min/max="
                f"{int(np.min(inner_min[active]))}/{int(np.max(inner_max[active]))}  "
                f"final k_outer max={int(np.max(outer_max[active]))}"
            )
            print(
                "    first-inner max(Fnorm/DnuNorm) min/max="
                f"{np.min(first_residual[active]):.8e}/"
                f"{np.max(first_residual[active]):.8e}  "
                f"update-fraction min/max="
                f"{np.min(update_fraction[active]):.3f}/"
                f"{np.max(update_fraction[active]):.3f}"
            )

    if "PL_theta1" in c:
        i_energy = int(np.argmax(np.abs(c["PL_dEnergy_erg"])))
        print("  final-stage radiation limiter:")
        limiter_diagnostics_available = not (
            np.all(c["PL_theta1"] == 0.0)
            and np.all(c["PL_theta2"] == 0.0)
            and np.all(c["PL_dEnergy_erg"] == 0.0)
        )
        if limiter_diagnostics_available:
            i_theta1 = int(np.argmin(c["PL_theta1"]))
            i_theta2 = int(np.argmin(c["PL_theta2"]))
            print(
                f"    min(theta1)={c['PL_theta1'][i_theta1]:.8e} "
                f"at r={radius[i_theta1] * 1.0e-5:.8e} km"
            )
            print(
                f"    min(theta2)={c['PL_theta2'][i_theta2]:.8e} "
                f"at r={radius[i_theta2] * 1.0e-5:.8e} km"
            )
        else:
            print("    limiter factors unavailable (limiters disabled or not called)")
        print(
            f"    max|Delta E_PL|={abs(c['PL_dEnergy_erg'][i_energy]):.8e} erg "
            f"at r={radius[i_energy] * 1.0e-5:.8e} km"
        )

    emab_number_key = emab_ledger_name("Nnet", 1, "cm3_s")
    if emab_number_key in c:
        print("  fresh post-limiter EmAb collision ledger:")
        qh_species_sum = np.zeros_like(radius)
        for species in range(1, data.n_species + 1):
            nin_key = emab_ledger_name("Nin", species, "cm3_s")
            if nin_key not in c:
                continue

            nout = c[emab_ledger_name("Nout", species, "cm3_s")]
            nin = c[nin_key]
            nnet = c[emab_ledger_name("Nnet", species, "cm3_s")]
            ein = c[emab_ledger_name("Ein", species, "erg_cm3_s")]
            eout = c[emab_ledger_name("Eout", species, "erg_cm3_s")]
            qh = c[emab_ledger_name("QH", species, "erg_cm3_s")]
            qh_species_sum += qh

            number_scale = np.maximum(np.abs(nin) + np.abs(nout), np.finfo(float).tiny)
            energy_scale = np.abs(ein) + np.abs(eout)
            active = energy_scale > 100.0 * np.finfo(float).eps * np.max(energy_scale)
            number_rel = np.max(np.abs(nnet - (nin - nout)) / number_scale)
            scalar_mismatch_rel = np.max(
                np.abs(qh - (eout - ein))
                / np.maximum(energy_scale, np.finfo(float).tiny)
            )
            if np.any(active):
                cancellation_ratio = np.abs(eout[active] - ein[active]) / energy_scale[active]
                minimum_index_active = int(np.argmin(cancellation_ratio))
                active_indices = np.flatnonzero(active)
                minimum_index = int(active_indices[minimum_index_active])
                minimum_cancellation = float(cancellation_ratio[minimum_index_active])
                minimum_radius_km = radius[minimum_index] * 1.0e-5
            else:
                minimum_cancellation = float("nan")
                minimum_radius_km = float("nan")

            global_cancellation = weighted_l1(qh, volume) / max(
                weighted_l1(ein, volume) + weighted_l1(eout, volume),
                np.finfo(float).tiny,
            )
            print(
                f"    s={species}: "
                f"max|Nnet-(Nin-Nout)|/(|Nin|+|Nout|)={number_rel:.3e}  "
                f"max|QH-(Eout-Ein)|/(|Ein|+|Eout|)={scalar_mismatch_rel:.3e} "
                "[includes O(v^2) frame term]\n"
                f"         min active |Eout-Ein|/(|Ein|+|Eout|)="
                f"{minimum_cancellation:.3e} at r={minimum_radius_km:.3e} km  "
                f"global L1(QH)/(L1(Ein)+L1(Eout))={global_cancellation:.3e}"
            )

        qh_process = c.get("QH_process_EmAb")
        if qh_process is not None:
            scale = np.maximum(np.abs(qh_process), np.finfo(float).tiny)
            closure = np.max(np.abs(qh_species_sum - qh_process) / scale)
            print(f"    species sum vs QH_process_EmAb: max relative={closure:.3e}")

    nns_number_key = nns_ledger_name("Nnet", 1, "cm3_s")
    if nns_number_key in c:
        print("  NNS collision ledger:")
        qh_species_sum = np.zeros_like(radius)
        for species in range(1, data.n_species + 1):
            nin_key = nns_ledger_name("Nin", species, "cm3_s")
            if nin_key not in c:
                continue

            nout = c[nns_ledger_name("Nout", species, "cm3_s")]
            nin = c[nin_key]
            nnet = c[nns_ledger_name("Nnet", species, "cm3_s")]
            ein = c[nns_ledger_name("Ein", species, "erg_cm3_s")]
            eout = c[nns_ledger_name("Eout", species, "erg_cm3_s")]
            qh = c[nns_ledger_name("QH", species, "erg_cm3_s")]
            qh_species_sum += qh

            number_scale = np.maximum(np.abs(nin) + np.abs(nout), np.finfo(float).tiny)
            energy_scale = np.maximum(np.abs(ein) + np.abs(eout), np.finfo(float).tiny)
            number_rel = np.max(np.abs(nnet) / number_scale)
            scalar_mismatch_rel = np.max(np.abs(qh - (eout - ein)) / energy_scale)
            cancellation = np.min(np.abs(qh) / energy_scale)
            print(
                f"    s={species}: max|Nnet|/(|Nin|+|Nout|)={number_rel:.3e}  "
                f"max|QH-(Eout-Ein)|/(|Ein|+|Eout|)={scalar_mismatch_rel:.3e}  "
                f"min|QH|/(|Ein|+|Eout|)={cancellation:.3e}"
            )

        qh_process = c.get("QH_process_NNS")
        if qh_process is not None:
            scale = np.maximum(np.abs(qh_process), np.finfo(float).tiny)
            closure = np.max(np.abs(qh_species_sum - qh_process) / scale)
            print(f"    species sum vs QH_process_NNS: max relative={closure:.3e}")

    pairing_names = [
        process for process in GENERAL_LEPTON_PROCESS_NAMES
        if f"pairing_abs_max_{process}_cm3_s" in c
    ]
    if pairing_names:
        print("  paired general-lepton collision identities:")
        for process in pairing_names:
            absolute = c[f"pairing_abs_max_{process}_cm3_s"]
            relative = c[f"pairing_rel_max_{process}"]
            print(
                f"    {process:4s}: max absolute={np.max(absolute):.8e} cm^-3 s^-1; "
                f"max relative={np.max(relative):.8e}"
            )

    print("  inverse-MFP minima [cm^-1]:")
    for species in range(1, data.n_species + 1):
        names = [
            f"inv_mfp_{process}_s{species}_cm_inv" for process in PROCESS_NAMES
            if f"inv_mfp_{process}_s{species}_cm_inv" in c
        ]
        if not names:
            continue
        minimum = min(float(np.min(c[name])) for name in names)
        print(f"    species {species}: {minimum: .8e}")

    print("  tau=2/3 surfaces and temperature checks:")
    for species in range(1, data.n_species + 1):
        required = [f"inv_mfp_{process}_s{species}_cm_inv" for process in PROCESS_NAMES]
        if any(name not in c for name in required):
            continue

        noniso = sum(process_inverse_mfp(data, species, process) for process in NONISO_PROCESS_NAMES)
        iso = process_inverse_mfp(data, species, "Iso")
        transport = noniso + iso
        effective_squared = noniso * transport

        if np.any(effective_squared < 0.0):
            energy_radius = float("nan")
        else:
            energy_radius = radius_at_tau(radius, width, np.sqrt(effective_squared))

        if np.any(transport < 0.0):
            transport_radius = float("nan")
        else:
            transport_radius = radius_at_tau(radius, width, transport)

        t_fluid = c["T_fluid_MeV"]
        t_spec = c[f"T_spec_s{species}"]
        t_flux = c[f"T_flux_s{species}"]
        flux_factor = c[f"flux_factor_s{species}"]

        at_energy = nearest(t_spec / t_fluid, radius, energy_radius)
        free_stream = (radius >= transport_radius) & (flux_factor >= 0.9) & (t_spec != 0.0)
        if np.any(free_stream):
            flux_ratio_median = float(np.median(t_flux[free_stream] / t_spec[free_stream]))
        else:
            flux_ratio_median = float("nan")

        print(
            f"    species {species}: r_energy={energy_radius: .8e} cm  "
            f"r_transport={transport_radius: .8e} cm  "
            f"Tspec/Tfluid@r_energy={at_energy: .8e}  "
            f"median(Tflux/Tspec; F>=0.9)={flux_ratio_median: .8e}"
        )


def comparison_columns(columns: Iterable[str]) -> Iterable[str]:
    for name in columns:
        if (
            name.startswith("QH_")
            or name.startswith("EmAb_")
            or name.startswith("NNS_")
            or name.startswith("inv_mfp_")
        ):
            yield name


def report_successive_convergence(series: list[DiagnosticFile]) -> None:
    if len(series) < 2:
        return

    series = sorted(series, key=lambda item: item.n_energy_elements)
    print("\nSuccessive-resolution volume-weighted differences:")

    for coarse, fine in zip(series[:-1], series[1:], strict=True):
        if coarse.n_species != fine.n_species:
            raise ValueError(
                "successive-resolution comparison requires the same species count; "
                f"{coarse.path} has {coarse.n_species}, while "
                f"{fine.path} has {fine.n_species}"
            )
        r_fine = fine.columns["r_cm"]
        dV_fine = fine.columns["dV_cm3"]
        print(f"  nE {coarse.n_energy_elements} -> {fine.n_energy_elements}")

        common = set(comparison_columns(coarse.columns)) & set(comparison_columns(fine.columns))
        for name in sorted(common):
            coarse_on_fine = np.interp(r_fine, coarse.columns["r_cm"], coarse.columns[name])
            difference = weighted_l1(fine.columns[name] - coarse_on_fine, dV_fine)
            scale = weighted_l1(fine.columns[name], dV_fine)
            relative = difference / scale if scale != 0.0 else float("nan")
            print(f"    {name:28s} abs={difference: .8e}  rel={relative: .8e}")

        nns_key = nns_ledger_name("Ein", 1, "erg_cm3_s")
        if nns_key not in coarse.columns or nns_key not in fine.columns:
            continue

        print("  NNS weak heating comparison:")
        for species in range(1, fine.n_species + 1):
            names = {
                quantity: nns_ledger_name(quantity, species, "erg_cm3_s")
                for quantity in ("Ein", "Eout", "QH")
            }
            if any(name not in coarse.columns or name not in fine.columns for name in names.values()):
                continue

            ein_coarse = np.interp(r_fine, coarse.columns["r_cm"], coarse.columns[names["Ein"]])
            eout_coarse = np.interp(r_fine, coarse.columns["r_cm"], coarse.columns[names["Eout"]])
            qh_coarse = np.interp(r_fine, coarse.columns["r_cm"], coarse.columns[names["QH"]])
            ein_fine = fine.columns[names["Ein"]]
            eout_fine = fine.columns[names["Eout"]]
            qh_fine = fine.columns[names["QH"]]

            transfer_scale = weighted_integral(ein_fine + eout_fine, dV_fine)
            ein_difference = weighted_l1(ein_fine - ein_coarse, dV_fine)
            eout_difference = weighted_l1(eout_fine - eout_coarse, dV_fine)
            cumulative_difference = np.max(np.abs(
                cumulative_volume_integral(qh_fine - qh_coarse, dV_fine)
            ))
            qh_l1 = weighted_l1(qh_fine, dV_fine)
            signed_difference = abs(
                weighted_integral(qh_fine - qh_coarse, dV_fine)
            )

            if transfer_scale != 0.0:
                ein_relative = ein_difference / transfer_scale
                eout_relative = eout_difference / transfer_scale
                cumulative_transfer_relative = (
                    cumulative_difference / transfer_scale
                )
                signed_relative = signed_difference / transfer_scale
            else:
                ein_relative = float("nan")
                eout_relative = float("nan")
                cumulative_transfer_relative = float("nan")
                signed_relative = float("nan")

            if qh_l1 != 0.0:
                cumulative_qh_relative = cumulative_difference / qh_l1
            else:
                cumulative_qh_relative = float("nan")

            print(
                f"    s={species}: "
                f"L1(Ein)/transfer={ein_relative:.3e}  "
                f"L1(Eout)/transfer={eout_relative:.3e}  "
                f"max|Delta integral QH|/transfer={cumulative_transfer_relative:.3e}  "
                f"/L1(QH)={cumulative_qh_relative:.3e}  "
                f"|Delta total QH|/transfer={signed_relative:.3e}"
            )


def plot_file(
    data: DiagnosticFile,
    output: Path,
    include_nns_ledger: bool = False,
    include_nns_ledger_terms: bool = False,
    include_emab_ledger: bool = False,
    include_emab_ledger_terms: bool = False,
    include_process_heating: bool = False,
    include_step_diagnostics: bool = False,
) -> None:
    """Write a compact radial overview for the electron-flavour diagnostics."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def place_legend(axis, *, fontsize="small") -> None:
        """Place an axis legend beside, rather than over, the plotted data."""

        handles, labels = axis.get_legend_handles_labels()
        if handles:
            axis.legend(
                handles,
                labels,
                loc="upper left",
                bbox_to_anchor=(1.02, 1.0),
                borderaxespad=0.0,
                fontsize=fontsize,
            )

    c = data.columns
    radius_km = c["r_cm"] * 1.0e-5
    width = c["dr_cm"]
    electron_species_info = plotted_species(data, maximum_species=2)
    process_heating_species_info = plotted_species(data, maximum_species=3)
    electron_species = tuple(info[0] for info in electron_species_info)

    if include_nns_ledger or include_nns_ledger_terms:
        required_ledger_columns = [
            nns_ledger_name(quantity, i_s, "erg_cm3_s")
            for i_s, _, _ in electron_species_info
            for quantity in ("Ein", "Eout")
        ]
        missing = [name for name in required_ledger_columns if name not in c]
        if missing:
            raise ValueError(
                "the NNS-ledger plot options require an NNS-ledger diagnostic file; "
                f"missing {', '.join(missing)}"
            )

    if include_emab_ledger or include_emab_ledger_terms:
        required_ledger_columns = [
            emab_ledger_name(quantity, i_s, "erg_cm3_s")
            for i_s, _, _ in electron_species_info
            for quantity in ("Ein", "Eout")
        ]
        missing = [name for name in required_ledger_columns if name not in c]
        if missing:
            raise ValueError(
                "the EmAb-ledger plot options require an EmAb-ledger diagnostic file; "
                f"missing {', '.join(missing)}"
            )

    if include_process_heating:
        required_process_columns = [
            f"QH_{process}_s{species}"
            for species, _, _ in process_heating_species_info
            for process in HEATING_PROCESS_NAMES
        ]
        required_process_columns.append("QH_applied_last_step_erg_cm3_s")
        missing = [name for name in required_process_columns if name not in c]
        if missing:
            raise ValueError(
                "--plot-process-heating requires process-by-species heating columns; "
                f"missing {', '.join(missing)}"
            )
        plotted_heating_processes = tuple(
            process for process in HEATING_PROCESS_NAMES
            if (
                any(
                    np.any(c[f"QH_{process}_s{i_s}"] != 0.0)
                    for i_s, _, _ in process_heating_species_info
                )
                or np.any(
                    c[f"QH_IMEX_source_process_{process}_erg_cm3_s"] != 0.0
                )
            )
        )
        process_heating_rows = (len(plotted_heating_processes) + 1) // 2
    else:
        plotted_heating_processes = ()
        process_heating_rows = 0

    if include_step_diagnostics:
        required_step_columns = [
            "QH_total_erg_cm3_s",
            "QH_applied_last_step_erg_cm3_s",
            "QH_IMEX_source_total_erg_cm3_s",
            "QH_IMEX_increment_erg_cm3_s",
            "fractional_eint_change_last_step",
            "PL_theta1",
            "PL_theta2",
            "PL_dEnergy_erg",
        ]
        missing = [name for name in required_step_columns if name not in c]
        if missing:
            raise ValueError(
                "--plot-step-diagnostics requires the IMEX/limiter diagnostic format; "
                f"missing {', '.join(missing)}"
            )

    if (
        include_nns_ledger or include_nns_ledger_terms
        or include_emab_ledger or include_emab_ledger_terms
        or include_process_heating or include_step_diagnostics
    ):
        n_rows = (
            2 + int(include_nns_ledger) + int(include_nns_ledger_terms)
            + int(include_emab_ledger) + int(include_emab_ledger_terms)
            + process_heating_rows
            + int(include_step_diagnostics)
        )
        fig = plt.figure(
            figsize=(18, 11 + 2.75 * (n_rows - 2)),
            constrained_layout=True,
        )
        grid = fig.add_gridspec(n_rows, 2)
        axes = np.empty((2, 2), dtype=object)
        axes[0, 0] = fig.add_subplot(grid[0, 0])
        axes[0, 1] = fig.add_subplot(grid[0, 1])
        axes[1, 0] = fig.add_subplot(grid[1, 0])
        axes[1, 1] = fig.add_subplot(grid[1, 1])
        next_row = 2
        if include_nns_ledger:
            ledger_axis = fig.add_subplot(grid[next_row, :])
            next_row += 1
        else:
            ledger_axis = None
        if include_nns_ledger_terms:
            ledger_terms_axis = fig.add_subplot(grid[next_row, :])
            next_row += 1
        else:
            ledger_terms_axis = None
        if include_emab_ledger:
            emab_ledger_axis = fig.add_subplot(grid[next_row, :])
            next_row += 1
        else:
            emab_ledger_axis = None
        if include_emab_ledger_terms:
            emab_ledger_terms_axis = fig.add_subplot(grid[next_row, :])
            next_row += 1
        else:
            emab_ledger_terms_axis = None
        if include_process_heating:
            process_heating_axes = {}
            for i_process, process in enumerate(plotted_heating_processes):
                process_heating_axes[process] = fig.add_subplot(
                    grid[next_row + i_process // 2, i_process % 2]
                )
            next_row += process_heating_rows
        else:
            process_heating_axes = None
        if include_step_diagnostics:
            step_heating_axis = fig.add_subplot(grid[next_row, 0])
            limiter_axis = fig.add_subplot(grid[next_row, 1])
        else:
            step_heating_axis = None
            limiter_axis = None
    else:
        fig, axes = plt.subplots(2, 2, figsize=(18, 11), constrained_layout=True)
        ledger_axis = None
        ledger_terms_axis = None
        emab_ledger_axis = None
        emab_ledger_terms_axis = None
        process_heating_axes = None
        step_heating_axis = None
        limiter_axis = None
    fig.suptitle(
        f"{data.path.name}: nE={data.n_energy_elements}, "
        f"nDOFE={data.n_energy_nodes}"
    )

    axis = axes[0, 0]
    for i_s, label, color in electron_species_info:
        values = c[f"QH_species_s{i_s}"]
        if np.any(values != 0.0):
            axis.plot(radius_km, values, color=color, label=label)
    axis.set_xscale("log")
    axis.set_yscale("symlog", linthresh=1.0e30)
    axis.set_xlabel("radius [km]")
    axis.set_ylabel(r"$Q_H$ [erg cm$^{-3}$ s$^{-1}$]")
    axis.set_title("instantaneous post-limiter collision heating by species")
    place_legend(axis)

    axis = axes[0, 1]
    for i_s, style in zip(electron_species, ("-", "--"), strict=True):
        for process, color in (("EmAb", "C0"), ("NNS", "C3")):
            values = process_inverse_mfp(data, i_s, process)
            positive = values > 0.0
            if np.any(positive):
                axis.plot(
                    radius_km[positive],
                    values[positive],
                    linestyle=style,
                    color=color,
                    label=f"{process}, s={i_s}",
                )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("radius [km]")
    axis.set_ylabel(r"energy-averaged inverse MFP [cm$^{-1}$]")
    axis.set_title("electron-flavour inverse MFP")
    place_legend(axis)

    if ledger_axis is not None:
        for i_s, label, color in electron_species_info:
            ein = c[nns_ledger_name("Ein", i_s, "erg_cm3_s")]
            eout = c[nns_ledger_name("Eout", i_s, "erg_cm3_s")]
            denominator = np.abs(ein) + np.abs(eout)
            imbalance = np.divide(
                eout - ein,
                denominator,
                out=np.zeros_like(denominator),
                where=denominator > 0.0,
            )
            ledger_axis.plot(radius_km, imbalance, color=color, label=label)
        ledger_axis.axhline(0.0, color="k", linewidth=0.8)
        ledger_axis.set_xscale("log")
        ledger_axis.set_yscale("symlog", linthresh=1.0e-6)
        ledger_axis.set_xlabel("radius [km]")
        ledger_axis.set_ylabel(r"$(E_{\rm out}-E_{\rm in})/(|E_{\rm out}|+|E_{\rm in}|)$")
        ledger_axis.set_title("NNS normalized energy imbalance")
        place_legend(ledger_axis)

    if ledger_terms_axis is not None:
        for i_s, label, color in electron_species_info:
            ein = c[nns_ledger_name("Ein", i_s, "erg_cm3_s")]
            eout = c[nns_ledger_name("Eout", i_s, "erg_cm3_s")]
            positive_ein = ein > 0.0
            positive_eout = eout > 0.0
            if np.any(positive_ein):
                ledger_terms_axis.plot(
                    radius_km[positive_ein], ein[positive_ein], color=color,
                    linestyle="--", label=fr"$E_{{\rm in}}$, {label}"
                )
            if np.any(positive_eout):
                ledger_terms_axis.plot(
                    radius_km[positive_eout], eout[positive_eout], color=color,
                    linestyle="-", label=fr"$E_{{\rm out}}$, {label}"
                )
        ledger_terms_axis.set_xscale("log")
        ledger_terms_axis.set_yscale("log")
        ledger_terms_axis.set_xlabel("radius [km]")
        ledger_terms_axis.set_ylabel(r"energy transfer [erg cm$^{-3}$ s$^{-1}$]")
        ledger_terms_axis.set_title("NNS energy-transfer terms (solid: out; dashed: in)")
        place_legend(ledger_terms_axis)

    if emab_ledger_axis is not None:
        for i_s, label, color in electron_species_info:
            ein = c[emab_ledger_name("Ein", i_s, "erg_cm3_s")]
            eout = c[emab_ledger_name("Eout", i_s, "erg_cm3_s")]
            denominator = np.abs(ein) + np.abs(eout)
            imbalance = np.divide(
                eout - ein,
                denominator,
                out=np.zeros_like(denominator),
                where=denominator > 0.0,
            )
            emab_ledger_axis.plot(radius_km, imbalance, color=color, label=label)
        emab_ledger_axis.axhline(0.0, color="k", linewidth=0.8)
        emab_ledger_axis.set_xscale("log")
        emab_ledger_axis.set_yscale("symlog", linthresh=1.0e-12)
        emab_ledger_axis.set_xlabel("radius [km]")
        emab_ledger_axis.set_ylabel(
            r"$(E_{\rm out}-E_{\rm in})/(|E_{\rm out}|+|E_{\rm in}|)$"
        )
        emab_ledger_axis.set_title(
            "fresh post-limiter EmAb normalized energy imbalance"
        )
        place_legend(emab_ledger_axis)

    if emab_ledger_terms_axis is not None:
        for i_s, label, color in electron_species_info:
            ein = c[emab_ledger_name("Ein", i_s, "erg_cm3_s")]
            eout = c[emab_ledger_name("Eout", i_s, "erg_cm3_s")]
            positive_ein = ein > 0.0
            positive_eout = eout > 0.0
            if np.any(positive_ein):
                emab_ledger_terms_axis.plot(
                    radius_km[positive_ein], ein[positive_ein], color=color,
                    linestyle="--", label=fr"$E_{{\rm in}}$, {label}"
                )
            if np.any(positive_eout):
                emab_ledger_terms_axis.plot(
                    radius_km[positive_eout], eout[positive_eout], color=color,
                    linestyle="-", label=fr"$E_{{\rm out}}$, {label}"
                )
        emab_ledger_terms_axis.set_xscale("log")
        emab_ledger_terms_axis.set_yscale("log")
        emab_ledger_terms_axis.set_xlabel("radius [km]")
        emab_ledger_terms_axis.set_ylabel(
            r"energy transfer [erg cm$^{-3}$ s$^{-1}$]"
        )
        emab_ledger_terms_axis.set_title(
            "fresh post-limiter EmAb terms "
            "(solid: neutrino removal; dashed: neutrino production)"
        )
        place_legend(emab_ledger_terms_axis)

    if process_heating_axes is not None:
        for i_process, (process, process_heating_axis) in enumerate(
            process_heating_axes.items()
        ):
            values_by_species = [
                c[f"QH_{process}_s{i_s}"]
                for i_s, _, _ in process_heating_species_info
            ]
            maximum = max(
                float(np.max(np.abs(values))) for values in values_by_species
            )

            # Show the instantaneous total once for scale.  Finite-step
            # quantities belong in the separate last-step panel below.
            if i_process == 0:
                total = c["QH_total_erg_cm3_s"]
                maximum = max(
                    maximum,
                    float(np.max(np.abs(total))),
                )

            for values, (_, label, color) in zip(
                values_by_species, process_heating_species_info, strict=True
            ):
                process_heating_axis.plot(
                    radius_km, values, color=color, linestyle="-",
                    label=fr"{label}, {process}",
                )
            if i_process == 0:
                process_heating_axis.plot(
                    radius_km, total, color="k", linestyle="--", linewidth=1.2,
                    label="instantaneous final-state total RHS (all processes)",
                )

            process_heating_axis.axhline(0.0, color="k", linewidth=0.8)
            process_heating_axis.set_xscale("log")
            if maximum > 0.0:
                process_heating_axis.set_yscale("symlog", linthresh=maximum * 1.0e-8)
            process_heating_axis.set_xlabel("radius [km]")
            process_heating_axis.set_ylabel(r"$Q_H$ [erg cm$^{-3}$ s$^{-1}$]")
            process_heating_axis.set_title(
                f"instantaneous final-state {process} heating"
            )
            place_legend(process_heating_axis, fontsize="x-small")
    if step_heating_axis is not None:
        heating_fields = [
            ("QH_IMEX_source_total_erg_cm3_s", "IMEX-weighted stage source", "C4", "-."),
            ("QH_IMEX_increment_erg_cm3_s", "IMEX-weighted solver increment", "C2", "-"),
            ("QH_applied_last_step_erg_cm3_s", "actual fluid change / dt", "0.35", ":"),
        ]
        for name, label, color, style in heating_fields:
            step_heating_axis.plot(
                radius_km, c[name], color=color, linestyle=style,
                linewidth=1.5, label=label,
            )
        maximum = max(float(np.max(np.abs(c[name]))) for name, _, _, _ in heating_fields)
        step_heating_axis.axhline(0.0, color="k", linewidth=0.8)
        step_heating_axis.set_xscale("log")
        if maximum > 0.0:
            step_heating_axis.set_yscale("symlog", linthresh=maximum * 1.0e-8)
        step_heating_axis.set_xlabel("radius [km]")
        step_heating_axis.set_ylabel(r"$Q_H$ [erg cm$^{-3}$ s$^{-1}$]")
        step_heating_axis.set_title("last-step heating comparison")
        place_legend(step_heating_axis, fontsize="x-small")

        limiter_diagnostics_available = not (
            np.all(c["PL_theta1"] == 0.0)
            and np.all(c["PL_theta2"] == 0.0)
            and np.all(c["PL_dEnergy_erg"] == 0.0)
        )
        limiter_axis.set_xscale("log")
        limiter_axis.set_xlabel("radius [km]")
        fractional_change = np.abs(c["fractional_eint_change_last_step"])
        positive_fractional_change = fractional_change[fractional_change > 0.0]

        if limiter_diagnostics_available:
            limiter_axis.plot(radius_km, c["PL_theta1"], label=r"$\theta_1$")
            limiter_axis.plot(radius_km, c["PL_theta2"], label=r"$\theta_2$")
            limiter_axis.set_ylabel(r"limiter factor $\theta$")
            theta_min = min(
                float(np.min(c["PL_theta1"])),
                float(np.min(c["PL_theta2"])),
            )
            theta_range = max(1.0 - theta_min, 0.05)
            limiter_axis.set_ylim(max(0.0, theta_min - 0.1 * theta_range), 1.02)
            fractional_axis = limiter_axis.twinx()
            limiter_axis.set_title("limiter and fractional fluid update")
        else:
            fractional_axis = limiter_axis
            limiter_axis.set_title("fractional fluid internal-energy update")

        fractional_axis.plot(
            radius_km,
            fractional_change,
            color="C2",
            label=r"$|\Delta e_{\rm int}/e_{\rm int}|$",
        )
        fractional_axis.set_yscale("log")
        fractional_axis.set_ylabel(r"$|\Delta e_{\rm int}/e_{\rm int}|$")
        if positive_fractional_change.size:
            lower_decade = 10.0 ** np.floor(
                np.log10(float(np.min(positive_fractional_change)))
            )
            upper_decade = 10.0 ** np.ceil(
                np.log10(float(np.max(positive_fractional_change)))
            )
            if upper_decade <= lower_decade:
                upper_decade = 10.0 * lower_decade
            fractional_axis.set_ylim(lower_decade, upper_decade)

        limiter_handles, limiter_labels = limiter_axis.get_legend_handles_labels()
        if fractional_axis is limiter_axis:
            fractional_handles, fractional_labels = [], []
        else:
            fractional_handles, fractional_labels = (
                fractional_axis.get_legend_handles_labels()
            )
        limiter_axis.legend(
            limiter_handles + fractional_handles,
            limiter_labels + fractional_labels,
            loc="upper left",
            bbox_to_anchor=(1.18, 1.0),
            borderaxespad=0.0,
            fontsize="x-small",
        )

    axis = axes[1, 0]
    for i_s, _, color in electron_species_info:
        noniso = sum(process_inverse_mfp(data, i_s, process) for process in NONISO_PROCESS_NAMES)
        iso = process_inverse_mfp(data, i_s, "Iso")
        transport = noniso + iso
        effective = np.sqrt(np.maximum(noniso * transport, 0.0))
        axis.plot(radius_km, optical_depth_at_centers(width, transport), color=color,
                  label=f"transport, s={i_s}")
        axis.plot(radius_km, optical_depth_at_centers(width, effective), color=color,
                  linestyle="--", label=f"effective, s={i_s}")
    axis.axhline(TAU_SURFACE, color="k", linewidth=0.8)
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("radius [km]")
    axis.set_ylabel(r"$\tau(r)$")
    axis.set_title(r"gray optical depths")
    place_legend(axis)

    axis = axes[1, 1]
    axis.plot(radius_km, c["T_fluid_MeV"], color="k", linewidth=1.5,
              label=r"$T_{\rm fluid}$")
    for i_s, _, color in electron_species_info:
        t_spec = c[f"T_spec_s{i_s}"]
        t_flux = c[f"T_flux_s{i_s}"]
        noniso = sum(process_inverse_mfp(data, i_s, process) for process in NONISO_PROCESS_NAMES)
        transport = noniso + process_inverse_mfp(data, i_s, "Iso")
        effective_squared = noniso * transport
        if np.all(effective_squared >= 0.0):
            energy_radius = radius_at_tau(c["r_cm"], width, np.sqrt(effective_squared))
            if np.isfinite(energy_radius):
                axis.axvline(energy_radius * 1.0e-5, color=color, linestyle=":", linewidth=1.2,
                             label=fr"$r_{{\rm energy}}$, s={i_s}")
        if np.all(transport >= 0.0):
            transport_radius = radius_at_tau(c["r_cm"], width, transport)
            if np.isfinite(transport_radius):
                axis.axvline(transport_radius * 1.0e-5, color=color, linestyle="-.", linewidth=1.2,
                             label=fr"$r_{{\rm transport}}$, s={i_s}")
        valid_spec = t_spec > 0.0
        valid_flux = t_flux > 0.0
        axis.plot(radius_km[valid_spec], t_spec[valid_spec], color=color,
                  label=fr"$T_\nu$, s={i_s}")
        axis.plot(radius_km[valid_flux], t_flux[valid_flux], color=color,
                  linestyle="--", label=fr"$T_{{\rm flux}}$, s={i_s}")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("radius [km]")
    axis.set_ylabel("temperature [MeV]")
    axis.set_title("fluid and neutrino temperatures")
    place_legend(axis)

    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument(
        "--compare",
        action="store_true",
        help="report successive spectral-resolution differences after per-file summaries",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        metavar="FILE.png",
        help="write a radial overview; requires exactly one diagnostic input",
    )
    parser.add_argument(
        "--plot-nns-ledger",
        action="store_true",
        help=(
            "add a panel of the NNS normalized energy imbalance; requires "
            "--plot and an NNS-ledger diagnostic file"
        ),
    )
    parser.add_argument(
        "--plot-nns-ledger-terms",
        action="store_true",
        help=(
            "add a panel of the positive NNS ledger energy-transfer terms; "
            "requires --plot and an NNS-ledger diagnostic file"
        ),
    )
    parser.add_argument(
        "--plot-emab-ledger",
        action="store_true",
        help=(
            "add a panel of the fresh post-limiter EmAb normalized energy "
            "imbalance; requires --plot and an EmAb-ledger diagnostic file"
        ),
    )
    parser.add_argument(
        "--plot-emab-ledger-terms",
        action="store_true",
        help=(
            "add a panel of the fresh post-limiter positive EmAb "
            "neutrino-production/removal terms; requires --plot and an "
            "EmAb-ledger diagnostic file"
        ),
    )
    parser.add_argument(
        "--plot-process-heating",
        action="store_true",
        help=(
            "add separate process heating panels for the available nue, nuebar, "
            "and numu species; "
            "requires --plot and the newer standalone diagnostic format"
        ),
    )
    parser.add_argument(
        "--plot-step-diagnostics",
        action="store_true",
        help=(
            "add IMEX-weighted/applied-heating and limiter panels; requires "
            "--plot and the newest standalone diagnostic format"
        ),
    )
    args = parser.parse_args()

    diagnostics = [read_diagnostic(path) for path in args.files]
    for diagnostic in diagnostics:
        report_file(diagnostic)
    if args.compare:
        report_successive_convergence(diagnostics)
    if args.plot is not None:
        if len(diagnostics) != 1:
            parser.error("--plot requires exactly one diagnostic input")
        plot_file(
            diagnostics[0], args.plot, args.plot_nns_ledger,
            args.plot_nns_ledger_terms, args.plot_emab_ledger,
            args.plot_emab_ledger_terms, args.plot_process_heating,
            args.plot_step_diagnostics,
        )
    elif (
        args.plot_nns_ledger or args.plot_nns_ledger_terms
        or args.plot_emab_ledger or args.plot_emab_ledger_terms
        or args.plot_process_heating or args.plot_step_diagnostics
    ):
        parser.error("--plot options require --plot")


if __name__ == "__main__":
    main()
