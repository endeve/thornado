#!/usr/bin/env python3

"""Compare self-consistent standard-LLF and source-aware-AP trajectories."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import analyze_matter_exchange_attribution as attribution


@dataclass
class RunData:
    label: str
    color: str
    prefix: Path
    elements: np.ndarray
    nodes: np.ndarray
    main: dict[str, np.ndarray]
    density: dict[tuple[int, int], float]
    histories: dict[tuple[str, str], tuple[np.ndarray, ...]]
    inner_histories: dict[tuple[str, str], tuple[np.ndarray, ...]]


def save_figure(fig: plt.Figure, output: Path) -> None:
    """Write a high-resolution preview and a vector PDF."""

    fig.savefig(output, dpi=300, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "directory",
        type=Path,
        help="directory produced with --case both_no_projection",
    )
    parser.add_argument("--cycle", type=int, default=None)
    parser.add_argument("--r-max-km", type=float, default=200.0)
    parser.add_argument("--inner-r-max-km", type=float, default=25.0)
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=None,
        help="default: <directory>/self_consistent_llf_ap",
    )
    return parser.parse_args()


def diagnostic_prefix(directory: Path, stem: str) -> Path:
    return directory / f"matter_diagnostics_{stem}.dat"


def verify_actual_only(path: Path) -> None:
    mode = None
    component_index: int | None = None
    components: set[str] = set()
    with path.open() as stream:
        for line in stream:
            if line.startswith("# cycle"):
                component_index = line[1:].split().index("component")
            if line.startswith("# attribution_mode"):
                mode = line.split()[-1]
            elif line.startswith("#") or not line.strip():
                continue
            else:
                if component_index is None:
                    raise ValueError(f"{path}: missing column header")
                components.add(line.split()[component_index])
                if len(components) > 1:
                    break
    if mode != "actual_only" or components != {"actual"}:
        raise ValueError(
            f"{path} is not an actual-only attribution file: "
            f"mode={mode}, components={sorted(components)}"
        )


def load_run(
    prefix: Path,
    label: str,
    color: str,
    cycle: int,
    r_max_cm: float,
    inner_r_max_cm: float,
) -> RunData:
    attribution_prefix = Path(str(prefix) + ".exchange_attribution")
    element_path = Path(str(attribution_prefix) + ".elements")
    node_path = Path(str(attribution_prefix) + ".nodes")
    verify_actual_only(element_path)
    verify_actual_only(node_path)
    elements, histories, inner_histories = (
        attribution.load_element_cycle_and_histories(
            element_path, cycle, np.inf, inner_r_max_cm
        )
    )
    return RunData(
        label=label,
        color=color,
        prefix=prefix,
        elements=elements,
        nodes=attribution.load_cycle_table(node_path, cycle),
        main=attribution.load_main_diagnostic(prefix),
        density=attribution.load_nodal_density(
            Path(str(prefix) + ".frozen_transient")
        ),
        histories=histories,
        inner_histories=inner_histories,
    )


def element_specific_rate(run: RunData, quantity: str) -> tuple[np.ndarray, np.ndarray]:
    rows = attribution.select(
        run.elements, int(run.elements["cycle"][0]), quantity, "actual", np.inf
    )
    radius = rows["r_cm"] / 1.0e5
    mask = run.main["r_cm"] <= float(np.max(rows["r_cm"])) + 1.0e-12
    density = run.main["D_g_cm3"][mask]
    if density.size != rows.size:
        raise ValueError(f"{run.prefix}: incompatible matter and element grids")
    if quantity == "QH":
        rate = rows["step_rate_mean"] / density / attribution.SPECIFIC_HEATING_SCALE
    else:
        rate = (
            attribution.ATOMIC_MASS_UNIT_G
            * rows["step_rate_mean"]
            / density
        )
    return radius, rate


def local_timescale(run: RunData, quantity: str) -> tuple[np.ndarray, np.ndarray]:
    cycle = int(run.elements["cycle"][0])
    rows = attribution.select(run.elements, cycle, quantity, "actual", np.inf)
    reference = attribution.element_reference(run.main, rows, quantity)
    magnitude = np.abs(rows["step_rate_mean"])
    timescale = np.full(magnitude.shape, np.nan)
    active = magnitude > 0.0
    timescale[active] = reference[active] / magnitude[active]
    return rows["r_cm"] / 1.0e5, timescale


def plot_developer_summary(
    runs: tuple[RunData, RunData],
    cycle: int,
    r_max_km: float,
    output: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.5))
    for row, quantity in enumerate(("QH", "QL")):
        for run in runs:
            radius, rate = element_specific_rate(run, quantity)
            within_plot = radius <= r_max_km
            axes[row, 0].plot(
                radius[within_plot],
                rate[within_plot],
                color=run.color,
                linewidth=1.5,
                label=run.label,
            )
            radius, timescale = local_timescale(run, quantity)
            within_plot = radius <= r_max_km
            axes[row, 1].plot(
                radius[within_plot],
                timescale[within_plot],
                color=run.color,
                linewidth=1.5,
                label=run.label,
            )
        for column in range(2):
            axis = axes[row, column]
            axis.set_xscale("log")
            axis.set_xlim(max(0.5, float(np.min(radius))), r_max_km)
            axis.set_xlabel(r"radius [$\mathrm{km}$]")
            axis.grid(alpha=0.25, which="both")
            axis.legend(fontsize=8.0, loc="best")
        axes[row, 0].axhline(0.0, color="0.45", linewidth=0.65)

        if quantity == "QH":
            axes[row, 0].set_ylabel(
                r"$\mathrm{d}\epsilon/\mathrm{d}t$ "
                r"[$10^{21}\ \mathrm{erg\,g^{-1}\,s^{-1}}$]"
            )
            axes[row, 0].set_title("Element-average specific heating")
            axes[row, 1].set_title(r"$e_{\rm int}/|Q_H|$")
        else:
            axes[row, 0].set_ylabel(r"$\mathrm{d}Y_e/\mathrm{d}t$ [$\mathrm{s^{-1}}$]")
            axes[row, 0].set_title("Element-average electron-fraction change")
            axes[row, 1].set_title(r"$n_e/|Q_L|$")
        axes[row, 1].set_ylabel("local matter-change timescale [s]")
        axes[row, 1].set_yscale("log")

    time_s = float(runs[0].elements["time_s"][0])
    fig.suptitle(
        "Self-consistent standard and source-aware AP LLF trajectories\n"
        rf"cycle {cycle}, $t={time_s:.3e}\,\mathrm{{s}}$",
        fontsize=15,
    )
    fig.tight_layout()
    save_figure(fig, output)
    plt.close(fig)


def element_intervals(
    run: RunData,
    quantity: str,
    cycle: int,
    values: np.ndarray,
    r_max_km: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return physical element edges paired with element-average values."""

    nodal_rows = attribution.select(
        run.nodes, cycle, quantity, "actual", np.inf
    )
    left: list[float] = []
    right: list[float] = []
    elements = np.unique(nodal_rows["element"].astype(int))
    selected_elements = []
    for element in elements:
        block = nodal_rows[nodal_rows["element"].astype(int) == element]
        if float(np.mean(block["r_cm"])) <= r_max_km * 1.0e5:
            selected_elements.append(element)
    if len(selected_elements) != values.size:
        raise ValueError(f"{run.prefix}: incompatible nodal and element grids")
    for element in selected_elements:
        block = nodal_rows[nodal_rows["element"].astype(int) == element]
        block = block[np.argsort(block["node"])]
        gauss_nodes = np.polynomial.legendre.leggauss(block.size)[0]
        radius_nodes = np.asarray(block["r_cm"], dtype=float) / 1.0e5
        center = float(np.mean(radius_nodes))
        half_width = float(
            np.dot(gauss_nodes, radius_nodes - center)
            / np.dot(gauss_nodes, gauss_nodes)
        )
        left.append(center - half_width)
        right.append(center + half_width)
    return np.asarray(left), np.asarray(right), values


def plot_horizontal_intervals(
    axis: plt.Axes,
    left: np.ndarray,
    right: np.ndarray,
    values: np.ndarray,
    *,
    color: str,
    linestyle: str,
    linewidth: float,
    label: str,
) -> None:
    for index in range(values.size):
        axis.plot(
            [left[index], right[index]],
            [values[index], values[index]],
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            label=label if index == 0 else None,
        )


def plot_flash_finite_volume_projection(
    runs: tuple[RunData, RunData],
    cycle: int,
    r_max_km: float,
    output: Path,
) -> None:
    """Compare Thornado element and projected FLASH-X cell averages."""

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(13.0, 8.5),
        sharex="col",
        sharey="row",
        gridspec_kw={"hspace": 0.10, "wspace": 0.08},
    )
    minimum_radius = 0.5 * min(
        float(np.min(run.nodes["r_cm"][run.nodes["r_cm"] > 0.0])) / 1.0e5
        for run in runs
    )
    for column, run in enumerate(runs):
        axes[0, column].set_title(run.label)
        for row, quantity in enumerate(("QH", "QL")):
            radius, element_rate = element_specific_rate(run, quantity)
            within_plot = radius <= r_max_km
            element_left, element_right, element_rate = element_intervals(
                run,
                quantity,
                cycle,
                element_rate[within_plot],
                r_max_km,
            )
            cell_left, cell_right, cell_rate = (
                attribution.project_matter_exchange_to_flash_cell_intervals(
                    run.nodes,
                    run.density,
                    cycle,
                    quantity,
                    "actual",
                    r_max_km * 1.0e5,
                )
            )
            axis = axes[row, column]
            plot_horizontal_intervals(
                axis,
                cell_left,
                cell_right,
                cell_rate,
                color="#0072B2",
                linestyle="-",
                linewidth=1.2,
                label="FLASH-X finite-volume cell average",
            )
            plot_horizontal_intervals(
                axis,
                element_left,
                element_right,
                element_rate,
                color="#D55E00",
                linestyle="--",
                linewidth=1.0,
                label="Thornado element average",
            )
            axis.set_xscale("log")
            axis.axhline(0.0, color="0.45", linewidth=0.65)
            axis.grid(alpha=0.22, which="both")
            axis.legend(fontsize=7.5, loc="best")

    for axis in axes.flat:
        axis.set_xlim(minimum_radius, r_max_km)
    axes[0, 0].set_ylabel(
        r"$\mathrm{d}\epsilon/\mathrm{d}t$ "
        r"[$10^{21}\ \mathrm{erg\,g^{-1}\,s^{-1}}$]"
    )
    axes[1, 0].set_ylabel(r"$\mathrm{d}Y_e/\mathrm{d}t$ [$\mathrm{s^{-1}}$]")
    axes[1, 0].set_xlabel(r"radius [$\mathrm{km}$]")
    axes[1, 1].set_xlabel(r"radius [$\mathrm{km}$]")
    fig.suptitle(
        "Matter exchange after projection to the FLASH-X finite-volume grid\n"
        rf"cycle {cycle}, $r\leq {r_max_km:g}\,\mathrm{{km}}$",
        fontsize=14,
    )
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.09, top=0.90)
    save_figure(fig, output)
    plt.close(fig)


def plot_final_fluid_state_separation(
    runs: tuple[RunData, RunData],
    r_max_km: float,
    output: Path,
) -> None:
    """Plot the final AP-minus-standard matter-state separation."""

    standard, ap = runs
    radius = standard.main["r_cm"] / 1.0e5
    mask = radius <= r_max_km
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(8.0, 7.0),
        sharex=True,
        gridspec_kw={"hspace": 0.0},
    )
    axes[0].plot(
        radius[mask],
        (ap.main["T_fluid_MeV"] - standard.main["T_fluid_MeV"])[mask],
        color="#0072B2",
        linewidth=1.5,
    )
    axes[1].plot(
        radius[mask],
        (ap.main["Ye"] - standard.main["Ye"])[mask],
        color="#0072B2",
        linewidth=1.5,
    )
    for axis in axes:
        axis.set_xscale("log")
        axis.axhline(0.0, color="0.45", linewidth=0.65)
        axis.grid(alpha=0.22, which="both")
    axes[0].set_ylabel(r"$T_{\rm AP}-T_{\rm LLF}$ [MeV]")
    axes[1].set_ylabel(r"$Y_{e,{\rm AP}}-Y_{e,{\rm LLF}}$")
    axes[1].set_xlabel(r"radius [$\mathrm{km}$]")
    fig.suptitle("Final fluid-state separation between the two trajectories")
    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.10, top=0.92)
    save_figure(fig, output)
    plt.close(fig)


def plot_specific_exchange_polynomials(
    runs: tuple[RunData, RunData],
    cycle: int,
    r_max_km: float,
    output: Path,
) -> None:
    """Show the nodal polynomial and element average without joining faces."""

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(12.0, 8.0),
        sharex=True,
        gridspec_kw={"hspace": 0.0},
    )
    for run in runs:
        for axis, quantity in zip(axes, ("QH", "QL"), strict=True):
            rows = attribution.select(
                run.nodes, cycle, quantity, "actual", r_max_km * 1.0e5
            )
            mean_radius, mean_rate = element_specific_rate(run, quantity)
            mean_by_element = {
                int(element): float(value)
                for element, value in zip(
                    np.unique(rows["element"].astype(int)),
                    mean_rate[mean_radius <= r_max_km],
                    strict=True,
                )
            }
            first_polynomial = True
            first_average = True
            for element_value in np.unique(rows["element"]):
                element = int(element_value)
                block = rows[rows["element"] == element_value]
                block = block[np.argsort(block["node"])]
                n_nodes = block.size
                gauss_nodes = np.polynomial.legendre.leggauss(n_nodes)[0]
                radius_nodes = np.asarray(block["r_cm"], dtype=float) / 1.0e5
                center = float(np.mean(radius_nodes))
                half_width = float(
                    np.dot(gauss_nodes, radius_nodes - center)
                    / np.dot(gauss_nodes, gauss_nodes)
                )
                density = np.asarray(
                    [
                        run.density[(element, int(node))]
                        for node in block["node"]
                    ]
                )
                if quantity == "QH":
                    nodal_rate = (
                        block["step_rate"]
                        / density
                        / attribution.SPECIFIC_HEATING_SCALE
                    )
                else:
                    nodal_rate = (
                        attribution.ATOMIC_MASS_UNIT_G
                        * block["step_rate"]
                        / density
                    )
                coefficients = np.polynomial.polynomial.polyfit(
                    gauss_nodes, nodal_rate, n_nodes - 1
                )
                xi = np.linspace(-1.0, 1.0, 24)
                radius_polynomial = center + half_width * xi
                rate_polynomial = np.polynomial.polynomial.polyval(
                    xi, coefficients
                )
                axis.plot(
                    radius_polynomial,
                    rate_polynomial,
                    color=run.color,
                    linewidth=0.75,
                    alpha=0.82,
                    label=(
                        f"{run.label}: DG polynomial"
                        if first_polynomial
                        else None
                    ),
                )
                axis.plot(
                    radius_nodes,
                    nodal_rate,
                    linestyle="none",
                    marker="o",
                    markersize=1.8,
                    color=run.color,
                    alpha=0.90,
                )
                axis.plot(
                    [center - half_width, center + half_width],
                    [mean_by_element[element], mean_by_element[element]],
                    color=run.color,
                    linewidth=1.0,
                    linestyle="--",
                    alpha=0.95,
                    label=(
                        f"{run.label}: element average"
                        if first_average
                        else None
                    ),
                )
                first_polynomial = False
                first_average = False

    axes[0].set_ylabel(
        r"$\mathrm{d}\epsilon/\mathrm{d}t$ "
        r"[$10^{21}\ \mathrm{erg\,g^{-1}\,s^{-1}}$]"
    )
    axes[1].set_ylabel(r"$\mathrm{d}Y_e/\mathrm{d}t$ [$\mathrm{s^{-1}}$]")
    positive_radii = np.concatenate(
        [
            run.nodes["r_cm"][run.nodes["r_cm"] > 0.0] / 1.0e5
            for run in runs
        ]
    )
    minimum_radius = 0.5 * float(np.min(positive_radii))
    axes[1].set_xlabel(r"radius [$\mathrm{km}$]")
    for axis in axes:
        axis.set_xscale("log")
        axis.axhline(0.0, color="0.45", linewidth=0.65)
        axis.set_xlim(minimum_radius, r_max_km)
        axis.grid(alpha=0.22, which="both")
        axis.margins(x=0.0)
    axes[0].tick_params(labelbottom=False)
    axes[0].legend(loc="best", fontsize=8.0, ncol=2)
    fig.suptitle(
        "Matter exchange represented within each Thornado element\n"
        f"cycle {cycle}",
        fontsize=14,
    )
    fig.align_ylabels(axes)
    fig.subplots_adjust(
        left=0.11,
        right=0.985,
        bottom=0.09,
        top=0.91,
        hspace=0.0,
    )
    save_figure(fig, output)
    plt.close(fig)


def plot_cumulative(
    runs: tuple[RunData, RunData],
    inner_r_max_km: float,
    output: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), sharex="col")
    for row, quantity in enumerate(("QH", "QL")):
        for column, region in enumerate(("full", "inner")):
            for run in runs:
                source = run.histories if region == "full" else run.inner_histories
                time, _, _, cumulative_signed, _ = (
                    source[(quantity, "actual")]
                )
                radius_mask = np.ones(run.main["r_cm"].shape, dtype=bool)
                if region == "inner":
                    radius_mask = (
                        run.main["r_cm"] <= inner_r_max_km * 1.0e5
                    )
                final_name = (
                    "eint_after_erg_cm3"
                    if quantity == "QH"
                    else "ne_after_cm3"
                )
                final_total = float(
                    np.sum(
                        run.main[final_name][radius_mask]
                        * run.main["dV_cm3"][radius_mask]
                    )
                )
                initial_total = final_total - float(cumulative_signed[-1])
                axes[row, column].plot(
                    time * 1.0e3,
                    cumulative_signed / initial_total,
                    color=run.color,
                    linewidth=1.5,
                    label=run.label,
                )
            axes[row, column].axhline(0.0, color="0.45", linewidth=0.65)
            axes[row, column].grid(alpha=0.25)
            axes[row, column].set_title(
                "full domain"
                if region == "full"
                else rf"inner ${inner_r_max_km:g}\,\mathrm{{km}}$"
            )
            axes[row, column].legend(fontsize=7.0, loc="best")
        axes[row, 0].set_ylabel(
            r"$\Delta E_{\rm int}/E_{{\rm int},0}$"
            if quantity == "QH"
            else r"$\Delta N_e/N_{e,0}$"
        )
    axes[1, 0].set_xlabel("time [ms]")
    axes[1, 1].set_xlabel("time [ms]")
    fig.suptitle(
        "Signed cumulative fluid change relative to the initial regional total"
    )
    fig.tight_layout()
    save_figure(fig, output)
    plt.close(fig)


def plot_radial_cumulative(
    runs: tuple[RunData, RunData],
    r_max_km: float,
    output: Path,
) -> None:
    """Show where the signed cumulative matter change is accumulated."""

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(8.0, 7.0),
        sharex=True,
        gridspec_kw={"hspace": 0.0},
    )
    for row, quantity in enumerate(("QH", "QL")):
        for run in runs:
            cycle = int(run.elements["cycle"][0])
            rows = attribution.select(
                run.elements, cycle, quantity, "actual", np.inf
            )
            shell_change = (
                rows["cumulative_signed_mean"] * rows["volume_cm3"]
            )
            cumulative = np.cumsum(shell_change)
            final_name = (
                "eint_after_erg_cm3"
                if quantity == "QH"
                else "ne_after_cm3"
            )
            final_total = float(
                np.sum(run.main[final_name] * run.main["dV_cm3"])
            )
            initial_total = final_total - float(cumulative[-1])
            axes[row].plot(
                rows["r_cm"] / 1.0e5,
                cumulative / initial_total,
                color=run.color,
                linewidth=1.5,
                label=run.label,
            )
        axes[row].axhline(0.0, color="0.45", linewidth=0.65)
        axes[row].grid(alpha=0.22, which="both")
        axes[row].legend(fontsize=8.0, loc="best")
        axes[row].set_xlim(0.0, r_max_km)
    axes[0].set_ylabel(r"$\Delta E_{\rm int}(<r)/E_{{\rm int},0}$")
    axes[1].set_ylabel(r"$\Delta N_e(<r)/N_{e,0}$")
    axes[1].set_xlabel(r"radius [$\mathrm{km}$]")
    fig.suptitle("Radial accumulation of the signed fluid change")
    fig.subplots_adjust(left=0.15, right=0.98, bottom=0.10, top=0.92)
    save_figure(fig, output)
    plt.close(fig)


def write_summary(
    runs: tuple[RunData, RunData],
    cycle: int,
    r_max_cm: float,
    output: Path,
) -> None:
    lines = [
        "Self-consistent standard-LLF/source-aware-AP comparison",
        f"cycle: {cycle}",
        "",
    ]
    for run in runs:
        lines.append(run.label)
        for quantity in ("QH", "QL"):
            rows = attribution.select(
                run.elements, cycle, quantity, "actual", r_max_cm
            )
            reference = attribution.element_reference(run.main, rows, quantity)
            magnitude = np.abs(rows["step_rate_mean"])
            active = magnitude > 0.0
            timescale = np.full(magnitude.shape, np.inf)
            timescale[active] = reference[active] / magnitude[active]
            minimum_index = int(np.argmin(timescale))
            minimum_timescale = float(timescale[minimum_index])
            _, element_rate = element_specific_rate(run, quantity)
            _, projected_rate = (
                attribution.project_matter_exchange_to_flash_cells(
                    run.nodes,
                    run.density,
                    cycle,
                    quantity,
                    "actual",
                    r_max_cm,
                )
            )
            history = run.histories[(quantity, "actual")]
            inner_history = run.inner_histories[(quantity, "actual")]
            signed = float(history[3][-1])
            inner_signed = float(inner_history[3][-1])
            lines.extend(
                [
                    f"  {quantity}:",
                    f"    peak element-specific rate: {np.max(np.abs(element_rate)):.8e}",
                    f"    peak projected specific rate: {np.max(np.abs(projected_rate)):.8e}",
                    f"    minimum local timescale [s]: {minimum_timescale:.8e}",
                    f"    radius of minimum timescale [km]: {rows['r_cm'][minimum_index] / 1.0e5:.8e}",
                    f"    cumulative signed, full: {signed:.8e}",
                    f"    cumulative signed, inner: {inner_signed:.8e}",
                ]
            )
        lines.append("")

    standard, ap = runs
    volume = standard.main["dV_cm3"]
    lines.extend(["AP-minus-standard trajectory separation", ""])
    for quantity, after_name in (
        ("QH", "eint_after_erg_cm3"),
        ("QL", "ne_after_cm3"),
    ):
        standard_history = standard.histories[(quantity, "actual")]
        ap_history = ap.histories[(quantity, "actual")]
        cumulative_difference = float(
            ap_history[3][-1] - standard_history[3][-1]
        )
        final_state_difference = float(
            np.sum(
                (ap.main[after_name] - standard.main[after_name]) * volume
            )
        )
        standard_final_total = float(
            np.sum(standard.main[after_name] * volume)
        )
        initial_total = standard_final_total - float(
            standard_history[3][-1]
        )
        closure_scale = max(
            abs(cumulative_difference), abs(final_state_difference), 1.0
        )
        closure = abs(cumulative_difference - final_state_difference) / closure_scale
        lines.extend(
            [
                f"  {quantity}:",
                f"    cumulative difference, full: {cumulative_difference:.8e}",
                f"    integrated final-state difference: {final_state_difference:.8e}",
                f"    cumulative/final-state closure: {closure:.8e}",
                f"    difference / initial fluid total: {final_state_difference / abs(initial_total):.8e}",
            ]
        )
    temperature_difference = (
        ap.main["T_fluid_MeV"] - standard.main["T_fluid_MeV"]
    )
    ye_difference = ap.main["Ye"] - standard.main["Ye"]
    temperature_index = int(np.argmax(np.abs(temperature_difference)))
    ye_index = int(np.argmax(np.abs(ye_difference)))
    lines.extend(
        [
            "",
            "  maximum local final-state differences:",
            f"    |delta T|: {abs(temperature_difference[temperature_index]):.8e} MeV at {standard.main['r_cm'][temperature_index] / 1.0e5:.8e} km",
            f"    |delta Ye|: {abs(ye_difference[ye_index]):.8e} at {standard.main['r_cm'][ye_index] / 1.0e5:.8e} km",
            "",
        ]
    )
    output.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    directory = args.directory.expanduser().resolve()
    standard_prefix = diagnostic_prefix(directory, "standard_no_projection")
    ap_prefix = diagnostic_prefix(directory, "source_aware_no_projection")
    standard_elements = Path(
        str(standard_prefix) + ".exchange_attribution.elements"
    )
    ap_elements = Path(str(ap_prefix) + ".exchange_attribution.elements")
    final_cycle = min(
        attribution.last_data_cycle(standard_elements),
        attribution.last_data_cycle(ap_elements),
    )
    cycle = final_cycle if args.cycle is None else args.cycle
    if cycle > final_cycle:
        raise ValueError(f"cycle {cycle} exceeds common final cycle {final_cycle}")
    r_max_cm = args.r_max_km * 1.0e5
    inner_r_max_cm = min(args.inner_r_max_km, args.r_max_km) * 1.0e5
    runs = (
        load_run(
            standard_prefix,
            "standard LLF",
            "#202020",
            cycle,
            r_max_cm,
            inner_r_max_cm,
        ),
        load_run(
            ap_prefix,
            "source-aware AP LLF",
            "#0072B2",
            cycle,
            r_max_cm,
            inner_r_max_cm,
        ),
    )
    output_prefix = args.output_prefix
    if output_prefix is None:
        output_prefix = directory / "self_consistent_llf_ap"
    output_prefix = output_prefix.expanduser().resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    developer = Path(str(output_prefix) + ".developer_summary.png")
    projection = Path(str(output_prefix) + ".flash_finite_volume.png")
    projection_core = Path(
        str(output_prefix) + ".flash_finite_volume_inner.png"
    )
    cumulative = Path(str(output_prefix) + ".cumulative.png")
    polynomials = Path(str(output_prefix) + ".specific_exchange_polynomials.png")
    state_separation = Path(str(output_prefix) + ".state_separation.png")
    radial_cumulative = Path(str(output_prefix) + ".radial_cumulative.png")
    summary = Path(str(output_prefix) + ".summary.txt")
    plot_developer_summary(runs, cycle, args.r_max_km, developer)
    plot_flash_finite_volume_projection(
        runs, cycle, args.r_max_km, projection
    )
    plot_flash_finite_volume_projection(
        runs,
        cycle,
        min(args.inner_r_max_km, args.r_max_km),
        projection_core,
    )
    plot_specific_exchange_polynomials(
        runs, cycle, args.r_max_km, polynomials
    )
    plot_cumulative(
        runs, min(args.inner_r_max_km, args.r_max_km), cumulative
    )
    plot_final_fluid_state_separation(
        runs, args.r_max_km, state_separation
    )
    plot_radial_cumulative(runs, args.r_max_km, radial_cumulative)
    write_summary(runs, cycle, r_max_cm, summary)
    print(f"wrote {developer} and {developer.with_suffix('.pdf')}")
    print(f"wrote {projection} and {projection.with_suffix('.pdf')}")
    print(
        f"wrote {projection_core} and {projection_core.with_suffix('.pdf')}"
    )
    print(f"wrote {polynomials} and {polynomials.with_suffix('.pdf')}")
    print(f"wrote {cumulative} and {cumulative.with_suffix('.pdf')}")
    print(
        f"wrote {state_separation} and {state_separation.with_suffix('.pdf')}"
    )
    print(
        f"wrote {radial_cumulative} and "
        f"{radial_cumulative.with_suffix('.pdf')}"
    )
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()
