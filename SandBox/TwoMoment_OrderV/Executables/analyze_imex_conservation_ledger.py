#!/usr/bin/env python3

"""Validate and plot the per-timestep IMEX conservation decomposition."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


COMPONENTS = (
    "collision_pair_residual",
    "transport_boundary_residual",
    "matter_remainder",
    "radiation_remainder",
    "conservation_residual",
)

COLORS = {
    "collision_pair_residual": "tab:red",
    "transport_boundary_residual": "tab:blue",
    "matter_remainder": "tab:orange",
    "radiation_remainder": "tab:green",
    "conservation_residual": "black",
}

LABELS = {
    "collision_pair_residual": r"implicit matter+radiation",
    "transport_boundary_residual": r"explicit update$-$boundary tally",
    "matter_remainder": r"matter assembly/limiter remainder",
    "radiation_remainder": r"radiation assembly/limiter remainder",
    "conservation_residual": r"measured full-step residual",
}

METHOD_COLORS = ("0.15", "#0072B2", "#D55E00", "#009E73")


def read_ledger(path: Path) -> np.ndarray:
    data = np.atleast_1d(
        np.genfromtxt(
            path,
            names=True,
            dtype=None,
            encoding="utf-8",
        )
    )
    required = {
        "cycle",
        "time_s",
        "dt_s",
        "quantity",
        "component",
        "step_value",
        "cumulative_value",
    }
    missing = required.difference(data.dtype.names or ())
    if missing:
        raise RuntimeError(f"{path}: missing columns: {', '.join(sorted(missing))}")
    if data.size == 0:
        raise RuntimeError(f"{path}: empty IMEX conservation ledger")
    return data


def series(data: np.ndarray, quantity: str, component: str) -> np.ndarray:
    selected = data[
        (data["quantity"] == quantity) & (data["component"] == component)
    ]
    if selected.size == 0:
        raise RuntimeError(f"missing {quantity}/{component} ledger rows")
    return selected


def has_component(data: np.ndarray, component: str) -> bool:
    return bool(np.any(data["component"] == component))


def validate_ledger(path: Path, data: np.ndarray) -> None:
    for quantity in ("energy", "lepton"):
        closure = series(data, quantity, "algebraic_closure")["step_value"]
        scales = np.zeros_like(closure)
        for component in COMPONENTS:
            scales += np.abs(series(data, quantity, component)["step_value"])
        relative = np.max(np.abs(closure) / np.maximum(scales, np.finfo(float).tiny))
        if relative > 5.0e-12:
            raise RuntimeError(
                f"{path}: {quantity} algebraic closure failed: {relative:.3e}"
            )


def validate_inventory(path: Path, data: np.ndarray, history_path: Path) -> None:
    history = np.atleast_1d(np.genfromtxt(history_path, names=True))[1:]
    cycles = series(data, "energy", "conservation_residual")["cycle"]
    if history.size != cycles.size or not np.array_equal(history["cycle"], cycles):
        raise RuntimeError(f"{path}: ledger and inventory histories have different cycles")

    checks = (
        (
            "energy",
            "step_energy_residual_erg",
            "cumulative_energy_residual_erg",
            "total_energy_erg",
        ),
        (
            "lepton",
            "step_lepton_residual",
            "cumulative_lepton_residual",
            "total_lepton",
        ),
    )
    for quantity, step_field, cumulative_field, total_field in checks:
        residual = series(data, quantity, "conservation_residual")
        inventory_scale = np.max(np.abs(history[total_field]))
        step_error = np.max(
            np.abs(residual["step_value"] - history[step_field])
        ) / inventory_scale
        cumulative_error = np.max(
            np.abs(residual["cumulative_value"] - history[cumulative_field])
        ) / inventory_scale
        if step_error > 32.0 * np.finfo(float).eps or cumulative_error > 32.0 * np.finfo(float).eps:
            raise RuntimeError(
                f"{path}: {quantity} ledger disagrees with inventory history: "
                f"step={step_error:.3e}, cumulative={cumulative_error:.3e}"
            )


def read_inventory_history(path: Path) -> np.ndarray:
    history = np.atleast_1d(np.genfromtxt(path, names=True))
    if history.size < 2:
        raise RuntimeError(f"{path}: inventory history has no evolved state")
    return history


def safe_ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    result = np.full_like(numerator, np.nan, dtype=float)
    np.divide(
        np.abs(numerator),
        np.abs(denominator),
        out=result,
        where=np.abs(denominator) > 0.0,
    )
    return result


def endpoint_label(axis, x: np.ndarray, y: np.ndarray, color: str) -> None:
    finite = np.flatnonzero(np.isfinite(y) & (y > 0.0))
    if finite.size == 0:
        return
    index = finite[-1]
    axis.annotate(
        rf"${y[index]:.2e}$",
        xy=(x[index], y[index]),
        xytext=(-5, 6),
        textcoords="offset points",
        ha="right",
        va="bottom",
        color=color,
        fontsize=8,
    )


def plot_conservation_summary(
    ledgers: list[np.ndarray],
    histories: list[np.ndarray],
    labels: list[str],
) -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4), sharex=True)
    quantity_metadata = (
        ("energy", "total_energy_erg", "Energy"),
        ("lepton", "total_lepton", "Electron lepton number"),
    )

    for row, (quantity, total_field, title) in enumerate(quantity_metadata):
        for index, (data, history, label) in enumerate(
            zip(ledgers, histories, labels)
        ):
            color = METHOD_COLORS[index % len(METHOD_COLORS)]
            residual = series(data, quantity, "conservation_residual")
            time_ms = residual["time_s"] * 1.0e3
            initial_total = abs(history[total_field][0])
            relative_residual = np.abs(residual["cumulative_value"]) / initial_total

            implicit_matter = series(
                data, quantity, "implicit_matter"
            )["cumulative_value"]
            implicit_radiation = series(
                data, quantity, "implicit_radiation"
            )["cumulative_value"]
            pairing_fraction = safe_ratio(
                implicit_matter + implicit_radiation,
                np.abs(implicit_matter) + np.abs(implicit_radiation),
            )

            axes[row, 0].plot(
                time_ms,
                relative_residual,
                color=color,
                linewidth=1.8,
                label=label,
            )
            axes[row, 1].plot(
                time_ms,
                pairing_fraction,
                color=color,
                linewidth=1.8,
                label=label,
            )
            endpoint_label(axes[row, 0], time_ms, relative_residual, color)
            endpoint_label(axes[row, 1], time_ms, pairing_fraction, color)

        axes[row, 0].set_yscale("log")
        axes[row, 1].set_yscale("log")
        axes[row, 0].set_title(f"{title}: total conservation error")
        axes[row, 1].set_title(f"{title}: implicit-pair mismatch")
        axes[row, 0].set_ylabel("fraction of initial total")
        axes[row, 1].set_ylabel("fraction of implicit exchange")

    axes[1, 0].set_xlabel(r"$t\,[\mathrm{ms}]$")
    axes[1, 1].set_xlabel(r"$t\,[\mathrm{ms}]$")
    axes[0, 0].legend(frameon=False, loc="upper left")
    fig.suptitle(
        "Conservation history and cancellation of the implicit "
        "matter--radiation exchange",
        y=0.985,
    )
    fig.text(
        0.5,
        0.012,
        r"Implicit-pair mismatch "
        r"$=|\Delta Q_{\rm matter}^{\rm I}+\Delta Q_{\rm rad}^{\rm I}|/"
        r"(|\Delta Q_{\rm matter}^{\rm I}|+|\Delta Q_{\rm rad}^{\rm I}|)$; "
        "exact pairing is zero.",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0.0, 0.045, 1.0, 0.95))
    return fig


def plot_implicit_pairing(
    ledgers: list[np.ndarray], labels: list[str]
) -> plt.Figure:
    n_methods = len(ledgers)
    fig, axes = plt.subplots(
        2,
        n_methods,
        figsize=(5.0 * n_methods, 7.2),
        sharex=True,
        squeeze=False,
    )
    quantity_metadata = (("energy", "Energy", "erg"), ("lepton", "Electron lepton number", ""))

    for column, (data, method_label) in enumerate(zip(ledgers, labels)):
        for row, (quantity, quantity_label, unit) in enumerate(quantity_metadata):
            matter_rows = series(data, quantity, "implicit_matter")
            radiation_rows = series(data, quantity, "implicit_radiation")
            time_ms = matter_rows["time_s"] * 1.0e3
            matter = matter_rows["cumulative_value"]
            opposite_radiation = -radiation_rows["cumulative_value"]
            scale = max(
                np.max(np.abs(matter)),
                np.max(np.abs(opposite_radiation)),
                np.finfo(float).tiny,
            )

            axis = axes[row, column]
            axis.plot(
                time_ms,
                matter / scale,
                color="#D55E00",
                linewidth=1.7,
                label="matter",
            )
            axis.plot(
                time_ms,
                opposite_radiation / scale,
                color="#0072B2",
                linewidth=1.7,
                linestyle="--",
                label="opposite radiation",
            )
            axis.axhline(0.0, color="0.6", linewidth=0.7)
            axis.set_title(f"{method_label}: {quantity_label}")
            axis.set_ylabel("cumulative implicit change / panel scale")
            scale_text = rf"scale $={scale:.3e}$"
            if unit:
                scale_text += rf" $\mathrm{{{unit}}}$"
            axis.text(
                0.03,
                0.05,
                scale_text,
                transform=axis.transAxes,
                fontsize=8,
                ha="left",
                va="bottom",
            )
            if row == 1:
                axis.set_xlabel(r"$t\,[\mathrm{ms}]$")

    axes[0, 0].legend(frameon=False, loc="best")
    fig.suptitle(
        "Cumulative implicit exchange: conservative pairing requires "
        "the two curves to coincide",
        y=0.985,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    return fig


def plot_residual_attribution(
    ledgers: list[np.ndarray], labels: list[str]
) -> plt.Figure:
    n_methods = len(ledgers)
    fig, axes = plt.subplots(
        2,
        n_methods,
        figsize=(5.0 * n_methods, 7.2),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    plotted_components = (
        "collision_pair_residual",
        "transport_boundary_residual",
        "matter_remainder",
        "radiation_remainder",
    )

    for column, (data, method_label) in enumerate(zip(ledgers, labels)):
        for row, (quantity, quantity_label) in enumerate(
            (("energy", "Energy"), ("lepton", "Electron lepton number"))
        ):
            full = series(data, quantity, "conservation_residual")
            time_ms = full["time_s"] * 1.0e3
            denominator = full["cumulative_value"]
            axis = axes[row, column]
            for component in plotted_components:
                values = series(data, quantity, component)["cumulative_value"]
                ratio = safe_ratio(values, denominator)
                if not np.any(np.isfinite(ratio) & (ratio > 0.0)):
                    continue
                axis.plot(
                    time_ms,
                    ratio,
                    color=COLORS[component],
                    linewidth=1.5,
                    label=LABELS[component],
                )
            axis.axhline(1.0, color="0.5", linewidth=0.7, linestyle=":")
            axis.set_yscale("log")
            axis.set_ylim(1.0e-13, 3.0)
            axis.set_title(f"{method_label}: {quantity_label}")
            axis.set_ylabel("absolute fraction of measured residual")
            if row == 1:
                axis.set_xlabel(r"$t\,[\mathrm{ms}]$")

    handles, legend_labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            legend_labels,
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.94),
            ncol=2,
            fontsize=8,
        )
    fig.suptitle(
        "Origin of the measured conservation residual",
        y=0.985,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.88))
    return fig


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("ledger", type=Path, nargs="+")
    parser.add_argument("--labels", nargs="+")
    parser.add_argument("--inventory-history", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.labels is not None and len(args.labels) != len(args.ledger):
        raise RuntimeError("--labels must provide one label per ledger")
    if args.inventory_history is not None and len(args.inventory_history) != len(
        args.ledger
    ):
        raise RuntimeError(
            "--inventory-history must provide one history per ledger"
        )

    labels = args.labels or [path.stem for path in args.ledger]
    ledgers = [read_ledger(path) for path in args.ledger]
    histories = None
    if args.inventory_history is not None:
        histories = [read_inventory_history(path) for path in args.inventory_history]
    for index, (path, data) in enumerate(zip(args.ledger, ledgers)):
        validate_ledger(path, data)
        if args.inventory_history is not None:
            validate_inventory(path, data, args.inventory_history[index])

        print(labels[index])
        for quantity in ("energy", "lepton"):
            print(f"  {quantity} final cumulative decomposition:")
            for component in COMPONENTS:
                value = series(data, quantity, component)["cumulative_value"][-1]
                unit = "erg" if quantity == "energy" else ""
                print(f"    {component:29s} {value:+.8e} {unit}")
            if has_component(data, "direct_implicit_pair_residual"):
                direct = series(
                    data, quantity, "direct_implicit_pair_residual"
                )["cumulative_value"][-1]
                collision = series(
                    data, quantity, "collision_pair_residual"
                )["cumulative_value"][-1]
                scale = max(abs(direct), abs(collision), np.finfo(float).tiny)
                print(
                    "    direct accepted-pair check    "
                    f"{direct:+.8e} {unit}; relative mismatch="
                    f"{abs(direct - collision) / scale:.3e}"
                )

    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "cm",
            "axes.grid": True,
            "grid.alpha": 0.25,
        }
    )
    figures = []
    if histories is not None:
        figures.append(plot_conservation_summary(ledgers, histories, labels))
    figures.append(plot_implicit_pairing(ledgers, labels))
    figures.append(plot_residual_attribution(ledgers, labels))

    if args.output.suffix.lower() == ".pdf":
        with PdfPages(args.output) as pdf:
            for figure in figures:
                pdf.savefig(figure)
        page_names = (
            ("summary", "implicit-pairing", "residual-attribution")
            if histories is not None
            else ("implicit-pairing", "residual-attribution")
        )
        for page_name, figure in zip(page_names, figures):
            page_path = args.output.with_name(
                f"{args.output.stem}.{page_name}.pdf"
            )
            figure.savefig(page_path)
    else:
        figures[0].savefig(args.output, dpi=250)
        for page, figure in enumerate(figures[1:], start=2):
            companion = args.output.with_name(
                f"{args.output.stem}.page{page}{args.output.suffix}"
            )
            figure.savefig(companion, dpi=250)
    for figure in figures:
        plt.close(figure)
    print(f"figure: {args.output}")


if __name__ == "__main__":
    main()
