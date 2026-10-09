#!/usr/bin/env python3
"""Quantify matter exchange caused by the standard and source-aware AP fluxes.

The analysis uses the paired shadow steps written by the standalone driver.  It
also decomposes the stage-weighted EmAb exchange into the element-constant and
nonconstant spatial parts of D-D0.  For a p=1 spatial discretization, the latter
is the complete linear spatial mode; it does not by itself distinguish a
backward slope from a physically oriented slope.
"""

from __future__ import annotations

import argparse
import io
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


matplotlib.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "mathtext.fontset": "cm",
    }
)


ATOMIC_MASS_UNIT_G = 1.66053906660e-24
SPECIFIC_HEATING_SCALE = 1.0e21
SPEED_OF_LIGHT_MKS = 2.99792458e8
GRAVITATIONAL_CONSTANT_MKS = 6.673e-11
CENTIMETER_CODE = 1.0e-2
GRAM_CODE = (
    1.0e-3 * GRAVITATIONAL_CONSTANT_MKS / SPEED_OF_LIGHT_MKS**2
)
MASS_DENSITY_CGS_UNIT_CODE = GRAM_CODE / CENTIMETER_CODE**3


COMPONENTS = (
    "actual",
    "no_stream_shadow",
    "ap_minus_no_stream",
    "standard_minus_ap",
    "standard_total",
    "ap_total",
    "emab_full",
    "emab_mean",
    "emab_high",
)

REPORT_COMPONENTS = (
    "standard_total",
    "ap_total",
    "no_stream_shadow",
    "ap_minus_no_stream",
    "emab_full",
    "emab_mean",
    "emab_high",
)

LABELS = {
    "actual": "production step",
    "standard_total": "standard LLF",
    "ap_total": "source-aware AP LLF",
    "no_stream_shadow": "no-streaming collision response",
    "ap_minus_no_stream": "AP transport response",
    "standard_minus_ap": r"standard $-$ AP shadow-step difference",
    "emab_full": "full EmAb exchange",
    "emab_mean": r"EmAb: element-constant mode of $D-D_0$",
    "emab_high": r"EmAb nonconstant spatial-mode ($p=1$) response",
}

COLORS = {
    "standard_total": "#202020",
    "ap_total": "#0072B2",
    "no_stream_shadow": "#009E73",
    "ap_minus_no_stream": "#56B4E9",
    "standard_minus_ap": "#D55E00",
    "emab_full": "#202020",
    "emab_mean": "#009E73",
    "emab_high": "#CC79A7",
}

QUANTITY = {
    "QH": {
        "rate_unit": r"erg cm$^{-3}$ s$^{-1}$",
        "global_rate_unit": r"erg s$^{-1}$",
        "increment_unit": "erg",
        "reference": "matter internal energy",
    },
    "QL": {
        "rate_unit": r"cm$^{-3}$ s$^{-1}$",
        "global_rate_unit": r"s$^{-1}$",
        "increment_unit": "leptons",
        "reference": "matter electron number",
    },
}


@dataclass
class Statistics:
    peak_rate: float
    signed_global_rate: float
    element_l1_rate: float
    nodal_l1_rate: float
    min_local_timescale: float
    global_element_l1_timescale: float
    global_nodal_l1_timescale: float
    nodal_cancellation_fraction: float
    cumulative_signed: float
    cumulative_nodal_absolute: float
    cumulative_absolute_fraction: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "prefix",
        type=Path,
        help=(
            "MatterDiagnosticsFileName, the *.exchange_attribution prefix, "
            "or either complete attribution filename"
        ),
    )
    parser.add_argument(
        "--cycle",
        type=int,
        default=None,
        help="cycle to analyze (default: final recorded cycle)",
    )
    parser.add_argument(
        "--r-max-km",
        type=float,
        default=200.0,
        help="outer radius of plots and integrals (default: 200 km)",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=None,
        help="output prefix (default: next to the attribution files)",
    )
    return parser.parse_args()


def attribution_prefix(path: Path) -> Path:
    text = str(path)
    for suffix in (
        ".exchange_attribution.elements",
        ".exchange_attribution.nodes",
    ):
        if text.endswith(suffix):
            return Path(text[: -len(suffix)] + ".exchange_attribution")
    if text.endswith(".exchange_attribution"):
        return path
    return Path(text + ".exchange_attribution")


def main_diagnostic_path(prefix: Path) -> Path:
    text = str(prefix)
    suffix = ".exchange_attribution"
    if not text.endswith(suffix):
        raise ValueError(f"cannot infer the main diagnostic path from {prefix}")
    return Path(text[: -len(suffix)])


def load_table(path: Path) -> np.ndarray:
    data = np.genfromtxt(path, names=True, dtype=None, encoding="utf-8")
    if data.size == 0:
        raise RuntimeError(f"no rows found in {path}")
    return np.atleast_1d(data)


def last_data_cycle(path: Path) -> int:
    """Read the cycle number from the last data record without scanning."""

    with path.open("rb") as stream:
        stream.seek(0, 2)
        end = stream.tell()
        stream.seek(max(0, end - 65536))
        lines = stream.read().splitlines()
    for line in reversed(lines):
        stripped = line.strip()
        if stripped and not stripped.startswith(b"#"):
            return int(stripped.split(None, 1)[0])
    raise RuntimeError(f"no data rows found in {path}")


def load_element_cycle_and_histories(
    path: Path,
    selected_cycle: int,
    r_max_cm: float,
    inner_r_max_cm: float,
) -> tuple[
    np.ndarray,
    dict[tuple[str, str], tuple[np.ndarray, ...]],
    dict[tuple[str, str], tuple[np.ndarray, ...]],
]:
    """Stream element history once, retaining full and inner histories."""

    names: list[str] | None = None
    indices: dict[str, int] = {}
    selected_rows: list[str] = []
    records: dict[tuple[str, str], list[tuple[float, ...]]] = {}
    inner_records: dict[tuple[str, str], list[tuple[float, ...]]] = {}
    current_cycle: int | None = None
    current: dict[tuple[str, str], list[float]] = {}
    inner_current: dict[tuple[str, str], list[float]] = {}

    def flush() -> None:
        if current_cycle is None:
            return
        for key, values in current.items():
            records.setdefault(key, []).append(tuple(values))
        for key, values in inner_current.items():
            inner_records.setdefault(key, []).append(tuple(values))

    def accumulate(
        collection: dict[tuple[str, str], list[float]],
        key: tuple[str, str],
        words: list[str],
    ) -> None:
        if key not in collection:
            collection[key] = [
                float(words[indices["time_s"]]),
                0.0,
                0.0,
                0.0,
                0.0,
            ]
        volume = float(words[indices["volume_cm3"]])
        dt = float(words[indices["dt_s"]])
        values = collection[key]
        values[1] += float(words[indices["step_rate_mean"]]) * volume
        values[2] += (
            float(words[indices["step_abs_node_mean"]]) / dt * volume
        )
        values[3] += (
            float(words[indices["cumulative_signed_mean"]]) * volume
        )
        values[4] += (
            float(words[indices["cumulative_absolute_node_mean"]])
            * volume
        )

    with path.open() as stream:
        for line in stream:
            if line.startswith("# cycle"):
                names = line[1:].split()
                indices = {name: index for index, name in enumerate(names)}
                continue
            if line.startswith("#") or not line.strip():
                continue
            if names is None:
                raise ValueError(f"{path}: data precede the column header")
            words = line.split()
            cycle = int(words[indices["cycle"]])
            if cycle > selected_cycle:
                break
            if current_cycle != cycle:
                flush()
                current_cycle = cycle
                current = {}
                inner_current = {}
            if cycle == selected_cycle:
                selected_rows.append(line)
            radius = float(words[indices["r_cm"]])
            if radius > r_max_cm:
                continue
            quantity = words[indices["quantity"]]
            component = words[indices["component"]]
            key = (quantity, component)
            accumulate(current, key, words)
            if radius <= inner_r_max_cm:
                accumulate(inner_current, key, words)
    flush()
    if names is None:
        raise ValueError(f"{path}: missing column header")
    if not selected_rows:
        raise ValueError(f"{path}: no rows for cycle {selected_cycle}")
    elements = np.atleast_1d(
        np.genfromtxt(
            io.StringIO("".join(selected_rows)),
            names=names,
            dtype=None,
            encoding="utf-8",
        )
    )
    histories = {
        key: tuple(np.asarray(column) for column in zip(*values, strict=True))
        for key, values in records.items()
    }
    inner_histories = {
        key: tuple(np.asarray(column) for column in zip(*values, strict=True))
        for key, values in inner_records.items()
    }
    return elements, histories, inner_histories


def load_cycle_table(path: Path, cycle: int) -> np.ndarray:
    """Read one cycle without retaining the full nodal history in memory."""

    names: list[str] | None = None
    rows: list[str] = []
    with path.open() as stream:
        for line in stream:
            if line.startswith("# cycle"):
                names = line[1:].split()
                continue
            if line.startswith("#") or not line.strip():
                continue
            row_cycle = int(line.split(None, 1)[0])
            if row_cycle == cycle:
                rows.append(line)
            elif rows and row_cycle > cycle:
                break
    if names is None:
        raise ValueError(f"{path}: missing column header")
    if not rows:
        raise ValueError(f"{path}: no rows for cycle {cycle}")
    return np.atleast_1d(
        np.genfromtxt(
            io.StringIO("".join(rows)),
            names=names,
            dtype=None,
            encoding="utf-8",
        )
    )


def load_main_diagnostic(path: Path) -> dict[str, np.ndarray]:
    names: list[str] | None = None
    with path.open() as stream:
        for line in stream:
            if line.startswith("# r_cm"):
                names = line[1:].split()
                break
    if names is None:
        raise ValueError(f"{path}: missing column header")
    raw = np.atleast_2d(np.loadtxt(path, comments="#"))
    if raw.shape[1] != len(names):
        raise ValueError(
            f"{path}: found {raw.shape[1]} columns, expected {len(names)}"
        )
    return dict(zip(names, raw.T, strict=True))


def load_nodal_density(
    frozen_path: Path,
) -> dict[tuple[int, int], float]:
    """Recover nodal mass density from the saved Thornado state."""

    names: list[str] | None = None
    indices: dict[str, int] = {}
    density_code: dict[tuple[int, int], float] = {}
    density_name: str | None = None
    with frozen_path.open() as stream:
        for line in stream:
            if line.startswith("# element spatial_node") or line.startswith(
                "# stage state element spatial_node"
            ):
                names = line[1:].split()
                indices = {name: index for index, name in enumerate(names)}
                density_name = (
                    "fluid_D" if "fluid_D" in indices else "rho_code"
                )
                continue
            if line.startswith("#") or not line.strip():
                continue
            if names is None:
                raise ValueError(f"{frozen_path}: data precede the column header")
            words = line.split()
            if "state" in indices and int(words[indices["state"]]) != 4:
                continue
            if (
                int(words[indices["energy_element"]]) != 1
                or int(words[indices["energy_node"]]) != 1
                or int(words[indices["species"]]) != 1
            ):
                continue
            key = (
                int(words[indices["element"]]),
                int(words[indices["spatial_node"]]),
            )
            if density_name is None:
                raise ValueError(f"{frozen_path}: missing density column")
            density_code[key] = float(words[indices[density_name]])
    if not density_code:
        raise RuntimeError(f"no nodal density found in {frozen_path}")

    return {
        key: value / MASS_DENSITY_CGS_UNIT_CODE
        for key, value in density_code.items()
    }


def select(
    data: np.ndarray,
    cycle: int,
    quantity: str,
    component: str,
    r_max_cm: float,
) -> np.ndarray:
    mask = (
        (data["cycle"] == cycle)
        & (data["quantity"] == quantity)
        & (data["component"] == component)
        & (data["r_cm"] <= r_max_cm)
    )
    rows = data[mask]
    if rows.size == 0:
        raise RuntimeError(
            f"missing cycle={cycle}, quantity={quantity}, component={component}"
        )
    return rows[np.argsort(rows["r_cm"])]


def element_reference(
    main: dict[str, np.ndarray], rows: np.ndarray, quantity: str
) -> np.ndarray:
    mask = main["r_cm"] <= float(np.max(rows["r_cm"])) + 1.0e-12
    radius = main["r_cm"][mask]
    if radius.size != rows.size or not np.allclose(
        radius, rows["r_cm"], rtol=0.0, atol=1.0e-10
    ):
        raise ValueError("matter reference and attribution element grids differ")
    if quantity == "QH":
        return 0.5 * (
            main["eint_before_erg_cm3"][mask]
            + main["eint_after_erg_cm3"][mask]
        )
    return 0.5 * (
        main["ne_before_cm3"][mask] + main["ne_after_cm3"][mask]
    )


def positive(values: np.ndarray) -> np.ndarray:
    return np.where(values > 0.0, values, np.nan)


def project_matter_exchange_to_flash_cell_intervals(
    nodes: np.ndarray,
    nodal_density: dict[tuple[int, int], float],
    cycle: int,
    quantity: str,
    component: str,
    r_max_cm: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project a volumetric matter-exchange rate to spherical FLASH cells."""

    if quantity not in ("QH", "QL"):
        raise ValueError(f"unsupported projected quantity: {quantity}")

    # Select complete elements before applying the radial plotting limit.  A
    # cutoff through an element would otherwise leave an incomplete nodal set
    # and corrupt the element-to-FV projection used only for visualization.
    rows = select(nodes, cycle, quantity, component, np.inf)
    cell_left: list[float] = []
    cell_right: list[float] = []
    specific_rate: list[float] = []
    for element_value in np.unique(rows["element"]):
        element = int(element_value)
        block = rows[rows["element"] == element_value]
        block = block[np.argsort(block["node"])]
        if float(np.mean(block["r_cm"])) > r_max_cm:
            continue
        n_nodes = block.size
        gauss_nodes, _ = np.polynomial.legendre.leggauss(n_nodes)
        radius_nodes = np.asarray(block["r_cm"], dtype=float)
        center = float(np.mean(radius_nodes))
        half_width = float(
            np.dot(gauss_nodes, radius_nodes - center)
            / np.dot(gauss_nodes, gauss_nodes)
        )
        if half_width <= 0.0:
            raise ValueError(f"element {element}: invalid spatial width")

        basis: list[np.polynomial.Polynomial] = []
        for j, x_j in enumerate(gauss_nodes):
            polynomial = np.polynomial.Polynomial([1.0])
            denominator = 1.0
            for k, x_k in enumerate(gauss_nodes):
                if k == j:
                    continue
                polynomial *= np.polynomial.Polynomial([-x_k, 1.0])
                denominator *= x_j - x_k
            basis.append(polynomial / denominator)

        radius_polynomial = np.polynomial.Polynomial(
            [center, half_width]
        )
        volume_weight = radius_polynomial * radius_polynomial
        cell_edges = np.linspace(-1.0, 1.0, n_nodes + 1)
        projection = np.empty((n_nodes, n_nodes))
        for i_cell in range(n_nodes):
            x_left = cell_edges[i_cell]
            x_right = cell_edges[i_cell + 1]
            weight_integral = volume_weight.integ()
            denominator = (
                weight_integral(x_right) - weight_integral(x_left)
            )
            for j, polynomial in enumerate(basis):
                numerator_integral = (polynomial * volume_weight).integ()
                projection[i_cell, j] = (
                    numerator_integral(x_right)
                    - numerator_integral(x_left)
                ) / denominator

        density_nodes = np.asarray(
            [
                nodal_density[(element, int(node))]
                for node in block["node"]
            ]
        )
        projected_rate = projection @ np.asarray(
            block["step_rate"], dtype=float
        )
        projected_density = projection @ density_nodes
        gauss_weights = np.polynomial.legendre.leggauss(n_nodes)[1]
        expected_nodal_volume = (
            4.0
            * np.pi
            * half_width
            * gauss_weights
            * radius_nodes**2
        )
        recorded_nodal_volume = np.asarray(
            block["volume_weight_cm3"], dtype=float
        )
        geometry_scale = max(
            float(np.max(np.abs(expected_nodal_volume))), 1.0e-300
        )
        if (
            np.max(
                np.abs(recorded_nodal_volume - expected_nodal_volume)
            )
            > 1.0e-11 * geometry_scale
        ):
            raise RuntimeError(
                "the attribution file contains incorrectly indexed spatial "
                "quadrature weights; rebuild with the corrected "
                "WriteMatterExchangeAttributionStep interface and rerun"
            )
        cell_volume = (
            4.0
            * np.pi
            / 3.0
            * (
                (center + half_width * cell_edges[1:]) ** 3
                - (center + half_width * cell_edges[:-1]) ** 3
            )
        )
        for nodal_values, projected_values, name in (
            (
                np.asarray(block["step_rate"], dtype=float),
                projected_rate,
                quantity,
            ),
            (density_nodes, projected_density, "density"),
        ):
            nodal_integral = float(
                np.sum(nodal_values * expected_nodal_volume)
            )
            projected_integral = float(
                np.sum(projected_values * cell_volume)
            )
            scale = max(
                abs(nodal_integral),
                abs(projected_integral),
                1.0e-300,
            )
            if abs(nodal_integral - projected_integral) > 1.0e-11 * scale:
                raise RuntimeError(
                    f"element {element}: nonconservative {name} projection"
                )
        for i_cell in range(n_nodes):
            cell_left.append(
                center + half_width * cell_edges[i_cell]
            )
            cell_right.append(
                center + half_width * cell_edges[i_cell + 1]
            )
            if quantity == "QH":
                specific_rate.append(
                    projected_rate[i_cell]
                    / projected_density[i_cell]
                    / SPECIFIC_HEATING_SCALE
                )
            else:
                specific_rate.append(
                    ATOMIC_MASS_UNIT_G
                    * projected_rate[i_cell]
                    / projected_density[i_cell]
                )

    order = np.argsort(cell_left)
    return (
        np.asarray(cell_left)[order] / 1.0e5,
        np.asarray(cell_right)[order] / 1.0e5,
        np.asarray(specific_rate)[order],
    )


def project_matter_exchange_to_flash_cells(
    nodes: np.ndarray,
    nodal_density: dict[tuple[int, int], float],
    cycle: int,
    quantity: str,
    component: str,
    r_max_cm: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return centers and averages of spherical FLASH finite-volume cells."""

    cell_left, cell_right, specific_rate = (
        project_matter_exchange_to_flash_cell_intervals(
            nodes,
            nodal_density,
            cycle,
            quantity,
            component,
            r_max_cm,
        )
    )
    return 0.5 * (cell_left + cell_right), specific_rate


def relative_closure(lhs: np.ndarray, rhs: np.ndarray) -> float:
    scale = max(float(np.max(np.abs(lhs))), float(np.max(np.abs(rhs))), 1.0e-300)
    return float(np.max(np.abs(lhs - rhs))) / scale


def correlation(lhs: np.ndarray, rhs: np.ndarray) -> float:
    lhs = lhs - np.mean(lhs)
    rhs = rhs - np.mean(rhs)
    denominator = np.linalg.norm(lhs) * np.linalg.norm(rhs)
    if denominator == 0.0:
        return float("nan")
    return float(np.dot(lhs, rhs) / denominator)


def component_statistics(rows: np.ndarray, reference: np.ndarray) -> Statistics:
    rate = rows["step_rate_mean"]
    volume = rows["volume_cm3"]
    dt = rows["dt_s"]
    nodal_absolute_rate = rows["step_abs_node_mean"] / dt
    element_l1_rate = float(np.sum(np.abs(rate) * volume))
    nodal_l1_rate = float(np.sum(nodal_absolute_rate * volume))
    reference_total = float(np.sum(reference * volume))
    local = np.full(rate.shape, np.inf)
    active = np.abs(rate) > 0.0
    local[active] = reference[active] / np.abs(rate[active])
    cancellation = 0.0
    if nodal_l1_rate > 0.0:
        cancellation = 1.0 - element_l1_rate / nodal_l1_rate
    cumulative_signed = float(
        np.sum(rows["cumulative_signed_mean"] * volume)
    )
    cumulative_absolute = float(
        np.sum(rows["cumulative_absolute_node_mean"] * volume)
    )
    return Statistics(
        peak_rate=float(np.max(np.abs(rate))),
        signed_global_rate=float(np.sum(rate * volume)),
        element_l1_rate=element_l1_rate,
        nodal_l1_rate=nodal_l1_rate,
        min_local_timescale=float(np.min(local)),
        global_element_l1_timescale=(
            reference_total / element_l1_rate
            if element_l1_rate > 0.0
            else float("inf")
        ),
        global_nodal_l1_timescale=(
            reference_total / nodal_l1_rate
            if nodal_l1_rate > 0.0
            else float("inf")
        ),
        nodal_cancellation_fraction=cancellation,
        cumulative_signed=cumulative_signed,
        cumulative_nodal_absolute=cumulative_absolute,
        cumulative_absolute_fraction=(
            cumulative_absolute / reference_total
            if reference_total > 0.0
            else float("nan")
        ),
    )


def report(
    elements: np.ndarray,
    main: dict[str, np.ndarray],
    cycle: int,
    r_max_cm: float,
) -> str:
    lines = [
        f"Matter-exchange attribution: cycle {cycle}, r <= {r_max_cm / 1.0e5:g} km",
        "",
        "Exact decomposition:",
        "  standard total = no-streaming response + AP transport response",
        "                   + (standard - AP) shadow-step difference",
        "  full EmAb       = element-constant mode + nonconstant spatial modes",
        "",
    ]
    for quantity in ("QH", "QL"):
        base = select(
            elements, cycle, quantity, "standard_total", r_max_cm
        )
        reference = element_reference(main, base, quantity)
        arrays = {
            component: select(
                elements, cycle, quantity, component, r_max_cm
            )["step_rate_mean"]
            for component in COMPONENTS
        }
        standard_sum = (
            arrays["no_stream_shadow"]
            + arrays["ap_minus_no_stream"]
            + arrays["standard_minus_ap"]
        )
        ap_sum = arrays["no_stream_shadow"] + arrays["ap_minus_no_stream"]
        modal_sum = arrays["emab_mean"] + arrays["emab_high"]
        lines.extend(
            (
                f"{quantity} ({QUANTITY[quantity]['reference']} reference):",
                "  algebraic closures:",
                "    standard shadow = "
                f"{relative_closure(arrays['standard_total'], standard_sum):.6e}",
                "    AP shadow       = "
                f"{relative_closure(arrays['ap_total'], ap_sum):.6e}",
                "    EmAb modes      = "
                f"{relative_closure(arrays['emab_full'], modal_sum):.6e}",
                "  attribution:",
                "    corr(standard-AP, full EmAb) = "
                f"{correlation(arrays['standard_minus_ap'], arrays['emab_full']):.6e}",
                "    corr(standard-AP, nonconstant EmAb modes) = "
                f"{correlation(arrays['standard_minus_ap'], arrays['emab_high']):.6e}",
            )
        )
        statistics: dict[str, Statistics] = {}
        for component in REPORT_COMPONENTS:
            rows = select(elements, cycle, quantity, component, r_max_cm)
            statistics[component] = component_statistics(rows, reference)
        ap_transport = statistics["ap_minus_no_stream"]
        removed = component_statistics(
            select(
                elements, cycle, quantity, "standard_minus_ap", r_max_cm
            ),
            reference,
        )
        ap_total = statistics["ap_total"]
        lines.extend(
            (
                "    peak |standard-AP| / peak |AP transport| = "
                f"{removed.peak_rate / max(ap_transport.peak_rate, 1.0e-300):.6e}",
                "    element-L1 |standard-AP| / |AP transport| = "
                f"{removed.element_l1_rate / max(ap_transport.element_l1_rate, 1.0e-300):.6e}",
                "    peak |standard-AP| / peak |AP total| = "
                f"{removed.peak_rate / max(ap_total.peak_rate, 1.0e-300):.6e}",
                "    element-L1 |standard-AP| / |AP total| = "
                f"{removed.element_l1_rate / max(ap_total.element_l1_rate, 1.0e-300):.6e}",
                "    |global signed standard-AP| / |AP total| = "
                f"{abs(removed.signed_global_rate) / max(abs(ap_total.signed_global_rate), 1.0e-300):.6e}",
                "    cumulative nodal-absolute standard-AP / AP total = "
                f"{removed.cumulative_nodal_absolute / max(ap_total.cumulative_nodal_absolute, 1.0e-300):.6e}",
                "",
                "  final-step rates and timescales:",
            )
        )
        for component in REPORT_COMPONENTS:
            stats = statistics[component]
            lines.extend(
                (
                    f"    {LABELS[component]}:",
                    f"      peak local |rate|       = {stats.peak_rate:.6e} "
                    f"{QUANTITY[quantity]['rate_unit']}",
                    f"      volume-integrated rate = {stats.signed_global_rate:.6e} "
                    f"{QUANTITY[quantity]['global_rate_unit']}",
                    f"      element-mean L1 rate   = {stats.element_l1_rate:.6e} "
                    f"{QUANTITY[quantity]['global_rate_unit']}",
                    f"      nodal L1 rate          = {stats.nodal_l1_rate:.6e} "
                    f"{QUANTITY[quantity]['global_rate_unit']}",
                    f"      min local timescale    = {stats.min_local_timescale:.6e} s",
                    "      within-element sign cancellation = "
                    f"{100.0 * stats.nodal_cancellation_fraction:.3f} %",
                )
            )
        lines.append("")
        lines.append("  cumulative through the selected cycle:")
        for component in REPORT_COMPONENTS:
            stats = statistics[component]
            lines.extend(
                (
                    f"    {LABELS[component]}:",
                    f"      signed = {stats.cumulative_signed:.6e} "
                    f"{QUANTITY[quantity]['increment_unit']}",
                    f"      nodal absolute = {stats.cumulative_nodal_absolute:.6e} "
                    f"{QUANTITY[quantity]['increment_unit']}",
                    "      nodal absolute/reference = "
                    f"{stats.cumulative_absolute_fraction:.6e}",
                )
            )
        lines.append("")
    lines.extend(
        (
            "Interpretation:",
            "  The (standard - AP) quantity is the direct difference between",
            "  otherwise identical standard-LLF and source-aware-AP shadow steps.",
            "  It is not an analytic reference solution.  The nonconstant EmAb",
            "  contribution contains every nonconstant spatial mode; for nDOFX=2",
            "  this is the complete linear mode, including both physically aligned",
            "  and backward within-element slopes.",
            "  Cumulative shadow quantities sum one-step counterfactual responses",
            "  evaluated along the production trajectory; they are not independent",
            "  source-aware-AP or no-streaming evolutions.",
        )
    )
    return "\n".join(lines) + "\n"


def headline_report(
    elements: np.ndarray,
    main: dict[str, np.ndarray],
    cycle: int,
    r_max_cm: float,
) -> str:
    lines = [
        f"cycle {cycle}, r <= {r_max_cm / 1.0e5:g} km",
        "  EmAb nonconstant is the complete p=1 spatial response, not only its backward branch",
    ]
    for quantity in ("QH", "QL"):
        base = select(elements, cycle, quantity, "standard_total", r_max_cm)
        reference = element_reference(main, base, quantity)
        stats = {
            component: component_statistics(
                select(elements, cycle, quantity, component, r_max_cm),
                reference,
            )
            for component in (
                "standard_total",
                "ap_total",
                "emab_high",
            )
        }
        lines.append(f"  {quantity}:")
        for component in (
            "standard_total",
            "ap_total",
            "emab_high",
        ):
            item = stats[component]
            lines.append(
                f"    {LABELS[component]}: peak={item.peak_rate:.3e} "
                f"{QUANTITY[quantity]['rate_unit']}; "
                f"min local time={item.min_local_timescale:.3e} s; "
                "within-element sign cancellation="
                f"{100.0 * item.nodal_cancellation_fraction:.1f}%"
            )
        lines.append(
            "    cumulative standard/AP/nonconstant absolute = "
            + "/".join(
                f"{stats[component].cumulative_nodal_absolute:.3e}"
                for component in (
                    "standard_total",
                    "ap_total",
                    "emab_high",
                )
            )
            + f" {QUANTITY[quantity]['increment_unit']}"
        )
    return "\n".join(lines) + "\n"


def plot_developer_summary(
    elements: np.ndarray,
    nodes: np.ndarray,
    main: dict[str, np.ndarray],
    frozen_path: Path,
    cycle: int,
    r_max_cm: float,
    output: Path,
) -> None:
    fig = plt.figure(figsize=(17.0, 9.5))
    grid = fig.add_gridspec(2, 3)
    axes = np.asarray(
        [
            [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])],
            [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])],
        ]
    )
    projection_axes = np.asarray(
        [fig.add_subplot(grid[0, 2]), fig.add_subplot(grid[1, 2])]
    )
    components = (
        "standard_total",
        "ap_total",
    )
    for row, quantity in enumerate(("QH", "QL")):
        base = select(elements, cycle, quantity, "standard_total", r_max_cm)
        radius = base["r_cm"] / 1.0e5
        reference = element_reference(main, base, quantity)
        main_mask = main["r_cm"] <= float(np.max(base["r_cm"])) + 1.0e-12
        density = main["D_g_cm3"][main_mask]

        axis = axes[row, 0]
        for component in components:
            rows = select(elements, cycle, quantity, component, r_max_cm)
            if quantity == "QH":
                values = (
                    rows["step_rate_mean"]
                    / density
                    / SPECIFIC_HEATING_SCALE
                )
            else:
                values = (
                    ATOMIC_MASS_UNIT_G
                    * rows["step_rate_mean"]
                    / density
                )
            axis.plot(
                radius,
                values,
                color=COLORS[component],
                linewidth=1.5,
                label=LABELS[component],
            )
        axis.set_xscale("log")
        axis.set_xlim(float(np.min(radius)), r_max_cm / 1.0e5)
        axis.axhline(0.0, color="0.45", linewidth=0.65)
        axis.set_xlabel(r"radius [$\mathrm{km}$]")
        axis.set_ylabel(
            r"$\mathrm{d}\epsilon/\mathrm{d}t$ "
            r"[$10^{21}\ \mathrm{erg\,g^{-1}\,s^{-1}}$]"
            if quantity == "QH"
            else r"$\mathrm{d}Y_e/\mathrm{d}t$ [$\mathrm{s^{-1}}$]"
        )
        axis.set_title(
            "Element-average specific heating"
            if quantity == "QH"
            else "Element-average electron-fraction change"
        )
        axis.grid(alpha=0.25, which="both")
        axis.legend(fontsize=7.5, loc="best")

        axis = axes[row, 1]
        for component in components:
            rows = select(elements, cycle, quantity, component, r_max_cm)
            magnitude = np.abs(rows["step_rate_mean"])
            timescale = np.full(magnitude.shape, np.nan)
            active = magnitude > 0.0
            timescale[active] = reference[active] / magnitude[active]
            axis.plot(
                radius,
                timescale,
                color=COLORS[component],
                linewidth=1.5,
                label=LABELS[component],
            )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlim(float(np.min(radius)), r_max_cm / 1.0e5)
        axis.set_xlabel("radius [km]")
        axis.set_ylabel("local matter-change timescale [s]")
        axis.set_title(
            r"$e_{\rm int}/|Q_H|$" if quantity == "QH"
            else r"$n_e/|Q_L|$"
        )
        axis.grid(alpha=0.25, which="both")
        axis.legend(fontsize=7.5, loc="best")

    nodal_density = load_nodal_density(frozen_path)
    for row, quantity in enumerate(("QH", "QL")):
        projection_axis = projection_axes[row]
        for component in components:
            projected_radius, projected_rate = (
                project_matter_exchange_to_flash_cells(
                    nodes,
                    nodal_density,
                    cycle,
                    quantity,
                    component,
                    r_max_cm,
                )
            )
            projection_axis.plot(
                projected_radius,
                projected_rate,
                color=COLORS[component],
                linewidth=1.0,
                drawstyle="steps-mid",
                label=LABELS[component],
            )
        projection_axis.set_xscale("log")
        projection_axis.set_xlim(
            float(np.min(projected_radius)), r_max_cm / 1.0e5
        )
        projection_axis.axhline(0.0, color="0.45", linewidth=0.65)
        projection_axis.set_xlabel(r"radius [$\mathrm{km}$]")
        projection_axis.set_ylabel(
            r"projected $\mathrm{d}\epsilon/\mathrm{d}t$ "
            r"[$10^{21}\ \mathrm{erg\,g^{-1}\,s^{-1}}$]"
            if quantity == "QH"
            else r"projected $\mathrm{d}Y_e/\mathrm{d}t$ "
                 r"[$\mathrm{s^{-1}}$]"
        )
        projection_axis.set_title(
            "Projected FLASH-X finite-volume heating"
            if quantity == "QH"
            else "Projected FLASH-X finite-volume lepton exchange"
        )
        projection_axis.grid(alpha=0.25, which="both")
        projection_axis.legend(fontsize=8.0, loc="best")

    time_s = float(base["time_s"][0])
    fig.suptitle(
        "Standard and source-aware AP LLF matter-exchange attribution\n"
        f"cycle {cycle}, $t={time_s:.3e}\\,\\mathrm{{s}}$, "
        f"$r\\leq {r_max_cm / 1.0e5:g}\\,\\mathrm{{km}}$",
        fontsize=15,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(output, dpi=200)
    plt.close(fig)


def plot_cumulative_exchange(
    elements: np.ndarray,
    histories: dict[tuple[str, str], tuple[np.ndarray, ...]],
    inner_histories: dict[tuple[str, str], tuple[np.ndarray, ...]],
    main: dict[str, np.ndarray],
    cycle: int,
    r_max_cm: float,
    inner_r_max_cm: float,
    output_prefix: Path,
) -> None:
    """Compare full-domain and inner-core cumulative shadow responses."""

    components = ("standard_total", "ap_total")
    regions = (
        (r_max_cm, histories),
        (inner_r_max_cm, inner_histories),
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), sharex=True)
    for row, quantity in enumerate(("QH", "QL")):
        for column, (radius_limit, region_histories) in enumerate(regions):
            axis = axes[row, column]
            base = select(
                elements,
                cycle,
                quantity,
                "standard_total",
                radius_limit,
            )
            reference = element_reference(main, base, quantity)
            reference_total = float(
                np.sum(reference * base["volume_cm3"])
            )
            for component in components:
                time, _, _, _, cumulative_absolute = history(
                    region_histories, quantity, component
                )
                axis.plot(
                    time / 1.0e-3,
                    positive(cumulative_absolute / reference_total),
                    color=COLORS[component],
                    linewidth=1.5,
                    label=LABELS[component],
                )
            axis.set_yscale("log")
            axis.set_title(
                rf"$r\leq{radius_limit / 1.0e5:g}\,\mathrm{{km}}$"
            )
            axis.set_ylabel(
                "cumulative nodal-absolute exchange / "
                + (
                    "internal energy"
                    if quantity == "QH"
                    else "electron number"
                )
            )
            axis.grid(alpha=0.25, which="both")
            axis.legend(fontsize=8.0, loc="best")
    for axis in axes[-1, :]:
        axis.set_xlabel(r"time along production trajectory [$\mathrm{ms}$]")
    fig.suptitle(
        "Cumulative counterfactual matter exchange",
        fontsize=15,
    )
    fig.text(
        0.5,
        0.012,
        "Each curve sums one-step shadow responses along the production trajectory; it is not a separate evolution.",
        ha="center",
        fontsize=8.5,
    )
    fig.tight_layout(rect=(0.0, 0.05, 1.0, 0.95))
    for suffix in ("png", "pdf"):
        fig.savefig(
            Path(str(output_prefix) + f".cumulative_exchange.{suffix}"),
            dpi=200,
        )
    plt.close(fig)


def plot_spatial_structure(
    nodes: np.ndarray,
    cycle: int,
    r_max_cm: float,
    output: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.0), sharex=True)
    structure = ("standard_total", "ap_total")
    for row, quantity in enumerate(("QH", "QL")):
        for column, component in enumerate(structure):
            rows = select(
                nodes, cycle, quantity, component, r_max_cm
            )
            axis = axes[row, column]
            peak = float(np.max(np.abs(rows["step_rate"])))
            first_segment = True
            first_mean = True
            for element in np.unique(rows["element"]):
                block = rows[rows["element"] == element]
                block = block[np.argsort(block["node"])]
                rate = block["step_rate"]
                scale = float(np.max(np.abs(rate)))
                if scale <= 1.0e-10 * peak:
                    continue
                radius = block["r_cm"] / 1.0e5
                normalized = rate / scale
                weight = block["volume_weight_cm3"]
                mean_radius = float(np.sum(weight * radius) / np.sum(weight))
                mean_rate = float(np.sum(weight * rate) / np.sum(weight)) / scale
                axis.plot(
                    radius,
                    normalized,
                    color=COLORS[component],
                    marker="o",
                    markersize=2.5,
                    linewidth=0.8,
                    alpha=0.75,
                    label="two DG nodes in each element" if first_segment else None,
                )
                axis.plot(
                    mean_radius,
                    mean_rate,
                    marker=".",
                    color="black",
                    markersize=3.0,
                    linestyle="none",
                    label="signed element average" if first_mean else None,
                )
                first_segment = False
                first_mean = False
            axis.set_xscale("log")
            axis.set_xlim(
                float(np.min(rows["r_cm"])) / 1.0e5,
                r_max_cm / 1.0e5,
            )
            axis.set_ylim(-1.05, 1.05)
            axis.axhline(0.0, color="#777777", linewidth=0.6)
            axis.set_title(LABELS[component])
            axis.set_ylabel(
                rf"signed nodal {quantity} / element max $|{quantity}|$"
            )
            axis.set_xlabel("radius [km]")
            axis.grid(alpha=0.25, which="both")
            axis.legend(fontsize=8, loc="best")

    fig.suptitle(
        "Signed nodal structure at the final recorded cycle\n"
        "each colored segment is the p=1 polynomial in one spatial element; "
        "normalization is local to that element"
    )
    fig.tight_layout()
    fig.savefig(output, dpi=190)
    plt.close(fig)


def plot_specific_exchange(
    elements: np.ndarray,
    main: dict[str, np.ndarray],
    cycle: int,
    r_max_cm: float,
    output_prefix: Path,
) -> None:
    """Plot signed specific heating and electron-fraction rates."""

    components = (
        "standard_total",
        "ap_total",
    )
    base = select(elements, cycle, "QH", "standard_total", r_max_cm)
    main_mask = main["r_cm"] <= float(np.max(base["r_cm"])) + 1.0e-12
    radius = base["r_cm"] / 1.0e5
    density = main["D_g_cm3"][main_mask]
    if density.size != radius.size or not np.allclose(
        main["r_cm"][main_mask], base["r_cm"], rtol=0.0, atol=1.0e-10
    ):
        raise ValueError("matter density and attribution element grids differ")

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10.5, 7.5),
        sharex=True,
        gridspec_kw={"hspace": 0.0},
    )
    for component in components:
        heating = select(
            elements, cycle, "QH", component, r_max_cm
        )["step_rate_mean"] / density
        dye_dt = (
            ATOMIC_MASS_UNIT_G
            * select(elements, cycle, "QL", component, r_max_cm)[
                "step_rate_mean"
            ]
            / density
        )
        axes[0].plot(
            radius,
            heating / SPECIFIC_HEATING_SCALE,
            color=COLORS[component],
            linewidth=1.35,
            label=LABELS[component],
        )
        axes[1].plot(
            radius,
            dye_dt,
            color=COLORS[component],
            linewidth=1.35,
            label=LABELS[component],
        )

    axes[0].set_ylabel(
        r"$\mathrm{d}\epsilon/\mathrm{d}t$ "
        r"[$10^{21}\ \mathrm{erg\,g^{-1}\,s^{-1}}$]"
    )
    axes[1].set_ylabel(r"$dY_e/dt$ [$\mathrm{s^{-1}}$]")
    axes[1].set_xlabel("radius [km]")
    for axis in axes:
        axis.axhline(0.0, color="0.45", linewidth=0.65)
        axis.set_xlim(0.0, r_max_cm / 1.0e5)
        axis.grid(alpha=0.22)
        axis.margins(x=0.0)
    axes[0].tick_params(labelbottom=False)
    axes[0].legend(loc="best", fontsize=8.5)
    fig.align_ylabels(axes)
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.10, top=0.98, hspace=0.0)
    for suffix in ("png", "pdf"):
        fig.savefig(Path(str(output_prefix) + f".specific_exchange.{suffix}"), dpi=200)
    plt.close(fig)


def plot_nodal_specific_exchange(
    nodes: np.ndarray,
    elements: np.ndarray,
    main: dict[str, np.ndarray],
    frozen_path: Path,
    cycle: int,
    r_max_cm: float,
    output_prefix: Path,
) -> None:
    """Plot signed nodal specific exchange without joining element faces."""

    components = ("standard_total", "ap_total")
    nodal_density = load_nodal_density(frozen_path)
    element_reference = select(
        elements, cycle, "QH", "standard_total", r_max_cm
    )
    element_mask = (
        main["r_cm"]
        <= float(np.max(element_reference["r_cm"])) + 1.0e-12
    )
    element_radius = element_reference["r_cm"] / 1.0e5
    element_density = main["D_g_cm3"][element_mask]
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10.5, 7.5),
        sharex=True,
        gridspec_kw={"hspace": 0.0},
    )

    for component in components:
        heating_rows = select(nodes, cycle, "QH", component, r_max_cm)
        lepton_rows = select(nodes, cycle, "QL", component, r_max_cm)
        for axis, rows, quantity in (
            (axes[0], heating_rows, "QH"),
            (axes[1], lepton_rows, "QL"),
        ):
            first = True
            for element in np.unique(rows["element"].astype(int)):
                block = rows[rows["element"].astype(int) == element]
                block = block[np.argsort(block["node"])]
                density = np.array(
                    [
                        nodal_density[(element, int(node))]
                        for node in block["node"]
                    ]
                )
                if quantity == "QH":
                    values = (
                        block["step_rate"]
                        / density
                        / SPECIFIC_HEATING_SCALE
                    )
                else:
                    values = (
                        ATOMIC_MASS_UNIT_G
                        * block["step_rate"]
                        / density
                    )
                axis.plot(
                    block["r_cm"] / 1.0e5,
                    values,
                    color=COLORS[component],
                    linewidth=0.8,
                    marker="o",
                    markersize=2.2,
                    alpha=0.80,
                    label=(
                        f"{LABELS[component]}: nodal"
                        if first
                        else None
                    ),
                )
                first = False

        heating_mean = select(
            elements, cycle, "QH", component, r_max_cm
        )["step_rate_mean"] / element_density / SPECIFIC_HEATING_SCALE
        dye_dt_mean = (
            ATOMIC_MASS_UNIT_G
            * select(elements, cycle, "QL", component, r_max_cm)[
                "step_rate_mean"
            ]
            / element_density
        )
        for axis, values in (
            (axes[0], heating_mean),
            (axes[1], dye_dt_mean),
        ):
            axis.plot(
                element_radius,
                values,
                color=COLORS[component],
                linewidth=0.85,
                linestyle="--",
                alpha=0.95,
                label=f"{LABELS[component]}: element average",
            )

    axes[0].set_ylabel(
        r"$\mathrm{d}\epsilon/\mathrm{d}t$ "
        r"[$10^{21}\ \mathrm{erg\,g^{-1}\,s^{-1}}$]"
    )
    axes[1].set_ylabel(r"$dY_e/dt$ [$\mathrm{s^{-1}}$]")
    axes[1].set_xlabel(r"radius [$\mathrm{km}$]")
    for axis in axes:
        axis.axhline(0.0, color="0.45", linewidth=0.65)
        axis.set_xlim(0.0, r_max_cm / 1.0e5)
        axis.grid(alpha=0.22)
        axis.margins(x=0.0)
    axes[0].tick_params(labelbottom=False)
    axes[0].legend(loc="best", fontsize=8.5)
    fig.align_ylabels(axes)
    fig.subplots_adjust(
        left=0.12,
        right=0.98,
        bottom=0.10,
        top=0.98,
        hspace=0.0,
    )
    for suffix in ("png", "pdf"):
        fig.savefig(
            Path(str(output_prefix) + f".specific_exchange_nodes.{suffix}"),
            dpi=200,
        )
    plt.close(fig)


def history(
    histories: dict[tuple[str, str], tuple[np.ndarray, ...]],
    quantity: str,
    component: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    key = (quantity, component)
    if key not in histories:
        raise RuntimeError(
            f"missing history for quantity={quantity}, component={component}"
        )
    return histories[key]


def main() -> None:
    args = parse_args()
    prefix = attribution_prefix(args.prefix)
    element_path = Path(str(prefix) + ".elements")
    cycle = (
        last_data_cycle(element_path)
        if args.cycle is None
        else args.cycle
    )
    r_max_cm = args.r_max_km * 1.0e5
    inner_r_max_cm = min(25.0e5, r_max_cm)
    elements, histories, inner_histories = load_element_cycle_and_histories(
        element_path,
        cycle,
        r_max_cm,
        inner_r_max_cm,
    )
    main_path = main_diagnostic_path(prefix)
    main_data = load_main_diagnostic(main_path)
    nodes = load_cycle_table(Path(str(prefix) + ".nodes"), cycle)
    output_prefix = args.output_prefix
    if output_prefix is None:
        output_prefix = Path(str(prefix) + f".cycle_{cycle:08d}")

    summary = report(elements, main_data, cycle, r_max_cm)
    print(headline_report(elements, main_data, cycle, r_max_cm), end="")
    summary_path = Path(str(output_prefix) + ".summary.txt")
    summary_path.write_text(summary)

    developer_path = Path(str(output_prefix) + ".developer_summary.png")
    structure_path = Path(str(output_prefix) + ".nodal_spatial_structure.png")
    frozen_path = Path(str(main_path) + ".frozen_transient")
    plot_developer_summary(
        elements,
        nodes,
        main_data,
        frozen_path,
        cycle,
        r_max_cm,
        developer_path,
    )
    plot_cumulative_exchange(
        elements,
        histories,
        inner_histories,
        main_data,
        cycle,
        r_max_cm,
        inner_r_max_cm,
        output_prefix,
    )
    cumulative_path = Path(
        str(output_prefix) + ".cumulative_exchange.png"
    )
    plot_spatial_structure(
        nodes, cycle, min(r_max_cm, 30.0e5), structure_path
    )
    plot_specific_exchange(
        elements, main_data, cycle, r_max_cm, output_prefix
    )
    specific_path = Path(str(output_prefix) + ".specific_exchange.png")
    plot_nodal_specific_exchange(
        nodes,
        elements,
        main_data,
        frozen_path,
        cycle,
        r_max_cm,
        output_prefix,
    )
    nodal_specific_path = Path(
        str(output_prefix) + ".specific_exchange_nodes.png"
    )
    print(f"wrote {summary_path}")
    print(f"wrote {developer_path}")
    print(f"wrote {cumulative_path}")
    print(f"wrote {structure_path}")
    print(f"wrote {specific_path}")
    print(f"wrote {nodal_specific_path}")


if __name__ == "__main__":
    main()
