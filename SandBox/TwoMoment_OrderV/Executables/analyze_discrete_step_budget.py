#!/usr/bin/env python3

"""Validate and analyze the one-step discrete transport/collision budget."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ATOMIC_MASS_UNIT_G = 1.66053906892e-24


COLUMNS = (
    "stage",
    "component",
    "element",
    "r_cm",
    "volume_cm3",
    "species",
    "wEX",
    "wIM",
    "QL",
    "QE",
    "QH",
    "QM1",
    "QH_exact",
)

REQUIRED_STAGE_COMPONENTS = {
    "COLLISION_SOURCE_PRE",
    "COLLISION_SOURCE_POST",
    "COLLISION_SOURCE_POST_LIMITER",
    "COLLISION_MATTER_ACTUAL",
    "COLLISION_RADIATION_ACTUAL",
    "COLLISION_CONSERVATION",
    "STREAMING_TOTAL",
    "STREAMING_CENTERED",
    "STREAMING_SURFACE",
    "STREAMING_VOLUME",
    "STREAMING_LLF",
    "STREAMING_OBSERVER",
    "STREAMING_DG_CLOSURE",
    "STREAMING_TOTAL_CLOSURE",
    "SEMIDISCRETE_BALANCE_POST",
    "PRE_STAGE_LIMITER_MATTER",
    "PRE_STAGE_LIMITER_RADIATION",
    "PRE_STAGE_LIMITER_CONSERVATION",
    "POST_STAGE_LIMITER_MATTER",
    "POST_STAGE_LIMITER_RADIATION",
    "POST_STAGE_LIMITER_CONSERVATION",
}

REQUIRED_FULL_COMPONENTS = {
    "FULL_STEP_MATTER_ACTUAL",
    "FULL_STEP_RADIATION_ACTUAL",
    "FULL_WEIGHTED_STREAMING",
    "FULL_WEIGHTED_COLLISION_MATTER",
    "FULL_WEIGHTED_COLLISION_RADIATION",
    "FULL_REMAINDER_MATTER",
    "FULL_REMAINDER_RADIATION",
    "FULL_ALGEBRAIC_CLOSURE",
}


@dataclass(frozen=True)
class Ledger:
    stage: np.ndarray
    component: np.ndarray
    element: np.ndarray
    r_cm: np.ndarray
    volume_cm3: np.ndarray
    species: np.ndarray
    w_ex: np.ndarray
    w_im: np.ndarray
    values: dict[str, np.ndarray]


def read_ledger(path: Path) -> Ledger:
    rows: list[list[str]] = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            words = line.split()
            if not words or words[0].startswith("#"):
                continue
            if len(words) != len(COLUMNS):
                raise ValueError(
                    f"{path}: expected {len(COLUMNS)} columns, got "
                    f"{len(words)} in: {line.rstrip()}"
                )
            rows.append(words)
    if not rows:
        raise ValueError(f"{path}: no discrete-budget rows")

    array = np.asarray(rows, dtype=object)
    numeric = lambda index, kind=float: np.asarray(array[:, index], dtype=kind)
    values = {
        name: numeric(index)
        for index, name in enumerate(COLUMNS[8:], start=8)
    }
    ledger = Ledger(
        stage=numeric(0, int),
        component=np.asarray(array[:, 1], dtype=str),
        element=numeric(2, int),
        r_cm=numeric(3),
        volume_cm3=numeric(4),
        species=numeric(5, int),
        w_ex=numeric(6),
        w_im=numeric(7),
        values=values,
    )
    for name, values_ in ledger.values.items():
        if not np.all(np.isfinite(values_)):
            raise ValueError(f"{path}: non-finite values in {name}")
    return ledger


def rows(
    ledger: Ledger, component: str, *, stage: int | None = None, species: int = 0
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mask = (ledger.component == component) & (ledger.species == species)
    if stage is not None:
        mask &= ledger.stage == stage
    indices = np.flatnonzero(mask)
    if not len(indices):
        raise ValueError(
            f"missing component={component}, stage={stage}, species={species}"
        )
    order = np.argsort(ledger.element[indices])
    indices = indices[order]
    return indices, ledger.r_cm[indices], ledger.volume_cm3[indices]


def field(
    ledger: Ledger,
    component: str,
    quantity: str,
    *,
    stage: int | None = None,
    species: int = 0,
) -> np.ndarray:
    indices, _, _ = rows(ledger, component, stage=stage, species=species)
    return ledger.values[quantity][indices]


def relative_l1(error: np.ndarray, *scales: np.ndarray) -> float:
    denominator = max(
        *(float(np.sum(np.abs(scale))) for scale in scales),
        np.finfo(float).tiny,
    )
    return float(np.sum(np.abs(error)) / denominator)


def validate(ledger: Ledger, tolerance: float) -> list[int]:
    stages = sorted(set(ledger.stage[ledger.stage > 0]))
    if not stages:
        raise ValueError("no IMEX stage rows")
    for stage in stages:
        available = set(ledger.component[ledger.stage == stage])
        missing = REQUIRED_STAGE_COMPONENTS - available
        if missing:
            raise ValueError(f"stage {stage} is incomplete: missing {sorted(missing)}")
    available = set(ledger.component[ledger.stage == 0])
    missing = REQUIRED_FULL_COMPONENTS - available
    if missing:
        raise ValueError(f"full-step ledger is incomplete: missing {sorted(missing)}")

    failures: list[str] = []
    quantities = ("QL", "QE", "QH", "QM1")
    for stage in stages:
        for quantity in quantities:
            total = field(ledger, "STREAMING_TOTAL", quantity, stage=stage)
            centered = field(ledger, "STREAMING_CENTERED", quantity, stage=stage)
            surface = field(ledger, "STREAMING_SURFACE", quantity, stage=stage)
            volume = field(ledger, "STREAMING_VOLUME", quantity, stage=stage)
            llf = field(ledger, "STREAMING_LLF", quantity, stage=stage)
            observer = field(ledger, "STREAMING_OBSERVER", quantity, stage=stage)
            dg_closure = field(
                ledger, "STREAMING_DG_CLOSURE", quantity, stage=stage
            )
            total_closure = field(
                ledger, "STREAMING_TOTAL_CLOSURE", quantity, stage=stage
            )
            error = total - centered - llf - observer
            rel = relative_l1(error, total, centered, llf, observer)
            if rel > tolerance:
                failures.append(
                    f"stage {stage} {quantity}: streaming split rel-L1={rel:.3e}"
                )
            error = centered - surface - volume - dg_closure
            rel = relative_l1(error, centered, surface, volume, dg_closure)
            if rel > tolerance:
                failures.append(
                    f"stage {stage} {quantity}: DG split rel-L1={rel:.3e}"
                )
            rel = relative_l1(total_closure, total)
            if rel > tolerance:
                failures.append(
                    f"stage {stage} {quantity}: stored total closure={rel:.3e}"
                )

            matter = field(
                ledger, "COLLISION_MATTER_ACTUAL", quantity, stage=stage
            )
            radiation = field(
                ledger, "COLLISION_RADIATION_ACTUAL", quantity, stage=stage
            )
            conservation = field(
                ledger, "COLLISION_CONSERVATION", quantity, stage=stage
            )
            error = conservation - matter - radiation
            rel = relative_l1(error, conservation, matter, radiation)
            if rel > tolerance:
                failures.append(
                    f"stage {stage} {quantity}: collision ledger identity={rel:.3e}"
                )

    for quantity in quantities:
        full_matter = field(ledger, "FULL_STEP_MATTER_ACTUAL", quantity)
        full_radiation = field(ledger, "FULL_STEP_RADIATION_ACTUAL", quantity)
        streaming = field(ledger, "FULL_WEIGHTED_STREAMING", quantity)
        collision_matter = field(
            ledger, "FULL_WEIGHTED_COLLISION_MATTER", quantity
        )
        collision_radiation = field(
            ledger, "FULL_WEIGHTED_COLLISION_RADIATION", quantity
        )
        remainder_matter = field(ledger, "FULL_REMAINDER_MATTER", quantity)
        remainder_radiation = field(
            ledger, "FULL_REMAINDER_RADIATION", quantity
        )
        error_matter = full_matter - collision_matter - remainder_matter
        error_radiation = (
            full_radiation - streaming - collision_radiation - remainder_radiation
        )
        rel_m = relative_l1(
            error_matter, full_matter, collision_matter, remainder_matter
        )
        rel_r = relative_l1(
            error_radiation,
            full_radiation,
            streaming,
            collision_radiation,
            remainder_radiation,
        )
        if rel_m > tolerance:
            failures.append(f"full {quantity}: matter identity={rel_m:.3e}")
        if rel_r > tolerance:
            failures.append(f"full {quantity}: radiation identity={rel_r:.3e}")

    exact_heating = field(ledger, "FULL_STEP_MATTER_ACTUAL", "QH_exact")
    weighted_heating = field(
        ledger, "FULL_WEIGHTED_COLLISION_MATTER", "QH"
    )
    exact_remainder = field(ledger, "FULL_REMAINDER_MATTER", "QH_exact")
    rel = relative_l1(
        exact_heating - weighted_heating - exact_remainder,
        exact_heating,
        weighted_heating,
        exact_remainder,
    )
    if rel > tolerance:
        failures.append(f"full exact internal-energy identity={rel:.3e}")

    if failures:
        raise ValueError("inconsistent discrete ledger:\n  " + "\n  ".join(failures))
    return stages


def high_pass(values: np.ndarray) -> np.ndarray:
    result = np.zeros_like(values)
    if len(values) < 3:
        return result
    result[1:-1] = values[1:-1] - 0.5 * (values[:-2] + values[2:])
    return result


def high_frequency_fraction(values: np.ndarray) -> float:
    centered = values - np.mean(values)
    power = np.abs(np.fft.rfft(centered)) ** 2
    if len(power) <= 1 or np.sum(power[1:]) == 0.0:
        return 0.0
    start = max(1, int(np.ceil(0.75 * (len(power) - 1))))
    return float(np.sum(power[start:]) / np.sum(power[1:]))


def correlation(left: np.ndarray, right: np.ndarray) -> float:
    left = left - np.mean(left)
    right = right - np.mean(right)
    norm = np.linalg.norm(left) * np.linalg.norm(right)
    return 0.0 if norm == 0.0 else float(np.dot(left, right) / norm)


def stage_weight(ledger: Ledger, stage: int, name: str) -> float:
    mask = ledger.stage == stage
    values = ledger.w_ex[mask] if name == "explicit" else ledger.w_im[mask]
    return float(values[0])


def weighted_stage_sum(
    ledger: Ledger, stages: list[int], component: str, quantity: str, weight: str
) -> np.ndarray:
    result = np.zeros_like(field(ledger, component, quantity, stage=stages[0]))
    for stage in stages:
        result += stage_weight(ledger, stage, weight) * field(
            ledger, component, quantity, stage=stage
        )
    return result


def norm_l1(values: np.ndarray, volume: np.ndarray) -> float:
    return float(np.sum(np.abs(values) * volume))


def read_reference_state(path: Path) -> dict[str, np.ndarray]:
    """Read the element-averaged matter state written beside the budget."""

    header = None
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("# r_cm "):
                header = line[2:].split()
                break
    if header is None:
        raise ValueError(f"{path}: no element-state column header")
    values = np.loadtxt(path, comments="#", ndmin=2)
    if values.shape[1] != len(header):
        raise ValueError(
            f"{path}: header has {len(header)} columns, data have {values.shape[1]}"
        )
    state = {name: values[:, index] for index, name in enumerate(header)}
    required = {"r_cm", "D_g_cm3", "Ye", "eint_after_erg_cm3"}
    missing = required - state.keys()
    if missing:
        raise ValueError(f"{path}: missing reference fields {sorted(missing)}")
    if "ne_after_cm3" not in state:
        state["ne_after_cm3"] = (
            state["D_g_cm3"] * state["Ye"] / ATOMIC_MASS_UNIT_G
        )
    return state


def infer_reference_state_path(ledger_path: Path) -> Path:
    suffix = ".discrete_budget"
    if str(ledger_path).endswith(suffix):
        return Path(str(ledger_path)[: -len(suffix)])
    raise ValueError(
        "cannot infer the reference-state file; pass --reference-state"
    )


def reference_reservoirs(
    state: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return rho and step-centered internal-energy/electron reservoirs."""

    rho = state["D_g_cm3"]
    eint = state["eint_after_erg_cm3"]
    ne = state["ne_after_cm3"]
    if "eint_before_erg_cm3" in state:
        eint = 0.5 * (state["eint_before_erg_cm3"] + eint)
    if "ne_before_cm3" in state:
        ne = 0.5 * (state["ne_before_cm3"] + ne)
    return rho, eint, ne


def origin_rates(
    ledger: Ledger, stages: list[int], quantity: str
) -> dict[str, np.ndarray]:
    """Return numerical origins in matter-gain sign for one quantity."""

    exact_quantity = "QH_exact" if quantity == "QH" else quantity
    pre_stage_limiters = np.zeros_like(
        field(ledger, "FULL_STEP_MATTER_ACTUAL", exact_quantity)
    )
    post_stage_limiters = np.zeros_like(pre_stage_limiters)
    for stage in stages:
        pre_stage_limiters += field(
            ledger, "PRE_STAGE_LIMITER_MATTER", exact_quantity, stage=stage
        )
        post_stage_limiters += field(
            ledger, "POST_STAGE_LIMITER_MATTER", exact_quantity, stage=stage
        )

    return {
        "applied matter change": field(
            ledger, "FULL_STEP_MATTER_ACTUAL", exact_quantity
        ),
        "implicit collisions": field(
            ledger, "FULL_WEIGHTED_COLLISION_MATTER", quantity
        ),
        "matter remainder": field(
            ledger, "FULL_REMAINDER_MATTER", exact_quantity
        ),
        "pre-stage fluid limiters": pre_stage_limiters,
        "post-stage fluid limiters": post_stage_limiters,
        "streaming-equivalent": field(
            ledger, "FULL_WEIGHTED_STREAMING", quantity
        ),
        "LLF streaming": weighted_stage_sum(
            ledger, stages, "STREAMING_LLF", quantity, "explicit"
        ),
        "observer corrections": weighted_stage_sum(
            ledger, stages, "STREAMING_OBSERVER", quantity, "explicit"
        ),
        "centered streaming": weighted_stage_sum(
            ledger, stages, "STREAMING_CENTERED", quantity, "explicit"
        ),
        "centered surface": weighted_stage_sum(
            ledger, stages, "STREAMING_SURFACE", quantity, "explicit"
        ),
        "volume derivative": weighted_stage_sum(
            ledger, stages, "STREAMING_VOLUME", quantity, "explicit"
        ),
        "DG split closure": weighted_stage_sum(
            ledger, stages, "STREAMING_DG_CLOSURE", quantity, "explicit"
        ),
    }


def report_exchange_impact(
    ledger: Ledger,
    stages: list[int],
    state: dict[str, np.ndarray],
) -> None:
    """Normalize the grid-scale exchange by the local matter reservoirs."""

    _, radius, volume = rows(ledger, "FULL_STEP_MATTER_ACTUAL")
    if len(state["r_cm"]) != len(radius) or not np.allclose(
        state["r_cm"], radius, rtol=1.0e-12, atol=0.0
    ):
        raise ValueError("reference state and discrete budget use different elements")

    rho, eint, ne = reference_reservoirs(state)
    if np.any(rho <= 0.0) or np.any(eint <= 0.0) or np.any(ne <= 0.0):
        raise ValueError("exchange normalization requires positive rho, eint, and ne")

    print("  grid-scale exchange relative to the step-centered matter state:")
    for quantity, reservoir, specific_factor, unit in (
        ("QH", eint, 1.0 / rho, "erg g^-1 s^-1"),
        ("QL", ne, ATOMIC_MASS_UNIT_G / rho, "s^-1 in Ye"),
    ):
        print(f"    {quantity}:")
        for name, source in origin_rates(ledger, stages, quantity).items():
            grid_scale = high_pass(source)
            fractional_rate = grid_scale / reservoir
            specific_rate = grid_scale * specific_factor
            global_rate = norm_l1(grid_scale, volume) / max(
                float(np.sum(reservoir * volume)), np.finfo(float).tiny
            )
            global_timescale = (
                np.inf if global_rate == 0.0 else 1.0 / global_rate
            )
            finite_local = np.abs(fractional_rate) > 0.0
            local_timescale = (
                float(np.min(1.0 / np.abs(fractional_rate[finite_local])))
                if np.any(finite_local)
                else np.inf
            )
            print(
                f"      {name:23s} "
                f"max|grid-scale specific rate|={np.max(np.abs(specific_rate)):.6e} {unit}; "
                f"max|fractional rate|={np.max(np.abs(fractional_rate)):.6e} s^-1; "
                f"min local tau={local_timescale:.6e} s; "
                f"global L1 tau={global_timescale:.6e} s"
            )


def plot_exchange_impact(
    ledger: Ledger,
    stages: list[int],
    state: dict[str, np.ndarray],
    output: Path,
    rmax_km: float | None = None,
) -> None:
    """Plot grid-scale heating and lepton exchange relative to matter."""

    _, radius_cm, _ = rows(ledger, "FULL_STEP_MATTER_ACTUAL")
    radius_km = radius_cm / 1.0e5
    rho, eint, ne = reference_reservoirs(state)
    heating = origin_rates(ledger, stages, "QH")
    lepton = origin_rates(ledger, stages, "QL")

    if rmax_km is None:
        rmax_km = float(radius_km[-1])

    selected_origins = (
        "applied matter change",
        "implicit collisions",
        "matter remainder",
        "pre-stage fluid limiters",
        "post-stage fluid limiters",
        "streaming-equivalent",
        "LLF streaming",
        "centered surface",
        "volume derivative",
    )

    def positive_magnitude(values: np.ndarray) -> np.ndarray:
        magnitude = np.abs(values)
        return np.where(magnitude > 0.0, magnitude, np.nan)

    fig, axes = plt.subplots(2, 3, figsize=(19, 10), sharex=True)
    for name in selected_origins:
        qh = high_pass(heating[name])
        ql = high_pass(lepton[name])
        axes[0, 0].plot(
            radius_km, positive_magnitude(qh / rho), label=name
        )
        axes[0, 1].plot(
            radius_km, positive_magnitude(qh / eint), label=name
        )
        tau_h = np.full_like(qh, np.nan)
        nonzero_h = qh != 0.0
        tau_h[nonzero_h] = eint[nonzero_h] / np.abs(qh[nonzero_h])
        axes[0, 2].plot(radius_km, tau_h, label=name)

        axes[1, 0].plot(
            radius_km,
            positive_magnitude(ATOMIC_MASS_UNIT_G * ql / rho),
            label=name,
        )
        axes[1, 1].plot(
            radius_km, positive_magnitude(ql / ne), label=name
        )
        tau_l = np.full_like(ql, np.nan)
        nonzero_l = ql != 0.0
        tau_l[nonzero_l] = ne[nonzero_l] / np.abs(ql[nonzero_l])
        axes[1, 2].plot(radius_km, tau_l, label=name)

    axes[0, 0].set_ylabel(r"grid-scale $|\dot e|$ [erg g$^{-1}$ s$^{-1}$]")
    axes[0, 1].set_ylabel(r"grid-scale $|Q_H|/e_{\rm int}$ [s$^{-1}$]")
    axes[0, 2].set_ylabel(r"$\tau_H=e_{\rm int}/|Q_H|$ [s]")
    axes[1, 0].set_ylabel(r"grid-scale $|\dot Y_e|$ [s$^{-1}$]")
    axes[1, 1].set_ylabel(r"grid-scale $|Q_L|/n_e$ [s$^{-1}$]")
    axes[1, 2].set_ylabel(r"$\tau_L=n_e/|Q_L|$ [s]")
    axes[0, 0].set_title("A. Specific heating by numerical origin")
    axes[0, 1].set_title("B. Heating relative to internal energy")
    axes[0, 2].set_title("C. Heating timescale")
    axes[1, 0].set_title("D. Electron-fraction change by numerical origin")
    axes[1, 1].set_title("E. Lepton exchange relative to electron density")
    axes[1, 2].set_title("F. Lepton-exchange timescale")
    for axis in axes.flat:
        axis.set_yscale("log")
        axis.set_xlim(radius_km[0], rmax_km)
        axis.set_xlabel("radius [km]")
        axis.grid(alpha=0.2)
        axis.legend(fontsize=7, loc="best")
    fig.suptitle(
        "Grid-scale matter exchange by numerical origin, normalized by the local matter state",
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def read_face_ledger(path: Path) -> dict[str, np.ndarray] | None:
    if not path.exists():
        return None
    header = None
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("# stage limiter_applied "):
                header = line[2:].split()
                break
    if header is None:
        raise ValueError(f"No face-ledger header found in {path}")
    values = np.loadtxt(path, comments="#", ndmin=2)
    if values.shape[1] != len(header):
        raise ValueError(
            f"{path}: header has {len(header)} columns, data have {values.shape[1]}"
        )
    return {name: values[:, index] for index, name in enumerate(header)}


def report_faces(faces: dict[str, np.ndarray], stages: list[int]) -> None:
    print("  face-level origin of the LLF term:")
    for stage in stages:
        if stage_weight_from_faces(faces, stage) == 0:
            continue
        for species in sorted(set(faces["species"].astype(int))):
            mask = (
                (faces["stage"].astype(int) == stage)
                & (faces["limiter_applied"].astype(int) == 1)
                & (faces["species"].astype(int) == species)
                & (faces["equilibrium_split_available"].astype(int) == 1)
            )
            indices = np.flatnonzero(mask)
            if not len(indices):
                continue
            area = 4.0 * np.pi * faces["r_cm"][indices] ** 2
            centered = area * faces["centered_flux_erg_cm2_s"][indices]
            llf = area * faces["LLF_penalty_flux_erg_cm2_s"][indices]
            llf_eq = area * faces["equilibrium_LLF_penalty_flux_erg_cm2_s"][indices]
            llf_neq = area * faces["nonequilibrium_LLF_penalty_flux_erg_cm2_s"][indices]
            strongest = indices[np.argmax(np.abs(llf))]
            l2_llf = np.linalg.norm(llf)
            centered_ratio = l2_llf / max(np.linalg.norm(centered), np.finfo(float).tiny)
            eq_ratio = np.linalg.norm(llf_eq) / max(l2_llf, np.finfo(float).tiny)
            neq_ratio = np.linalg.norm(llf_neq) / max(l2_llf, np.finfo(float).tiny)
            tau = np.sqrt(
                abs(
                    faces["tau_total_left"][strongest]
                    * faces["tau_total_right"][strongest]
                )
            )
            print(
                f"    stage {stage}, species {species}: "
                f"||LLF||/||centered||={centered_ratio:.4e}, "
                f"equilibrium/LLF={eq_ratio:.4e}, "
                f"nonequilibrium/LLF={neq_ratio:.4e}; "
                f"strongest at r={faces['r_cm'][strongest]/1.0e5:.3f} km, "
                f"tau={tau:.4e}"
            )


def stage_weight_from_faces(faces: dict[str, np.ndarray], stage: int) -> int:
    mask = faces["stage"].astype(int) == stage
    return int(np.max(faces["explicit_rhs_used"][mask]))


def report(
    ledger: Ledger,
    stages: list[int],
    faces: dict[str, np.ndarray] | None = None,
) -> None:
    target = field(ledger, "FULL_STEP_MATTER_ACTUAL", "QH_exact")
    _, radius, volume = rows(ledger, "FULL_STEP_MATTER_ACTUAL")
    scale = max(norm_l1(target, volume), np.finfo(float).tiny)

    weighted_collision = field(
        ledger, "FULL_WEIGHTED_COLLISION_MATTER", "QH"
    )
    matter_remainder = field(ledger, "FULL_REMAINDER_MATTER", "QH_exact")
    weighted_streaming = field(ledger, "FULL_WEIGHTED_STREAMING", "QH")
    weighted_llf = weighted_stage_sum(
        ledger, stages, "STREAMING_LLF", "QH", "explicit"
    )
    weighted_centered = weighted_stage_sum(
        ledger, stages, "STREAMING_CENTERED", "QH", "explicit"
    )
    weighted_surface = weighted_stage_sum(
        ledger, stages, "STREAMING_SURFACE", "QH", "explicit"
    )
    weighted_volume = weighted_stage_sum(
        ledger, stages, "STREAMING_VOLUME", "QH", "explicit"
    )
    collision_conservation = field(ledger, "FULL_ALGEBRAIC_CLOSURE", "QH")

    explicit_norm = max(norm_l1(weighted_streaming, volume), np.finfo(float).tiny)
    target_hp = high_pass(target)
    stream_hp = high_pass(weighted_streaming)
    llf_hp = high_pass(weighted_llf)
    centered_hp = high_pass(weighted_centered)
    surface_hp = high_pass(weighted_surface)
    volume_hp = high_pass(weighted_volume)
    target_hp_norm = max(norm_l1(target_hp, volume), np.finfo(float).tiny)

    print("Discrete one-step budget: all algebraic identities passed")
    print(
        "  full matter heating: "
        f"integral={np.sum(target * volume): .8e} erg/s  "
        f"L1={scale:.8e} erg/s  "
        f"high-frequency power fraction={high_frequency_fraction(target):.6e}"
    )

    print("  exact full-step attribution:")
    print(
        "    actual matter vs weighted collision increment: "
        f"relative L1={norm_l1(target-weighted_collision, volume)/scale:.6e}"
    )
    print(
        "    matter limiter/method remainder: "
        f"relative L1={norm_l1(matter_remainder, volume)/scale:.6e}"
    )
    print(
        "    weighted collision conservation residual: "
        f"relative to applied heating={norm_l1(collision_conservation, volume)/scale:.6e}"
    )
    print("  origin of the grid-scale structure:")
    print(
        "    weighted streaming vs applied heating: "
        f"high-pass correlation={correlation(stream_hp,target_hp):.6f}, "
        f"high-pass L1 ratio={norm_l1(stream_hp,volume)/target_hp_norm:.6e}"
    )
    print(
        "    LLF contribution (including any source-aware scaling): "
        f"L1={norm_l1(weighted_llf,volume)/explicit_norm:.6e}, "
        f"high-pass/applied={norm_l1(llf_hp,volume)/target_hp_norm:.6e}, "
        f"high-pass correlation with applied heating={correlation(llf_hp,target_hp):.6f}"
    )
    print(
        "    centered share of weighted streaming: "
        f"L1={norm_l1(weighted_centered,volume)/explicit_norm:.6e}, "
        f"high-pass/applied={norm_l1(centered_hp,volume)/target_hp_norm:.6e}, "
        f"high-pass correlation with applied heating={correlation(centered_hp,target_hp):.6f}"
    )
    print(
        "      centered surface contribution: "
        f"L1={norm_l1(weighted_surface,volume)/explicit_norm:.6e}, "
        f"high-pass/applied={norm_l1(surface_hp,volume)/target_hp_norm:.6e}, "
        f"high-pass correlation with applied heating={correlation(surface_hp,target_hp):.6f}"
    )
    print(
        "      volume contribution: "
        f"L1={norm_l1(weighted_volume,volume)/explicit_norm:.6e}, "
        f"high-pass/applied={norm_l1(volume_hp,volume)/target_hp_norm:.6e}, "
        f"high-pass correlation with applied heating={correlation(volume_hp,target_hp):.6f}"
    )

    for stage in stages:
        matter = field(ledger, "COLLISION_MATTER_ACTUAL", "QH", stage=stage)
        radiation = field(
            ledger, "COLLISION_RADIATION_ACTUAL", "QH", stage=stage
        )
        defect = matter + radiation
        denominator = max(
            float(np.sum(np.abs(matter) * volume)),
            float(np.sum(np.abs(radiation) * volume)),
            np.finfo(float).tiny,
        )
        balance = field(
            ledger, "SEMIDISCRETE_BALANCE_POST", "QH", stage=stage
        )
        if stage_weight(ledger, stage, "implicit") == 0.0:
            continue
        print(
            f"  stage {stage}: collision conservation rel-L1="
            f"{np.sum(np.abs(defect) * volume) / denominator:.6e}; "
            f"semidiscrete-balance HF-power={high_frequency_fraction(balance):.6e}"
        )

    full_matter = field(ledger, "FULL_STEP_MATTER_ACTUAL", "QH_exact")
    remainder_matter = field(ledger, "FULL_REMAINDER_MATTER", "QH_exact")
    full_radiation = field(ledger, "FULL_STEP_RADIATION_ACTUAL", "QH")
    remainder_radiation = field(ledger, "FULL_REMAINDER_RADIATION", "QH")
    print(
        "  full-step non-tableau remainder relative L1: "
        f"matter={np.sum(np.abs(remainder_matter)*volume) / max(np.sum(np.abs(full_matter)*volume), np.finfo(float).tiny):.6e} "
        f"radiation={np.sum(np.abs(remainder_radiation)*volume) / max(np.sum(np.abs(full_radiation)*volume), np.finfo(float).tiny):.6e}"
    )
    if faces is not None:
        report_faces(faces, stages)

    if (
        correlation(llf_hp, target_hp) > 0.9
        and norm_l1(weighted_llf, volume) / explicit_norm > 0.5
        and norm_l1(matter_remainder, volume) / scale < 1.0e-3
    ):
        print(
            "  ATTRIBUTION: the grid-scale applied heating is generated by the "
            "LLF part of the explicit streaming operator and transferred to matter "
            "by the subsequent implicit collisions. It is not generated by a limiter "
            "or by the final matter mapping in this step."
        )
    elif (
        correlation(surface_hp, target_hp) > 0.9
        and norm_l1(surface_hp, volume) / target_hp_norm > 0.5
        and norm_l1(matter_remainder, volume) / scale < 1.0e-3
    ):
        print(
            "  ATTRIBUTION: after the LLF contribution is reduced, the "
            "grid-scale applied heating follows the centered DG surface "
            "contribution and is transferred to matter by the subsequent "
            "implicit collisions. The element-integrated volume contribution "
            "is negligible in this step."
        )


def plot(
    ledger: Ledger,
    stages: list[int],
    output: Path,
    rmax_km: float | None = None,
) -> None:
    _, radius_cm, volume = rows(ledger, "FULL_STEP_MATTER_ACTUAL")
    radius_km = radius_cm / 1.0e5
    target = field(ledger, "FULL_STEP_MATTER_ACTUAL", "QH_exact")
    weighted_collision = field(
        ledger, "FULL_WEIGHTED_COLLISION_MATTER", "QH"
    )
    remainder_matter = field(ledger, "FULL_REMAINDER_MATTER", "QH_exact")
    weighted_streaming = field(ledger, "FULL_WEIGHTED_STREAMING", "QH")
    weighted_llf = weighted_stage_sum(
        ledger, stages, "STREAMING_LLF", "QH", "explicit"
    )
    weighted_centered = weighted_stage_sum(
        ledger, stages, "STREAMING_CENTERED", "QH", "explicit"
    )
    weighted_surface = weighted_stage_sum(
        ledger, stages, "STREAMING_SURFACE", "QH", "explicit"
    )
    weighted_volume = weighted_stage_sum(
        ledger, stages, "STREAMING_VOLUME", "QH", "explicit"
    )
    collision_conservation = field(ledger, "FULL_ALGEBRAIC_CLOSURE", "QH")

    if rmax_km is None:
        activity = np.abs(high_pass(target)) * volume
        if np.sum(activity) > 0.0:
            cumulative = np.cumsum(activity) / np.sum(activity)
            last = min(np.searchsorted(cumulative, 0.999) + 3, len(radius_km) - 1)
            rmax_km = float(radius_km[last])
        else:
            rmax_km = float(radius_km[-1])

    fig, axes = plt.subplots(2, 3, figsize=(19, 10), sharex=False)
    axes[0, 0].plot(
        radius_km, target, label=r"actual fluid $\Delta e_{int}/\Delta t$"
    )
    axes[0, 0].plot(radius_km, weighted_collision, label="weighted collision matter")
    axes[0, 0].plot(radius_km, remainder_matter, label="matter remainder")
    axes[0, 0].set_title("A. Exact full-step matter attribution")

    axes[0, 1].plot(radius_km, target, label="actual applied heating")
    axes[0, 1].plot(
        radius_km, weighted_streaming, label="weighted streaming RHS"
    )
    axes[0, 1].plot(radius_km, weighted_llf, label="LLF streaming contribution")
    axes[0, 1].set_title("B. Streaming pattern inherited by matter")

    axes[0, 2].plot(radius_km, weighted_surface, label="centered surface")
    axes[0, 2].plot(radius_km, weighted_volume, label="volume")
    axes[0, 2].plot(radius_km, weighted_llf, label="LLF (after scaling)")
    axes[0, 2].plot(
        radius_km,
        weighted_centered,
        "k--",
        lw=1.0,
        label="surface + volume",
    )
    axes[0, 2].set_title("C. Explicit streaming decomposition")

    colors = plt.cm.tab10(np.linspace(0.0, 1.0, max(len(stages), 2)))
    for color, stage in zip(colors, stages):
        weight = stage_weight(ledger, stage, "implicit")
        if weight == 0.0:
            continue
        axes[1, 0].plot(
            radius_km,
            weight * field(
                ledger, "COLLISION_MATTER_ACTUAL", "QH_exact", stage=stage
            ),
            color=color,
            label=f"stage {stage}, weight={weight:g}",
        )
    axes[1, 0].plot(radius_km, target, "k--", lw=1.2, label="stage sum")
    axes[1, 0].set_title("D. Actual implicit heating by stage")

    target_hp = high_pass(target)
    normalization = max(np.max(np.abs(target_hp)), np.finfo(float).tiny)
    axes[1, 1].plot(
        radius_km, target_hp / normalization, label="applied heating"
    )
    axes[1, 1].plot(
        radius_km,
        high_pass(weighted_llf) / normalization,
        label="LLF (after scaling)",
    )
    axes[1, 1].plot(
        radius_km,
        high_pass(weighted_surface) / normalization,
        label="centered surface",
    )
    axes[1, 1].plot(
        radius_km,
        high_pass(weighted_volume) / normalization,
        label="volume",
    )
    axes[1, 1].set_ylabel("high-pass / max |applied high-pass|")
    axes[1, 1].set_title("E. Grid-scale components")

    target_norm = max(norm_l1(target, volume), np.finfo(float).tiny)
    stream_norm = max(norm_l1(weighted_streaming, volume), np.finfo(float).tiny)
    limiter_norm = 0.0
    for stage in stages:
        limiter_norm += norm_l1(
            field(ledger, "PRE_STAGE_LIMITER_MATTER", "QH_exact", stage=stage),
            volume,
        )
        limiter_norm += norm_l1(
            field(ledger, "POST_STAGE_LIMITER_MATTER", "QH_exact", stage=stage),
            volume,
        )
    labels = [
        "LLF/\nstream",
        "surface/\nstream",
        "volume/\nstream",
        "collision\nnonconservation",
        "matter\nremainder",
        "fluid\nlimiters",
    ]
    fractions = np.asarray(
        [
            norm_l1(weighted_llf, volume) / stream_norm,
            norm_l1(weighted_surface, volume) / stream_norm,
            norm_l1(weighted_volume, volume) / stream_norm,
            norm_l1(collision_conservation, volume) / target_norm,
            norm_l1(remainder_matter, volume) / target_norm,
            limiter_norm / target_norm,
        ]
    )
    plot_floor = 1.0e-16
    plot_fractions = np.maximum(fractions, plot_floor)
    axes[1, 2].bar(np.arange(len(labels)), plot_fractions)
    axes[1, 2].set_xticks(np.arange(len(labels)), labels, rotation=20, ha="right")
    axes[1, 2].set_yscale("log")
    axes[1, 2].set_ylim(0.5 * plot_floor, 2.0)
    axes[1, 2].set_ylabel("relative volume-weighted $L^1$")
    axes[1, 2].set_title("F. Quantitative attribution")
    for index, (fraction, plotted) in enumerate(zip(fractions, plot_fractions)):
        label = "0" if fraction == 0.0 else f"{fraction:.1e}"
        axes[1, 2].text(index, plotted * 1.4, label, ha="center", fontsize=8)

    radial_axes = (
        axes[0, 0],
        axes[0, 1],
        axes[0, 2],
        axes[1, 0],
        axes[1, 1],
    )
    for axis in radial_axes:
        axis.axhline(0.0, color="0.5", lw=0.8)
        axis.set_xlim(radius_km[0], rmax_km)
        axis.set_xlabel("radius [km]")
        if axis is not axes[1, 1]:
            axis.set_ylabel(r"rate [erg cm$^{-3}$ s$^{-1}$]")
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8, loc="best")
    axes[1, 2].grid(axis="y", alpha=0.2)
    fig.suptitle(
        "Discrete attribution of grid-scale fluid-heating oscillations\n"
        "explicit DG streaming components and implicit collision response",
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_stage_decomposition(
    ledger: Ledger,
    stages: list[int],
    output: Path,
    rmax_km: float | None = None,
) -> None:
    """Plot the transport split and the following matter response at each stage."""

    _, radius_cm, _ = rows(ledger, "FULL_STEP_MATTER_ACTUAL")
    radius_km = radius_cm / 1.0e5
    active_stages = [
        stage
        for stage in stages
        if stage_weight(ledger, stage, "explicit") != 0.0
        or stage_weight(ledger, stage, "implicit") != 0.0
    ]
    if not active_stages:
        raise ValueError("no active explicit or implicit IMEX stages")

    if rmax_km is None:
        target = field(ledger, "FULL_STEP_MATTER_ACTUAL", "QH_exact")
        activity = np.abs(high_pass(target))
        if np.sum(activity) > 0.0:
            cumulative = np.cumsum(activity) / np.sum(activity)
            last = min(np.searchsorted(cumulative, 0.999) + 3, len(radius_km) - 1)
            rmax_km = float(radius_km[last])
        else:
            rmax_km = float(radius_km[-1])

    fig, axes = plt.subplots(
        len(active_stages),
        2,
        figsize=(15, 4.2 * len(active_stages)),
        squeeze=False,
        sharex=True,
    )
    for row, stage in enumerate(active_stages):
        w_ex = stage_weight(ledger, stage, "explicit")
        w_im = stage_weight(ledger, stage, "implicit")
        surface = field(ledger, "STREAMING_SURFACE", "QH", stage=stage)
        volume = field(ledger, "STREAMING_VOLUME", "QH", stage=stage)
        llf = field(ledger, "STREAMING_LLF", "QH", stage=stage)
        total = field(ledger, "STREAMING_TOTAL", "QH", stage=stage)
        matter = field(
            ledger, "COLLISION_MATTER_ACTUAL", "QH_exact", stage=stage
        )

        left = axes[row, 0]
        left.plot(radius_km, surface, label="centered surface")
        left.plot(radius_km, volume, label="volume")
        left.plot(radius_km, llf, label="LLF (after scaling)")
        left.plot(radius_km, total, "k--", lw=1.1, label="total streaming")
        left.set_title(f"Stage {stage}: explicit RHS (tableau weight {w_ex:g})")

        right = axes[row, 1]
        right.plot(radius_km, total, label="streaming RHS")
        right.plot(
            radius_km,
            matter,
            label="actual implicit matter exchange",
        )
        right.set_title(f"Stage {stage}: collision response (tableau weight {w_im:g})")

        for axis in (left, right):
            axis.axhline(0.0, color="0.5", lw=0.8)
            axis.set_xlim(radius_km[0], rmax_km)
            axis.set_ylabel(r"rate [erg cm$^{-3}$ s$^{-1}$]")
            axis.grid(alpha=0.2)
            axis.legend(fontsize=8, loc="best")

    for axis in axes[-1, :]:
        axis.set_xlabel("radius [km]")
    fig.suptitle(
        "Stage-resolved DG streaming decomposition and implicit matter response",
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("ledger", type=Path)
    parser.add_argument("--plot", type=Path)
    parser.add_argument(
        "--stage-plot",
        type=Path,
        help="write a stage-by-stage surface/volume/scaled-LLF decomposition",
    )
    parser.add_argument(
        "--reference-state",
        type=Path,
        help=(
            "element-averaged matter-diagnostics file used to normalize "
            "heating and lepton exchange; inferred from *.discrete_budget"
        ),
    )
    parser.add_argument(
        "--impact-plot",
        type=Path,
        help="plot grid-scale heating/lepton impact and timescales by origin",
    )
    parser.add_argument(
        "--rmax-km",
        type=float,
        help="right plot boundary; default encloses 99.9 percent of high-pass activity",
    )
    parser.add_argument(
        "--face-ledger",
        type=Path,
        help="optional .imex_faces file; inferred from the budget name by default",
    )
    parser.add_argument(
        "--identity-tolerance",
        type=float,
        default=5.0e-12,
        help="relative L1 tolerance for identities assembled from stored columns",
    )
    args = parser.parse_args()

    ledger = read_ledger(args.ledger)
    stages = validate(ledger, args.identity_tolerance)
    face_path = args.face_ledger
    if face_path is None and str(args.ledger).endswith(".discrete_budget"):
        face_path = Path(
            str(args.ledger)[: -len(".discrete_budget")] + ".imex_faces"
        )
    faces = read_face_ledger(face_path) if face_path is not None else None
    report(ledger, stages, faces)
    state_path = args.reference_state
    if state_path is None:
        inferred_state_path = infer_reference_state_path(args.ledger)
        if inferred_state_path.exists():
            state_path = inferred_state_path
    state = read_reference_state(state_path) if state_path is not None else None
    if state is not None:
        report_exchange_impact(ledger, stages, state)
    if args.plot is not None:
        plot(ledger, stages, args.plot, args.rmax_km)
    if args.stage_plot is not None:
        plot_stage_decomposition(
            ledger, stages, args.stage_plot, args.rmax_km
        )
    if args.impact_plot is not None:
        if state is None:
            raise ValueError("--impact-plot requires --reference-state")
        plot_exchange_impact(
            ledger, stages, state, args.impact_plot, args.rmax_km
        )


if __name__ == "__main__":
    main()
