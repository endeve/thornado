#!/usr/bin/env python3

"""Audit DG face jumps and the dissipative Lax--Friedrichs flux.

In addition to plotting the independently reconstructed face states, this
script can compare the divergence of their LLF dissipative flux with the LLF
contribution recorded inside the streaming operator.  That comparison checks
the face ordering, spherical face areas, divergence sign, and the jump-flux
formula without relying on the gray collision diagnostic.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_ledger(path: Path) -> np.ndarray:
    header = None
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("# stage limiter_applied"):
                header = line[2:].split()
                break
    if header is None:
        raise RuntimeError(f"Could not find column header in {path}")
    return np.genfromtxt(path, names=header, comments="#")


def read_streaming_ledger(path: Path) -> np.ndarray:
    header = None
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("# stage explicit_rhs_used"):
                header = line[2:].split()
                break
    if header is None:
        raise RuntimeError(f"Could not find streaming column header in {path}")
    return np.genfromtxt(path, names=header, comments="#")


def symlog(axis: plt.Axes, arrays: tuple[np.ndarray, ...]) -> None:
    maximum = max((np.nanmax(np.abs(array)) for array in arrays), default=0.0)
    if np.isfinite(maximum) and maximum > 0.0:
        axis.set_yscale("symlog", linthresh=max(maximum * 1.0e-8, 1.0e-300))
    axis.axhline(0.0, color="0.7", linewidth=0.7)
    axis.grid(alpha=0.2)


def select(data: np.ndarray, stage: int, species: int, limited: int) -> np.ndarray:
    selected = data[
        (data["stage"] == stage)
        & (data["species"] == species)
        & (data["limiter_applied"] == limited)
    ]
    return np.sort(selected, order="r_cm")


def plot(data: np.ndarray, species: int, output: Path) -> None:
    stages = np.unique(data["stage"]).astype(int)
    figure, axes = plt.subplots(
        7,
        len(stages),
        figsize=(5.2 * len(stages), 22.0),
        sharex="col",
        squeeze=False,
    )

    for column, stage in enumerate(stages):
        pre = select(data, stage, species, 0)
        post = select(data, stage, species, 1)
        radius = post["r_cm"] / 1.0e5
        used = bool(np.max(post["explicit_rhs_used"]))

        axis = axes[0, column]
        axis.plot(radius, pre["energy_jump_relative_L1"], label="pre-limiter")
        axis.plot(radius, post["energy_jump_relative_L1"], label="post-limiter")
        positive = np.concatenate(
            (
                pre["energy_jump_relative_L1"][pre["energy_jump_relative_L1"] > 0.0],
                post["energy_jump_relative_L1"][post["energy_jump_relative_L1"] > 0.0],
            )
        )
        if positive.size:
            axis.set_yscale("log")
        axis.grid(alpha=0.2)
        axis.set_title(f"stage {stage}; explicit RHS used={used}")
        axis.set_ylabel("relative spectral L1 jump")
        axis.legend(loc="best", fontsize="small")

        available = post["equilibrium_split_available"] == 1

        axis = axes[1, column]
        actual_relative = post["energy_jump_relative_L1"]
        equilibrium_relative = post["equilibrium_jump_relative_to_actual_L1"]
        nonequilibrium_relative = post[
            "nonequilibrium_jump_relative_to_actual_L1"
        ]
        axis.plot(radius[available], actual_relative[available], label="actual")
        axis.plot(
            radius[available],
            equilibrium_relative[available],
            label="limited D_0",
        )
        axis.plot(
            radius[available],
            nonequilibrium_relative[available],
            label="actual - limited D_0",
        )
        positive = np.concatenate(
            (
                actual_relative[available][actual_relative[available] > 0.0],
                equilibrium_relative[available][
                    equilibrium_relative[available] > 0.0
                ],
                nonequilibrium_relative[available][
                    nonequilibrium_relative[available] > 0.0
                ],
            )
        )
        if positive.size:
            axis.set_yscale("log")
        axis.grid(alpha=0.2)
        axis.set_ylabel("post-limiter jump / actual energy")
        axis.legend(loc="best", fontsize="small")

        axis = axes[2, column]
        area = 4.0 * np.pi * post["r_cm"] ** 2
        penalty_pre = area * pre["LLF_penalty_flux_erg_cm2_s"]
        penalty_post = area * post["LLF_penalty_flux_erg_cm2_s"]
        axis.plot(radius, penalty_pre, label="pre-limiter")
        axis.plot(radius, penalty_post, label="post-limiter")
        symlog(axis, (penalty_pre, penalty_post))
        axis.set_ylabel(r"LLF penalty luminosity [erg s$^{-1}$]")

        axis = axes[3, column]
        actual_penalty = area * post["LLF_penalty_flux_erg_cm2_s"]
        equilibrium_penalty = area * post[
            "equilibrium_LLF_penalty_flux_erg_cm2_s"
        ]
        nonequilibrium_penalty = area * post[
            "nonequilibrium_LLF_penalty_flux_erg_cm2_s"
        ]
        axis.plot(
            radius[available], actual_penalty[available], label="actual"
        )
        axis.plot(
            radius[available],
            equilibrium_penalty[available],
            label="limited D_0",
        )
        axis.plot(
            radius[available],
            nonequilibrium_penalty[available],
            label="actual - limited D_0",
        )
        symlog(
            axis,
            (
                actual_penalty[available],
                equilibrium_penalty[available],
                nonequilibrium_penalty[available],
            ),
        )
        axis.set_ylabel(r"post-limiter LLF decomposition [erg s$^{-1}$]")
        axis.legend(loc="best", fontsize="small")

        axis = axes[4, column]
        centered = area * post["centered_flux_erg_cm2_s"]
        penalty = area * post["LLF_penalty_flux_erg_cm2_s"]
        total = area * post["total_flux_erg_cm2_s"]
        axis.plot(radius, centered, label="centered")
        axis.plot(radius, penalty, label="LLF penalty")
        axis.plot(radius, total, label="total", color="black", linewidth=1.2)
        symlog(axis, (centered, penalty, total))
        axis.set_ylabel(r"post-limiter luminosity [erg s$^{-1}$]")
        axis.legend(loc="best", fontsize="small")

        axis = axes[5, column]
        signed_actual = post["energy_jump_signed_erg_cm3"]
        signed_equilibrium = post["equilibrium_energy_jump_signed_erg_cm3"]
        signed_nonequilibrium = post[
            "nonequilibrium_energy_jump_signed_erg_cm3"
        ]
        axis.plot(radius[available], signed_actual[available], label="actual")
        axis.plot(
            radius[available],
            signed_equilibrium[available],
            label="limited D_0",
        )
        axis.plot(
            radius[available],
            signed_nonequilibrium[available],
            label="actual - limited D_0",
        )
        symlog(
            axis,
            (
                signed_actual[available],
                signed_equilibrium[available],
                signed_nonequilibrium[available],
            ),
        )
        axis.set_ylabel(r"post-limiter signed jump [erg cm$^{-3}$]")

        axis = axes[6, column]
        tau_iso = np.sqrt(
            np.abs(post["tau_iso_left"] * post["tau_iso_right"])
        )
        tau_noniso = np.sqrt(
            np.abs(post["tau_noniso_left"] * post["tau_noniso_right"])
        )
        tau_total = np.sqrt(
            np.abs(post["tau_total_left"] * post["tau_total_right"])
        )
        axis.plot(radius, tau_iso, label="Iso")
        axis.plot(radius, tau_noniso, label="non-Iso")
        axis.plot(radius, tau_total, label="total", color="black")
        axis.axhline(1.0, color="0.5", linestyle="--", linewidth=0.8)
        if np.any(tau_total > 0.0):
            axis.set_yscale("log")
        axis.grid(alpha=0.2)
        axis.set_ylabel(r"final-state element optical depth")
        axis.set_xlabel("face radius [km]")
        axis.legend(loc="best", fontsize="small")

    figure.suptitle(f"DG radial-face diagnostic; species={species}")
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def streaming_selection(
    data: np.ndarray, stage: int, species: int
) -> np.ndarray:
    selected = data[
        (data["stage"] == stage) & (data["species"] == species)
    ]
    return np.sort(selected, order="element")


def relative_l1(error: np.ndarray, reference: np.ndarray, volume: np.ndarray) -> float:
    denominator = np.sum(np.abs(reference) * volume)
    if denominator == 0.0:
        return 0.0 if np.all(error == 0.0) else np.inf
    return float(np.sum(np.abs(error) * volume) / denominator)


def report_streaming_closure(
    faces: np.ndarray,
    streaming: np.ndarray,
    stage: int,
    species: int,
) -> None:
    post = select(faces, stage, species, 1)
    cells = streaming_selection(streaming, stage, species)

    if post.size != cells.size + 1:
        raise RuntimeError(
            f"stage={stage}, species={species}: expected one more face than "
            f"elements, found {post.size} faces and {cells.size} elements"
        )

    expected_faces = np.arange(post.size)
    expected_elements = np.arange(1, cells.size + 1)
    if not np.array_equal(post["face"].astype(int), expected_faces):
        raise RuntimeError(
            f"stage={stage}, species={species}: noncontiguous face ordering"
        )
    if not np.array_equal(cells["element"].astype(int), expected_elements):
        raise RuntimeError(
            f"stage={stage}, species={species}: noncontiguous element ordering"
        )

    luminosity = (
        4.0
        * np.pi
        * post["r_cm"] ** 2
        * post["LLF_penalty_flux_erg_cm2_s"]
    )
    reconstructed_rhs = -(
        luminosity[1:] - luminosity[:-1]
    ) / cells["volume_cm3"]
    recorded_rhs = cells["QH_LLF_penalty_erg_cm3_s"]
    difference = reconstructed_rhs - recorded_rhs

    l1_closure = relative_l1(difference, recorded_rhs, cells["volume_cm3"])
    maximum_scale = max(float(np.max(np.abs(recorded_rhs))), np.finfo(float).tiny)
    linf_closure = float(np.max(np.abs(difference)) / maximum_scale)

    print(
        "    face-flux divergence vs streaming LLF RHS: "
        f"relative L1={l1_closure:.6e}, relative Linf={linf_closure:.6e}"
    )


def report(
    data: np.ndarray,
    species: int,
    streaming: np.ndarray | None = None,
) -> None:
    for stage in np.unique(data["stage"]).astype(int):
        pre = select(data, stage, species, 0)
        post = select(data, stage, species, 1)
        pre_jump = pre["energy_jump_relative_L1"]
        post_jump = post["energy_jump_relative_L1"]
        area = 4.0 * np.pi * post["r_cm"] ** 2
        pre_penalty = area * pre["LLF_penalty_flux_erg_cm2_s"]
        post_penalty = area * post["LLF_penalty_flux_erg_cm2_s"]
        jump_scale = np.maximum(pre_jump, np.finfo(float).tiny)
        penalty_change = post_penalty - pre_penalty
        index = np.argmax(np.abs(penalty_change))
        available = post["equilibrium_split_available"] == 1
        equilibrium_penalty = area * post[
            "equilibrium_LLF_penalty_flux_erg_cm2_s"
        ]
        nonequilibrium_penalty = area * post[
            "nonequilibrium_LLF_penalty_flux_erg_cm2_s"
        ]
        decomposition_closure = area * post[
            "LLF_decomposition_closure_erg_cm2_s"
        ]
        actual_norm = np.linalg.norm(post_penalty[available])
        if actual_norm > 0.0:
            equilibrium_ratio = (
                np.linalg.norm(equilibrium_penalty[available]) / actual_norm
            )
            nonequilibrium_ratio = (
                np.linalg.norm(nonequilibrium_penalty[available]) / actual_norm
            )
            closure_ratio = (
                np.linalg.norm(decomposition_closure[available]) / actual_norm
            )
        else:
            equilibrium_ratio = 0.0
            nonequilibrium_ratio = 0.0
            closure_ratio = 0.0

        interior = available
        interior_indices = np.flatnonzero(interior)
        if interior_indices.size:
            local_index = int(
                np.argmax(np.abs(post_penalty[interior_indices]))
            )
            maximum_index = int(interior_indices[local_index])
            tau_left = post["tau_total_left"][maximum_index]
            tau_right = post["tau_total_right"][maximum_index]
            tau_face = np.sqrt(abs(tau_left * tau_right))
            actual_at_maximum = post_penalty[maximum_index]
            if actual_at_maximum != 0.0:
                equilibrium_fraction = (
                    equilibrium_penalty[maximum_index] / actual_at_maximum
                )
                nonequilibrium_fraction = (
                    nonequilibrium_penalty[maximum_index] / actual_at_maximum
                )
            else:
                equilibrium_fraction = 0.0
                nonequilibrium_fraction = 0.0
        else:
            maximum_index = 0
            tau_face = 0.0
            equilibrium_fraction = 0.0
            nonequilibrium_fraction = 0.0

        centered_luminosity = area * post["centered_flux_erg_cm2_s"]
        centered_norm = np.linalg.norm(centered_luminosity[interior])
        dissipative_norm = np.linalg.norm(post_penalty[interior])
        if centered_norm > 0.0:
            dissipative_to_centered = dissipative_norm / centered_norm
        else:
            dissipative_to_centered = np.inf if dissipative_norm > 0.0 else 0.0

        print(
            f"stage={stage} used={bool(np.max(post['explicit_rhs_used']))} "
            f"species={species}: max relative jump pre={np.max(pre_jump):.6e} "
            f"post={np.max(post_jump):.6e}; "
            f"max post/pre={np.max(post_jump / jump_scale):.6e}; "
            f"max limiter-induced |delta LLF luminosity|="
            f"{np.abs(penalty_change[index]):.6e} at "
            f"r={post['r_cm'][index] / 1.0e5:.6e} km; "
            f"post-limiter L2 ratios limited-D0/actual="
            f"{equilibrium_ratio:.6e}, nonequilibrium/actual="
            f"{nonequilibrium_ratio:.6e}, closure/actual="
            f"{closure_ratio:.6e}"
        )
        print(
            "    strongest interior LLF dissipative luminosity: "
            f"r={post['r_cm'][maximum_index] / 1.0e5:.6e} km, "
            f"tau_face={tau_face:.6e}, "
            f"relative_jump={post['energy_jump_relative_L1'][maximum_index]:.6e}, "
            f"equilibrium/actual={equilibrium_fraction:.6e}, "
            f"nonequilibrium/actual={nonequilibrium_fraction:.6e}"
        )
        print(
            "    face L2(LLF dissipative luminosity) / "
            f"L2(centered luminosity)={dissipative_to_centered:.6e}"
        )
        if streaming is not None:
            report_streaming_closure(data, streaming, stage, species)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "diagnostic",
        type=Path,
        help="main matter-diagnostic file or its .imex_faces ledger",
    )
    parser.add_argument("--species", type=int, default=1)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--streaming-ledger",
        type=Path,
        help=(
            "optional .imex_streaming ledger; by default the corresponding "
            "file is used when it exists"
        ),
    )
    args = parser.parse_args()

    ledger = args.diagnostic
    if not str(ledger).endswith(".imex_faces"):
        ledger = Path(str(ledger) + ".imex_faces")
    data = read_ledger(ledger)

    streaming_path = args.streaming_ledger
    if streaming_path is None:
        streaming_path = Path(
            str(ledger)[: -len(".imex_faces")] + ".imex_streaming"
        )
    streaming = (
        read_streaming_ledger(streaming_path)
        if streaming_path.exists()
        else None
    )

    species_values = np.unique(data["species"]).astype(int)
    if args.species not in species_values:
        raise ValueError(
            f"species={args.species} is unavailable; choices are "
            f"{species_values.tolist()}"
        )

    output = args.output
    if output is None:
        output = Path(str(ledger) + f"_s{args.species}.png")

    plot(data, args.species, output)
    report(data, args.species, streaming)


if __name__ == "__main__":
    main()
