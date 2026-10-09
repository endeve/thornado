#!/usr/bin/env python3
"""Audit the nodal DG balance in the absorption--diffusion test."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.legendre import leggauss, legvander


D_FIELDS = (
    "streaming_D",
    "streaming_central_D",
    "streaming_penalty_D",
    "streaming_volume_D",
    "streaming_face_dn_D",
    "streaming_face_up_D",
    "absorption_D",
    "residual_D",
)
I_FIELDS = (
    "streaming_I1",
    "streaming_central_I1",
    "streaming_penalty_I1",
    "streaming_volume_I1",
    "streaming_face_dn_I1",
    "streaming_face_up_I1",
    "collision_I1",
    "residual_I1",
)


def relative_max(value: np.ndarray, scale: np.ndarray) -> float:
    return float(np.max(np.abs(value)) / max(np.max(np.abs(scale)), np.finfo(float).tiny))


def parity_error(values: np.ndarray, sign: float) -> float:
    reflected = values[::-1, ::-1]
    return relative_max(values - sign * reflected, values)


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(values**2)))


def correlation(first: np.ndarray, second: np.ndarray) -> float:
    denominator = np.linalg.norm(first.ravel()) * np.linalg.norm(second.ravel())
    if denominator == 0.0:
        return 0.0
    return float(np.dot(first.ravel(), second.ravel()) / denominator)


def analyze(path: Path, output_dir: Path) -> None:
    data = np.genfromtxt(path, names=True)
    elements = np.asarray(data["element"], dtype=int)
    nodes = np.asarray(data["node"], dtype=int)
    n_elements = int(elements.max())
    n_nodes = int(nodes.max())
    if data.size != n_elements * n_nodes:
        raise ValueError(f"{path}: incomplete element/node tensor product")

    order = np.lexsort((nodes, elements))
    x = np.asarray(data["x_km"])[order].reshape(n_elements, n_nodes)
    gl_nodes, gl_weights = leggauss(n_nodes)
    normalized_weights = 0.5 * gl_weights
    vandermonde = legvander(gl_nodes, n_nodes - 1)

    values: dict[str, np.ndarray] = {}
    modes: dict[str, np.ndarray] = {}
    for name in D_FIELDS + I_FIELDS:
        values[name] = np.asarray(data[name])[order].reshape(n_elements, n_nodes)
        modes[name] = np.linalg.solve(vandermonde, values[name].T).T

    decomposition_d = np.asarray(data["decomposition_D"])[order].reshape(
        n_elements, n_nodes
    )
    decomposition_i = np.asarray(data["decomposition_I1"])[order].reshape(
        n_elements, n_nodes
    )

    # On the uniform periodic test mesh, this is the physical element width.
    centers = x @ normalized_weights
    dx = float(np.median(np.diff(centers))) if n_elements > 1 else 1.0
    integral = lambda a: float(dx * np.sum(a @ normalized_weights))

    face_dn_integral_d = dx * (values["streaming_face_dn_D"] @ normalized_weights)
    face_up_integral_d = dx * (values["streaming_face_up_D"] @ normalized_weights)
    face_dn_integral_i = dx * (values["streaming_face_dn_I1"] @ normalized_weights)
    face_up_integral_i = dx * (values["streaming_face_up_I1"] @ normalized_weights)
    interface_d = face_up_integral_d + np.roll(face_dn_integral_d, -1)
    interface_i = face_up_integral_i + np.roll(face_dn_integral_i, -1)

    def interface_relative(mismatch: np.ndarray, down: np.ndarray, up: np.ndarray) -> float:
        scale = max(np.max(np.abs(down)), np.max(np.abs(up)), np.finfo(float).tiny)
        return float(np.max(np.abs(mismatch)) / scale)

    modal_rows = []
    for i_element in range(n_elements):
        row = [i_element + 1, centers[i_element]]
        for name in D_FIELDS + I_FIELDS:
            row.extend(modes[name][i_element])
        modal_rows.append(row)
    modal_header = ["element", "x_center_km"]
    for name in D_FIELDS + I_FIELDS:
        modal_header.extend(f"{name}_mode_{mode}" for mode in range(n_nodes))
    modal_path = output_dir / f"{path.stem}.modal.dat"
    np.savetxt(modal_path, np.asarray(modal_rows), header=" ".join(modal_header))

    report_path = output_dir / f"{path.stem}.audit.txt"
    with report_path.open("w") as report:
        def emit(line: str = "") -> None:
            print(line)
            print(line, file=report)

        emit(f"Operator audit: {path}")
        emit(f"  elements/nodes                  = {n_elements}/{n_nodes}")
        emit(
            "  D decomposition max abs./rel.   = "
            f"{np.max(np.abs(decomposition_d)):.6e} / "
            f"{relative_max(decomposition_d, values['streaming_D']):.6e}"
        )
        emit(
            "  I1 decomposition max abs./rel.  = "
            f"{np.max(np.abs(decomposition_i)):.6e} / "
            f"{relative_max(decomposition_i, values['streaming_I1']):.6e}"
        )
        emit(
            "  D interface conservation rel.   = "
            f"{interface_relative(interface_d, face_dn_integral_d, face_up_integral_d):.6e}"
        )
        emit(
            "  I1 interface conservation rel.  = "
            f"{interface_relative(interface_i, face_dn_integral_i, face_up_integral_i):.6e}"
        )
        emit(f"  integral streaming D            = {integral(values['streaming_D']):.16e}")
        emit(f"  integral streaming I1           = {integral(values['streaming_I1']):.16e}")
        emit(f"  integral absorption D           = {integral(values['absorption_D']):.16e}")
        emit(f"  integral residual D             = {integral(values['residual_D']):.16e}")
        emit(f"  integral collision I1           = {integral(values['collision_I1']):.16e}")
        emit(f"  integral residual I1            = {integral(values['residual_I1']):.16e}")
        emit(
            "  RMS D residual / LLF term       = "
            f"{rms(values['residual_D']):.6e} / "
            f"{rms(values['streaming_penalty_D']):.6e}"
        )
        emit(
            "  corr(D residual, LLF term)      = "
            f"{correlation(values['residual_D'], values['streaming_penalty_D']):.12f}"
        )
        nonpenalty_d = values["streaming_central_D"] + values["absorption_D"]
        emit(
            "  RMS central+absorption D        = "
            f"{rms(nonpenalty_d):.6e}"
        )
        emit()
        emit("  Reflection parity errors (even / odd):")
        for name in D_FIELDS + I_FIELDS:
            emit(
                f"    {name:24s} {parity_error(values[name], 1.0):.6e} / "
                f"{parity_error(values[name], -1.0):.6e}"
            )
        emit("  Reflected lower/upper face pairing (same / opposite):")
        for suffix in ("D", "I1"):
            lower = values[f"streaming_face_dn_{suffix}"]
            reflected_upper = values[f"streaming_face_up_{suffix}"][::-1, ::-1]
            emit(
                f"    {suffix:24s} "
                f"{relative_max(lower-reflected_upper, lower):.6e} / "
                f"{relative_max(lower+reflected_upper, lower):.6e}"
            )
        emit()
        emit("  RMS Legendre coefficients by local mode:")
        for name in D_FIELDS + I_FIELDS:
            modal_rms = np.sqrt(np.mean(modes[name] ** 2, axis=0))
            emit(
                f"    {name:24s} "
                + " ".join(f"{value:.6e}" for value in modal_rms)
            )

    fig, axes = plt.subplots(3, 2, figsize=(13.0, 10.0), constrained_layout=True)
    for name in ("streaming_D", "streaming_central_D", "streaming_penalty_D",
                 "streaming_volume_D", "streaming_face_dn_D",
                 "streaming_face_up_D", "absorption_D", "residual_D"):
        axes[0, 0].plot(x.ravel(), values[name].ravel(), marker=".", ms=2, lw=0.8,
                        label=name)
    axes[0, 0].set_title("Density equation: nodal terms")
    axes[0, 0].legend(fontsize=7, ncol=2)

    for name in ("streaming_I1", "streaming_central_I1", "streaming_penalty_I1",
                 "streaming_volume_I1", "streaming_face_dn_I1",
                 "streaming_face_up_I1", "collision_I1", "residual_I1"):
        axes[0, 1].plot(x.ravel(), values[name].ravel(), marker=".", ms=2, lw=0.8,
                        label=name)
    axes[0, 1].set_title("Flux equation: nodal terms")
    axes[0, 1].legend(fontsize=7, ncol=2)

    for mode in range(n_nodes):
        axes[1, 0].plot(centers, modes["residual_D"][:, mode], marker="o", ms=3,
                        label=f"mode {mode}")
        axes[1, 1].plot(centers, modes["residual_I1"][:, mode], marker="o", ms=3,
                        label=f"mode {mode}")
    axes[1, 0].set_title("Density residual: element Legendre modes")
    axes[1, 1].set_title("Flux residual: element Legendre modes")
    axes[1, 0].legend(fontsize=8)
    axes[1, 1].legend(fontsize=8)

    axes[2, 0].plot(np.arange(1, n_elements + 1), interface_d, marker="o", ms=3)
    axes[2, 0].set_title("D face cancellation at each interface")
    axes[2, 0].set_xlabel("interface")
    axes[2, 0].set_ylabel("upper(left) + lower(right)")
    axes[2, 1].plot(np.arange(1, n_elements + 1), interface_i, marker="o", ms=3)
    axes[2, 1].set_title("I1 face cancellation at each interface")
    axes[2, 1].set_xlabel("interface")
    axes[2, 1].set_ylabel("upper(left) + lower(right)")

    for axis in axes.flat:
        axis.grid(alpha=0.2)
    for axis in axes[:2, :].flat:
        axis.set_xlabel("x [km]")
    fig.suptitle(path.parent.name + ": " + path.stem)
    fig.savefig(output_dir / f"{path.stem}.audit.pdf", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    for path in args.files:
        output_dir = args.output_dir or path.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        analyze(path, output_dir)


if __name__ == "__main__":
    main()
