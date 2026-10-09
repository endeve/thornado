#!/usr/bin/env python3
"""Fourier and frozen-profile audit of paired one-sided DG traces.

The periodic Cartesian calculation compares Thornado's centered physical
fluxes with the two alternating choices

    I_hat = I^- , K_hat = K^+,
    I_hat = I^+ , K_hat = K^-,

while retaining Thornado's LLF jump terms.  The frozen-profile calculation
changes only those physical face fluxes in a saved zero-velocity
DeleptonizationWave state.  Its pressure trace uses the same Minerbo closure
as the executable, so the comparison remains defined through the
free-streaming region; it does not impose K=D/3 there.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from analyze_frozen_transient_operator import (
    column_map,
    lagrange_values,
    read_table,
    select_block,
)
from centered_p1_modal_audit import modal_symbol
from frozen_matter_schur_audit import build_operators


TRACE_PAIRS = ("centered", "left-right", "right-left")
TRACE_LABELS = {
    "centered": "centered",
    "left-right": r"$\widehat I=I^{-}$, $\widehat K=K^{+}$",
    "right-left": r"$\widehat I=I^{+}$, $\widehat K=K^{-}$",
}


def minerbo_eddington_factor(flux_factor: np.ndarray) -> np.ndarray:
    """MOMENT_CLOSURE_MINERBO specialization used by this sandbox."""

    f = np.clip(flux_factor, 0.0, 1.0)
    closure_polynomial = 0.6 * f**2 * (1.0 - f / 3.0 + f**2)
    return 1.0 / 3.0 + (2.0 / 3.0) * closure_polynomial


def minerbo_eddington_derivative(flux_factor: np.ndarray) -> np.ndarray:
    """Derivative of the sandbox Minerbo Eddington factor with respect to f."""

    f = np.clip(flux_factor, 0.0, 1.0)
    return 0.4 * (2.0 * f - f**2 + 4.0 * f**3)


def fourier_audit(n_cells: int) -> dict[str, object]:
    """Return exact block-circulant p=1 symbols for each trace pair."""

    width = 1.0 / n_cells
    # modal_symbol applies the operator on a finite periodic grid, so only
    # Fourier phases compatible with that grid are exact symbols.
    theta = 2.0 * np.pi * np.arange(n_cells // 2 + 1) / n_cells
    operators = {
        name: build_operators(
            n_cells,
            2,
            mean_optical_depth=1.0,
            opacity_contrast=0.0,
            source_aware=False,
            trace_pair=name,
        )
        for name in TRACE_PAIRS
    }
    result: dict[str, object] = {
        "width": width,
        "theta": theta,
        "operators": operators,
    }
    background_flux_factor = np.linspace(0.0, 1.0, 101)
    result["background_flux_factor"] = background_flux_factor
    for name, op in operators.items():
        diffusion_eigenvalues = np.empty((theta.size, 2), dtype=complex)
        hyperbolic_eigenvalues = np.empty((theta.size, 4), dtype=complex)
        symbols: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []
        for index, angle in enumerate(theta):
            a = modal_symbol(op.density_from_flux, angle, n_cells)
            b = modal_symbol(op.flux_from_density, angle, n_cells)
            jd = modal_symbol(op.density_penalty, angle, n_cells)
            ji = modal_symbol(op.flux_penalty, angle, n_cells)
            symbols.append((a, b, jd, ji))
            diffusion_eigenvalues[index] = np.linalg.eigvals(a @ b)
            hyperbolic_eigenvalues[index] = np.linalg.eigvals(
                np.block([[jd, a], [b, ji]])
            )
        # Sort the diffusion branches by magnitude at every wavenumber.  The
        # first is the physical branch near theta=0 and the second is the
        # element-internal branch.
        order = np.argsort(np.abs(diffusion_eigenvalues), axis=1)
        diffusion_eigenvalues = np.take_along_axis(
            diffusion_eigenvalues, order, axis=1
        )
        m1_maximum_real = np.empty(
            (background_flux_factor.size, theta.size)
        )
        for i_flux, flux_factor in enumerate(background_flux_factor):
            eddington = float(minerbo_eddington_factor(flux_factor))
            derivative = float(minerbo_eddington_derivative(flux_factor))
            d_k_d_d = eddington - flux_factor * derivative
            d_k_d_i = derivative
            for i_theta, (a, b, jd, ji) in enumerate(symbols):
                # b is the paired derivative of D/3.  Multiplication by
                # three therefore gives the derivative used by the
                # linearized pressure flux
                #   delta K = (dK/dD) delta D + (dK/dI) delta I.
                linearized_m1 = np.block(
                    [
                        [jd, a],
                        [
                            3.0 * d_k_d_d * b,
                            ji + 3.0 * d_k_d_i * b,
                        ],
                    ]
                )
                m1_maximum_real[i_flux, i_theta] = np.max(
                    np.linalg.eigvals(linearized_m1).real
                )
        result[name] = {
            "diffusion": diffusion_eigenvalues,
            "hyperbolic": hyperbolic_eigenvalues,
            "m1_maximum_real": m1_maximum_real,
        }
    return result


def frozen_profile_audit(
    path: Path,
    species: int,
    energy: int,
) -> dict[str, np.ndarray | dict[str, np.ndarray] | float]:
    """Apply paired-trace face corrections to one saved production RHS."""

    metadata, columns, table = read_table(path)
    data = column_map(columns, table)
    block = select_block(data, species, energy)
    order = np.lexsort((block["spatial_node"], block["element"]))
    block = {name: field[order] for name, field in block.items()}

    n_cells = int(block["element"].max())
    n_nodes = int(block["spatial_node"].max())
    size = n_cells * n_nodes
    if size != order.size:
        raise ValueError("selected block is not a complete spatial grid")

    def shaped(name: str) -> np.ndarray:
        return block[name].reshape(n_cells, n_nodes)

    radius = shaped("r_code")
    widths = shaped("dr_code")[:, 0]
    weights = shaped("wx")[0]
    centers = np.mean(radius, axis=1)
    nodes = (radius[0] - centers[0]) / widths[0]
    lower = lagrange_values(nodes, -0.5)
    upper = lagrange_values(nodes, +0.5)

    h1 = shaped("h1")
    h2 = shaped("h2")
    h3 = shaped("h3")
    sqrtgm = shaped("sqrtgm")
    sqrtgm_face = np.empty(n_cells + 1)
    sqrtgm_face[0] = (
        (h1[0] @ lower) * (h2[0] @ lower) * (h3[0] @ lower)
    )
    sqrtgm_face[-1] = (
        (h1[-1] @ upper) * (h2[-1] @ upper) * (h3[-1] @ upper)
    )
    for face in range(1, n_cells):
        h1_face = 0.5 * (h1[face - 1] @ upper + h1[face] @ lower)
        h2_face = 0.5 * (h2[face - 1] @ upper + h2[face] @ lower)
        h3_face = 0.5 * (h3[face - 1] @ upper + h3[face] @ lower)
        sqrtgm_face[face] = h1_face * h2_face * h3_face

    density = shaped("D_before")
    first = shaped("I1_before")
    gm11 = shaped("gm11")
    flux_factor = np.sqrt(gm11) * np.abs(first) / np.maximum(
        density, np.finfo(float).tiny
    )
    flux_factor = np.clip(flux_factor, 0.0, 1.0)
    d_lower, d_upper = density @ lower, density @ upper
    i_lower, i_upper = first @ lower, first @ upper

    # Thornado traces the conserved state first and evaluates the nonlinear
    # physical flux at that one-sided face state.  Do not interpolate nodal
    # pressure values: those operations do not commute outside diffusion.
    f_lower = np.clip(
        np.abs(i_lower) / np.maximum(d_lower, np.finfo(float).tiny),
        0.0,
        1.0,
    )
    f_upper = np.clip(
        np.abs(i_upper) / np.maximum(d_upper, np.finfo(float).tiny),
        0.0,
        1.0,
    )
    k_lower = minerbo_eddington_factor(f_lower) * d_lower
    k_upper = minerbo_eddington_factor(f_upper) * d_upper

    def surface_response(face_delta: np.ndarray) -> np.ndarray:
        response = np.empty((n_cells, n_nodes))
        for cell in range(n_cells):
            raw = (
                lower * sqrtgm_face[cell] * face_delta[cell]
                - upper * sqrtgm_face[cell + 1] * face_delta[cell + 1]
            )
            response[cell] = raw / (
                weights * widths[cell] * sqrtgm[cell]
            )
        return response

    production_density = shaped("stream_total_N")
    production_first = shaped("stream_total_G1")
    variants: dict[str, np.ndarray] = {}
    for name in TRACE_PAIRS[1:]:
        delta_density_flux = np.zeros(n_cells + 1)
        delta_first_flux = np.zeros(n_cells + 1)
        for face in range(1, n_cells):
            centered_i = 0.5 * (i_upper[face - 1] + i_lower[face])
            centered_k = 0.5 * (k_upper[face - 1] + k_lower[face])
            if name == "left-right":
                paired_i = i_upper[face - 1]
                paired_k = k_lower[face]
            else:
                paired_i = i_lower[face]
                paired_k = k_upper[face - 1]
            delta_density_flux[face] = paired_i - centered_i
            delta_first_flux[face] = paired_k - centered_k
        variants[f"density_{name}"] = (
            production_density + surface_response(delta_density_flux)
        )
        variants[f"first_{name}"] = (
            production_first + surface_response(delta_first_flux)
        )

    opacity = shaped("kappa_transport")
    optical_depth = np.maximum(opacity * widths[:, None], 0.0)
    centimeter_code = metadata["centimeter_code"]
    radius_km = radius / centimeter_code / 1.0e5
    cell_weights = weights[None, :] * widths[:, None] * sqrtgm

    def region_norm(values: np.ndarray, mask: np.ndarray) -> float:
        if not np.any(mask):
            return float("nan")
        return float(np.sqrt(np.sum(cell_weights[mask] * values[mask] ** 2)))

    regions = {
        "opaque (tau>=10)": optical_depth >= 10.0,
        "transition (0.1<tau<10)": (
            (optical_depth > 0.1) & (optical_depth < 10.0)
        ),
        "transparent (tau<=0.1)": optical_depth <= 0.1,
    }
    print(
        f"Frozen DeleptonizationWave trace audit: species={species}, "
        f"flat energy={energy}, r={radius_km.min():.3f}--"
        f"{radius_km.max():.3f} km"
    )
    for label, mask in regions.items():
        density_scale = max(
            region_norm(production_density, mask), np.finfo(float).tiny
        )
        first_scale = max(
            region_norm(production_first, mask), np.finfo(float).tiny
        )
        print(f"  {label}")
        for name in TRACE_PAIRS[1:]:
            density_change = variants[f"density_{name}"] - production_density
            first_change = variants[f"first_{name}"] - production_first
            print(
                f"    {name:11s}: "
                f"||delta R_D||/||R_D||={region_norm(density_change, mask) / density_scale:.6e}, "
                f"||delta R_I||/||R_I||={region_norm(first_change, mask) / first_scale:.6e}"
            )

    free_streaming_tail = (flux_factor >= 0.9) & (optical_depth <= 0.1)
    all_nodes = np.ones_like(free_streaming_tail, dtype=bool)
    if np.any(free_streaming_tail):
        density_change = (
            variants["density_left-right"] - production_density
        )
        first_change = variants["first_left-right"] - production_first
        print(
            "  free-streaming tail correction / global production RHS: "
            f"density={region_norm(density_change, free_streaming_tail) / max(region_norm(production_density, all_nodes), np.finfo(float).tiny):.6e}, "
            f"first={region_norm(first_change, free_streaming_tail) / max(region_norm(production_first, all_nodes), np.finfo(float).tiny):.6e}"
        )

    density_closure = (
        variants["density_left-right"]
        + variants["density_right-left"]
        - 2.0 * production_density
    )
    first_closure = (
        variants["first_left-right"]
        + variants["first_right-left"]
        - 2.0 * production_first
    )
    print(
        "  orientation-average closure (max abs): "
        f"density={np.max(np.abs(density_closure)):.6e}, "
        f"first={np.max(np.abs(first_closure)):.6e}"
    )

    return {
        "radius_km": radius_km,
        "optical_depth": optical_depth,
        "flux_factor": flux_factor,
        "production_density": production_density,
        "production_first": production_first,
        "variants": variants,
        "cell_weights": cell_weights,
    }


def write_figure(
    output: Path,
    fourier: dict[str, object],
    profile: dict[str, object] | None,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    theta = np.asarray(fourier["theta"])
    width = float(fourier["width"])
    fig, axes = plt.subplots(
        3, 2, figsize=(14, 12), constrained_layout=True
    )

    for name in TRACE_PAIRS:
        diffusion = np.asarray(fourier[name]["diffusion"])
        axes[0, 0].plot(
            theta[1:],
            -width**2 * diffusion[1:, 0].real,
            label=TRACE_LABELS[name],
        )
        axes[0, 1].plot(
            theta[1:],
            -width**2 * diffusion[1:, 1].real,
            label=TRACE_LABELS[name],
        )
    axes[0, 0].plot(theta, theta**2 / 3.0, "k--", label=r"$\varphi^2/3$")
    axes[0, 0].set_title("(a) Physical branch of paired mixed derivative")
    axes[0, 1].set_title("(b) Element-internal branch")
    axes[0, 0].set_ylabel(r"$-h^2\,\mathrm{Re}(\lambda)$")
    axes[0, 1].set_ylabel(r"$-h^2\,\mathrm{Re}(\lambda)$")

    background_flux_factor = np.asarray(fourier["background_flux_factor"])
    for name in TRACE_PAIRS:
        m1_maximum_real = width * np.asarray(
            fourier[name]["m1_maximum_real"]
        )
        axes[1, 0].plot(
            theta,
            np.max(m1_maximum_real, axis=0),
            label=TRACE_LABELS[name],
        )
        axes[1, 1].plot(
            background_flux_factor,
            np.max(m1_maximum_real, axis=1),
            label=TRACE_LABELS[name],
        )
    axes[1, 0].axhline(0.0, color="k", linestyle="--")
    axes[1, 1].axhline(0.0, color="k", linestyle="--")
    axes[1, 0].set_title("(c) Worst Minerbo-M1 growth over flux factor")
    axes[1, 1].set_title("(d) Worst Minerbo-M1 growth over wavenumber")
    axes[1, 0].set_ylabel(r"$h\max\mathrm{Re}(\lambda)$")
    axes[1, 1].set_ylabel(r"$h\max\mathrm{Re}(\lambda)$")
    axes[1, 1].set_xlabel("background flux factor")

    if profile is not None:
        radius = np.asarray(profile["radius_km"]).reshape(-1)
        tau = np.asarray(profile["optical_depth"]).reshape(-1)
        flux_factor = np.asarray(profile["flux_factor"]).reshape(-1)
        order = np.argsort(radius)
        axes[2, 0].plot(radius[order], flux_factor[order], label="flux factor")
        axes[2, 0].plot(radius[order], tau[order], label=r"$\kappa_{\rm tr}\Delta r$")
        axes[2, 0].set_yscale("log")
        axes[2, 0].set_ylim(bottom=1.0e-8)
        axes[2, 0].set_title("(e) Saved initial profile: transport regimes")
        axes[2, 0].set_ylabel("flux factor or element optical depth")

        production_density = np.asarray(profile["production_density"])
        production_first = np.asarray(profile["production_first"])
        variants = profile["variants"]
        cell_weights = np.asarray(profile["cell_weights"])

        def cell_norm(values: np.ndarray) -> np.ndarray:
            return np.sqrt(np.sum(cell_weights * values**2, axis=1))

        delta_density = (
            np.asarray(variants["density_left-right"])
            - production_density
        )
        delta_first = (
            np.asarray(variants["first_left-right"])
            - production_first
        )
        density_scale = cell_norm(production_density)
        first_scale = cell_norm(production_first)
        density_floor = 1.0e-14 * max(
            np.max(density_scale), np.finfo(float).tiny
        )
        first_floor = 1.0e-14 * max(
            np.max(first_scale), np.finfo(float).tiny
        )
        cell_radius = np.mean(np.asarray(profile["radius_km"]), axis=1)
        axes[2, 1].plot(
            cell_radius,
            cell_norm(delta_density) / np.maximum(density_scale, density_floor),
            label=r"density equation $\|\delta R_D\|/\|R_D\|$",
        )
        axes[2, 1].plot(
            cell_radius,
            cell_norm(delta_first) / np.maximum(first_scale, first_floor),
            label=r"first-moment equation $\|\delta R_I\|/\|R_I\|$",
        )
        axes[2, 1].set_yscale("log")
        axes[2, 1].set_ylim(bottom=1.0e-10)
        axes[2, 1].set_title("(f) Effect of either paired orientation by element")
        axes[2, 1].set_ylabel("relative change from centered production RHS")
    else:
        axes[2, 0].axis("off")
        axes[2, 1].axis("off")

    axes[0, 0].set_xlabel(r"element wavenumber $\varphi$")
    axes[0, 1].set_xlabel(r"element wavenumber $\varphi$")
    axes[1, 0].set_xlabel(r"element wavenumber $\varphi$")
    for axis in axes[2, :].flat:
        if axis.has_data():
            axis.set_xlabel("radius [km]")
    for axis in axes.flat:
        if axis.has_data():
            axis.grid(alpha=0.2)
            axis.legend(loc="best", fontsize=8)
    fig.savefig(output, dpi=180)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", type=int, default=64)
    parser.add_argument("--frozen-file", type=Path)
    parser.add_argument("--species", type=int, default=1)
    parser.add_argument("--energy", type=int, default=25)
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args()

    fourier = fourier_audit(args.cells)
    print("Paired-trace p=1 Fourier audit")
    for name in TRACE_PAIRS:
        hyperbolic = np.asarray(fourier[name]["hyperbolic"])
        diffusion = np.asarray(fourier[name]["diffusion"])
        m1_maximum_real = float(fourier["width"]) * np.asarray(
            fourier[name]["m1_maximum_real"]
        )
        growth_by_flux_factor = np.max(m1_maximum_real, axis=1)
        unstable = np.flatnonzero(growth_by_flux_factor > 1.0e-12)
        onset = (
            float(np.asarray(fourier["background_flux_factor"])[unstable[0]])
            if unstable.size
            else float("nan")
        )
        print(
            f"  {name:11s}: max h Re(lambda_hyp)="
            f"{float(fourier['width']) * np.max(hyperbolic.real):.6e}, "
            f"max h Re(lambda_M1)={np.max(m1_maximum_real):.6e}, "
            f"instability onset f={onset:.3f}, "
            f"theta=0 diffusion eigenvalues h^2="
            f"{float(fourier['width'])**2 * diffusion[0]}"
        )

    profile = None
    if args.frozen_file is not None:
        profile = frozen_profile_audit(
            args.frozen_file, args.species, args.energy
        )
    if args.plot is not None:
        write_figure(args.plot, fourier, profile)


if __name__ == "__main__":
    main()
