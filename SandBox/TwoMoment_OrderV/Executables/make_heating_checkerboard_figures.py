#!/usr/bin/env python3
"""Generate the figures used by DiffusionLimitNumericalFlux.tex."""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_named(path: Path) -> dict[str, np.ndarray]:
    with path.open() as stream:
        for line in stream:
            if line.startswith("# r_cm") or line.startswith("# x_km"):
                names = line[1:].split()
                break
        else:
            raise ValueError(f"missing named header in {path}")
    values = np.atleast_2d(np.loadtxt(path, comments="#"))
    return dict(zip(names, values.T, strict=True))


def signed_axis(
    axis: plt.Axes, linear_threshold: float, ticks: tuple[float, ...] | None = None
) -> None:
    axis.axhline(0.0, color="0.45", lw=0.7)
    axis.set_yscale("symlog", linthresh=linear_threshold, linscale=0.7)
    if ticks is not None:
        axis.set_yticks(ticks)
    axis.grid(alpha=0.2)


def heating_diagnostics(data_root: Path, output: Path) -> None:
    pdars = read_named(
        data_root / "md_PDARS_All_Iso_lim_tol_1e-8_1e-8_no_cons2_wb_all.dat"
    )

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.2), constrained_layout=True)

    r = pdars["r_cm"] / 1.0e5
    mask = r <= 40.0
    fields = (
        ("QH_total_erg_cm3_s", "instantaneous post-step collision RHS", "C0"),
        ("QH_IMEX_source_total_erg_cm3_s", "IMEX-weighted collision RHS", "C1"),
        ("QH_IMEX_increment_erg_cm3_s", "IMEX-weighted implicit increment", "C2"),
        ("QH_applied_last_step_erg_cm3_s", "actual fluid increment", "k"),
    )
    for field, label, color in fields:
        axes[0].plot(r[mask], pdars[field][mask], lw=1.1, label=label, color=color)
    signed_axis(
        axes[0], 1.0e27,
        (-1.0e36, -1.0e33, -1.0e30, 0.0, 1.0e30, 1.0e33, 1.0e36),
    )
    axes[0].set_title("A. PDARS: completed-step diagnostics")
    axes[0].legend(fontsize=7)

    species_labels = (
        (1, r"$\nu_e$"), (2, r"$\bar\nu_e$"),
        (3, r"$\nu_\mu$"), (4, r"$\bar\nu_\mu$"),
        (5, r"$\nu_\tau$"), (6, r"$\bar\nu_\tau$"),
    )
    for species, label in species_labels:
        axes[1].plot(
            r[mask], pdars[f"QH_species_s{species}"][mask], lw=1.0, label=label
        )
    signed_axis(
        axes[1], 1.0e27,
        (-1.0e36, -1.0e33, -1.0e30, 0.0, 1.0e30, 1.0e33, 1.0e36),
    )
    axes[1].set_title("B. Instantaneous heating by species")
    axes[1].legend(fontsize=7, ncol=2)

    for process in ("EmAb", "NNS", "NES", "Pair", "NuPair", "Brem"):
        axes[2].plot(
            r[mask], pdars[f"QH_process_{process}"][mask], lw=1.0, label=process
        )
    signed_axis(
        axes[2], 1.0e27,
        (-1.0e36, -1.0e33, -1.0e30, 0.0, 1.0e30, 1.0e33, 1.0e36),
    )
    axes[2].set_title("C. Instantaneous heating by process")
    axes[2].legend(fontsize=7, ncol=2)

    for axis in axes.flat:
        axis.set_xlim(0.0, axis.get_xlim()[1])
        axis.set_xlabel("radius [km]")
        axis.set_ylabel(r"$Q_H$ [erg cm$^{-3}$ s$^{-1}$]")

    fig.savefig(output, dpi=200)
    plt.close(fig)


def prescribed_absorption(data_root: Path, output: Path) -> None:
    sweep = data_root / "absorption_diffusion_optical_depth_sweep"
    cases = (
        (
            "absorption_diffusion_balance_tau_1p000000e05.dat",
            "A. Standard LLF",
            2.5e-9,
        ),
        (
            "absorption_diffusion_balance_scaled_tau_1p000000e05.dat",
            r"B. LLF factor $\min(1,1/\bar\tau_h)$",
            4.0e-11,
        ),
    )

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.2), constrained_layout=True)
    plotted_data = []
    for axis, (name, title, limit) in zip(axes[:2], cases, strict=True):
        data = read_named(sweep / name)
        plotted_data.append(data)
        axis.plot(
            data["x_km"], data["source_exact"],
            lw=2.0, color="k", label=r"prescribed $Q_*$", zorder=2,
        )
        axis.plot(
            data["x_km"], data["source_num"],
            marker="o", ms=2.8, lw=0.8, color="C3",
            label=r"nodal $Q_h$", zorder=3,
        )
        # This controlled test uses two equal-weight Gauss--Legendre nodes
        # per spatial element.
        x_element = data["x_km"].reshape(-1, 2).mean(axis=1)
        q_element = data["source_num"].reshape(-1, 2).mean(axis=1)
        q_exact_element = data["source_exact"].reshape(-1, 2).mean(axis=1)
        axis.plot(
            x_element, q_exact_element,
            ls="--", lw=1.2, color="0.35",
            label=r"element-average $Q_*$", zorder=2,
        )
        axis.plot(
            x_element, q_element,
            marker="s", ms=3.0, lw=1.2, color="C0",
            label=r"element-average $Q_h$", zorder=4,
        )
        axis.axhline(0.0, color="0.55", lw=0.7)
        axis.set_ylim(-limit, limit)
        axis.set_title(title)
        axis.set_xlabel("x [km]")
        axis.set_ylabel(r"absorption exchange rate $Q$ [code units]")
        axis.set_xlim(0.0, 1.0)
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8)

    standard, scaled = plotted_data
    axes[2].semilogy(
        standard["x_km"], np.abs(standard["source_exact"]),
        lw=2.0, color="k", label=r"prescribed $|Q_*|$",
    )
    axes[2].semilogy(
        standard["x_km"], np.abs(standard["source_num"]),
        marker="o", ms=2.8, lw=0.8, color="C3",
        label=r"standard LLF $|Q_h|$",
    )
    axes[2].semilogy(
        scaled["x_km"], np.abs(scaled["source_num"]),
        marker="s", ms=2.8, lw=0.8, color="C0",
        label=r"scaled LLF $|Q_h|$",
    )
    axes[2].set_title("C. Exchange-rate magnitudes")
    axes[2].set_xlabel("x [km]")
    axes[2].set_ylabel(r"$|Q|$ [code units]")
    axes[2].set_xlim(0.0, 1.0)
    axes[2].grid(which="both", alpha=0.2)
    axes[2].legend(fontsize=8)
    fig.suptitle(r"Prescribed absorption--diffusion balance, $\bar\tau=10^5$")
    fig.savefig(output, dpi=200)
    plt.close(fig)


def polynomial_order_comparison(data_root: Path, output: Path) -> None:
    local_root = Path(__file__).resolve().parent
    cases = (
        (
            "16 elements, 2 nodes",
            data_root / "absorption_diffusion_optical_depth_sweep"
                / "absorption_diffusion_balance_scaled_tau_1p000000e04.dat",
            2,
            "C3",
        ),
        (
            "24 elements, 2 nodes",
            local_root / "absorption_diffusion_optical_depth_sweep_nX24_nNodes2"
                / "absorption_diffusion_balance_scaled_tau_1p000000e04.dat",
            2,
            "C1",
        ),
        (
            "16 elements, 3 nodes",
            local_root / "absorption_diffusion_optical_depth_sweep_nNodes3"
                / "absorption_diffusion_balance_scaled_tau_1p000000e04.dat",
            3,
            "C0",
        ),
    )
    weights = {
        2: np.asarray((0.5, 0.5)),
        3: np.asarray((5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0)),
    }
    loaded = [(label, read_named(path), nodes, color)
              for label, path, nodes, color in cases]

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.2), constrained_layout=True)
    for axis, index, title in (
        (axes[0], 0, "A. 16 elements, 2 nodes"),
        (axes[1], 2, "B. 16 elements, 3 nodes"),
    ):
        _, data, nodes, _ = loaded[index]
        x_element = data["x_km"].reshape(-1, nodes) @ weights[nodes]
        q_element = data["source_num"].reshape(-1, nodes) @ weights[nodes]
        q_exact_element = data["source_exact"].reshape(-1, nodes) @ weights[nodes]
        axis.plot(
            data["x_km"], data["source_exact"],
            lw=2.0, color="k", label=r"prescribed $Q_*$",
        )
        axis.plot(
            data["x_km"], data["source_num"],
            marker="o", ms=2.7, lw=0.8, color="C3", label=r"nodal $Q_h$",
        )
        axis.plot(
            x_element, q_exact_element,
            ls="--", lw=1.2, color="0.35", label=r"element-average $Q_*$",
        )
        axis.plot(
            x_element, q_element,
            marker="s", ms=3.0, lw=1.2, color="C0",
            label=r"element-average $Q_h$",
        )
        axis.axhline(0.0, color="0.55", lw=0.7)
        axis.set_title(title)
        axis.set_xlabel("x [km]")
        axis.set_ylabel(r"absorption exchange rate $Q$ [code units]")
        axis.set_xlim(0.0, 1.0)
        axis.set_ylim(-3.5e-10, 1.0e-10)
        axis.grid(alpha=0.2)
        axis.legend(fontsize=7)

    for label, data, _, color in loaded:
        axes[2].semilogy(
            data["x_km"], np.abs(data["source_error"]),
            marker="o", ms=2.5, lw=0.8, color=color, label=label,
        )
    axes[2].set_title(r"C. Scaled-LLF error, $\bar\tau=10^4$")
    axes[2].set_xlabel("x [km]")
    axes[2].set_ylabel(r"$|Q_h-Q_*|$ [code units]")
    axes[2].set_xlim(0.0, 1.0)
    axes[2].grid(which="both", alpha=0.2)
    axes[2].legend(fontsize=8)

    fig.suptitle(r"Spatial polynomial-order comparison, $\bar\tau=10^4$")
    fig.savefig(output, dpi=200)
    plt.close(fig)


def standard_llf_polynomial_order_comparison(
    data_root: Path, output: Path,
) -> None:
    local_root = Path(__file__).resolve().parent
    cases = (
        (
            "16 elements, 2 nodes",
            data_root / "absorption_diffusion_optical_depth_sweep"
                / "absorption_diffusion_balance_tau_1p000000e05.dat",
            2,
            "C3",
        ),
        (
            "24 elements, 2 nodes",
            local_root / "absorption_diffusion_optical_depth_sweep_nX24_nNodes2"
                / "absorption_diffusion_balance_standard_tau_1p000000e05.dat",
            2,
            "C1",
        ),
        (
            "16 elements, 3 nodes",
            local_root / "absorption_diffusion_optical_depth_sweep_nNodes3"
                / "absorption_diffusion_balance_standard_tau_1p000000e05.dat",
            3,
            "C0",
        ),
    )
    weights = {
        2: np.asarray((0.5, 0.5)),
        3: np.asarray((5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0)),
    }
    loaded = [(label, read_named(path), nodes, color)
              for label, path, nodes, color in cases]

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.2), constrained_layout=True)
    for axis, index, title in (
        (axes[0], 0, "A. 16 elements, 2 nodes"),
        (axes[1], 2, "B. 16 elements, 3 nodes"),
    ):
        _, data, nodes, _ = loaded[index]
        x_element = data["x_km"].reshape(-1, nodes) @ weights[nodes]
        q_element = data["source_num"].reshape(-1, nodes) @ weights[nodes]
        q_exact_element = data["source_exact"].reshape(-1, nodes) @ weights[nodes]
        axis.plot(
            data["x_km"], data["source_exact"],
            lw=2.0, color="k", label=r"prescribed $Q_*$",
        )
        axis.plot(
            data["x_km"], data["source_num"],
            marker="o", ms=2.7, lw=0.8, color="C3", label=r"nodal $Q_h$",
        )
        axis.plot(
            x_element, q_exact_element,
            ls="--", lw=1.2, color="0.35", label=r"element-average $Q_*$",
        )
        axis.plot(
            x_element, q_element,
            marker="s", ms=3.0, lw=1.2, color="C0",
            label=r"element-average $Q_h$",
        )
        axis.axhline(0.0, color="0.55", lw=0.7)
        axis.set_title(title)
        axis.set_xlabel("x [km]")
        axis.set_ylabel(r"absorption exchange rate $Q$ [code units]")
        axis.set_xlim(0.0, 1.0)
        axis.set_ylim(-2.5e-9, 2.5e-9)
        axis.grid(alpha=0.2)
        axis.legend(fontsize=7)

    for label, data, _, color in loaded:
        axes[2].semilogy(
            data["x_km"], np.abs(data["source_error"]),
            marker="o", ms=2.5, lw=0.8, color=color, label=label,
        )
    axes[2].set_title(r"C. Standard-LLF error, $\bar\tau=10^5$")
    axes[2].set_xlabel("x [km]")
    axes[2].set_ylabel(r"$|Q_h-Q_*|$ [code units]")
    axes[2].set_xlim(0.0, 1.0)
    axes[2].grid(which="both", alpha=0.2)
    axes[2].legend(fontsize=8)

    fig.suptitle(r"Standard LLF: spatial polynomial-order comparison, $\bar\tau=10^5$")
    fig.savefig(output, dpi=200)
    plt.close(fig)


def operator_decomposition_audit(output: Path) -> None:
    local_root = Path(__file__).resolve().parent
    paths = {
        (2, "standard"): local_root
            / "absorption_diffusion_operator_audit_nX16_nNodes2"
            / "absorption_diffusion_operator_balance_standard_tau_1p000000e05.dat",
        (2, "scaled"): local_root
            / "absorption_diffusion_operator_audit_nX16_nNodes2"
            / "absorption_diffusion_operator_balance_scaled_tau_1p000000e05.dat",
        (3, "standard"): local_root
            / "absorption_diffusion_operator_audit_nX16_nNodes3"
            / "absorption_diffusion_operator_balance_standard_tau_1p000000e05.dat",
        (3, "scaled"): local_root
            / "absorption_diffusion_operator_audit_nX16_nNodes3"
            / "absorption_diffusion_operator_balance_scaled_tau_1p000000e05.dat",
    }
    data = {key: read_named(path) for key, path in paths.items()}

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.2), constrained_layout=True)
    for axis, nodes, title in (
        (axes[0, 0], 2, "A. Original LLF, 2 nodes"),
        (axes[0, 1], 3, "B. Original LLF, 3 nodes"),
    ):
        case = data[(nodes, "standard")]
        central_source = case["streaming_central_D"] + case["absorption_D"]
        axis.plot(
            case["x_km"], case["residual_D"], marker="o", ms=2.8,
            lw=0.8, color="C3", label=r"complete residual $R_D$",
        )
        axis.plot(
            case["x_km"], case["streaming_penalty_D"], marker=".", ms=2.5,
            lw=1.0, color="C0", label=r"LLF density-jump contribution",
        )
        axis.plot(
            case["x_km"], central_source, lw=1.5, color="k",
            label=r"central flux + absorption",
        )
        axis.axhline(0.0, color="0.6", lw=0.7)
        axis.set_title(title)
        axis.set_xlabel("x [km]")
        axis.set_ylabel("density-equation rate [code units]")
        axis.set_xlim(0.0, 1.0)
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8)

    labels = ("2 std.", "3 std.", "2 scaled", "3 scaled")
    keys = ((2, "standard"), (3, "standard"),
            (2, "scaled"), (3, "scaled"))
    residual_rms = []
    llf_rms = []
    central_rms = []
    for key in keys:
        case = data[key]
        residual_rms.append(np.sqrt(np.mean(case["residual_D"] ** 2)))
        llf_rms.append(np.sqrt(np.mean(case["streaming_penalty_D"] ** 2)))
        central_source = case["streaming_central_D"] + case["absorption_D"]
        central_rms.append(np.sqrt(np.mean(central_source ** 2)))
    xbar = np.arange(len(labels), dtype=float)
    width = 0.25
    axes[1, 0].bar(xbar - width, residual_rms, width, color="C3", label=r"$R_D$")
    axes[1, 0].bar(xbar, llf_rms, width, color="C0", label="LLF contribution")
    axes[1, 0].bar(
        xbar + width, central_rms, width, color="0.25",
        label="central + absorption",
    )
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_xticks(xbar, labels)
    axes[1, 0].set_ylabel("nodal RMS rate [code units]")
    axes[1, 0].set_title(r"C. Operator decomposition, $\bar\tau=10^5$")
    axes[1, 0].grid(axis="y", which="both", alpha=0.2)
    axes[1, 0].legend(fontsize=8)

    modal_series = []
    for nodes in (2, 3):
        case = data[(nodes, "standard")]
        element = case["element"].astype(int)
        node = case["node"].astype(int)
        order = np.lexsort((node, element))
        n_elements = int(element.max())
        gl_nodes, _ = np.polynomial.legendre.leggauss(nodes)
        vandermonde = np.polynomial.legendre.legvander(gl_nodes, nodes - 1)
        for field in ("residual_D", "streaming_penalty_D"):
            nodal = case[field][order].reshape(n_elements, nodes)
            coefficients = np.linalg.solve(vandermonde, nodal.T).T
            modal_series.append(
                (nodes, field, np.sqrt(np.mean(coefficients**2, axis=0)))
            )
    mode_x = np.arange(3, dtype=float)
    offsets = (-0.27, -0.09, 0.09, 0.27)
    colors = ("C3", "C0", "C1", "C2")
    names = (
        r"2 nodes: $R_D$", "2 nodes: LLF",
        r"3 nodes: $R_D$", "3 nodes: LLF",
    )
    for offset, color, name, (nodes, _, amplitudes) in zip(
        offsets, colors, names, modal_series, strict=True,
    ):
        axes[1, 1].bar(
            mode_x[:nodes] + offset, amplitudes, 0.18,
            color=color, label=name,
        )
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xticks(mode_x, ("mode 0", "mode 1", "mode 2"))
    axes[1, 1].set_ylabel("RMS Legendre coefficient [code units]")
    axes[1, 1].set_title("D. Local DG-mode amplitudes, original LLF")
    axes[1, 1].grid(axis="y", which="both", alpha=0.2)
    axes[1, 1].legend(fontsize=8, ncol=2)

    fig.suptitle(r"Semidiscrete operator audit at $\bar\tau=10^5$")
    fig.savefig(output, dpi=200)
    plt.close(fig)


def extract_log_value(path: Path, pattern: str) -> float:
    text = path.read_text()
    match = re.search(pattern, text)
    if match is None:
        raise ValueError(f"{pattern!r} not found in {path}")
    return float(match.group(1).replace("D", "E"))


def reduced_problem_sweeps(data_root: Path, output: Path) -> None:
    absorption = data_root / "absorption_diffusion_optical_depth_sweep"
    taus = np.asarray([1.0e2, 1.0e3, 1.0e4, 1.0e5])
    tags = ("1p000000e02", "1p000000e03", "1p000000e04", "1p000000e05")
    standard = []
    scaled = []
    for tag in tags:
        for prefix, target in (("tau", standard), ("scaled_tau", scaled)):
            path = absorption / f"{prefix}_{tag}.log"
            target.append(
                extract_log_value(
                    path,
                    r"source rel\. L1/L2/Linf\s*=\s*([+\-0-9.EeDd]+)",
                )
            )

    scattering = np.genfromtxt(
        data_root / "scattering_diffusion_limit_sweep" / "summary.csv",
        delimiter=",", names=True,
    )
    coupled = np.genfromtxt(
        data_root / "coupled_absorption_diffusion_sweep" / "summary.csv",
        delimiter=",", names=True, dtype=None, encoding="utf-8",
    )

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.0), constrained_layout=True)
    axes[0].loglog(taus, standard, "o-", label="standard LLF")
    axes[0].loglog(taus, scaled, "s-", label=r"LLF factor $\min(1,1/\bar\tau_h)$")
    axes[0].set_title("A. Prescribed absorption--diffusion")
    axes[0].set_ylabel("one-step source relative $L^1$")

    for nx, marker in ((16, "o"), (32, "s")):
        mask = scattering["spatial_elements"] == nx
        axes[1].loglog(
            scattering["optical_depth"][mask], scattering["density_l1"][mask],
            marker + "-", label=f"{nx} elements",
        )
    axes[1].set_title("B. Variable-opacity scattering diffusion")
    axes[1].set_ylabel("density relative $L^1$")

    for scheme, line in (("IMEX_ARS_111", "-"), ("IMEX_PDARS", "--")):
        for nx, marker in ((16, "o"), (32, "s")):
            mask = (coupled["scheme"] == scheme) & (coupled["spatial_elements"] == nx)
            axes[2].loglog(
                coupled["optical_depth"][mask], coupled["exchange_l1"][mask],
                marker=marker, ls=line, label=f"{scheme.removeprefix('IMEX_')}, {nx}",
            )
    axes[2].set_title("C. Coupled absorption--diffusion")
    axes[2].set_ylabel("exchange relative $L^1$")

    for axis in axes:
        axis.set_xlabel("domain optical depth")
        axis.grid(which="both", alpha=0.2)
        axis.legend(fontsize=7)
    fig.savefig(output, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-root", type=Path,
        default=Path("/Users/vv4/repos/github/thornado/SandBox/TwoMoment_OrderV/Executables"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    heating_diagnostics(args.data_root, args.output_dir / "heating_diagnostics.pdf")
    prescribed_absorption(args.data_root, args.output_dir / "prescribed_absorption.pdf")
    polynomial_order_comparison(
        args.data_root, args.output_dir / "polynomial_order_comparison.pdf"
    )
    standard_llf_polynomial_order_comparison(
        args.data_root,
        args.output_dir / "standard_llf_polynomial_order_comparison.pdf",
    )
    operator_decomposition_audit(
        args.output_dir / "operator_decomposition_audit.pdf"
    )
    reduced_problem_sweeps(args.data_root, args.output_dir / "reduced_problem_sweeps.pdf")


if __name__ == "__main__":
    main()
