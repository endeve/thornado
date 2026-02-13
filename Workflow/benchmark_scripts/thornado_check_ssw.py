#!/usr/bin/env python3
import sys
import re
import argparse
from typing import Dict, Tuple, List

# ----------------------------
# Baseline (from your example)
# ----------------------------
BASELINE_ERRORS = {
    # species: (N, G1, G2, G3)
    1: (5.6484e-02, 5.9192e-02, 1.7667e-16, 1.7667e-16),
    2: (5.6484e-02, 5.9192e-02, 1.7667e-16, 1.7667e-16),
    3: (5.6484e-02, 5.9192e-02, 1.7667e-16, 1.7667e-16),
    4: (5.6484e-02, 5.9192e-02, 1.7667e-16, 1.7667e-16),
    5: (5.6484e-02, 5.9192e-02, 1.7667e-16, 1.7667e-16),
    6: (5.6484e-02, 5.9192e-02, 1.7667e-16, 1.7667e-16),
}

BASELINE_TIMERS = {
    "Timer_IMEX": 2.824543e+01,
    "Timer_Streaming": 2.683871e+01,
    "Timer_Streaming_Divergence": 1.952667e+01,
    "Timer_Streaming_ObserverCorrections": 7.067150e+00,
    "Timer_Streaming_Derivatives": 6.064119e-01,
    "Timer_Streaming_Eigenvalues": 2.174729e-02,
    "Timer_Streaming_NumericalFlux": 7.089531e-01,
    "Timer_Streaming_NumericalFlux_InOut": 2.323853e+00,
    "Timer_Streaming_NumericalFlux_RHS": 1.050864e+01,
    "Timer_Streaming_NumericalFlux_LS": 2.120300e+00,
    "Timer_Streaming_NumericalFlux_Update": 6.042609e+00,
    "Timer_Streaming_PrimitiveTwoMoment": 2.167529e+01,
    "Timer_Streaming_Sources": 1.160514e-01,
    "Timer_Streaming_LinearAlgebra": 1.502844e+00,

    "Timer_PL": 1.275936e+00,
    "Timer_PL_Permute": 2.541254e-01,
    "Timer_PL_CellAverage": 6.680557e-02,
    "Timer_PL_PointValues": 1.720351e-01,
    "Timer_PL_Theta_1": 4.616201e-02,
    "Timer_PL_Theta_2": 2.566116e-01,
    "Timer_PL_EnergyLimiter": 2.815912e-01,

    "Timer_TimeStepper": 2.752848e-01,
}

# ----------------------------
# Regex helpers
# ----------------------------
# Scientific float: captures 1.23, 1.23E+04, -4.5e-6, etc.
FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"

# Error table line:
#  01  5.6484E-02  5.9192E-02  1.7667E-16  1.7667E-16
ERROR_LINE_RE = re.compile(
    rf"^\s*(\d+)\s+({FLOAT_RE})\s+({FLOAT_RE})\s+({FLOAT_RE})\s+({FLOAT_RE})\s*$"
)

# Timer line (may have extra trailing numbers; we capture the first seconds value):
# Timer_Streaming_Divergence : 1.952667E+01 s 7.275562E-01
TIMER_LINE_RE = re.compile(
    rf"^\s*(Timer_[A-Za-z0-9_]+)\s*:\s*({FLOAT_RE})\s*s\b"
)

# ----------------------------
# Parsing
# ----------------------------
def parse_error_check(text: str) -> Dict[int, Tuple[float, float, float, float]]:
    """
    Parses the SineWaveStreaming error table rows into {species: (N,G1,G2,G3)}.
    We don’t require the header lines; we just match the numeric rows.
    """
    out: Dict[int, Tuple[float, float, float, float]] = {}
    for line in text.splitlines():
        m = ERROR_LINE_RE.match(line)
        if not m:
            continue
        sp = int(m.group(1))
        N  = float(m.group(2))
        G1 = float(m.group(3))
        G2 = float(m.group(4))
        G3 = float(m.group(5))
        out[sp] = (N, G1, G2, G3)
    return out

def parse_timers(text: str) -> Dict[str, float]:
    """
    Extracts Timer_* lines of the form:
      Timer_Name :  1.234567E+00 s
    """
    timers: Dict[str, float] = {}
    for line in text.splitlines():
        m = TIMER_LINE_RE.match(line)
        if m:
            timers[m.group(1)] = float(m.group(2))
    return timers

def select_timer_subset(all_timers: Dict[str, float]) -> Dict[str, float]:
    """
    Keep only:
      Timer_IMEX, Timer_Streaming, Timer_TimeStepper,
      all Timer_Streaming_*, all Timer_PL_*
    """
    keep: Dict[str, float] = {}
    for k, v in all_timers.items():
        if k in ("Timer_IMEX", "Timer_Streaming", "Timer_PL", "Timer_TimeStepper"):
            keep[k] = v
        elif k.startswith("Timer_Streaming_"):
            keep[k] = v
        elif k.startswith("Timer_PL_"):
            keep[k] = v
    return keep

# ----------------------------
# Checks + reporting
# ----------------------------
def eval_errors(
    found: Dict[int, Tuple[float, float, float, float]],
    abs_tol: float,
    rel_tol: float,
    require_all_species: bool) -> Tuple[bool, List[str], float]:
    """
    PASS if for each species considered:
      err <= max(abs_tol, rel_tol * baseline_err)
    Also supports requiring all baseline species to be present.
    """
    msgs: List[str] = []
    ok = True

    if require_all_species:
        missing = sorted(set(BASELINE_ERRORS.keys()) - set(found.keys()))
        if missing:
            ok = False
            msgs.append(f"Missing species lines: {missing}")

    species_to_check = sorted(set(found.keys()) & set(BASELINE_ERRORS.keys()))
    if not species_to_check:
        return False, ["No SineWaveStreaming error rows found."]

    worst = -1.0
    for sp in species_to_check:
        (N, G1, G2, G3) = found[sp]
        (bN, bG1, bG2, bG3) = BASELINE_ERRORS[sp]

        for label, err, base_err in (("N", N, bN), ("G1", G1, bG1), ("G2", G2, bG2), ("G3", G3, bG3)):
            worst = max(worst, err)
            this_ok, msg = check_value(f"Species {sp:02d} {label:2} Error", err, base_err, abs_tol, rel_tol)
            msgs.append(msg)
            if not this_ok:
                ok = False

    return ok, msgs, worst

def check_value(name: str, err: float, base_err: float, abs_tol: float, rel_tol: float) -> Tuple[bool, str]:
    """
    PASS if |err-base_err| <= max(abs_tol, rel_tol*|base_err|)
    """
    diff = abs(err - base_err)
    thresh = max(abs_tol, rel_tol * abs(base_err))
    ok = diff <= thresh
    status = "OK" if ok else "BAD"
    msg = f"{name}: err={err:.6e} base_err={base_err:.6e} |Δ|={diff:.6e} thresh={thresh:.6e} => {status}"
    return ok, msg

def format_timer_table(
    new: Dict[str, float],
    baseline: Dict[str, float]) -> str:
    """
    Pretty fixed-width table with delta and percent change.
    """
    keys = sorted(new.keys())
    # widths
    name_w = max(28, max((len(k) for k in keys), default=28))
    header = (
        f"{'TIMER'.ljust(name_w)}  {'NEW(s)'.rjust(12)}  {'BASE(s)'.rjust(12)}  "
        f"{'Δ(s)'.rjust(12)}  {'Δ%'.rjust(9)}"
    )
    lines = [header, "-" * len(header)]
    for k in keys:
        nv = new[k]
        bv = baseline.get(k)
        if bv is None:
            lines.append(f"{k.ljust(name_w)}  {nv:12.6e}  {'(none)'.rjust(12)}  {'':>12}  {'':>9}")
            continue
        dv = nv - bv
        dp = (dv / bv * 100.0) if bv != 0 else float("inf")
        lines.append(f"{k.ljust(name_w)}  {nv:12.6e}  {bv:12.6e}  {dv:12.6e}  {dp:9.2f}")
    return "\n".join(lines)

def eval_perf_regression(
    new: Dict[str, float],
    baseline: Dict[str, float],
    rel_regress: float,
    abs_regress_s: float) -> Tuple[bool, List[str]]:
    """
    Regression if (new - base) > max(abs_regress_s, rel_regress*base)
    """
    ok = True
    msgs: List[str] = []
    for k, bv in baseline.items():
        if k not in new:
            continue
        nv = new[k]
        dv = nv - bv
        thresh = max(abs_regress_s, rel_regress * bv)
        if dv > thresh:
            ok = False
            msgs.append(f"{k}: +{dv:.3f}s exceeds thresh={thresh:.3f}s (base={bv:.3f}s)")
    return ok, msgs

# ----------------------------
# Main
# ----------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description="Check Thornado SineWaveStreaming benchmark stdout (errors + timers).")
    ap.add_argument("logfile", nargs="?", default="-", help="stdout log path, or '-' for stdin")

    # Error tolerances
    ap.add_argument("--abs-tol", type=float, default=1e-8,
                    help="Absolute tolerance on each error (default: 1e-8).")
    ap.add_argument("--rel-tol", type=float, default=0.05,
                    help="Relative tolerance vs baseline on each error (default: 0.05 = 5%).")
    ap.add_argument("--require-all-species", action="store_true",
                    help="Fail if any baseline species line is missing.")

    # Perf
    ap.add_argument("--perf-rel-regress", type=float, default=0.10,
                    help="Relative perf regression threshold (default: 0.10).")
    ap.add_argument("--perf-abs-regress-s", type=float, default=1.0,
                    help="Absolute perf regression threshold in seconds (default: 1.0s).")
    ap.add_argument("--fail-on-perf-regression", action="store_true",
                    help="Exit nonzero if performance regression is detected.")
    args = ap.parse_args()

    if args.logfile == "-" or args.logfile == "":
        text = sys.stdin.read()
        src = "stdin"
    else:
        with open(args.logfile, "r", encoding="utf-8", errors="replace") as f:
            text = f.read()
        src = args.logfile

    found_errors = parse_error_check(text)
    all_timers = parse_timers(text)
    picked_timers = select_timer_subset(all_timers)

    # ----------------------------
    # Correctness report (extractable)
    # ----------------------------
    ok_err, err_msgs, worst = eval_errors(found_errors, args.abs_tol, args.rel_tol, args.require_all_species)

    print("=== THORNADO_TEST_RESULT_BEGIN ===")
    print("THORNADO_BENCHMARK=SineWaveStreaming")
    print(f"THORNADO_LOG_SOURCE={src}")
    print(f"THORNADO_ERROR_STATUS={'PASS' if ok_err else 'FAIL'}")
    print(f"THORNADO_ERROR_ABS_TOL={args.abs_tol:.3e}")
    print(f"THORNADO_ERROR_REL_TOL={args.rel_tol:.3e}")
    print(f"THORNADO_ERROR_FOUND_SPECIES={len(found_errors)}")
    if found_errors:
        max_err = max(err for tup in found_errors.values() for err in tup)
        print(f"THORNADO_ERROR_MAX_ERROR={max_err:.6e}")
    for m in err_msgs:
        print(f"THORNADO_ERROR_DETAIL={m}")
    print("=== THORNADO_TEST_RESULT_END ===")

    # ----------------------------
    # Timers report (extractable + pretty)
    # ----------------------------
    print("=== THORNADO_TIMERS_SUMMARY_BEGIN ===")
    if not picked_timers:
        print("THORNADO_TIMERS_STATUS=NO_TIMERS_FOUND")
        ok_perf = True
        perf_msgs: List[str] = []
    else:
        print("THORNADO_TIMERS_STATUS=OK")
        print(format_timer_table(picked_timers, BASELINE_TIMERS))

        ok_perf, perf_msgs = eval_perf_regression(
            picked_timers, BASELINE_TIMERS, args.perf_rel_regress, args.perf_abs_regress_s
        )
        print(f"THORNADO_PERF_STATUS={'PASS' if ok_perf else 'REGRESSION'}")
        print(f"THORNADO_PERF_REL_THRESHOLD={args.perf_rel_regress:.2f}")
        print(f"THORNADO_PERF_ABS_THRESHOLD_S={args.perf_abs_regress_s:.2f}")
        for m in perf_msgs:
            print(f"THORNADO_PERF_DETAIL={m}")
    print("=== THORNADO_TIMERS_SUMMARY_END ===")

    # Exit codes:
    # 0 = all ok (and perf ok or not enforced)
    # 2 = correctness failure
    # 3 = perf regression (only if enforced)
    if not ok_err:
        return 2
    if args.fail_on_perf_regression and picked_timers and not ok_perf:
        return 3
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
