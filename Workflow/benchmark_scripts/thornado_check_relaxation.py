#!/usr/bin/env python3
import sys
import re
import argparse
from typing import Dict, Tuple, List

# ----------------------------
# Baselines (from your example)
# ----------------------------
BASELINE_ERRORS = {
    1: 6.7280e-11,
    2: 8.6806e-11,
    3: 2.4620e-07,
    4: 2.4629e-07,
    5: 2.4620e-07,
    6: 2.4629e-07,
}

BASELINE_TIMERS = {
    "Timer_IMEX": 4.592535e+01,
    "Timer_Collisions": 4.560897e+01,
    "Timer_Collisions_PrimitiveFluid": 2.184375e-02,
    "Timer_Collisions_PrimitiveTwoMoment": 6.812166e+00,
    "Timer_Collisions_Solve": 3.749735e+01,
    "Timer_Collisions_OuterLoop": 3.555637e+01,
    "Timer_Collisions_InnerLoop": 2.757357e+01,
    "Timer_Collisions_ComputeOpacity": 6.672805e+00,
    "Timer_Collisions_ComputeRates": 9.139744e+00,
    "Timer_Collisions_InitializeRHS": 3.392771e-01,
    "Timer_Collisions_NeutrinoRHS": 7.851303e+00,
    "Timer_Collisions_MatterRHS": 9.714965e-01,
    "Timer_Collisions_SolveLS": 4.546268e+00,
    "Timer_Collisions_UpdateFP": 7.458959e+00,
    "Timer_Collisions_CheckOuter": 2.312379e-01,
    "Timer_Collisions_CheckInner": 2.958580e+00,
    "Timer_Opacity_D0": 2.569843e-01,
    "Timer_Opacity_LimitD0": 2.673732e-01,
    "Timer_Opacity_EC": 1.276625e+00,
    "Timer_Opacity_ES": 2.936564e-01,
    "Timer_Opacity_NES": 1.056617e+00,
    "Timer_Opacity_Pair": 1.039805e+00,
    "Timer_Opacity_Brem": 1.262378e+00,
    "Timer_OpacityRate_NES": 3.164959e+00,
    "Timer_OpacityRate_Pair": 2.268004e+00,
    "Timer_OpacityRate_Brem": 1.435534e+00,
    "Timer_TimeStepper": 3.232194e-01,
}

# ----------------------------
# Regex helpers
# ----------------------------
# Scientific float: captures 1.23, 1.23E+04, -4.5e-6, etc.
FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"

# Error table line:
#  Species = 01, Inf Error = 6.7280E-11
ERROR_LINE_RE = re.compile(
    rf"^\s*Species\s*=\s*(\d+)\s*,\s*Inf\s+Error\s*=\s*({FLOAT_RE})\s*$"
)

# Timer line (may have extra trailing numbers; we capture the first seconds value):
# Timer_Collisions : 4.560897E+01 s
TIMER_LINE_RE = re.compile(
    rf"^\s*(Timer_[A-Za-z0-9_]+)\s*:\s*({FLOAT_RE})\s*s\b"
)

# ----------------------------
# Parsing
# ----------------------------
def parse_error_check(text: str) -> Dict[int, float]:
    """
    Parses the Relaxation error table rows.
    We don’t require the header lines; we just match the numeric rows.
    """
    errors: Dict[int, float] = {}
    for line in text.splitlines():
        m = ERROR_LINE_RE.match(line)
        if m:
            species = int(m.group(1))
            val = float(m.group(2))
            errors[species] = val
    return errors

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
      Timer_IMEX, Timer_Collisions, Timer_TimeStepper,
      all Timer_Collisions_*, all Timer_Opacity_*
    """
    keep = {}
    for k, v in all_timers.items():
        if k in ("Timer_IMEX", "Timer_Collisions", "Timer_TimeStepper"):
            keep[k] = v
        elif k.startswith("Timer_Collisions_"):
            keep[k] = v
        elif k.startswith("Timer_Opacity_") or k.startswith("Timer_OpacityRate_"):
            keep[k] = v
    return keep

# ----------------------------
# Checks + reporting
# ----------------------------
def eval_errors(
    found: Dict[int, float],
    abs_tol: float,
    rel_tol: float,
    require_all_species: bool) -> Tuple[bool, List[str], float]:
    """
    PASS if for each species considered:
      err <= max(abs_tol, rel_tol * baseline_err)
    Also supports requiring all baseline species to be present.
    """
    msgs = []
    ok = True

    if require_all_species:
        missing = sorted(set(BASELINE_ERRORS.keys()) - set(found.keys()))
        if missing:
            ok = False
            msgs.append(f"Missing species lines: {missing}")

    species_to_check = sorted(set(found.keys()) & set(BASELINE_ERRORS.keys()))
    if not species_to_check:
        return False, ["No Species/Inf Error lines found."], float("inf")

    worst = -1.0
    for sp in species_to_check:
        errors = found[sp]
        base_errors = BASELINE_ERRORS[sp]

        for label, err, base_err in (("N", errors, base_errors),):
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
    msgs = []
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
    ap = argparse.ArgumentParser(description="Check Thornado Relaxation benchmark stdout (errors + timers).")
    ap.add_argument("logfile", nargs="?", default="-", help="stdout log path, or '-' for stdin")

    # Error tolerances
    ap.add_argument("--abs-tol", type=float, default=1.0e-6,
                    help="Absolute tolerance on each error (default: 1e-6).")
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
    print("THORNADO_BENCHMARK=Relaxation")
    print(f"THORNADO_LOG_SOURCE={src}")
    print(f"THORNADO_ERROR_STATUS={'PASS' if ok_err else 'FAIL'}")
    print(f"THORNADO_ERROR_ABS_TOL={args.abs_tol:.3e}")
    print(f"THORNADO_ERROR_REL_TOL={args.rel_tol:.3e}")
    print(f"THORNADO_ERROR_FOUND_SPECIES={len(found_errors)}")
    if found_errors:
        max_err = max(found_errors.values())
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
