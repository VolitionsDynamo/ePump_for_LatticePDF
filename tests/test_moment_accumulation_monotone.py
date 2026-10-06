"""
test_moment_accumulation_monotone.py

Runs MomentAccumulationScanner on CT18NNLO with max_n=3 and checks that
σ_after/σ_before is non-increasing as more moments are constrained.

Uses a temporary output directory so it does not pollute scan_results.
"""

import sys, os, shutil, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scan_window_moments import MomentAccumulationScanner

# ── config ────────────────────────────────────────────────────────────────────
TMP_DIR = tempfile.mkdtemp(prefix='macc_test_')

cfg = {
    "pdf":    "CT18NNLO",
    "flavor": "u-d",
    "Q2":     4.0,
    "nx":     50,           # fewer points for speed

    "mc2h_neig":    50,
    "mc2h_Q":       1.0,
    "mc2h_epsilon": 1000.0,
    "mc2h_max_nf":  3,

    "moment_scan": {
        "x0":      0.30,
        "w":       0.10,
        "weight":  "gaussian",
        "max_n":   3,
        "rel_unc": 0.10,
        "corr":    None,
    },

    "moment_xmin":  1e-4,
    "moment_xmax":  0.999,

    "charge_xmin": 1e-4,
    "charge_xmax": 0.999,
    "charge_observables": [
        {"label": "gT(u-d)", "flavor": "u-d", "moment": 0, "weight": "1"},
        {"label": "gT(u)",   "flavor": "u",   "moment": 0, "weight": "1"},
    ],

    "pdf_label":   "CT18NNLO",
    "output_dir":  TMP_DIR,
    "epump_path":  "./ePump_kp20221218/src/UpdatePDFs",
    "lhapdf_path": None,
}

# ── run ───────────────────────────────────────────────────────────────────────
print(f"Output dir: {TMP_DIR}")
scanner = MomentAccumulationScanner(cfg)
scanner.run()

# ── check monotonicity ────────────────────────────────────────────────────────
import numpy as np

TOLERANCE = 1e-6   # allow tiny floating-point noise

failures = []

def check_monotone(name, series):
    """series[i] = ratio after constraining i+1 moments."""
    for i in range(len(series) - 1):
        if series[i+1] > series[i] + TOLERANCE:
            failures.append(
                f"  {name}: ratio INCREASED from n={i+1} to n={i+2}: "
                f"{series[i]:.6f} → {series[i+1]:.6f}")

max_n = len(scanner.n_values)

print("\n── Window moments ──────────────────────────────────────────────────")
for k in range(max_n):
    series = scanner.ratio_window[k]
    label  = f"wm{k+1}"
    print(f"  {label}: {series}")
    check_monotone(label, series)

print("\n── Charge observables ──────────────────────────────────────────────")
for obs_idx, lbl in enumerate(scanner.charge_labels):
    series = scanner.ratio_charges[obs_idx]
    print(f"  {lbl}: {series}")
    check_monotone(lbl, series)

print("\n── Polynomial moments ──────────────────────────────────────────────")
for ni in range(scanner.ratio_poly_moments.shape[0]):
    series = scanner.ratio_poly_moments[ni]
    label  = f"poly_n={ni}"
    print(f"  {label}: {series}")
    check_monotone(label, series)

# ── result ────────────────────────────────────────────────────────────────────
print()
if failures:
    print("FAIL — non-monotone ratios detected:")
    for f in failures:
        print(f)
    sys.exit(1)
else:
    print("PASS — all ratios are non-increasing as moments are accumulated.")

# cleanup
shutil.rmtree(TMP_DIR)
print(f"Cleaned up {TMP_DIR}")

