"""
Verify run_moments() / load_moments() / plot_moments() /
       run_charges() / load_charges() / plot_charges() end-to-end.

Runs a tiny 2-point scan with CT18NNLO (Hessian, fast), then exercises
the window-to-moment and window-to-charge pipelines without additional ePump
invocations.
"""
import os, sys, tempfile, shutil
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import matplotlib
matplotlib.use('Agg')

from scan_window_moments import WindowMomentScanner

OUT = tempfile.mkdtemp(prefix='scan_moments_test_')
print(f"Output dir: {OUT}")

CFG = {
    "pdf":    "CT18NNLO",
    "flavor": "u-d",
    "Q2":     4.0,
    "nx":     50,
    "moment": 1,
    "weight": "gaussian",
    "mc2h_neig": 50,
    "mc2h_Q": 1.0,
    "mc2h_epsilon": 1000.0,
    "mc2h_max_nf": 3,
    "scan_points": [
        [0.25, 0.10, 0.10],
        [0.35, 0.10, 0.10],
    ],
    "moment_xmin": 1e-4,
    "moment_xmax": 0.999,
    "charge_xmin": 1e-4,
    "charge_xmax": 0.999,
    "charge_observables": [
        {"label": r"$g_T^{u-d}$",  "flavor": "u-d",           "moment": 0, "weight": "1"},
        {"label": r"$g_T^u$",      "flavor": "u",              "moment": 0, "weight": "1"},
        {"label": r"$g_T^d$",      "flavor": "d",              "moment": 0, "weight": "1"},
        {"label": r"$\langle x(\delta u^+ - \delta d^+)\rangle$",
                                   "flavor": "u+ubar-d-dbar",  "moment": 1, "weight": "1"},
    ],
    "output_dir": OUT,
    "epump_path": os.path.join(os.path.dirname(__file__), '..', 'ePump_kp20221218/src/UpdatePDFs'),
    "lhapdf_path": os.path.join(os.path.dirname(__file__), '..'),
}

# ── Step 1: run window-to-window scan ─────────────────────────────────────────
print("\n=== Step 1: run() ===")
scanner = WindowMomentScanner(CFG)
scanner.run()
assert scanner.ratio is not None and scanner.ratio.shape == (2, 1)
assert scanner.sigma_before_ww is not None and scanner.sigma_before_ww.shape == (2, 1)
assert scanner.central_ww is not None and scanner.central_ww.shape == (2, 1)
print(f"  window ratio = {scanner.ratio}")
print(f"  sigma_before_ww = {scanner.sigma_before_ww}")

# ── Step 2: run_moments() from same object ─────────────────────────────────────
print("\n=== Step 2: run_moments() ===")
scanner.run_moments()
assert scanner.ratio_moments is not None
assert scanner.ratio_moments.shape == (4, 2, 1), f"Bad shape: {scanner.ratio_moments.shape}"
assert scanner.sigma_before_moments is not None and scanner.sigma_before_moments.shape == (4,)
assert scanner.central_moments is not None and scanner.central_moments.shape == (4,)
print(f"  ratio_moments shape: {scanner.ratio_moments.shape}")
for ni in range(4):
    pct = 100 * scanner.sigma_before_moments[ni] / abs(scanner.central_moments[ni])
    print(f"  n={ni}: ratio={scanner.ratio_moments[ni].flatten()}  σ/μ={pct:.1f}%")

# ── Step 3: plot() and plot_moments() suffix = output_dir basename ─────────────
print("\n=== Step 3: plot() and plot_moments() ===")
_suffix = os.path.basename(OUT)
fig_ww = scanner.plot(save=True)
expected_ww = os.path.join(OUT, f'heatmap_{_suffix}.pdf')
assert os.path.exists(expected_ww), f"Missing: {expected_ww}"
print(f"  Window-to-window plot saved OK → {expected_ww}")

fig = scanner.plot_moments(save=True)
expected_plot = os.path.join(OUT, f'heatmap_moments_{_suffix}.pdf')
assert os.path.exists(expected_plot), f"Missing: {expected_plot}"
print(f"  Moments plot saved OK → {expected_plot}")

# ── Step 4: load_moments() in fresh object ─────────────────────────────────────
print("\n=== Step 4: load_moments() fresh object ===")
s2 = WindowMomentScanner(CFG)
s2.load_moments()
assert s2.ratio_moments is not None
assert s2.ratio_moments.shape == (4, 2, 1)
assert s2.sigma_before_moments is not None
assert all(
    abs(a - b) < 1e-10
    for a, b in zip(s2.ratio_moments.flatten(), scanner.ratio_moments.flatten())
    if not (a != a or b != b)
), "load_moments values differ from run_moments!"
print("  Loaded values match run_moments values.")

fig2 = s2.plot_moments(save=False)
print("  plot_moments() from loaded data: OK")

# ── Step 5: run_charges() from same object ────────────────────────────────────
print("\n=== Step 5: run_charges() ===")
scanner.run_charges()
assert scanner.ratio_charges is not None
assert scanner.ratio_charges.shape == (4, 2, 1), f"Bad shape: {scanner.ratio_charges.shape}"
assert scanner.sigma_before_charges is not None and scanner.sigma_before_charges.shape == (4,)
assert scanner.central_charges is not None and scanner.central_charges.shape == (4,)
assert scanner.charge_labels is not None and len(scanner.charge_labels) == 4
print(f"  ratio_charges shape: {scanner.ratio_charges.shape}")
for oi, lbl in enumerate(scanner.charge_labels):
    pct = 100 * scanner.sigma_before_charges[oi] / abs(scanner.central_charges[oi])
    print(f"  {lbl}: ratio={scanner.ratio_charges[oi].flatten()}  σ/μ={pct:.1f}%")

# ── Step 6: plot_charges() ────────────────────────────────────────────────────
print("\n=== Step 6: plot_charges() ===")
fig_ch = scanner.plot_charges(save=True)
expected_ch = os.path.join(OUT, f'heatmap_charges_{_suffix}.pdf')
assert os.path.exists(expected_ch), f"Missing: {expected_ch}"
print(f"  Charges plot saved OK → {expected_ch}")

# ── Step 7: load_charges() in fresh object ────────────────────────────────────
print("\n=== Step 7: load_charges() fresh object ===")
s3 = WindowMomentScanner(CFG)
s3.load_charges()
assert s3.ratio_charges is not None
assert s3.ratio_charges.shape == (4, 2, 1)
assert s3.sigma_before_charges is not None
assert all(
    abs(a - b) < 1e-10
    for a, b in zip(s3.ratio_charges.flatten(), scanner.ratio_charges.flatten())
    if not (a != a or b != b)
), "load_charges values differ from run_charges!"
print("  Loaded values match run_charges values.")

fig3 = s3.plot_charges(save=False)
print("  plot_charges() from loaded data: OK")

shutil.rmtree(OUT)
print(f"\nSUCCESS — all assertions passed.")

