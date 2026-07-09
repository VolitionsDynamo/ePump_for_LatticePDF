# runcard.py — configuration for scan_window_moments.py
# Run with: python scan_window_moments.py runcard.py

cfg = {

    # ── PDF settings ──────────────────────────────────────────────────────────
    "pdf":    "CT18NNLO",    # base PDF set (Hessian or MC replicas)
    "flavor": "u-d",         # algebraic combination: "u-d", "u", "2*u-d", …
    "Q2":     4.0,           # scale Q² in GeV²
    "nx":     100,           # x integration points
    "moment": 1,             # moment order n (x^n weighting inside the window)
    "weight": "gaussian",   # weight function: "gaussian" → g_n, "1" → a_n (flat)

    # ── MC-to-Hessian (ignored if pdf is already Hessian) ─────────────────────
    "mc2h_neig":    50,
    "mc2h_Q":       1.0,
    "mc2h_epsilon": 1000.0,
    "mc2h_max_nf":  3,

    # ── Scan grid: rows of [midpoint, window_width, relative_uncertainty] ─────
    # Each row is one independent ePump profiling run.
    # Unique (midpoint, width) pairs are pivoted into a 2D heat map.
    # Cells with x0 ± w/2 outside [0, 1] are automatically clipped and hatched.
    "scan_points": [
        # [x0,   w,    ε ]
        [0.20,  0.05, 0.10],
        [0.20,  0.10, 0.10],
        [0.20,  0.15, 0.10],
        [0.30,  0.05, 0.10],
        [0.30,  0.10, 0.10],
        [0.30,  0.15, 0.10],
        [0.40,  0.05, 0.10],
        [0.40,  0.10, 0.10],
        [0.40,  0.15, 0.10],
    ],

    # ── Output ────────────────────────────────────────────────────────────────
    "output_dir":  "scan_results",   # heatmap.pdf is saved here automatically
    "epump_path":  "./ePump_kp20221218/src/UpdatePDFs",
    "lhapdf_path": None,     # optional extra LHAPDF data path

}

# ── Helper: build a regular grid programmatically ─────────────────────────────
# Uncomment and adjust to replace the hand-written scan_points above.
#
# import numpy as np
# midpoints = np.linspace(0.10, 0.60, 6)
# widths    = np.linspace(0.04, 0.20, 5)
# epsilon   = 0.10   # uniform relative uncertainty
# cfg["scan_points"] = [
#     [float(x0), float(w), epsilon]
#     for x0 in midpoints
#     for w  in widths
# ]
