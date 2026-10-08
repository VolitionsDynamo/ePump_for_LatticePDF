# runcard_anchored.py — configuration for AnchoredWindowScanner
# Run with: python scan_window_moments.py runcard_anchored.py anchored_scan

cfg = {

    # ── PDF settings ──────────────────────────────────────────────────────────
    "pdf":    "CT18NNLO",    # base PDF set (Hessian or MC replicas)
    "flavor": "u-d",         # algebraic combination: "u-d", "u", "2*u-d", …
    "Q2":     4.0,           # scale Q² in GeV²
    "nx":     100,           # x integration points
    "moment": 1,             # moment order n
    "weight": "gaussian",    # "gaussian" (g_n) or "1" (a_n flat)

    # ── MC-to-Hessian (ignored if pdf is already Hessian) ─────────────────────
    "mc2h_neig":    50,
    "mc2h_Q":       1.0,
    "mc2h_epsilon": 1000.0,
    "mc2h_max_nf":  3,

    # ── Anchor: the fixed first measurement (best window from scan 1) ──────────
    # Width w is also used for all scan windows.
    "anchor": {"x0": 0.30, "w": 0.10},

    # ── Fractional uncertainty applied to both the anchor and scan measurements ─
    "rel_unc": 0.10,

    # ── Scan midpoints: x0 values for the second window ───────────────────────
    # Width is always anchor["w"].  Midpoints where |x0 - anchor_x0| < anchor_w
    # (overlapping windows) are automatically skipped.
    "scan_midpoints": [
        0.05, 0.10, 0.15, 0.20,
        # gap: 0.25, 0.30, 0.35 overlap anchor at x0=0.30, w=0.10
        0.40, 0.45, 0.50, 0.55, 0.60,
    ],

    # ── Tensor-charge observables (for future run_charges support) ────────────
    "charge_observables": [
        {"label": r"$g_T^{u-d}$",
         "flavor": "u-d",           "moment": 0, "weight": "1"},
        {"label": r"$g_T^u$",
         "flavor": "u",             "moment": 0, "weight": "1"},
        {"label": r"$g_T^d$",
         "flavor": "d",             "moment": 0, "weight": "1"},
        {"label": r"$\langle x(\delta u^+ - \delta d^+)\rangle$",
         "flavor": "u+ubar-d-dbar", "moment": 1, "weight": "1"},
    ],

    # ── Output ────────────────────────────────────────────────────────────────
    "pdf_label":   "CT18NNLO",
    "output_dir":  "scan_results_anchored",
    "epump_path":  "./ePump_kp20221218/src/UpdatePDFs",
    "lhapdf_path": None,

}
