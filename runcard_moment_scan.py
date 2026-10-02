# runcard_moment_scan.py — configuration for MomentAccumulationScanner
# Run with: python scan_window_moments.py runcard_moment_scan.py moment_scan

cfg = {

    # ── PDF settings ──────────────────────────────────────────────────────────
    "pdf":    "CT18NNLO",    # base PDF set (Hessian or MC replicas)
    "flavor": "u-d",         # algebraic combination: "u-d", "u", "2*u-d", …
    "Q2":     4.0,           # scale Q² in GeV²
    "nx":     100,           # x integration points

    # ── MC-to-Hessian (ignored if pdf is already Hessian) ─────────────────────
    "mc2h_neig":    50,
    "mc2h_Q":       1.0,
    "mc2h_epsilon": 1000.0,
    "mc2h_max_nf":  3,

    # ── Moment accumulation scan ───────────────────────────────────────────────
    # Profile a fixed Gaussian window with moments 1, then 1+2, then 1+2+3, …
    # tracking how σ_after/σ_before improves for window moments, tensor charges,
    # and full polynomial moments as simultaneous constraints are added.
    "moment_scan": {
        "x0":      0.30,        # window midpoint
        "w":       0.10,        # window width → xmin = x0 - w/2, xmax = x0 + w/2
        "weight":  "gaussian",  # "gaussian" (g_n) or "1" (a_n flat)
        "max_n":   3,           # accumulate moments 1 … max_n
        "rel_unc": 0.10,        # fractional pseudo-data uncertainty on each moment

        # corr: None → measurements are uncorrelated (stat error only).
        # To encode correlated uncertainties supply a max_n × max_n correlation
        # matrix as a nested list; the diagonal must be 1.  Example for max_n=3:
        #
        # "corr": [
        #     [1.00, 0.50, 0.30],
        #     [0.50, 1.00, 0.50],
        #     [0.30, 0.50, 1.00],
        # ],
        "corr":    None,
    },

    # ── Full-moment integration (tracked after profiling) ─────────────────────
    "moment_xmin":  1e-4,   # lower x bound for ∫ x^n f(x) dx  (n = 0,1,2,3)
    "moment_xmax":  0.999,  # upper x bound

    # ── Tensor-charge observables (tracked after profiling) ───────────────────
    "charge_xmin": 1e-4,
    "charge_xmax": 0.999,
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
    "output_dir":  "scan_results_moment_scan",
    "epump_path":  "./ePump_kp20221218/src/UpdatePDFs",
    "lhapdf_path": None,

}
