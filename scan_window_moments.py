#!/home/daniel/miniconda3/envs/apfelpp/bin/python3
"""
scan_window_moments.py — 2D scan of Gaussian window moment constraints via ePump.

Usage:
    python scan_window_moments.py runcard.py

The runcard is a Python file that defines a dict named `cfg`.
See runcard.py for the template.
"""

import os
import sys
import importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches

import lhapdf

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)
from e_profiler import (
    setup_lhapdf_path,
    detect_pdf_error_type,
    convert_mc_to_hessian,
    parse_flavor_expression,
    compute_integrated_moment,
    EProfiler,
)


def load_runcard(path):
    spec = importlib.util.spec_from_file_location("runcard", os.path.abspath(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.cfg


def convert_if_needed(cfg):
    """Convert MC replica set to asymmetric Hessian once.
    Returns (hessian_pdf_name, mc2h_dir | None).
    """
    pdf_name = cfg['pdf']
    err_type = detect_pdf_error_type(pdf_name)
    if err_type not in ('replicas', 'mc'):
        print(f"PDF '{pdf_name}' is {err_type!r} — no conversion needed.")
        return pdf_name, None

    mc2h_dir = os.path.join(os.path.abspath(cfg['output_dir']), '_hessian')
    os.makedirs(mc2h_dir, exist_ok=True)
    print(f"Converting '{pdf_name}' (MC replicas) → asymmetric Hessian in {mc2h_dir} …")
    hessian_name, _ = convert_mc_to_hessian(
        pdf_name,
        neig=cfg.get('mc2h_neig', 50),
        Q=float(cfg.get('mc2h_Q', 1.0)),
        epsilon=float(cfg.get('mc2h_epsilon', 1000.0)),
        output_dir=mc2h_dir,
        max_nf=int(cfg.get('mc2h_max_nf', 3)),
    )
    setup_lhapdf_path(custom_path=mc2h_dir)
    lhapdf.setPaths([mc2h_dir] + lhapdf.paths())
    print(f"  → '{hessian_name}'")
    return hessian_name, mc2h_dir


def bin_edges(arr):
    """Compute pcolormesh bin edges (len+1) for a sorted 1D array of cell centres."""
    arr = np.asarray(arr, dtype=float)
    if len(arr) == 1:
        return np.array([arr[0] - 0.005, arr[0] + 0.005])
    diffs = np.diff(arr)
    return np.concatenate([
        [arr[0] - diffs[0] / 2],
        (arr[:-1] + arr[1:]) / 2,
        [arr[-1] + diffs[-1] / 2],
    ])


def run_scan(cfg, pdf_name, mc2h_dir):
    scan_points = cfg['scan_points']
    flavor      = cfg['flavor']
    Q2          = float(cfg['Q2'])
    nx          = int(cfg['nx'])
    moment      = int(cfg.get('moment', 1))
    output_dir  = os.path.abspath(cfg['output_dir'])
    epump_path  = os.path.abspath(
        cfg.get('epump_path', './ePump_kp20221218/src/UpdatePDFs')
    )
    lhapdf_arg  = mc2h_dir  # None if no conversion was needed

    midpoints = sorted(set(float(p[0]) for p in scan_points))
    widths    = sorted(set(float(p[1]) for p in scan_points))
    mid_idx   = {v: i for i, v in enumerate(midpoints)}
    wid_idx   = {v: j for j, v in enumerate(widths)}

    ratio    = np.full((len(midpoints), len(widths)), np.nan)
    boundary = np.zeros((len(midpoints), len(widths)), dtype=bool)

    parsed_terms = parse_flavor_expression(flavor)

    # Load base set once; reused for all grid points to avoid per-iteration reloading
    base_set     = lhapdf.getPDFSet(pdf_name)
    base_central = base_set.mkPDF(0)
    base_members = base_set.mkPDFs()

    n_total = len(scan_points)
    for idx, row in enumerate(scan_points, 1):
        x0, w, rel_unc = float(row[0]), float(row[1]), float(row[2])
        xmin = max(1e-4, x0 - w / 2)
        xmax = min(0.999, x0 + w / 2)

        clipped = (x0 - w / 2 < 1e-4) or (x0 + w / 2 > 0.999)
        boundary[mid_idx[x0], wid_idx[w]] = clipped

        central_val = compute_integrated_moment(
            base_central, parsed_terms, xmin, xmax, nx, Q2,
            weight_type='gaussian', moment=moment,
        )
        stat_err = rel_unc * abs(central_val)

        label    = f"mid_{x0:.4f}_wid_{w:.4f}"
        run_name = os.path.join(output_dir, label, label)

        clip_tag = " [BOUNDARY-CLIPPED]" if clipped else ""
        print(f"\n({idx}/{n_total}) [{label}]  x=[{xmin:.4f}, {xmax:.4f}]"
              f"  central={central_val:.6g}  stat={stat_err:.6g}{clip_tag}")

        # Ensure mc2h_dir stays in LHAPDF's path list after repeated ep.run() prepends
        if mc2h_dir and mc2h_dir not in lhapdf.paths():
            lhapdf.setPaths([mc2h_dir] + lhapdf.paths())

        ep = EProfiler(pdf_name, run_name, epump_path=epump_path,
                       lhapdf_path=lhapdf_arg)
        # Reuse pre-loaded base set so EProfiler doesn't reload all members each iteration
        ep.pdf_set     = base_set
        ep.pdf_members = base_members

        ep.add_measurement(
            x=x0, Q2=Q2, value=central_val, stat=stat_err,
            obs_type='moment', flavor=flavor,
            xmin=xmin, xmax=xmax, nx=nx,
            weight='gaussian', moment=moment,
        )
        ep.generate_files()
        ep.run()
        print(f"  Profiling complete ({idx}/{n_total})")

        kwargs = dict(weight_type='gaussian', moment=moment)
        orig_vals = np.array([
            compute_integrated_moment(m, parsed_terms, xmin, xmax, nx, Q2, **kwargs)
            for m in ep.pdf_members
        ])
        prof_vals = np.array([
            compute_integrated_moment(m, parsed_terms, xmin, xmax, nx, Q2, **kwargs)
            for m in ep.profiled_members
        ])

        o = ep.pdf_set.uncertainty(orig_vals.tolist())
        p = ep.profiled_set.uncertainty(prof_vals.tolist())
        sigma_b = (o.errminus + o.errplus) / 2.0
        sigma_a = (p.errminus + p.errplus) / 2.0
        r = sigma_a / sigma_b if sigma_b > 0 else np.nan
        ratio[mid_idx[x0], wid_idx[w]] = r
        print(f"  σ_before={sigma_b:.5g}  σ_after={sigma_a:.5g}  ratio={r:.4f}")

    return midpoints, widths, ratio, boundary


def plot_heatmap(cfg, midpoints, widths, ratio, boundary):
    M   = np.array(midpoints)
    W   = np.array(widths)
    M_e = bin_edges(M)
    W_e = bin_edges(W)

    fig, ax = plt.subplots(figsize=(9, 6))
    masked = np.ma.masked_invalid(ratio)
    cm = ax.pcolormesh(W_e, M_e, masked, cmap='plasma_r', vmin=0, vmax=1)
    plt.colorbar(cm, ax=ax,
                 label=r'$\sigma_\mathrm{after}\ /\ \sigma_\mathrm{before}$')

    for i in range(len(midpoints)):
        for j in range(len(widths)):
            if boundary[i, j]:
                ax.add_patch(matplotlib.patches.Rectangle(
                    (W_e[j], M_e[i]),
                    W_e[j + 1] - W_e[j], M_e[i + 1] - M_e[i],
                    fill=False, hatch='///', edgecolor='white', linewidth=0.5,
                ))

    ax.set_xlabel('Window width  $w$')
    ax.set_ylabel('Window midpoint  $x_0$')
    ax.set_title(
        f"{cfg['pdf']}    {cfg['flavor']}    $Q^2 = {cfg['Q2']}$ GeV$^2$"
    )
    plt.tight_layout()

    out = cfg['output_plot']
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    plt.savefig(out, dpi=150)
    print(f"\nHeat map → {out}")


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} runcard.py", file=sys.stderr)
        sys.exit(1)

    cfg = load_runcard(sys.argv[1])
    os.makedirs(cfg['output_dir'], exist_ok=True)
    setup_lhapdf_path(cfg.get('lhapdf_path'))

    pdf_name, mc2h_dir = convert_if_needed(cfg)
    midpoints, widths, ratio, boundary = run_scan(cfg, pdf_name, mc2h_dir)
    plot_heatmap(cfg, midpoints, widths, ratio, boundary)


if __name__ == '__main__':
    main()
