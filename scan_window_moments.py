#!/home/daniel/miniconda3/envs/apfelpp/bin/python3
"""
scan_window_moments.py — 2D scan of Gaussian window moment constraints via ePump.

Terminal usage:
    python scan_window_moments.py runcard.py

Notebook usage:
    from scan_window_moments import WindowMomentScanner
    scanner = WindowMomentScanner('runcard.py')   # or pass a cfg dict directly
    scanner.run()
    fig = scanner.plot()
"""

import os
import sys
import importlib.util
import numpy as np
import matplotlib
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


def _bin_edges(arr):
    """pcolormesh bin edges (len+1) for a sorted 1D array of cell centres."""
    arr = np.asarray(arr, dtype=float)
    if len(arr) == 1:
        return np.array([arr[0] - 0.005, arr[0] + 0.005])
    diffs = np.diff(arr)
    return np.concatenate([
        [arr[0] - diffs[0] / 2],
        (arr[:-1] + arr[1:]) / 2,
        [arr[-1] + diffs[-1] / 2],
    ])


class WindowMomentScanner:
    """
    2D scan of Gaussian window moment constraints via ePump profiling.

    Parameters
    ----------
    cfg : dict or str
        Configuration dict, or path to a runcard.py file that defines ``cfg``.
    """

    def __init__(self, cfg):
        if isinstance(cfg, (str, os.PathLike)):
            cfg = load_runcard(str(cfg))
        self.cfg = cfg

        # Results — populated by run()
        self.midpoints = None
        self.widths    = None
        self.ratio     = None
        self.boundary  = None

        # Internal state — populated by setup()
        self._pdf_name = None
        self._mc2h_dir = None

    # ------------------------------------------------------------------
    def setup(self):
        """Configure LHAPDF paths and convert MC replicas to Hessian (once)."""
        cfg = self.cfg
        os.makedirs(cfg['output_dir'], exist_ok=True)
        setup_lhapdf_path(cfg.get('lhapdf_path'))

        pdf_name = cfg['pdf']
        err_type = detect_pdf_error_type(pdf_name)
        if err_type not in ('replicas', 'mc'):
            print(f"PDF '{pdf_name}' is {err_type!r} — no conversion needed.")
            self._pdf_name = pdf_name
            self._mc2h_dir = None
        else:
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
            self._pdf_name = hessian_name
            self._mc2h_dir = mc2h_dir

        return self

    # ------------------------------------------------------------------
    def run(self):
        """Run the full 2D scan. Calls setup() automatically if not already done."""
        if self._pdf_name is None:
            self.setup()

        cfg        = self.cfg
        scan_points = cfg['scan_points']
        flavor     = cfg['flavor']
        Q2         = float(cfg['Q2'])
        nx         = int(cfg['nx'])
        moment     = int(cfg.get('moment', 1))
        weight     = cfg.get('weight', 'gaussian')
        output_dir = os.path.abspath(cfg['output_dir'])
        epump_path = os.path.abspath(cfg.get('epump_path', './ePump_kp20221218/src/UpdatePDFs'))
        pdf_name   = self._pdf_name
        mc2h_dir   = self._mc2h_dir

        midpoints = sorted(set(float(p[0]) for p in scan_points))
        widths    = sorted(set(float(p[1]) for p in scan_points))
        mid_idx   = {v: i for i, v in enumerate(midpoints)}
        wid_idx   = {v: j for j, v in enumerate(widths)}

        ratio    = np.full((len(midpoints), len(widths)), np.nan)
        boundary = np.zeros((len(midpoints), len(widths)), dtype=bool)

        parsed_terms = parse_flavor_expression(flavor)

        # Load base set once; reused for all grid points
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
                weight_type=weight, moment=moment,
            )
            stat_err = rel_unc * abs(central_val)

            label    = f"mid_{x0:.4f}_wid_{w:.4f}"
            run_name = os.path.join(output_dir, label, label)

            clip_tag = " [BOUNDARY-CLIPPED]" if clipped else ""
            print(f"\n({idx}/{n_total}) [{label}]  x=[{xmin:.4f}, {xmax:.4f}]"
                  f"  central={central_val:.6g}  stat={stat_err:.6g}{clip_tag}")

            if mc2h_dir and mc2h_dir not in lhapdf.paths():
                lhapdf.setPaths([mc2h_dir] + lhapdf.paths())

            ep = EProfiler(pdf_name, run_name, epump_path=epump_path,
                           lhapdf_path=mc2h_dir)
            ep.pdf_set     = base_set
            ep.pdf_members = base_members

            ep.add_measurement(
                x=x0, Q2=Q2, value=central_val, stat=stat_err,
                obs_type='moment', flavor=flavor,
                xmin=xmin, xmax=xmax, nx=nx,
                weight=weight, moment=moment,
            )
            ep.generate_files()
            ep.run()
            print(f"  Profiling complete ({idx}/{n_total})")

            kwargs = dict(weight_type=weight, moment=moment)
            n_orig = len(ep.pdf_members)
            orig_vals = []
            for k, m in enumerate(ep.pdf_members, 1):
                print(f"  original  {k}/{n_orig}", end='\r', flush=True)
                orig_vals.append(compute_integrated_moment(
                    m, parsed_terms, xmin, xmax, nx, Q2, **kwargs))
            print()
            orig_vals = np.array(orig_vals)

            n_prof = len(ep.profiled_members)
            prof_vals = []
            for k, m in enumerate(ep.profiled_members, 1):
                print(f"  profiled  {k}/{n_prof}", end='\r', flush=True)
                prof_vals.append(compute_integrated_moment(
                    m, parsed_terms, xmin, xmax, nx, Q2, **kwargs))
            print()
            prof_vals = np.array(prof_vals)

            o = ep.pdf_set.uncertainty(orig_vals.tolist())
            p = ep.profiled_set.uncertainty(prof_vals.tolist())
            sigma_b = (o.errminus + o.errplus) / 2.0
            sigma_a = (p.errminus + p.errplus) / 2.0
            r = sigma_a / sigma_b if sigma_b > 0 else np.nan
            ratio[mid_idx[x0], wid_idx[w]] = r
            print(f"  σ_before={sigma_b:.5g}  σ_after={sigma_a:.5g}  ratio={r:.4f}")

        self.midpoints = midpoints
        self.widths    = widths
        self.ratio     = ratio
        self.boundary  = boundary

        results_path = os.path.join(output_dir, 'results.npz')
        np.savez(results_path,
                 midpoints=self.midpoints,
                 widths=self.widths,
                 ratio=self.ratio,
                 boundary=self.boundary,
                 pdf_name=np.array(self._pdf_name),
                 mc2h_dir=np.array(self._mc2h_dir or ''))
        print(f"Results saved → {results_path}")
        return self

    # ------------------------------------------------------------------
    def load(self):
        """Load ratio/boundary arrays from a previous run() call.  Enables plot()."""
        results_path = os.path.join(os.path.abspath(self.cfg['output_dir']), 'results.npz')
        if not os.path.exists(results_path):
            raise FileNotFoundError(
                f"No saved results at {results_path}. Run the scan first.")
        data = np.load(results_path, allow_pickle=True)
        self.midpoints = data['midpoints'].tolist()
        self.widths    = data['widths'].tolist()
        self.ratio     = data['ratio']
        self.boundary  = data['boundary']
        print(f"Results loaded ← {results_path}")
        return self

    # ------------------------------------------------------------------
    def recompute(self):
        """Recompute ratio from existing profiled PDFs on disk (no ePump re-run).

        Useful when you want to evaluate a different flavor or Q2 without
        re-running the full scan.  Requires a previous run() in the same output_dir.
        """
        output_dir   = os.path.abspath(self.cfg['output_dir'])
        results_path = os.path.join(output_dir, 'results.npz')
        if not os.path.exists(results_path):
            raise FileNotFoundError(
                f"No saved state at {results_path}. Run the scan first.")

        # Restore pdf_name / mc2h_dir — avoids re-running MC→Hessian conversion
        saved = np.load(results_path, allow_pickle=True)
        self._pdf_name = str(saved['pdf_name'])
        mc2h_str       = str(saved['mc2h_dir'])
        self._mc2h_dir = mc2h_str if mc2h_str else None

        setup_lhapdf_path(self.cfg.get('lhapdf_path'))
        if self._mc2h_dir:
            setup_lhapdf_path(custom_path=self._mc2h_dir)
            lhapdf.setPaths([self._mc2h_dir] + lhapdf.paths())

        cfg         = self.cfg
        scan_points = cfg['scan_points']
        flavor      = cfg['flavor']
        Q2          = float(cfg['Q2'])
        nx          = int(cfg['nx'])
        moment      = int(cfg.get('moment', 1))
        weight      = cfg.get('weight', 'gaussian')
        pdf_name    = self._pdf_name

        midpoints = sorted(set(float(p[0]) for p in scan_points))
        widths    = sorted(set(float(p[1]) for p in scan_points))
        mid_idx   = {v: i for i, v in enumerate(midpoints)}
        wid_idx   = {v: j for j, v in enumerate(widths)}

        ratio    = np.full((len(midpoints), len(widths)), np.nan)
        boundary = np.zeros((len(midpoints), len(widths)), dtype=bool)

        parsed_terms = parse_flavor_expression(flavor)
        base_set     = lhapdf.getPDFSet(pdf_name)
        base_members = base_set.mkPDFs()

        n_total = len(scan_points)
        for idx, row in enumerate(scan_points, 1):
            x0, w = float(row[0]), float(row[1])
            xmin = max(1e-4, x0 - w / 2)
            xmax = min(0.999, x0 + w / 2)
            boundary[mid_idx[x0], wid_idx[w]] = (x0 - w / 2 < 1e-4) or (x0 + w / 2 > 0.999)

            label   = f"mid_{x0:.4f}_wid_{w:.4f}"
            run_dir = os.path.join(output_dir, label)
            print(f"\n({idx}/{n_total}) [{label}] — loading profiled set …")

            if not os.path.isdir(os.path.join(run_dir, label)):
                print(f"  WARNING: profiled set not found at {run_dir}/{label}/ — skipping.")
                continue

            if run_dir not in lhapdf.paths():
                lhapdf.setPaths([run_dir] + lhapdf.paths())
            profiled_set     = lhapdf.getPDFSet(label)
            profiled_members = profiled_set.mkPDFs()

            kwargs    = dict(weight_type=weight, moment=moment)
            n_base    = len(base_members)
            n_prof    = len(profiled_members)
            orig_vals = []
            for k, m in enumerate(base_members, 1):
                print(f"  original  {k}/{n_base}", end='\r', flush=True)
                orig_vals.append(compute_integrated_moment(
                    m, parsed_terms, xmin, xmax, nx, Q2, **kwargs))
            print()
            prof_vals = []
            for k, m in enumerate(profiled_members, 1):
                print(f"  profiled  {k}/{n_prof}", end='\r', flush=True)
                prof_vals.append(compute_integrated_moment(
                    m, parsed_terms, xmin, xmax, nx, Q2, **kwargs))
            print()

            o = base_set.uncertainty(orig_vals)
            p = profiled_set.uncertainty(prof_vals)
            sigma_b = (o.errminus + o.errplus) / 2.0
            sigma_a = (p.errminus + p.errplus) / 2.0
            r = sigma_a / sigma_b if sigma_b > 0 else np.nan
            ratio[mid_idx[x0], wid_idx[w]] = r
            print(f"  σ_before={sigma_b:.5g}  σ_after={sigma_a:.5g}  ratio={r:.4f}")

        self.midpoints = midpoints
        self.widths    = widths
        self.ratio     = ratio
        self.boundary  = boundary

        np.savez(results_path,
                 midpoints=self.midpoints, widths=self.widths,
                 ratio=self.ratio, boundary=self.boundary,
                 pdf_name=np.array(self._pdf_name),
                 mc2h_dir=np.array(self._mc2h_dir or ''))
        print(f"Results saved → {results_path}")
        return self

    # ------------------------------------------------------------------
    def plot(self, save=True):
        """
        Generate the heat map.

        Parameters
        ----------
        save : bool
            Write the figure to cfg['output_plot'] (default True).

        Returns
        -------
        matplotlib.figure.Figure
            The figure, so Jupyter can display it inline.
        """
        if self.ratio is None:
            raise RuntimeError("Call run() before plot().")

        import matplotlib.pyplot as plt  # lazy import — respects notebook backend

        cfg = self.cfg
        M_e = _bin_edges(np.array(self.midpoints))
        W_e = _bin_edges(np.array(self.widths))

        moment = int(cfg.get('moment', 1))
        weight = cfg.get('weight', 'gaussian')
        obs_label = rf'$g_{{{moment}}}$' if weight == 'gaussian' else rf'$a_{{{moment}}}$'

        fig, ax = plt.subplots(figsize=(9, 6))
        masked = np.ma.masked_invalid(self.ratio)
        cm = ax.pcolormesh(W_e, M_e, masked, cmap='plasma_r', vmin=0, vmax=1)
        plt.colorbar(cm, ax=ax,
                     label=r'$\sigma_\mathrm{after}\ /\ \sigma_\mathrm{before}$')

        for i in range(len(self.midpoints)):
            for j in range(len(self.widths)):
                if self.boundary[i, j]:
                    ax.add_patch(matplotlib.patches.Rectangle(
                        (W_e[j], M_e[i]),
                        W_e[j + 1] - W_e[j], M_e[i + 1] - M_e[i],
                        fill=False, hatch='///', edgecolor='white', linewidth=0.5,
                    ))

        ax.set_xlabel('Window width  $w$')
        ax.set_ylabel('Window midpoint  $x_0$')
        ax.set_title(
            f"{cfg['pdf']}    {cfg['flavor']}    {obs_label}    $Q^2 = {cfg['Q2']}$ GeV$^2$"
        )
        plt.tight_layout()

        if save:
            out = os.path.join(os.path.abspath(cfg['output_dir']), 'heatmap.pdf')
            fig.savefig(out, dpi=150)
            print(f"Heat map → {out}")

        return fig


# ── Terminal entry point ───────────────────────────────────────────────────────

def main():
    matplotlib.use('Agg')   # non-interactive backend for terminal use

    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} runcard.py [load|recompute]", file=sys.stderr)
        sys.exit(1)

    scanner = WindowMomentScanner(sys.argv[1])
    cmd = sys.argv[2] if len(sys.argv) >= 3 else 'run'
    if cmd == 'load':
        scanner.load()
    elif cmd == 'recompute':
        scanner.recompute()
    else:
        scanner.run()
    scanner.plot(save=True)


if __name__ == '__main__':
    main()
