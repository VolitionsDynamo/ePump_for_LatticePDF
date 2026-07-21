#!/home/daniel/miniconda3/envs/apfelpp/bin/python3
"""
scan_window_moments.py — 2D scan of Gaussian window moment constraints via ePump.

Terminal usage:
    python scan_window_moments.py runcard.py
    python scan_window_moments.py runcard.py [run|load|recompute|moments|load_moments]

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
        self.midpoints       = None
        self.widths          = None
        self.ratio           = None
        self.boundary        = None
        self.sigma_before_ww = None  # shape (n_mid, n_wid)
        self.central_ww      = None  # shape (n_mid, n_wid)

        # Window-to-moment results — populated by run_moments()
        self.ratio_moments        = None  # shape (4, n_mid, n_wid)
        self.sigma_before_moments = None  # shape (4,)
        self.central_moments      = None  # shape (4,)

        # Tensor-charge results — populated by run_charges()
        self.ratio_charges        = None  # shape (n_obs, n_mid, n_wid)
        self.sigma_before_charges = None  # shape (n_obs,)
        self.central_charges      = None  # shape (n_obs,)
        self.charge_labels        = None  # list of label strings

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
    def run(self, force=False):
        """Run the full 2D scan. Calls setup() automatically if not already done.

        Parameters
        ----------
        force : bool
            If True, ignore any previously saved results and rerun ePump for
            every cell from scratch.  If False (default), resume from partial
            results saved by an earlier run().
        """
        if self._pdf_name is None:
            self.setup()

        cfg         = self.cfg
        scan_points = cfg['scan_points']
        flavor      = cfg['flavor']
        Q2          = float(cfg['Q2'])
        nx          = int(cfg['nx'])
        moment      = int(cfg.get('moment', 1))
        weight      = cfg.get('weight', 'gaussian')
        output_dir  = os.path.abspath(cfg['output_dir'])
        epump_path  = os.path.abspath(cfg.get('epump_path', './ePump_kp20221218/src/UpdatePDFs'))
        pdf_name    = self._pdf_name
        mc2h_dir    = self._mc2h_dir

        midpoints = sorted(set(float(p[0]) for p in scan_points))
        widths    = sorted(set(float(p[1]) for p in scan_points))
        mid_idx   = {v: i for i, v in enumerate(midpoints)}
        wid_idx   = {v: j for j, v in enumerate(widths)}

        ratio           = np.full((len(midpoints), len(widths)), np.nan)
        boundary        = np.zeros((len(midpoints), len(widths)), dtype=bool)
        sigma_before_ww = np.full((len(midpoints), len(widths)), np.nan)
        central_ww      = np.full((len(midpoints), len(widths)), np.nan)

        results_path = os.path.join(output_dir, 'results.npz')
        if not force and os.path.exists(results_path):
            saved = np.load(results_path, allow_pickle=True)
            if (np.array_equal(saved['midpoints'], midpoints) and
                    np.array_equal(saved['widths'], widths)):
                ratio[...]    = saved['ratio']
                boundary[...] = saved['boundary']
                if 'sigma_before_ww' in saved:
                    sigma_before_ww[...] = saved['sigma_before_ww']
                if 'central_ww' in saved:
                    central_ww[...]      = saved['central_ww']
                n_done = int(np.sum(np.isfinite(ratio)))
                if n_done:
                    print(f"Resuming: {n_done}/{ratio.size} cells already done.")

        parsed_terms = parse_flavor_expression(flavor)

        # Load base set once; reused for all grid points
        base_set     = lhapdf.getPDFSet(pdf_name)
        base_central = base_set.mkPDF(0)
        base_members = base_set.mkPDFs()

        n_total = len(scan_points)
        for idx, row in enumerate(scan_points, 1):
            x0, w, rel_unc = float(row[0]), float(row[1]), float(row[2])
            i, j = mid_idx[x0], wid_idx[w]

            if np.isfinite(ratio[i, j]):
                print(f"({idx}/{n_total}) [mid_{x0:.4f}_wid_{w:.4f}]  ← already done")
                continue

            xmin    = max(1e-4, x0 - w / 2)
            xmax    = min(0.999, x0 + w / 2)
            clipped = (x0 - w / 2 < 1e-4) or (x0 + w / 2 > 0.999)
            boundary[i, j] = clipped

            central_val = compute_integrated_moment(
                base_central, parsed_terms, xmin, xmax, nx, Q2,
                weight_type=weight, moment=moment,
            )
            stat_err = rel_unc * abs(central_val)

            label    = f"mid_{x0:.4f}_wid_{w:.4f}"
            run_name = os.path.join(output_dir, label, label)
            run_dir  = os.path.join(output_dir, label)

            if mc2h_dir and mc2h_dir not in lhapdf.paths():
                lhapdf.setPaths([mc2h_dir] + lhapdf.paths())

            already_profiled = os.path.isdir(os.path.join(run_dir, label))
            if already_profiled:
                if run_dir not in lhapdf.paths():
                    lhapdf.setPaths([run_dir] + lhapdf.paths())
                profiled_set = lhapdf.getPDFSet(label)
            else:
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
                profiled_set = ep.profiled_set

            kwargs    = dict(weight_type=weight, moment=moment)
            orig_vals = [compute_integrated_moment(
                             m, parsed_terms, xmin, xmax, nx, Q2, **kwargs)
                         for m in base_members]
            lhapdf.setVerbosity(0)
            prof_vals = [compute_integrated_moment(
                             profiled_set.mkPDF(k), parsed_terms, xmin, xmax, nx, Q2, **kwargs)
                         for k in range(profiled_set.size)]
            lhapdf.setVerbosity(1)

            o = base_set.uncertainty(orig_vals)
            p = profiled_set.uncertainty(prof_vals)
            sigma_b = (o.errminus + o.errplus) / 2.0
            sigma_a = (p.errminus + p.errplus) / 2.0
            r = sigma_a / sigma_b if sigma_b > 0 else np.nan
            ratio[i, j]           = r
            sigma_before_ww[i, j] = sigma_b
            central_ww[i, j]      = central_val

            clip_tag = "  [BOUNDARY-CLIPPED]" if clipped else ""
            print(f"({idx}/{n_total}) [{label}]  → ratio={r:.4f}  "
                  f"σ_before={sigma_b:.5g}  σ_after={sigma_a:.5g}{clip_tag}")
            np.savez(results_path,
                     midpoints=midpoints, widths=widths,
                     ratio=ratio, boundary=boundary,
                     sigma_before_ww=sigma_before_ww, central_ww=central_ww,
                     pdf_name=np.array(self._pdf_name),
                     mc2h_dir=np.array(self._mc2h_dir or ''))

        self.midpoints       = midpoints
        self.widths          = widths
        self.ratio           = ratio
        self.boundary        = boundary
        self.sigma_before_ww = sigma_before_ww
        self.central_ww      = central_ww
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
        self.midpoints       = data['midpoints'].tolist()
        self.widths          = data['widths'].tolist()
        self.ratio           = data['ratio']
        self.boundary        = data['boundary']
        self.sigma_before_ww = data['sigma_before_ww'] if 'sigma_before_ww' in data else None
        self.central_ww      = data['central_ww']      if 'central_ww'      in data else None
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

            if not os.path.isdir(os.path.join(run_dir, label)):
                print(f"({idx}/{n_total}) [{label}]  WARNING: profiled set not found — skipping.")
                continue

            if run_dir not in lhapdf.paths():
                lhapdf.setPaths([run_dir] + lhapdf.paths())
            profiled_set = lhapdf.getPDFSet(label)

            kwargs    = dict(weight_type=weight, moment=moment)
            orig_vals = [compute_integrated_moment(
                             m, parsed_terms, xmin, xmax, nx, Q2, **kwargs)
                         for m in base_members]
            lhapdf.setVerbosity(0)
            prof_vals = [compute_integrated_moment(
                             profiled_set.mkPDF(k), parsed_terms, xmin, xmax, nx, Q2, **kwargs)
                         for k in range(profiled_set.size)]
            lhapdf.setVerbosity(1)

            o = base_set.uncertainty(orig_vals)
            p = profiled_set.uncertainty(prof_vals)
            sigma_b = (o.errminus + o.errplus) / 2.0
            sigma_a = (p.errminus + p.errplus) / 2.0
            r = sigma_a / sigma_b if sigma_b > 0 else np.nan
            ratio[mid_idx[x0], wid_idx[w]] = r
            print(f"({idx}/{n_total}) [{label}]  → ratio={r:.4f}  σ_before={sigma_b:.5g}  σ_after={sigma_a:.5g}")

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
    def run_moments(self, force=False):
        """Compute window-to-moment ratios from existing profiled PDFs (no ePump re-run).

        For each scan point (x0, w), loads the profiled PDF set produced by run(),
        then evaluates the full polynomial moment ∫ x^n f(x) dx for n = 0, 1, 2, 3
        over [moment_xmin, moment_xmax] for both base and profiled PDFs.  Stores
        ratio_moments of shape (4, n_midpoints, n_widths).

        Parameters
        ----------
        force : bool
            If True, recompute all cells from scratch.  If False (default),
            resume from partial results saved by an earlier run_moments().
        """
        output_dir   = os.path.abspath(self.cfg['output_dir'])
        results_path = os.path.join(output_dir, 'results.npz')
        if not os.path.exists(results_path):
            raise FileNotFoundError(
                f"No saved state at {results_path}. Run the scan first.")

        saved = np.load(results_path, allow_pickle=True)
        self._pdf_name = str(saved['pdf_name'])
        mc2h_str       = str(saved['mc2h_dir'])
        self._mc2h_dir = mc2h_str if mc2h_str else None

        setup_lhapdf_path(self.cfg.get('lhapdf_path'))
        if self._mc2h_dir:
            setup_lhapdf_path(custom_path=self._mc2h_dir)
            lhapdf.setPaths([self._mc2h_dir] + lhapdf.paths())

        cfg          = self.cfg
        scan_points  = cfg['scan_points']
        flavor       = cfg['flavor']
        Q2           = float(cfg['Q2'])
        nx           = int(cfg['nx'])
        moment_xmin  = float(cfg.get('moment_xmin', 1e-4))
        moment_xmax  = float(cfg.get('moment_xmax', 0.999))
        pdf_name     = self._pdf_name

        midpoints = sorted(set(float(p[0]) for p in scan_points))
        widths    = sorted(set(float(p[1]) for p in scan_points))
        mid_idx   = {v: i for i, v in enumerate(midpoints)}
        wid_idx   = {v: j for j, v in enumerate(widths)}

        n_orders      = 4
        ratio_moments = np.full((n_orders, len(midpoints), len(widths)), np.nan)
        boundary      = saved['boundary'].copy()

        moments_path = os.path.join(output_dir, 'results_moments.npz')
        if not force and os.path.exists(moments_path):
            saved_m = np.load(moments_path, allow_pickle=True)
            if (np.array_equal(saved_m['midpoints'], midpoints) and
                    np.array_equal(saved_m['widths'], widths)):
                ratio_moments[...] = saved_m['ratio_moments']
                n_done = int(np.sum(np.all(np.isfinite(ratio_moments), axis=0)))
                if n_done:
                    print(f"Resuming moments: {n_done}/{len(midpoints)*len(widths)} cells already done.")

        parsed_terms = parse_flavor_expression(flavor)
        base_set     = lhapdf.getPDFSet(pdf_name)
        base_central = base_set.mkPDF(0)
        base_members = base_set.mkPDFs()

        # Compute base-PDF relative uncertainty for each n-moment once (constant across scan)
        sigma_before_moments = np.zeros(n_orders)
        central_moments_arr  = np.zeros(n_orders)
        for ni in range(n_orders):
            kwargs_base = dict(weight_type='1', moment=ni)
            base_vals = [compute_integrated_moment(
                             m, parsed_terms, moment_xmin, moment_xmax, nx, Q2, **kwargs_base)
                         for m in base_members]
            o_base = base_set.uncertainty(base_vals)
            sigma_before_moments[ni] = (o_base.errminus + o_base.errplus) / 2.0
            central_moments_arr[ni]  = compute_integrated_moment(
                base_central, parsed_terms, moment_xmin, moment_xmax, nx, Q2, **kwargs_base)

        n_total = len(scan_points)
        for idx, row in enumerate(scan_points, 1):
            x0, w = float(row[0]), float(row[1])
            i, j  = mid_idx[x0], wid_idx[w]

            if np.all(np.isfinite(ratio_moments[:, i, j])):
                print(f"({idx}/{n_total}) [mid_{x0:.4f}_wid_{w:.4f}]  ← already done")
                continue

            label   = f"mid_{x0:.4f}_wid_{w:.4f}"
            run_dir = os.path.join(output_dir, label)

            if not os.path.isdir(os.path.join(run_dir, label)):
                print(f"({idx}/{n_total}) [{label}]  WARNING: profiled set not found — skipping.")
                continue

            if run_dir not in lhapdf.paths():
                lhapdf.setPaths([run_dir] + lhapdf.paths())
            profiled_set = lhapdf.getPDFSet(label)

            ratio_strs = []
            lhapdf.setVerbosity(0)
            for ni in range(n_orders):
                kwargs    = dict(weight_type='1', moment=ni)
                orig_vals = [compute_integrated_moment(
                                 m, parsed_terms, moment_xmin, moment_xmax, nx, Q2, **kwargs)
                             for m in base_members]
                prof_vals = [compute_integrated_moment(
                                 profiled_set.mkPDF(k), parsed_terms,
                                 moment_xmin, moment_xmax, nx, Q2, **kwargs)
                             for k in range(profiled_set.size)]
                o = base_set.uncertainty(orig_vals)
                p = profiled_set.uncertainty(prof_vals)
                sigma_b = (o.errminus + o.errplus) / 2.0
                sigma_a = (p.errminus + p.errplus) / 2.0
                r = sigma_a / sigma_b if sigma_b > 0 else np.nan
                ratio_moments[ni, i, j] = r
                ratio_strs.append(f"n={ni}:{r:.3f}")
            lhapdf.setVerbosity(1)

            print(f"({idx}/{n_total}) [{label}]  → " + "  ".join(ratio_strs))
            np.savez(moments_path,
                     ratio_moments=ratio_moments,
                     sigma_before_moments=sigma_before_moments,
                     central_moments=central_moments_arr,
                     midpoints=np.array(midpoints),
                     widths=np.array(widths),
                     boundary=boundary,
                     pdf_name=np.array(self._pdf_name),
                     mc2h_dir=np.array(self._mc2h_dir or ''))

        self.midpoints            = midpoints
        self.widths               = widths
        self.boundary             = boundary
        self.ratio_moments        = ratio_moments
        self.sigma_before_moments = sigma_before_moments
        self.central_moments      = central_moments_arr
        print(f"Moment results saved → {moments_path}")
        return self

    # ------------------------------------------------------------------
    def load_moments(self):
        """Load window-to-moment results from a previous run_moments() call."""
        moments_path = os.path.join(os.path.abspath(self.cfg['output_dir']), 'results_moments.npz')
        if not os.path.exists(moments_path):
            raise FileNotFoundError(
                f"No saved moment results at {moments_path}. Run run_moments() first.")
        data = np.load(moments_path, allow_pickle=True)
        self.ratio_moments        = data['ratio_moments']
        self.midpoints            = data['midpoints'].tolist()
        self.widths               = data['widths'].tolist()
        self.boundary             = data['boundary']
        self.sigma_before_moments = data['sigma_before_moments'] if 'sigma_before_moments' in data else None
        self.central_moments      = data['central_moments']      if 'central_moments'      in data else None
        print(f"Moment results loaded ← {moments_path}")
        return self

    # ------------------------------------------------------------------
    def plot(self, save=True):
        """
        Generate the window-to-window heat map.

        Parameters
        ----------
        save : bool
            Write the figure to output_dir/heatmap.pdf (default True).

        Returns
        -------
        matplotlib.figure.Figure
        """
        if self.ratio is None:
            raise RuntimeError("Call run() before plot().")

        import matplotlib.pyplot as plt

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
                if (self.sigma_before_ww is not None and
                        self.central_ww is not None and
                        np.isfinite(self.sigma_before_ww[i, j]) and
                        abs(self.central_ww[i, j]) > 0):
                    pct = 100 * self.sigma_before_ww[i, j] / abs(self.central_ww[i, j])
                    cx = (W_e[j] + W_e[j + 1]) / 2
                    cy = (M_e[i] + M_e[i + 1]) / 2
                    ax.text(cx, cy, f"{pct:.0f}%",
                            ha='center', va='center', fontsize=7, color='white')

        ax.set_xlabel('Window width  $w$')
        ax.set_ylabel('Window midpoint  $x_0$')
        ax.set_title(
            f"{cfg['pdf']}    {cfg['flavor']}    {obs_label}    $Q^2 = {cfg['Q2']}$ GeV$^2$"
        )
        plt.tight_layout()

        if save:
            suffix = '_' + os.path.basename(os.path.abspath(cfg['output_dir']))
            out = os.path.join(os.path.abspath(cfg['output_dir']), f"heatmap{suffix}.pdf")
            fig.savefig(out, dpi=150)
            print(f"Heat map → {out}")

        return fig

    # ------------------------------------------------------------------
    def plot_moments(self, save=True):
        """
        Generate a 2×2 grid of window-to-moment heat maps (n = 0, 1, 2, 3).

        Each panel shows σ_after / σ_before for the full polynomial moment
        ∫ x^n f(x) dx, where the profiling was done with the window observable
        as input data.

        Parameters
        ----------
        save : bool
            Write the figure to output_dir/heatmap_moments[_runcard].pdf (default True).

        Returns
        -------
        matplotlib.figure.Figure
        """
        if self.ratio_moments is None:
            raise RuntimeError("Call run_moments() or load_moments() before plot_moments().")

        import matplotlib.pyplot as plt

        cfg = self.cfg
        M_e = _bin_edges(np.array(self.midpoints))
        W_e = _bin_edges(np.array(self.widths))

        moment = int(cfg.get('moment', 1))
        weight = cfg.get('weight', 'gaussian')
        obs_label = rf'$g_{{{moment}}}$' if weight == 'gaussian' else rf'$a_{{{moment}}}$'

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        axes_flat = axes.flatten()

        vmin, vmax = 0.0, 1.0
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

        for ni, ax in enumerate(axes_flat):
            masked = np.ma.masked_invalid(self.ratio_moments[ni])
            cm = ax.pcolormesh(W_e, M_e, masked, cmap='plasma_r', norm=norm)

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
            if (self.sigma_before_moments is not None and
                    self.central_moments is not None and
                    abs(self.central_moments[ni]) > 0):
                pct = 100 * self.sigma_before_moments[ni] / abs(self.central_moments[ni])
                ax.set_title(rf'$n = {ni}$   ($\sigma/\mu = {pct:.1f}\%$)')
            else:
                ax.set_title(rf'$n = {ni}$')

        fig.suptitle(
            f"{cfg['pdf']}    {cfg['flavor']}    window: {obs_label}"
            f"    $Q^2 = {cfg['Q2']}$ GeV$^2$",
            fontsize=13,
        )

        fig.subplots_adjust(top=0.92, right=0.87, hspace=0.38, wspace=0.30)
        cax = fig.add_axes([0.90, 0.12, 0.02, 0.74])
        fig.colorbar(
            matplotlib.cm.ScalarMappable(norm=norm, cmap='plasma_r'),
            cax=cax,
            label=r'$\sigma_\mathrm{after}\ /\ \sigma_\mathrm{before}$  (full moment)',
        )

        if save:
            suffix = '_' + os.path.basename(os.path.abspath(cfg['output_dir']))
            out = os.path.join(
                os.path.abspath(cfg['output_dir']),
                f"heatmap_moments{suffix}.pdf",
            )
            fig.savefig(out, dpi=150)
            print(f"Moment heat map → {out}")

        return fig


    # ------------------------------------------------------------------
    def run_charges(self, force=False):
        """Compute tensor-charge ratios from existing profiled PDFs (no ePump re-run).

        For each scan point (x0, w), loads the profiled PDF set produced by run(),
        then evaluates each observable in cfg['charge_observables'] over
        [charge_xmin, charge_xmax] for both base and profiled PDFs.  Stores
        ratio_charges of shape (n_obs, n_midpoints, n_widths).

        Parameters
        ----------
        force : bool
            If True, recompute all cells from scratch.  If False (default),
            resume from partial results saved by an earlier run_charges().
        """
        output_dir   = os.path.abspath(self.cfg['output_dir'])
        results_path = os.path.join(output_dir, 'results.npz')
        if not os.path.exists(results_path):
            raise FileNotFoundError(
                f"No saved state at {results_path}. Run the scan first.")

        saved = np.load(results_path, allow_pickle=True)
        self._pdf_name = str(saved['pdf_name'])
        mc2h_str       = str(saved['mc2h_dir'])
        self._mc2h_dir = mc2h_str if mc2h_str else None

        setup_lhapdf_path(self.cfg.get('lhapdf_path'))
        if self._mc2h_dir:
            setup_lhapdf_path(custom_path=self._mc2h_dir)
            lhapdf.setPaths([self._mc2h_dir] + lhapdf.paths())

        cfg                = self.cfg
        scan_points        = cfg['scan_points']
        Q2                 = float(cfg['Q2'])
        nx                 = int(cfg['nx'])
        charge_xmin        = float(cfg.get('charge_xmin', cfg.get('moment_xmin', 1e-4)))
        charge_xmax        = float(cfg.get('charge_xmax', cfg.get('moment_xmax', 0.999)))
        charge_observables = cfg.get('charge_observables', [])
        pdf_name           = self._pdf_name

        if not charge_observables:
            raise ValueError("cfg['charge_observables'] is empty or missing.")

        midpoints = sorted(set(float(p[0]) for p in scan_points))
        widths    = sorted(set(float(p[1]) for p in scan_points))
        mid_idx   = {v: i for i, v in enumerate(midpoints)}
        wid_idx   = {v: j for j, v in enumerate(widths)}

        n_obs                = len(charge_observables)
        ratio_charges        = np.full((n_obs, len(midpoints), len(widths)), np.nan)
        sigma_before_charges = np.zeros(n_obs)
        central_charges      = np.zeros(n_obs)
        boundary             = saved['boundary'].copy()

        charges_path = os.path.join(output_dir, 'results_charges.npz')
        if not force and os.path.exists(charges_path):
            saved_c = np.load(charges_path, allow_pickle=True)
            if (np.array_equal(saved_c['midpoints'], midpoints) and
                    np.array_equal(saved_c['widths'], widths) and
                    saved_c['ratio_charges'].shape[0] == n_obs):
                ratio_charges[...] = saved_c['ratio_charges']
                n_done = int(np.sum(np.all(np.isfinite(ratio_charges), axis=0)))
                if n_done:
                    print(f"Resuming charges: {n_done}/{len(midpoints)*len(widths)} cells already done.")

        base_set     = lhapdf.getPDFSet(pdf_name)
        base_central = base_set.mkPDF(0)
        base_members = base_set.mkPDFs()

        # Pre-compute base-PDF uncertainty for each charge observable
        parsed_charge_terms = []
        for obs_idx, obs in enumerate(charge_observables):
            parsed = parse_flavor_expression(obs['flavor'])
            parsed_charge_terms.append(parsed)
            kwargs = dict(weight_type='1', moment=obs.get('moment', 0))
            base_vals = [compute_integrated_moment(
                             m, parsed, charge_xmin, charge_xmax, nx, Q2, **kwargs)
                         for m in base_members]
            o = base_set.uncertainty(base_vals)
            sigma_before_charges[obs_idx] = (o.errminus + o.errplus) / 2.0
            central_charges[obs_idx]      = compute_integrated_moment(
                base_central, parsed, charge_xmin, charge_xmax, nx, Q2, **kwargs)

        n_total = len(scan_points)
        for idx, row in enumerate(scan_points, 1):
            x0, w = float(row[0]), float(row[1])
            i, j  = mid_idx[x0], wid_idx[w]

            if np.all(np.isfinite(ratio_charges[:, i, j])):
                print(f"({idx}/{n_total}) [mid_{x0:.4f}_wid_{w:.4f}]  ← already done")
                continue

            label   = f"mid_{x0:.4f}_wid_{w:.4f}"
            run_dir = os.path.join(output_dir, label)

            if not os.path.isdir(os.path.join(run_dir, label)):
                print(f"({idx}/{n_total}) [{label}]  WARNING: profiled set not found — skipping.")
                continue

            if run_dir not in lhapdf.paths():
                lhapdf.setPaths([run_dir] + lhapdf.paths())
            profiled_set = lhapdf.getPDFSet(label)

            ratio_strs = []
            lhapdf.setVerbosity(0)
            for obs_idx, obs in enumerate(charge_observables):
                parsed = parsed_charge_terms[obs_idx]
                kwargs = dict(weight_type='1', moment=obs.get('moment', 0))
                prof_vals = [compute_integrated_moment(
                                 profiled_set.mkPDF(k), parsed,
                                 charge_xmin, charge_xmax, nx, Q2, **kwargs)
                             for k in range(profiled_set.size)]
                p = profiled_set.uncertainty(prof_vals)
                sigma_a = (p.errminus + p.errplus) / 2.0
                r = sigma_a / sigma_before_charges[obs_idx] if sigma_before_charges[obs_idx] > 0 else np.nan
                ratio_charges[obs_idx, i, j] = r
                ratio_strs.append(f"obs{obs_idx}:{r:.3f}")
            lhapdf.setVerbosity(1)

            print(f"({idx}/{n_total}) [{label}]  → " + "  ".join(ratio_strs))
            np.savez(charges_path,
                     ratio_charges=ratio_charges,
                     sigma_before_charges=sigma_before_charges,
                     central_charges=central_charges,
                     charge_labels=np.array([obs.get('label', obs['flavor'])
                                             for obs in charge_observables]),
                     midpoints=np.array(midpoints),
                     widths=np.array(widths),
                     boundary=boundary,
                     pdf_name=np.array(self._pdf_name),
                     mc2h_dir=np.array(self._mc2h_dir or ''))

        self.midpoints            = midpoints
        self.widths               = widths
        self.boundary             = boundary
        self.ratio_charges        = ratio_charges
        self.sigma_before_charges = sigma_before_charges
        self.central_charges      = central_charges
        self.charge_labels        = [obs.get('label', obs['flavor']) for obs in charge_observables]
        print(f"Charge results saved → {charges_path}")
        return self

    # ------------------------------------------------------------------
    def load_charges(self):
        """Load tensor-charge results from a previous run_charges() call."""
        charges_path = os.path.join(os.path.abspath(self.cfg['output_dir']), 'results_charges.npz')
        if not os.path.exists(charges_path):
            raise FileNotFoundError(
                f"No saved charge results at {charges_path}. Run run_charges() first.")
        data = np.load(charges_path, allow_pickle=True)
        self.ratio_charges        = data['ratio_charges']
        self.sigma_before_charges = data['sigma_before_charges']
        self.central_charges      = data['central_charges']
        self.charge_labels        = data['charge_labels'].tolist()
        self.midpoints            = data['midpoints'].tolist()
        self.widths               = data['widths'].tolist()
        self.boundary             = data['boundary']
        print(f"Charge results loaded ← {charges_path}")
        return self

    # ------------------------------------------------------------------
    def plot_charges(self, save=True):
        """
        Generate a 2×2 grid of window-to-charge heat maps.

        Each panel shows σ_after / σ_before for one tensor-charge observable
        defined in cfg['charge_observables'], with the before-profiling relative
        uncertainty displayed in the subplot title.

        Parameters
        ----------
        save : bool
            Write figure to output_dir/heatmap_charges_{basename}.pdf (default True).

        Returns
        -------
        matplotlib.figure.Figure
        """
        if self.ratio_charges is None:
            raise RuntimeError("Call run_charges() or load_charges() before plot_charges().")

        import matplotlib.pyplot as plt

        cfg = self.cfg
        M_e = _bin_edges(np.array(self.midpoints))
        W_e = _bin_edges(np.array(self.widths))

        moment = int(cfg.get('moment', 1))
        weight = cfg.get('weight', 'gaussian')
        obs_label = rf'$g_{{{moment}}}$' if weight == 'gaussian' else rf'$a_{{{moment}}}$'

        n_obs     = self.ratio_charges.shape[0]
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        axes_flat = axes.flatten()

        vmin, vmax = 0.0, 1.0
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

        for obs_idx in range(min(n_obs, 4)):
            ax = axes_flat[obs_idx]
            masked = np.ma.masked_invalid(self.ratio_charges[obs_idx])
            ax.pcolormesh(W_e, M_e, masked, cmap='plasma_r', norm=norm)

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

            lbl = self.charge_labels[obs_idx] if self.charge_labels else f"obs {obs_idx}"
            if (self.sigma_before_charges is not None and
                    self.central_charges is not None and
                    abs(self.central_charges[obs_idx]) > 0):
                pct = 100 * self.sigma_before_charges[obs_idx] / abs(self.central_charges[obs_idx])
                ax.set_title(rf"{lbl}   ($\sigma/\mu = {pct:.1f}\%$)")
            else:
                ax.set_title(lbl)

        for obs_idx in range(n_obs, 4):
            axes_flat[obs_idx].set_visible(False)

        fig.suptitle(
            f"{cfg['pdf']}    window: {obs_label}"
            f"    $Q^2 = {cfg['Q2']}$ GeV$^2$",
            fontsize=13,
        )

        fig.subplots_adjust(top=0.92, right=0.87, hspace=0.38, wspace=0.30)
        cax = fig.add_axes([0.90, 0.12, 0.02, 0.74])
        fig.colorbar(
            matplotlib.cm.ScalarMappable(norm=norm, cmap='plasma_r'),
            cax=cax,
            label=r'$\sigma_\mathrm{after}\ /\ \sigma_\mathrm{before}$  (tensor charge)',
        )

        if save:
            suffix = '_' + os.path.basename(os.path.abspath(cfg['output_dir']))
            out = os.path.join(
                os.path.abspath(cfg['output_dir']),
                f"heatmap_charges{suffix}.pdf",
            )
            fig.savefig(out, dpi=150)
            print(f"Charge heat map → {out}")

        return fig


# ── Terminal entry point ───────────────────────────────────────────────────────

def main():
    matplotlib.use('Agg')   # non-interactive backend for terminal use

    if len(sys.argv) < 2:
        print(
            f"Usage: {sys.argv[0]} runcard.py [run|recompute|moments|charges]",
            file=sys.stderr,
        )
        sys.exit(1)

    scanner = WindowMomentScanner(sys.argv[1])
    cmd = sys.argv[2] if len(sys.argv) >= 3 else 'run'

    if cmd == 'recompute':
        scanner.recompute()
        scanner.plot(save=True)
    elif cmd == 'moments':
        scanner.run_moments()
        scanner.plot_moments(save=True)
    elif cmd == 'charges':
        scanner.run_charges()
        scanner.plot_charges(save=True)
    else:
        scanner.run()
        scanner.plot(save=True)


if __name__ == '__main__':
    main()
