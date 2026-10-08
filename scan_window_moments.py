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
    p = str(path)
    if not os.path.exists(p) and not p.endswith('.py'):
        p = p + '.py'
    spec = importlib.util.spec_from_file_location("runcard", os.path.abspath(p))
    if spec is None:
        raise FileNotFoundError(f"Runcard not found: {path!r}")
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


def _check_theory_file_format(theory_path, n_obs):
    """Assert that a theory file has Ncol=1 with one row per observable per member block.

    Raises AssertionError if the format is wrong (catches the multi-obs bug where
    all values were written on one line per member instead of one line each).
    """
    with open(theory_path) as f:
        lines = [l.strip() for l in f if l.strip()]
    for i, line in enumerate(lines):
        if line.startswith("PDF_0_Set"):
            for k in range(n_obs):
                val_line = lines[i + 1 + k]
                n_tokens = len(val_line.split())
                assert n_tokens == 1, (
                    f"Theory file format error at obs {k}: expected 1 value per line, "
                    f"got {n_tokens}. Line: {val_line!r}\n"
                    f"  File: {theory_path}")
            print(f"Theory file OK: {n_obs} obs/member, one row each  ({theory_path})")
            return
    print(f"Warning: PDF_0_Set marker not found in {theory_path} — format check skipped.")


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
        self._pdf_name       = str(saved['pdf_name'])
        mc2h_str             = str(saved['mc2h_dir'])
        self._mc2h_dir       = mc2h_str if mc2h_str else None
        self.sigma_before_ww = saved['sigma_before_ww'] if 'sigma_before_ww' in saved else None
        self.central_ww      = saved['central_ww']      if 'central_ww'      in saved else None

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
                     sigma_before_ww=self.sigma_before_ww if self.sigma_before_ww is not None else np.array([]),
                     central_ww=self.central_ww           if self.central_ww      is not None else np.array([]),
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
        ww = data['sigma_before_ww'] if 'sigma_before_ww' in data else None
        self.sigma_before_ww = ww if (ww is not None and ww.size > 0) else None
        cw = data['central_ww'] if 'central_ww' in data else None
        self.central_ww      = cw if (cw is not None and cw.size > 0) else None
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
                if np.isfinite(self.ratio[i, j]):
                    cx = (W_e[j] + W_e[j + 1]) / 2
                    cy = (M_e[i] + M_e[i + 1]) / 2
                    ax.text(cx, cy, f"{self.ratio[i, j]:.2f}",
                            ha='center', va='center', fontsize=7, color='white')

        ax.set_xlabel('Window width  $w$')
        ax.set_ylabel('Window midpoint  $x_0$')
        pdf_label = cfg.get('pdf_label', cfg['pdf'])
        fig.suptitle(
            f"{pdf_label} — Window moment profiling\n"
            f"{cfg['flavor']}    {obs_label}    $Q^2 = {cfg['Q2']}$ GeV$^2$",
            fontsize=12,
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
                    if np.isfinite(self.ratio_moments[ni, i, j]):
                        cx = (W_e[j] + W_e[j + 1]) / 2
                        cy = (M_e[i] + M_e[i + 1]) / 2
                        ax.text(cx, cy, f"{self.ratio_moments[ni, i, j]:.2f}",
                                ha='center', va='center', fontsize=7, color='white')

            ax.set_xlabel('Window width  $w$')
            ax.set_ylabel('Window midpoint  $x_0$')
            if (self.sigma_before_moments is not None and
                    self.central_moments is not None and
                    abs(self.central_moments[ni]) > 0):
                pct = 100 * self.sigma_before_moments[ni] / abs(self.central_moments[ni])
                ax.set_title(rf'$n = {ni}$   ($\sigma/\mu = {pct:.1f}\%$)')
            else:
                ax.set_title(rf'$n = {ni}$')

        pdf_label = cfg.get('pdf_label', cfg['pdf'])
        fig.suptitle(
            f"{pdf_label} — Full-moment profiling\n"
            f"{cfg['flavor']}    window: {obs_label}    $Q^2 = {cfg['Q2']}$ GeV$^2$",
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
        self._pdf_name       = str(saved['pdf_name'])
        mc2h_str             = str(saved['mc2h_dir'])
        self._mc2h_dir       = mc2h_str if mc2h_str else None
        self.sigma_before_ww = saved['sigma_before_ww'] if 'sigma_before_ww' in saved else None
        self.central_ww      = saved['central_ww']      if 'central_ww'      in saved else None

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
                     sigma_before_ww=self.sigma_before_ww if self.sigma_before_ww is not None else np.array([]),
                     central_ww=self.central_ww           if self.central_ww      is not None else np.array([]),
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
        ww = data['sigma_before_ww'] if 'sigma_before_ww' in data else None
        self.sigma_before_ww = ww if (ww is not None and ww.size > 0) else None
        cw = data['central_ww'] if 'central_ww' in data else None
        self.central_ww      = cw if (cw is not None and cw.size > 0) else None
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
                    if np.isfinite(self.ratio_charges[obs_idx, i, j]):
                        cx = (W_e[j] + W_e[j + 1]) / 2
                        cy = (M_e[i] + M_e[i + 1]) / 2
                        ax.text(cx, cy, f"{self.ratio_charges[obs_idx, i, j]:.2f}",
                                ha='center', va='center', fontsize=7, color='white')

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

        pdf_label = cfg.get('pdf_label', cfg['pdf'])
        fig.suptitle(
            f"{pdf_label} — Tensor charge profiling\n"
            f"window: {obs_label}    $Q^2 = {cfg['Q2']}$ GeV$^2$",
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


class AnchoredWindowScanner:
    """
    1D scan over window midpoints with a fixed anchor measurement.

    Each cell profiles from the base PDF using two simultaneous ePump measurements:
    (1) the fixed anchor window and (2) the current scan midpoint (same width).
    Midpoints that overlap the anchor window are automatically skipped.

    Ratios σ_after/σ_before use the base PDF as the denominator, directly
    comparable to WindowMomentScanner results for the same cells.

    Parameters
    ----------
    cfg : dict or str
        Configuration dict, or path to a runcard.py file that defines ``cfg``.
        Required keys beyond the shared PDF/flavor/Q2/nx/moment/weight keys:

        ``anchor``
            Dict with ``x0`` and ``w`` — the fixed first measurement.
        ``scan_midpoints``
            List of x0 values to scan.  Width is always ``anchor["w"]``.
        ``rel_unc``
            Fractional pseudo-data uncertainty applied to both measurements.
    """

    def __init__(self, cfg):
        if isinstance(cfg, (str, os.PathLike)):
            cfg = load_runcard(str(cfg))
        self.cfg = cfg

        self.midpoints    = None   # 1D list of scanned x0 values
        self.ratio        = None   # 1D ndarray, NaN for overlapping/skipped cells
        self.boundary     = None   # 1D bool array: True = overlaps anchor window
        self.sigma_before = None   # 1D ndarray: base PDF σ per midpoint
        self.central_vals = None   # 1D ndarray: central moment value per midpoint

        self.ratio_charges        = None   # (n_obs, n_midpoints)
        self.sigma_before_charges = None   # (n_obs,)
        self.charge_labels        = None   # list[str]
        self.ratio_poly_moments   = None   # (4, n_midpoints)
        self.sigma_before_poly    = None   # (4,)

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
        """Run the anchored 1D scan.  Calls setup() automatically if not already done.

        Parameters
        ----------
        force : bool
            Rerun ePump for every cell from scratch, ignoring any saved results.
        """
        if self._pdf_name is None:
            self.setup()

        cfg        = self.cfg
        anchor_cfg = cfg['anchor']
        ax0        = float(anchor_cfg['x0'])
        w          = float(anchor_cfg['w'])
        flavor     = cfg['flavor']
        Q2         = float(cfg['Q2'])
        nx         = int(cfg['nx'])
        moment     = int(cfg.get('moment', 1))
        weight     = cfg.get('weight', 'gaussian')
        rel_unc    = float(cfg['rel_unc'])
        output_dir = os.path.abspath(cfg['output_dir'])
        epump_path = os.path.abspath(cfg.get('epump_path', './ePump_kp20221218/src/UpdatePDFs'))
        pdf_name   = self._pdf_name
        mc2h_dir   = self._mc2h_dir

        charge_obs  = cfg.get('charge_observables', [])
        n_obs       = len(charge_obs)
        charge_xmin = float(cfg.get('charge_xmin', 1e-4))
        charge_xmax = float(cfg.get('charge_xmax', 0.999))
        n_orders    = 4
        moment_xmin = float(cfg.get('moment_xmin', 1e-4))
        moment_xmax = float(cfg.get('moment_xmax', 0.999))

        midpoints = sorted(float(x) for x in cfg['scan_midpoints'])
        n         = len(midpoints)

        ratio             = np.full(n, np.nan)
        boundary          = np.zeros(n, dtype=bool)
        sigma_before      = np.full(n, np.nan)
        central_vals      = np.full(n, np.nan)
        ratio_charges_buf = np.full((max(n_obs, 1), n), np.nan)
        sigma_before_ch   = np.full(max(n_obs, 1), np.nan)
        ratio_poly_buf    = np.full((n_orders, n), np.nan)
        sigma_before_poly = np.full(n_orders, np.nan)

        results_path = os.path.join(output_dir, 'results_anchored.npz')
        if not force and os.path.exists(results_path):
            saved = np.load(results_path, allow_pickle=True)
            if np.array_equal(saved['midpoints'], midpoints):
                ratio[...]        = saved['ratio']
                boundary[...]     = saved['boundary']
                if 'sigma_before' in saved:
                    sigma_before[...] = saved['sigma_before']
                if 'central_vals' in saved:
                    central_vals[...] = saved['central_vals']
                if 'ratio_charges' in saved and saved['ratio_charges'].shape == ratio_charges_buf.shape:
                    ratio_charges_buf[...] = saved['ratio_charges']
                if 'ratio_poly_moments' in saved:
                    ratio_poly_buf[...] = saved['ratio_poly_moments']
                n_done = int(np.sum(np.isfinite(ratio)))
                if n_done:
                    print(f"Resuming: {n_done}/{n} cells already done.")

        parsed_terms = parse_flavor_expression(flavor)
        poly_terms   = parse_flavor_expression(flavor)
        base_set     = lhapdf.getPDFSet(pdf_name)
        base_central = base_set.mkPDF(0)
        base_members = base_set.mkPDFs()

        # Pre-compute anchor measurement (constant across all cells)
        axmin       = max(1e-4, ax0 - w / 2)
        axmax       = min(0.999, ax0 + w / 2)
        anchor_val  = compute_integrated_moment(
            base_central, parsed_terms, axmin, axmax, nx, Q2,
            weight_type=weight, moment=moment,
        )
        anchor_stat = rel_unc * abs(anchor_val)

        orig_anchor_vals = [
            compute_integrated_moment(m, parsed_terms, axmin, axmax, nx, Q2,
                                      weight_type=weight, moment=moment)
            for m in base_members
        ]
        unc_anc = base_set.uncertainty(orig_anchor_vals)
        sigma_before_anchor = (unc_anc.errminus + unc_anc.errplus) / 2.0
        print(f"Anchor: x0={ax0}, w={w}  |  val={anchor_val:.5g}  σ_before={sigma_before_anchor:.5g}")

        # σ_before for charge observables (computed once from base PDF)
        for obs_idx, obs in enumerate(charge_obs):
            obs_terms = parse_flavor_expression(obs['flavor'])
            kw = dict(weight_type=obs['weight'], moment=int(obs['moment']))
            vals = [compute_integrated_moment(m, obs_terms, charge_xmin, charge_xmax,
                                              nx, Q2, **kw) for m in base_members]
            unc = base_set.uncertainty(vals)
            sigma_before_ch[obs_idx] = (unc.errminus + unc.errplus) / 2.0

        # σ_before for polynomial moments (computed once from base PDF)
        for ni in range(n_orders):
            kw = dict(weight_type='1', moment=ni)
            vals = [compute_integrated_moment(m, poly_terms, moment_xmin, moment_xmax,
                                              nx, Q2, **kw) for m in base_members]
            unc = base_set.uncertainty(vals)
            sigma_before_poly[ni] = (unc.errminus + unc.errplus) / 2.0

        theory_checked = False

        for i, x0 in enumerate(midpoints):
            if abs(x0 - ax0) < w - 1e-9:   # strict interior overlap; touching is fine
                boundary[i] = True
                print(f"({i+1}/{n}) [x0={x0:.4f}]  ← overlaps anchor — skipped")
                continue

            # Skip only if ratio AND charges AND moments are all computed
            cell_done = (np.isfinite(ratio[i]) and
                         (n_obs == 0 or np.all(np.isfinite(ratio_charges_buf[:n_obs, i]))) and
                         np.all(np.isfinite(ratio_poly_buf[:, i])))
            if cell_done and not force:
                print(f"({i+1}/{n}) [x0={x0:.4f}]  ← already done")
                continue

            xmin    = max(1e-4, x0 - w / 2)
            xmax    = min(0.999, x0 + w / 2)
            clipped = (x0 - w / 2 < 1e-4) or (x0 + w / 2 > 0.999)

            central_val = compute_integrated_moment(
                base_central, parsed_terms, xmin, xmax, nx, Q2,
                weight_type=weight, moment=moment,
            )
            stat_err = rel_unc * abs(central_val)

            label    = f"anc_{ax0:.4f}_scan_{x0:.4f}"
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
                    x=ax0, Q2=Q2, value=anchor_val, stat=anchor_stat,
                    obs_type='moment', flavor=flavor,
                    xmin=axmin, xmax=axmax, nx=nx,
                    weight=weight, moment=moment,
                )
                ep.add_measurement(
                    x=x0, Q2=Q2, value=central_val, stat=stat_err,
                    obs_type='moment', flavor=flavor,
                    xmin=xmin, xmax=xmax, nx=nx,
                    weight=weight, moment=moment,
                )
                ep.generate_files()

                if not theory_checked:
                    _check_theory_file_format(f"{run_name}.theory", n_obs=2)
                    theory_checked = True

                ep.run()
                profiled_set = ep.profiled_set

            # Evaluate σ on the SCAN window moment (consistent denominator with scan 1)
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

            ratio[i]        = r
            sigma_before[i] = sigma_b
            central_vals[i] = central_val

            # Evaluate tensor charges
            for obs_idx, obs in enumerate(charge_obs):
                obs_terms = parse_flavor_expression(obs['flavor'])
                kw = dict(weight_type=obs['weight'], moment=int(obs['moment']))
                lhapdf.setVerbosity(0)
                pv = [compute_integrated_moment(profiled_set.mkPDF(k), obs_terms,
                                                charge_xmin, charge_xmax, nx, Q2, **kw)
                      for k in range(profiled_set.size)]
                lhapdf.setVerbosity(1)
                unc = profiled_set.uncertainty(pv)
                ratio_charges_buf[obs_idx, i] = ((unc.errminus + unc.errplus) / 2.0 /
                                                  sigma_before_ch[obs_idx])

            # Evaluate polynomial moments (n = 0 … n_orders-1)
            for ni in range(n_orders):
                kw = dict(weight_type='1', moment=ni)
                lhapdf.setVerbosity(0)
                pv = [compute_integrated_moment(profiled_set.mkPDF(k), poly_terms,
                                                moment_xmin, moment_xmax, nx, Q2, **kw)
                      for k in range(profiled_set.size)]
                lhapdf.setVerbosity(1)
                unc = profiled_set.uncertainty(pv)
                ratio_poly_buf[ni, i] = ((unc.errminus + unc.errplus) / 2.0 /
                                          sigma_before_poly[ni])

            clip_tag = "  [BOUNDARY-CLIPPED]" if clipped else ""
            print(f"({i+1}/{n}) [x0={x0:.4f}]  → ratio={r:.4f}  "
                  f"σ_before={sigma_b:.5g}  σ_after={sigma_a:.5g}{clip_tag}")
            np.savez(results_path,
                     midpoints=midpoints, ratio=ratio, boundary=boundary,
                     sigma_before=sigma_before, central_vals=central_vals,
                     ratio_charges=ratio_charges_buf,
                     sigma_before_charges=sigma_before_ch,
                     charge_labels=np.array([o['label'] for o in charge_obs]),
                     ratio_poly_moments=ratio_poly_buf,
                     sigma_before_poly=sigma_before_poly,
                     anchor_x0=np.array(ax0), anchor_w=np.array(w),
                     pdf_name=np.array(self._pdf_name),
                     mc2h_dir=np.array(self._mc2h_dir or ''))

        self.midpoints            = midpoints
        self.ratio                = ratio
        self.boundary             = boundary
        self.sigma_before         = sigma_before
        self.central_vals         = central_vals
        self.ratio_charges        = ratio_charges_buf[:n_obs] if n_obs else np.empty((0, n))
        self.sigma_before_charges = sigma_before_ch[:n_obs]
        self.charge_labels        = [o['label'] for o in charge_obs]
        self.ratio_poly_moments   = ratio_poly_buf
        self.sigma_before_poly    = sigma_before_poly
        print(f"Results saved → {results_path}")
        return self

    # ------------------------------------------------------------------
    def load(self):
        """Load results from a previous run() call."""
        results_path = os.path.join(os.path.abspath(self.cfg['output_dir']),
                                    'results_anchored.npz')
        if not os.path.exists(results_path):
            raise FileNotFoundError(
                f"No saved results at {results_path}. Call run() first.")
        data = np.load(results_path, allow_pickle=True)
        self.midpoints            = data['midpoints'].tolist()
        self.ratio                = data['ratio']
        self.boundary             = data['boundary']
        self.sigma_before         = data['sigma_before']         if 'sigma_before'         in data else None
        self.central_vals         = data['central_vals']         if 'central_vals'         in data else None
        self.ratio_charges        = data['ratio_charges']        if 'ratio_charges'        in data else None
        self.sigma_before_charges = data['sigma_before_charges'] if 'sigma_before_charges' in data else None
        self.charge_labels        = (data['charge_labels'].tolist()
                                     if 'charge_labels' in data else [])
        self.ratio_poly_moments   = data['ratio_poly_moments']   if 'ratio_poly_moments'   in data else None
        self.sigma_before_poly    = data['sigma_before_poly']    if 'sigma_before_poly'    in data else None
        print(f"Results loaded ← {results_path}")
        return self

    # ------------------------------------------------------------------
    def plot(self, save=True):
        """Multi-panel plot of σ_after/σ_before vs scan midpoint.

        One panel per tensor-charge observable and polynomial moment order.
        Falls back to a single window-moment panel if charges/moments were not
        computed (old results file).

        Parameters
        ----------
        save : bool
            Save figure to output_dir/anchored_scan_<suffix>.pdf.

        Returns
        -------
        matplotlib.figure.Figure
        """
        if self.ratio is None:
            raise RuntimeError("Call run() or load() before plot().")

        import matplotlib.pyplot as plt

        cfg       = self.cfg
        ax0       = float(cfg['anchor']['x0'])
        w         = float(cfg['anchor']['w'])
        moment    = int(cfg.get('moment', 1))
        weight    = cfg.get('weight', 'gaussian')
        obs_label = rf'$g_{{{moment}}}$' if weight == 'gaussian' else rf'$a_{{{moment}}}$'
        pdf_label = cfg.get('pdf_label', cfg['pdf'])

        x        = np.array(self.midpoints)
        mask     = ~np.asarray(self.boundary) & np.isfinite(self.ratio)
        xlo      = max(0.0, min(x) - w / 2)
        xhi      = min(1.0, max(x) + w / 2)

        has_charges = (self.ratio_charges is not None and
                       self.ratio_charges.shape[0] > 0)
        has_poly    = self.ratio_poly_moments is not None

        def _draw_panel(ax, ydata, title):
            ax.plot(x[mask], ydata[mask], 'o-', color='C0')
            ax.axhline(1.0, color='gray', ls='--', lw=0.8)
            ax.axvline(ax0, color='C1', ls=':', lw=1.2)
            ax.axvspan(max(0, ax0 - w), min(1, ax0 + w), alpha=0.08, color='C1')
            ax.set_xlim(xlo, xhi)
            ax.set_ylim(0, 1.1)
            ax.set_xlabel('$x_0$')
            ax.set_ylabel(r'$\sigma_{\rm after}/\sigma_{\rm before}$')
            ax.set_title(title, fontsize=9)

        if not has_charges and not has_poly:
            fig, ax = plt.subplots(figsize=(8, 4))
            _draw_panel(ax, self.ratio, f'Window {obs_label}')
        else:
            n_charge_panels = self.ratio_charges.shape[0] if has_charges else 0
            n_poly_panels   = self.ratio_poly_moments.shape[0] if has_poly else 0
            n_panels        = n_charge_panels + n_poly_panels
            ncols = min(n_panels, 4)
            nrows = (n_panels + ncols - 1) // ncols
            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(4 * ncols, 3.5 * nrows),
                                     squeeze=False, sharey=True)
            axes_flat = axes.flatten()
            panel = 0

            if has_charges:
                labels = self.charge_labels or [f'obs {k}' for k in range(n_charge_panels)]
                for k in range(n_charge_panels):
                    _draw_panel(axes_flat[panel], self.ratio_charges[k], labels[k])
                    panel += 1

            if has_poly:
                for ni in range(n_poly_panels):
                    _draw_panel(axes_flat[panel], self.ratio_poly_moments[ni],
                                rf'Full $n={ni}$')
                    panel += 1

            for k in range(panel, len(axes_flat)):
                axes_flat[k].set_visible(False)

            for r in range(nrows):
                for c in range(1, ncols):
                    axes[r, c].tick_params(labelleft=False)
                    axes[r, c].set_ylabel('')

        fig.suptitle(
            f"{pdf_label} — Anchored window scan    {obs_label}    "
            f"$Q^2={cfg['Q2']}$ GeV$^2$    {cfg['flavor']}\n"
            f"Anchor: $x_0={ax0}$, $w={w}$",
            fontsize=11,
        )
        plt.tight_layout()
        plt.subplots_adjust(wspace=0)

        if save:
            suffix = '_' + os.path.basename(os.path.abspath(cfg['output_dir']))
            out = os.path.join(os.path.abspath(cfg['output_dir']),
                               f"anchored_scan{suffix}.pdf")
            fig.savefig(out, dpi=150)
            print(f"Plot → {out}")

        return fig


class MomentAccumulationScanner:
    """
    Profile a fixed Gaussian window with increasing simultaneous moment constraints.

    For n = 1 … max_n, runs ePump with window moments 1..n as simultaneous
    constraints, then records σ_after / σ_before for:
      - each window moment 1..max_n
      - every tensor-charge observable in cfg['charge_observables']
      - full polynomial moments ∫ x^n f(x) dx for n = 0, 1, 2, 3

    Parameters
    ----------
    cfg : dict or str
        Configuration dict (must contain a ``moment_scan`` sub-dict), or path
        to a runcard file.
    """

    def __init__(self, cfg):
        if isinstance(cfg, (str, os.PathLike)):
            cfg = load_runcard(str(cfg))
        self.cfg = cfg

        # Results — populated by run()
        self.n_values             = None  # [1, 2, ..., max_n]
        self.ratio_window         = None  # (max_n, max_n): window moment k vs n constraints
        self.ratio_charges        = None  # (n_obs, max_n)
        self.ratio_poly_moments   = None  # (4, max_n)
        self.sigma_before_window  = None  # (max_n,)
        self.sigma_before_charges = None  # (n_obs,)
        self.sigma_before_poly    = None  # (4,)
        self.central_window       = None  # (max_n,)
        self.central_charges      = None  # (n_obs,)
        self.central_poly         = None  # (4,)
        self.charge_labels        = None

        # Internal — populated by setup()
        self._pdf_name = None
        self._mc2h_dir = None

    # ------------------------------------------------------------------
    def setup(self):
        """Configure LHAPDF paths and convert MC replicas to Hessian (once)."""
        cfg        = self.cfg
        output_dir = os.path.abspath(cfg['output_dir'])
        ms_dir     = os.path.join(output_dir, 'moment_accumulation')
        os.makedirs(ms_dir, exist_ok=True)
        setup_lhapdf_path(cfg.get('lhapdf_path'))

        pdf_name = cfg['pdf']
        err_type = detect_pdf_error_type(pdf_name)
        if err_type not in ('replicas', 'mc'):
            print(f"PDF '{pdf_name}' is {err_type!r} — no conversion needed.")
            self._pdf_name = pdf_name
            self._mc2h_dir = None
        else:
            mc2h_dir = os.path.join(output_dir, '_hessian')
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
        """Run moment accumulation: profile moments 1..n for n = 1..max_n.

        For each level n, adds window moments 1..n as simultaneous ePump
        constraints, then evaluates uncertainties on all tracked quantities.

        Parameters
        ----------
        force : bool
            If True, ignore previously saved results and rerun all levels.
        """
        if self._pdf_name is None:
            self.setup()

        cfg    = self.cfg
        ms_cfg = cfg['moment_scan']

        x0      = float(ms_cfg['x0'])
        w       = float(ms_cfg['w'])
        weight  = ms_cfg.get('weight', 'gaussian')
        max_n   = int(ms_cfg['max_n'])
        rel_unc = float(ms_cfg.get('rel_unc', 0.10))
        corr    = ms_cfg.get('corr', None)

        start_n = ms_cfg.get('start_n', None)
        if start_n is None:
            start_n = 0 if weight == '1' else 1
        start_n = int(start_n)
        if weight == 'gaussian' and start_n == 0:
            raise ValueError("start_n=0 is invalid for gaussian weight (integral is zero).")
        n_levels = max_n - start_n + 1

        flavor             = cfg['flavor']
        Q2                 = float(cfg['Q2'])
        nx                 = int(cfg['nx'])
        charge_xmin        = float(cfg.get('charge_xmin', cfg.get('moment_xmin', 1e-4)))
        charge_xmax        = float(cfg.get('charge_xmax', cfg.get('moment_xmax', 0.999)))
        moment_xmin        = float(cfg.get('moment_xmin', 1e-4))
        moment_xmax        = float(cfg.get('moment_xmax', 0.999))
        charge_observables = cfg.get('charge_observables', [])

        output_dir = os.path.abspath(cfg['output_dir'])
        ms_dir     = os.path.join(output_dir, 'moment_accumulation')
        os.makedirs(ms_dir, exist_ok=True)
        epump_path = os.path.abspath(cfg.get('epump_path', './ePump_kp20221218/src/UpdatePDFs'))
        pdf_name   = self._pdf_name
        mc2h_dir   = self._mc2h_dir

        xmin = max(1e-4, x0 - w / 2)
        xmax = min(0.999, x0 + w / 2)

        parsed_flavor = parse_flavor_expression(flavor)

        base_set     = lhapdf.getPDFSet(pdf_name)
        base_central = base_set.mkPDF(0)
        base_members = base_set.mkPDFs()

        n_obs    = len(charge_observables)
        n_orders = 4

        # Central values and base uncertainties for window moments start_n..max_n
        window_central      = np.zeros(n_levels)
        sigma_before_window = np.zeros(n_levels)
        for k in range(n_levels):
            kw = dict(weight_type=weight, moment=k + start_n)
            window_central[k] = compute_integrated_moment(
                base_central, parsed_flavor, xmin, xmax, nx, Q2, **kw)
            base_vals = [compute_integrated_moment(
                             m, parsed_flavor, xmin, xmax, nx, Q2, **kw)
                         for m in base_members]
            o = base_set.uncertainty(base_vals)
            sigma_before_window[k] = (o.errminus + o.errplus) / 2.0

        # Base uncertainties for charge observables
        parsed_charge_terms  = []
        sigma_before_charges = np.zeros(n_obs)
        central_charges      = np.zeros(n_obs)
        for obs_idx, obs in enumerate(charge_observables):
            parsed = parse_flavor_expression(obs['flavor'])
            parsed_charge_terms.append(parsed)
            kw = dict(weight_type='1', moment=obs.get('moment', 0))
            base_vals = [compute_integrated_moment(
                             m, parsed, charge_xmin, charge_xmax, nx, Q2, **kw)
                         for m in base_members]
            o = base_set.uncertainty(base_vals)
            sigma_before_charges[obs_idx] = (o.errminus + o.errplus) / 2.0
            central_charges[obs_idx] = compute_integrated_moment(
                base_central, parsed, charge_xmin, charge_xmax, nx, Q2, **kw)

        # Base uncertainties for polynomial moments n = 0..3
        sigma_before_poly = np.zeros(n_orders)
        central_poly      = np.zeros(n_orders)
        for ni in range(n_orders):
            kw = dict(weight_type='1', moment=ni)
            base_vals = [compute_integrated_moment(
                             m, parsed_flavor, moment_xmin, moment_xmax, nx, Q2, **kw)
                         for m in base_members]
            o = base_set.uncertainty(base_vals)
            sigma_before_poly[ni] = (o.errminus + o.errplus) / 2.0
            central_poly[ni] = compute_integrated_moment(
                base_central, parsed_flavor, moment_xmin, moment_xmax, nx, Q2, **kw)

        # Result arrays — last axis is constraint level index (0 = only moment start_n)
        ratio_window       = np.full((n_levels, n_levels), np.nan)
        ratio_charges_buf  = np.full((max(n_obs, 1), n_levels), np.nan)
        ratio_poly_moments = np.full((n_orders, n_levels), np.nan)

        # Attempt to resume from a previous run
        results_path = os.path.join(ms_dir, 'results_moment_scan.npz')
        if not force and os.path.exists(results_path):
            saved = np.load(results_path, allow_pickle=True)
            if (int(saved['max_n']) == max_n and
                    int(saved.get('start_n', 1)) == start_n and
                    saved['ratio_charges'].shape[0] == max(n_obs, 1)):
                ratio_window[...]       = saved['ratio_window']
                ratio_charges_buf[...]  = saved['ratio_charges']
                ratio_poly_moments[...] = saved['ratio_poly_moments']
                n_done = int(np.sum(np.isfinite(ratio_window[0])))
                if n_done:
                    print(f"Resuming: {n_done}/{n_levels} levels already done.")

        # Main loop: n = start_n .. max_n
        for n in range(start_n, max_n + 1):
            lvl = n - start_n  # column index

            if (np.all(np.isfinite(ratio_window[:, lvl])) and
                    np.all(np.isfinite(ratio_poly_moments[:, lvl]))):
                print(f"n={n}: ← already done")
                continue

            vals   = window_central[: lvl + 1]
            sigmas = rel_unc * np.abs(vals)

            # Build Cholesky structure for correlated measurements
            L = None
            if corr is not None:
                rho = np.array(corr, dtype=float)[: lvl + 1, : lvl + 1]
                C   = np.outer(sigmas, sigmas) * rho
                jitter = 1e-12 * max(float(np.abs(C).max()), 1.0)
                C  += np.eye(lvl + 1) * jitter
                L   = np.linalg.cholesky(C)

            n_tag    = '_'.join(str(k + start_n) for k in range(lvl + 1))
            label    = f"macc_n{n_tag}"
            run_dir  = os.path.join(ms_dir, label)
            run_name = os.path.join(run_dir, label)

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

                for k in range(lvl + 1):
                    val_k = float(vals[k])
                    if L is not None:
                        stat_k    = 0.0
                        cor_sys_k = [L[k, j] / val_k * 100.0 for j in range(lvl + 1)]
                        uncor_k   = 0.0
                    else:
                        stat_k    = float(sigmas[k])
                        cor_sys_k = []
                        uncor_k   = 0.0

                    ep.add_measurement(
                        x=x0, Q2=Q2, value=val_k,
                        stat=stat_k, uncor_sys=uncor_k, cor_sys=cor_sys_k,
                        obs_type='moment', flavor=flavor,
                        xmin=xmin, xmax=xmax, nx=nx,
                        weight=weight, moment=k + start_n,
                    )

                ep.generate_files()
                ep.run()
                profiled_set = ep.profiled_set

            # Evaluate all tracked quantities on the profiled set
            lhapdf.setVerbosity(0)

            for k in range(n_levels):
                kw = dict(weight_type=weight, moment=k + start_n)
                prof_vals = [compute_integrated_moment(
                                 profiled_set.mkPDF(m), parsed_flavor,
                                 xmin, xmax, nx, Q2, **kw)
                             for m in range(profiled_set.size)]
                p = profiled_set.uncertainty(prof_vals)
                sigma_a = (p.errminus + p.errplus) / 2.0
                ratio_window[k, lvl] = (sigma_a / sigma_before_window[k]
                                        if sigma_before_window[k] > 0 else np.nan)

            for obs_idx, obs in enumerate(charge_observables):
                parsed = parsed_charge_terms[obs_idx]
                kw = dict(weight_type='1', moment=obs.get('moment', 0))
                prof_vals = [compute_integrated_moment(
                                 profiled_set.mkPDF(m), parsed,
                                 charge_xmin, charge_xmax, nx, Q2, **kw)
                             for m in range(profiled_set.size)]
                p = profiled_set.uncertainty(prof_vals)
                sigma_a = (p.errminus + p.errplus) / 2.0
                ratio_charges_buf[obs_idx, lvl] = (
                    sigma_a / sigma_before_charges[obs_idx]
                    if sigma_before_charges[obs_idx] > 0 else np.nan)

            for ni in range(n_orders):
                kw = dict(weight_type='1', moment=ni)
                prof_vals = [compute_integrated_moment(
                                 profiled_set.mkPDF(m), parsed_flavor,
                                 moment_xmin, moment_xmax, nx, Q2, **kw)
                             for m in range(profiled_set.size)]
                p = profiled_set.uncertainty(prof_vals)
                sigma_a = (p.errminus + p.errplus) / 2.0
                ratio_poly_moments[ni, lvl] = (sigma_a / sigma_before_poly[ni]
                                               if sigma_before_poly[ni] > 0 else np.nan)

            lhapdf.setVerbosity(1)

            win_str = '  '.join(
                f"wm{k+start_n}:{ratio_window[k, lvl]:.3f}" for k in range(n_levels))
            print(f"n={n}:  {win_str}")

            np.savez(results_path,
                     max_n=max_n, start_n=start_n, x0=x0, w=w, weight=weight,
                     ratio_window=ratio_window,
                     ratio_charges=ratio_charges_buf,
                     ratio_poly_moments=ratio_poly_moments,
                     sigma_before_window=sigma_before_window,
                     sigma_before_charges=sigma_before_charges,
                     sigma_before_poly=sigma_before_poly,
                     central_window=window_central,
                     central_charges=central_charges,
                     central_poly=central_poly,
                     charge_labels=np.array([obs.get('label', obs['flavor'])
                                             for obs in charge_observables]),
                     pdf_name=np.array(self._pdf_name),
                     mc2h_dir=np.array(self._mc2h_dir or ''))

        self.n_values             = list(range(start_n, max_n + 1))
        self.ratio_window         = ratio_window
        self.ratio_charges        = ratio_charges_buf[:n_obs] if n_obs else np.empty((0, n_levels))
        self.ratio_poly_moments   = ratio_poly_moments
        self.sigma_before_window  = sigma_before_window
        self.sigma_before_charges = sigma_before_charges
        self.sigma_before_poly    = sigma_before_poly
        self.central_window       = window_central
        self.central_charges      = central_charges
        self.central_poly         = central_poly
        self.charge_labels = [obs.get('label', obs['flavor'])
                               for obs in charge_observables]
        print(f"Results saved → {results_path}")
        return self

    # ------------------------------------------------------------------
    def load(self):
        """Load results saved by a previous run() call without re-running ePump."""
        output_dir   = os.path.abspath(self.cfg['output_dir'])
        results_path = os.path.join(output_dir, 'moment_accumulation',
                                    'results_moment_scan.npz')
        if not os.path.exists(results_path):
            raise FileNotFoundError(
                f"No saved results at {results_path}. Call run() first.")

        saved   = np.load(results_path, allow_pickle=True)
        max_n   = int(saved['max_n'])
        start_n = int(saved.get('start_n', 1))
        n_obs   = len(saved['charge_labels'])

        self.n_values             = list(range(start_n, max_n + 1))
        self.ratio_window         = saved['ratio_window']
        self.ratio_charges        = saved['ratio_charges'][:n_obs]
        self.ratio_poly_moments   = saved['ratio_poly_moments']
        self.sigma_before_window  = saved['sigma_before_window']
        self.sigma_before_charges = saved['sigma_before_charges']
        self.sigma_before_poly    = saved['sigma_before_poly']
        self.central_window       = saved['central_window']
        self.central_charges      = saved['central_charges']
        self.central_poly         = saved['central_poly']
        self.charge_labels        = saved['charge_labels'].tolist()
        self._pdf_name            = str(saved['pdf_name'])
        mc2h_str                  = str(saved['mc2h_dir'])
        self._mc2h_dir            = mc2h_str if mc2h_str else None
        print(f"Loaded ← {results_path}")
        return self

    # ------------------------------------------------------------------
    def plot(self, save=True):
        """Plot σ_after/σ_before vs number of simultaneous profiling moments.

        Panels (left to right, top to bottom):
          - one per window moment 1..max_n
          - one per tensor-charge observable
          - one per polynomial moment n = 0, 1, 2, 3

        Parameters
        ----------
        save : bool
            Write figure to output_dir/moment_accumulation/moment_accumulation.pdf.

        Returns
        -------
        matplotlib.figure.Figure
        """
        if self.ratio_window is None:
            raise RuntimeError("Call run() or load() before plot().")

        import matplotlib.pyplot as plt

        cfg    = self.cfg
        ms_cfg = cfg['moment_scan']
        n_vals = self.n_values
        max_n  = len(n_vals)
        x0     = ms_cfg['x0']
        w      = ms_cfg['w']
        weight = ms_cfg.get('weight', 'gaussian')

        n_obs    = self.ratio_charges.shape[0]
        n_orders = self.ratio_poly_moments.shape[0]
        n_panels = max_n + n_obs + n_orders

        ncols = min(n_panels, 4)
        nrows = (n_panels + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(4 * ncols, 3.5 * nrows),
                                 squeeze=False, sharey=True)
        axes_flat = axes.flatten()
        panel     = 0

        def _draw(ax, ydata, color, title):
            ax.plot(n_vals, ydata, 'o-', color=color)
            ax.axhline(1.0, color='gray', ls='--', lw=0.8)
            ax.set_ylim(0, 1.1)
            ax.set_xticks(n_vals)
            ax.set_xlabel('Moments constrained')
            ax.set_ylabel(r'$\sigma_{\rm after}/\sigma_{\rm before}$')
            ax.set_title(title)

        # Window moment panels
        for k in range(max_n):
            idx = n_vals[k]
            sym = (rf'$g_{{{idx}}}$' if weight == 'gaussian'
                   else rf'$a_{{{idx}}}$')
            if (self.sigma_before_window is not None and
                    self.central_window is not None and
                    abs(self.central_window[k]) > 0):
                pct = (100 * self.sigma_before_window[k]
                       / abs(self.central_window[k]))
                title = rf"Window {sym}   ($\sigma/\mu={pct:.1f}\%$)"
            else:
                title = rf"Window {sym}"
            _draw(axes_flat[panel], self.ratio_window[k], 'C2', title)
            panel += 1

        # Charge observable panels
        for obs_idx in range(n_obs):
            lbl = self.charge_labels[obs_idx] if self.charge_labels else f"obs {obs_idx}"
            if (self.sigma_before_charges is not None and
                    self.central_charges is not None and
                    abs(self.central_charges[obs_idx]) > 0):
                pct = (100 * self.sigma_before_charges[obs_idx]
                       / abs(self.central_charges[obs_idx]))
                title = rf"{lbl}   ($\sigma/\mu={pct:.1f}\%$)"
            else:
                title = lbl
            _draw(axes_flat[panel], self.ratio_charges[obs_idx], 'C0', title)
            panel += 1

        # Polynomial moment panels
        for ni in range(n_orders):
            if (self.sigma_before_poly is not None and
                    self.central_poly is not None and
                    abs(self.central_poly[ni]) > 0):
                pct = (100 * self.sigma_before_poly[ni]
                       / abs(self.central_poly[ni]))
                title = rf'Full $n={ni}$   ($\sigma/\mu={pct:.1f}\%$)'
            else:
                title = rf'Full $n={ni}$'
            _draw(axes_flat[panel], self.ratio_poly_moments[ni], 'C1', title)
            panel += 1

        for k in range(panel, len(axes_flat)):
            axes_flat[k].set_visible(False)

        for r in range(nrows):
            for c in range(1, ncols):
                axes[r, c].tick_params(labelleft=False)
                axes[r, c].set_ylabel('')

        pdf_label = cfg.get('pdf_label', cfg['pdf'])
        fig.suptitle(
            f"{pdf_label} — Moment accumulation  ({cfg['flavor']})\n"
            f"{weight}  $x_0={x0}$  $w={w}$    $Q^2={cfg['Q2']}$ GeV$^2$",
            fontsize=12,
        )
        plt.tight_layout()
        plt.subplots_adjust(wspace=0)

        if save:
            suffix = '_' + os.path.basename(os.path.abspath(cfg['output_dir']))
            out = os.path.join(os.path.abspath(cfg['output_dir']),
                               f"moment_accumulation{suffix}.pdf")
            fig.savefig(out, dpi=150)
            print(f"Plot → {out}")

        return fig

    # ------------------------------------------------------------------
    @classmethod
    def plot_compare(cls, scanners, labels=None, save=False, save_path=None):
        """Class-level comparison plot — overlay results from multiple scanners.

        Parameters
        ----------
        scanners : list[MomentAccumulationScanner]
            Two or more loaded scanners (run() or load() must have been called).
        labels : list[str], optional
            One legend label per scanner.  Auto-generated if None.
        save : bool
            If True, write to first scanner's output_dir.
        save_path : str, optional
            Override save destination.

        Returns
        -------
        matplotlib.figure.Figure

        Examples
        --------
        >>> s1 = MomentAccumulationScanner('runcard_a.py').load()
        >>> s2 = MomentAccumulationScanner('runcard_b.py').load()
        >>> fig = MomentAccumulationScanner.plot_compare([s1, s2])
        """
        return plot_moment_accumulation_comparison(
            scanners, labels=labels, colors=None,
            save=save, save_path=save_path)

    # ------------------------------------------------------------------
    def plot_comparison(self, other, labels=None, colors=None,
                        save=False, save_path=None):
        """Compare this scanner's results against one or more others on the same axes.

        Parameters
        ----------
        other : MomentAccumulationScanner or list[MomentAccumulationScanner]
        labels, colors, save, save_path : forwarded to plot_moment_accumulation_comparison
        """
        scanners = [self] + (other if isinstance(other, list) else [other])
        return plot_moment_accumulation_comparison(
            scanners, labels=labels, colors=colors,
            save=save, save_path=save_path)


# ── Comparison plot (module-level) ─────────────────────────────────────────────

def plot_moment_accumulation_comparison(scanners, labels=None, colors=None,
                                        save=False, save_path=None):
    """Overlay σ_after/σ_before curves from multiple MomentAccumulationScanner runs.

    Parameters
    ----------
    scanners : list[MomentAccumulationScanner]
        Two or more scanners with results loaded (via run() or load()).
    labels : list[str] or None
        Legend label for each scanner.  Defaults to each scanner's pdf_label/pdf.
    colors : list or None
        Matplotlib color specs, one per scanner.  Defaults to C0, C1, C2, …
    save : bool
        Save the figure.
    save_path : str or None
        File to save to.  Defaults to the first scanner's output_dir with name
        moment_accumulation_comparison_<basename>.pdf.

    Returns
    -------
    matplotlib.figure.Figure

    Example (notebook)
    ------------------
    >>> from scan_window_moments import MomentAccumulationScanner
    >>> s1 = MomentAccumulationScanner('runcard_a.py').load()
    >>> s2 = MomentAccumulationScanner('runcard_b.py').load()
    >>> fig = s1.plot_comparison(s2, labels=['CT18NNLO', 'MSHT20'])
    """
    import matplotlib.pyplot as plt

    if not scanners:
        raise ValueError("scanners list is empty")
    for i, s in enumerate(scanners):
        if s.ratio_window is None:
            raise RuntimeError(
                f"Scanner {i} has no results — call run() or load() first.")

    if labels is None:
        labels = [s.cfg.get('pdf_label', s.cfg.get('pdf', f'run {i}'))
                  for i, s in enumerate(scanners)]
    if colors is None:
        colors = [f'C{i}' for i in range(len(scanners))]

    ref    = scanners[0]
    cfg    = ref.cfg
    n_obs  = ref.ratio_charges.shape[0]
    n_orders = ref.ratio_poly_moments.shape[0]
    n_panels = n_obs + n_orders

    ncols = min(n_panels, 4)
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4 * ncols, 3.5 * nrows),
                             squeeze=False, sharey=True)
    axes_flat = axes.flatten()
    panel = 0

    def _draw_multi(ax, ydatas, title):
        for s, ydata, color, label in zip(scanners, ydatas, colors, labels):
            if ydata is None:
                continue
            ax.plot(s.n_values, ydata, 'o-', color=color, label=label)
        ax.axhline(1.0, color='gray', ls='--', lw=0.8)
        ax.set_ylim(0, 1.1)
        all_n = sorted({v for s in scanners for v in s.n_values})
        ax.set_xticks(all_n)
        ax.set_xlabel('Moments constrained')
        ax.set_ylabel(r'$\sigma_{\rm after}/\sigma_{\rm before}$')
        ax.set_title(title)
        ax.legend(fontsize=7)

    for obs_idx in range(n_obs):
        lbl = ref.charge_labels[obs_idx] if ref.charge_labels else f"obs {obs_idx}"
        _draw_multi(axes_flat[panel],
                    [s.ratio_charges[obs_idx] for s in scanners],
                    lbl)
        panel += 1

    for ni in range(n_orders):
        _draw_multi(axes_flat[panel],
                    [s.ratio_poly_moments[ni] for s in scanners],
                    rf'Full $n={ni}$')
        panel += 1

    for k in range(panel, len(axes_flat)):
        axes_flat[k].set_visible(False)

    for r in range(nrows):
        for c in range(1, ncols):
            axes[r, c].tick_params(labelleft=False)
            axes[r, c].set_ylabel('')

    fig.suptitle(
        f"{' vs '.join(labels)} — Moment accumulation ({cfg['flavor']})  "
        f"$Q^2={cfg['Q2']}$ GeV$^2$",
        fontsize=12,
    )
    plt.tight_layout()
    plt.subplots_adjust(wspace=0)

    if save:
        if save_path is None:
            suffix = '_' + os.path.basename(os.path.abspath(cfg['output_dir']))
            save_path = os.path.join(
                os.path.abspath(cfg['output_dir']),
                f"moment_accumulation_comparison{suffix}.pdf")
        fig.savefig(save_path, dpi=150)
        print(f"Plot → {save_path}")

    return fig


# ── Terminal entry point ───────────────────────────────────────────────────────

def main():
    matplotlib.use('Agg')   # non-interactive backend for terminal use

    if len(sys.argv) < 2:
        print(
            f"Usage: {sys.argv[0]} runcard.py [run|recompute|moments|charges|moment_scan]",
            file=sys.stderr,
        )
        sys.exit(1)

    cmd = sys.argv[2] if len(sys.argv) >= 3 else 'run'

    if cmd == 'moment_scan':
        scanner = MomentAccumulationScanner(sys.argv[1])
        scanner.run()
        scanner.plot(save=True)
        return

    scanner = WindowMomentScanner(sys.argv[1])

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
