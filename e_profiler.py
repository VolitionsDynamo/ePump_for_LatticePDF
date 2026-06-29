#!/home/daniel/miniconda3/envs/apfelpp/bin/python3
"""
e_profiler.py - Standalone, robust CLI script to interface Python observable computations with the ePump profiling engine.
"""

import sys
import os
import argparse
import re
import tempfile
import shutil
import atexit
import numpy as np
import lhapdf

# Find relative path of default LHAPDF examples in the repository
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LATTICE_TEST_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, 'ePump_kp20221218', 'lattice_moment_test1'))
MC2H_SRC = os.path.join(SCRIPT_DIR, 'mc2hessian', 'src')

def setup_lhapdf_path(custom_path=None):
    """
    Sets up the LHAPDF_DATA_PATH environment variable portably.
    """
    paths = []
    if custom_path:
        paths.append(os.path.abspath(custom_path))
    paths.append(LATTICE_TEST_DIR)
    paths.append(os.getcwd())
    
    # Append any pre-existing paths
    env_path = os.environ.get('LHAPDF_DATA_PATH', '')
    if env_path:
        paths.append(env_path)
        
    os.environ['LHAPDF_DATA_PATH'] = ':'.join(paths)

def find_pdf_dir(pdf_name):
    """
    Finds the absolute path of the starting PDF set directory.
    Checks:
    1. Current working directory.
    2. LHAPDF search paths via lhapdf.paths().
    3. LHAPDF_DATA_PATH environment variable.
    """
    # 1. Check current directory
    if os.path.exists(pdf_name):
        return os.path.abspath(pdf_name)
    # 2. Try importing lhapdf and using search paths
    try:
        import lhapdf
        for p in lhapdf.paths():
            full_path = os.path.join(p, pdf_name)
            if os.path.exists(full_path):
                return os.path.abspath(full_path)
    except Exception:
        pass
    # 3. Check environment LHAPDF_DATA_PATH
    env_path = os.environ.get("LHAPDF_DATA_PATH", "")
    for p in env_path.split(":"):
        if p:
            full_path = os.path.join(p, pdf_name)
            if os.path.exists(full_path):
                return os.path.abspath(full_path)
    return None

def ensure_pdf_symlink(pdf_name):
    """
    Ensures that a symlink or directory with name `pdf_name` exists in the current
    working directory, pointing to the actual PDF set directory in LHAPDF.
    """
    if os.path.exists(pdf_name):
        return
    
    pdf_dir = find_pdf_dir(pdf_name)
    if pdf_dir:
        abs_pdf_name = os.path.abspath(pdf_name)
        if pdf_dir != abs_pdf_name:
            print(f"Creating symlink for PDF set {pdf_name} -> {pdf_dir} ...")
            try:
                os.symlink(pdf_dir, pdf_name)
            except Exception as e:
                print(f"Warning: Failed to create symlink: {e}", file=sys.stderr)

def detect_pdf_error_type(pdf_name):
    """
    Reads the LHAPDF .info file for pdf_name and returns the ErrorType string.
    Returns one of: 'hessian', 'symmhessian', 'replicas', or 'unknown'.
    """
    pdf_dir = find_pdf_dir(pdf_name)
    if pdf_dir is None:
        print(f"Warning: Could not locate PDF set '{pdf_name}' to check ErrorType.", file=sys.stderr)
        return 'unknown'

    info_path = os.path.join(pdf_dir, f"{pdf_name}.info")
    if not os.path.exists(info_path):
        print(f"Warning: No .info file found at {info_path}", file=sys.stderr)
        return 'unknown'

    with open(info_path, 'r') as f:
        for line in f:
            if line.strip().startswith('ErrorType:'):
                return line.split(':', 1)[1].strip().lower()
    return 'unknown'


def _import_mc2hlib():
    """Imports mc2hlib, handling the fastcache dependency gracefully."""
    if MC2H_SRC not in sys.path:
        sys.path.insert(0, MC2H_SRC)
    try:
        import fastcache
    except ImportError:
        import functools
        import types
        fastcache_mod = types.ModuleType('fastcache')
        fastcache_mod.lru_cache = functools.lru_cache
        sys.modules['fastcache'] = fastcache_mod

    from mc2hlib.common import load_pdf, compress_X_abs, get_limits
    from mc2hlib.lh import hessian_from_lincomb, load_replica_2, write_replica
    return load_pdf, compress_X_abs, get_limits, hessian_from_lincomb, load_replica_2, write_replica


def expand_symm_to_asymm(set_dir, set_name, neig):
    """
    Expands a symmetric Hessian set (neig+1 members) to asymmetric +/- pairs
    (2*neig+1 members) as required by ePump.
    """
    _, _, _, _, load_replica_2, write_replica = _import_mc2hlib()

    base = os.path.join(set_dir, set_name)

    # Read central member grids
    central_header, central_grids = load_replica_2(0, base)

    # Process in reverse to avoid filename collisions
    for k in range(neig, 0, -1):
        _, grids_k = load_replica_2(k, base)

        # Compute minus direction: 2*central - plus
        grids_minus = central_grids.multiply(2).subtract(grids_k)

        # Rename original file k -> 2k-1 (plus direction)
        old_path = f"{base}_{k:04d}.dat"
        plus_path = f"{base}_{(2*k-1):04d}.dat"
        os.rename(old_path, plus_path)

        # Write minus direction as member 2k
        hess_header = b"PdfType: error\nFormat: lhagrid1\n"
        write_replica(2*k, base, hess_header, grids_minus)

    # Update .info file
    info_path = f"{base}.info"
    with open(info_path, 'r') as f:
        lines = f.readlines()
    with open(info_path, 'w') as f:
        for line in lines:
            if line.strip().startswith('NumMembers:'):
                f.write(f"NumMembers: {2*neig + 1}\n")
            elif line.strip().startswith('ErrorType:'):
                f.write("ErrorType: hessian\n")
            else:
                f.write(line)

    print(f"  Expanded symmetric Hessian ({neig+1} members) to asymmetric pairs ({2*neig+1} members).")


def convert_mc_to_hessian(pdf_name, neig=50, Q=1.0, epsilon=1000.0,
                          output_dir=None, max_nf=3):
    """
    Converts an MC replica PDF set to an asymmetric Hessian set compatible with ePump.
    Returns (converted_set_name, converted_set_dir).
    """
    load_pdf, compress_X_abs, get_limits, hessian_from_lincomb, _, _ = _import_mc2hlib()

    print(f"  MC-to-Hessian conversion at Q = {Q} GeV, neig = {neig}, epsilon = {epsilon}")

    # Ensure the PDF set's parent dir is first in LHAPDF_DATA_PATH so
    # hessian_from_lincomb (which reads from lhapdf.paths()[0]) can find it.
    pdf_dir = find_pdf_dir(pdf_name)
    if pdf_dir:
        parent_dir = os.path.dirname(pdf_dir)
        current = os.environ.get('LHAPDF_DATA_PATH', '')
        os.environ['LHAPDF_DATA_PATH'] = parent_dir + ':' + current
        lhapdf.setPaths([parent_dir] + lhapdf.paths())

    # Step 1: Load PDF and build covariance matrix
    pdf, fl, xgrid = load_pdf(pdf_name, Q, max_nf, False)
    nx, nf = xgrid.n, fl.n
    X = (pdf.xfxQ.reshape(pdf.n_rep, nx * nf) - pdf.f0.reshape(nx * nf)).T

    # Epsilon masking
    l = get_limits(X.T)
    diff = (l.up1s - l.low1s) / 2
    std = np.std(X, axis=1)
    mask = (np.abs((diff - std) / diff) < epsilon)
    print(f"  Keeping {np.count_nonzero(mask)} / {len(mask)} points (epsilon = {epsilon})")
    X = X[mask, :]

    # Step 2: SVD compression
    vec, cov = compress_X_abs(X, neig)

    # Step 3: Export symmetric Hessian set
    set_name = f"{pdf_name}_hessian_{neig}"
    hessian_from_lincomb(pdf, vec, set_name=set_name, folder=output_dir)
    print(f"  Symmetric Hessian set written to {os.path.join(output_dir or '', set_name)}/")

    # Step 4: Expand symmetric -> asymmetric pairs for ePump
    set_dir = os.path.join(output_dir or '', set_name)
    expand_symm_to_asymm(set_dir, set_name, neig)

    return set_name, os.path.abspath(set_dir)


def run_closure_test(mc_pdf_name, hessian_pdf_name, measurements, args):
    """
    Compares uncertainty bands between the original MC set and the converted
    Hessian set to validate the conversion.
    """
    print("\n============================================================")
    print("              MC-TO-HESSIAN CLOSURE TEST                     ")
    print("============================================================")

    # Load both sets
    mc_set = lhapdf.getPDFSet(mc_pdf_name)
    mc_members = mc_set.mkPDFs()
    hess_set = lhapdf.getPDFSet(hessian_pdf_name)
    hess_members = hess_set.mkPDFs()

    n_mc = len(mc_members)
    n_hess = len(hess_members)
    n_obs = len(measurements)

    mc_vals = np.zeros((n_mc, n_obs))
    hess_vals = np.zeros((n_hess, n_obs))

    parsed_flavors_cache = []
    for m in measurements:
        parsed_flavors_cache.append(parse_flavor_expression(m['flavor']))

    for i, member in enumerate(mc_members):
        print(f"\r  MC member {i+1}/{n_mc}", end="", flush=True)
        for j, m in enumerate(measurements):
            parsed_terms = parsed_flavors_cache[j]
            if m['obs_type'] == 'point':
                mc_vals[i, j] = evaluate_pdf_combination(member, parsed_terms, m['x'], m['Q2'])
            else:
                mc_vals[i, j] = compute_integrated_moment(
                    member, parsed_terms, m['xmin'], m['xmax'], m['nx'], m['Q2'],
                    weight_type=m['weight'], moment=m['moment'])
    print()

    for i, member in enumerate(hess_members):
        print(f"\r  Hessian member {i+1}/{n_hess}", end="", flush=True)
        for j, m in enumerate(measurements):
            parsed_terms = parsed_flavors_cache[j]
            if m['obs_type'] == 'point':
                hess_vals[i, j] = evaluate_pdf_combination(member, parsed_terms, m['x'], m['Q2'])
            else:
                hess_vals[i, j] = compute_integrated_moment(
                    member, parsed_terms, m['xmin'], m['xmax'], m['nx'], m['Q2'],
                    weight_type=m['weight'], moment=m['moment'])
    print()

    CENTRAL_TOL = 0.05
    SIGMA_TOL_LOW = 0.90
    SIGMA_TOL_HIGH = 1.10

    print(f"  {'Obs':<4} | {'MC Central':<12} {'MC Err':<12} | {'Hess Central':<12} {'Hess Err':<12} | {'C Ratio':<8} {'E Ratio':<8} | Result")
    print("-" * 100)

    all_pass = True
    for j in range(n_obs):
        # MC: central = mean of replicas (skip member 0), error = std
        mc_central = np.mean(mc_vals[1:, j])
        mc_sigma = np.std(mc_vals[1:, j], ddof=1)

        # Hessian: use LHAPDF uncertainty
        hess_col = [float(v) for v in hess_vals[:, j]]
        hess_unc = hess_set.uncertainty(hess_col)
        hess_central = hess_unc.central
        hess_sigma = (hess_unc.errplus + hess_unc.errminus) / 2.0

        central_ratio = abs(hess_central - mc_central) / (abs(mc_central) + 1e-30)
        sigma_ratio = hess_sigma / (mc_sigma + 1e-30)

        c_ok = central_ratio < CENTRAL_TOL
        e_ok = SIGMA_TOL_LOW <= sigma_ratio <= SIGMA_TOL_HIGH
        status = "PASS" if (c_ok and e_ok) else "FAIL"
        if status == "FAIL":
            all_pass = False

        print(f"  {j+1:<4} | {mc_central:<12.6f} {mc_sigma:<12.6f} | "
              f"{hess_central:<12.6f} {hess_sigma:<12.6f} | "
              f"{central_ratio:<8.4f} {sigma_ratio:<8.4f} | {status}")

    print("----")
    if all_pass:
        print("  Closure test PASSED. MC and Hessian uncertainties agree within tolerance.")
    else:
        print("  Closure test FAILED for some observables. Proceeding anyway — review results above.")
    print("============================================================\n")
    return all_pass


def generate_in_file(filepath, base_name, n_ev_pairs, n_obs, pdf_name):
    """
    Writes the .in control file for ePump.
    """
    ensure_pdf_symlink(pdf_name)
    pdf_in_path = f"./{pdf_name}/{pdf_name}"
    pdf_out_path = f"./{base_name}/{base_name}"
    
    with open(filepath, "w") as f:
        f.write("+++ N(EV pairs)                       N(Data Sets)   PDFtype(C/L/N)    DiagonalQuad(Y/N)    Dyn_Tol?(Y/N)  Tol_squared \n")
        f.write(f"        {n_ev_pairs:<38}{1:<17}L                  N                   N            1    \n")
        f.write("+++ ObservableFile                    N(Observables)  Data?(Y/N)      Error_type     Weight          \n")
        f.write(f"        {base_name:<36}{n_obs:<16}Y                4           1\n")
        f.write("+++     PDFname                       PDFout   \n")
        f.write(f"     {pdf_in_path}    {pdf_out_path}\n")
        f.write("# Generated dynamically by e_profiler.py\n")

def parse_flavor_expression(expression):
    """
    Safely parses an algebraic linear combination of flavor names or PIDs.
    Returns a list of tuples: (coefficient, PID)
    """
    clean_expr = expression.replace(' ', '')
    
    # Match terms: optional coefficient (e.g. 2, -1.5, +), operator, and flavor/PID
    pattern = r'([+-]?\s*\d*\.?\d*)\s*\*?\s*([a-zA-Z]+|-?\d+)'
    matches = re.findall(pattern, clean_expr)
    
    parsed_terms = []
    flavor_map = {
        'g': 21, 'gluon': 21,
        'd': 1, 'u': 2, 's': 3, 'c': 4, 'b': 5, 't': 6,
        'dbar': -1, 'ubar': -2, 'sbar': -3, 'cbar': -4, 'bbar': -5, 'tbar': -6
    }
    
    for coeff_str, name_str in matches:
        coeff_str = coeff_str.replace(' ', '')
        if not coeff_str or coeff_str == '+':
            coeff = 1.0
        elif coeff_str == '-':
            coeff = -1.0
        else:
            coeff = float(coeff_str)
            
        name_str = name_str.strip()
        if name_str.lower() in flavor_map:
            pid = flavor_map[name_str.lower()]
        else:
            try:
                pid = int(name_str)
            except ValueError:
                raise ValueError(f"Unknown PDF flavor or PID: '{name_str}'")
                
        parsed_terms.append((coeff, pid))
        
    if not parsed_terms:
        raise ValueError(f"Could not parse flavor expression: '{expression}'")
        
    return parsed_terms

def evaluate_pdf_combination(pdf_member, parsed_terms, x, Q2):
    """
    Evaluates any linear combination of PDFs for a given member, x, and Q^2.
    """
    # LHAPDF member xfxQ2 returns x * f(x, Q^2)
    return sum(coeff * pdf_member.xfxQ2(pid, x, Q2) for coeff, pid in parsed_terms)

def compute_integrated_moment(pdf_member, parsed_terms, xmin, xmax, nx, Q2, weight_type='1', moment=0):
    """
    Computes numerical moment integration over [xmin, xmax] using logarithmic x spacing and specific weight functions.
    """
    x_grid = 10**np.linspace(np.log10(xmin), np.log10(xmax), nx)
    
    # Evaluate flavor combination: returns x * f(x, Q^2)
    xf_values = np.array([evaluate_pdf_combination(pdf_member, parsed_terms, x, Q2) for x in x_grid])
    
    if weight_type == 'gaussian':
        mean = (xmin + xmax) / 2.0
        sigma = (xmax - xmin) / 2.0
        
        if sigma == 0:
            norm = 0.0
        else:
            norm = np.sqrt((moment**2) / (2.0 * np.pi * (sigma**2)))
            
        if sigma == 0:
            gaussian_weights = np.zeros_like(x_grid)
        else:
            gaussian_weights = norm * np.exp(- (moment**2) * ((x_grid - mean)**2) / (2.0 * (sigma**2)))
            
        # Replaces x^n completely: Integrand = Gaussian(x) * f(x) = Gaussian(x) * (x*f(x)) / x
        y_values = gaussian_weights * (xf_values / x_grid)
    else:
        # Default weight: Integrand = x^n * f(x) = x^{n-1} * (x*f(x))
        y_values = (x_grid**(moment - 1)) * xf_values
        
    # Trapezoidal integration
    integral = np.sum((y_values[1:] + y_values[:-1]) * np.diff(x_grid)) / 2.0
    return float(integral)

def compute_observable(pdf_set_members, parsed_terms, args):
    """
    Computes observable (point or moment) for all members of a PDF set.
    """
    results = []
    for member in pdf_set_members:
        if args.obs_type == 'point':
            val = evaluate_pdf_combination(member, parsed_terms, args.x, args.Q2)
        else:
            val = compute_integrated_moment(member, parsed_terms, args.xmin, args.xmax, args.nx, args.Q2, args.weight, args.moment)
        results.append(val)
    return np.array(results)

def generate_theory_file(filepath, pdf_members, measurements, args):
    """
    Computes theoretical predictions for all PDF members and all measurements,
    taking into account per-measurement flavor combinations, types, bounds, and weight functions.
    Writes them to a .theory file.
    """
    n_members = len(pdf_members)
    n_obs = len(measurements)
    
    Ms = np.zeros((n_members, n_obs))
    
    parsed_flavors_cache = []
    for m in measurements:
        parsed_terms = parse_flavor_expression(m['flavor'])
        parsed_flavors_cache.append(parsed_terms)
        
    for i, member in enumerate(pdf_members):
        for j, m in enumerate(measurements):
            x_val = m['x']
            Q2_val = m['Q2']
            parsed_terms = parsed_flavors_cache[j]
            
            if m['obs_type'] == 'point':
                val = evaluate_pdf_combination(member, parsed_terms, x_val, Q2_val)
            else:
                val = compute_integrated_moment(
                    member, parsed_terms, m['xmin'], m['xmax'], m['nx'], Q2_val,
                    weight_type=m['weight'], moment=m['moment']
                )
            Ms[i, j] = val
            
    with open(filepath, "w") as f:
        f.write("*\n*\n*\n")
        f.write("Theory column\n")
        f.write(f"{n_obs}\n")
        for i in range(n_members):
            f.write(f"PDF_{i}_Set\n")
            formatted_row = "".join(f"{Ms[i, j]:12.6f}" for j in range(n_obs))
            f.write(formatted_row + "\n")
            
    return Ms

def generate_data_file(filepath, measurements):
    """
    Writes measurements to a .data file in ePump format.
    """
    first_m = measurements[0]
    C = len(first_m['cor_sys'])
    
    with open(filepath, "w") as f:
        f.write("*\n* Generated dynamically by e_profiler.py\n*\n")
        f.write("*  x    Q2    value         stat        uncor_sys   cor_sys...\n")
        f.write("#corr_err  data stat uncor_sys cor_begin\n")
        f.write(f"{C} 3 4 5 6\n")
        f.write("x  Q2  value         error\n")
        
        for m in measurements:
            vals = [m['x'], m['Q2'], m['value'], m['stat'], m['uncor_sys']]
            cor_sys = list(m['cor_sys'])
            while len(cor_sys) < C:
                cor_sys.append(0.0)
            cor_sys = cor_sys[:C]
            vals.extend(cor_sys)
            
            row_str = " ".join(f"{val}" for val in vals)
            f.write(row_str + "\n")

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='ePump Hessian profiling pipeline for Point/Moment PDF observables.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Core/PDF Options
    parser.add_argument('--pdf', type=str, default='CT18NNLO',
                        help='Starting PDF set name.')
    
    # Observable Options
    parser.add_argument('--obs-type', type=str, choices=['point', 'moment'], default='moment',
                        help='Observable type to profile.')
    parser.add_argument('--x', type=float, default=0.5,
                        help='x value for point observable.')
    parser.add_argument('--xmin', type=float, default=0.1,
                        help='xmin for integrated moment.')
    parser.add_argument('--xmax', type=float, default=0.7,
                        help='xmax for integrated moment.')
    parser.add_argument('--nx', type=int, default=100,
                        help='Number of integration points.')
    parser.add_argument('--Q2', type=float, default=4.0,
                        help='Scale Q^2 for evaluation in GeV^2.')
    parser.add_argument('--flavor', type=str, default='u-d',
                        help='Algebraic linear combination of PDF flavors (e.g., "u-d", "2*u - 1.5*d", "ubar - dbar").')
    parser.add_argument('--weight', type=str, choices=['1', 'gaussian'], default='1',
                        help='Weight function definition to use. "1" for default x^n weight, "gaussian" for a normalized gaussian.')
    parser.add_argument('--moment', type=int, default=0,
                        help='Number of the moment (n) to calculate.')

    # Measurement/Data Options
    parser.add_argument('--measurement', action='append', nargs='+', type=float,
                        help='Repeatable measurement: x Q2 value stat [uncor_sys] [cor_sys...]')
    parser.add_argument('--data-file', type=str, default=None,
                        help='Path to a space- or comma-separated text file with columns: x Q2 value stat [uncor_sys] [cor_sys...]')

    # Output / ePump Options
    parser.add_argument('--name', type=str, default='CT18NNLO_mom_test1',
                        help='Output base name for ePump input/output files.')
    parser.add_argument('--epump-path', type=str, default='./ePump_kp20221218/src/UpdatePDFs',
                        help='Path to the compiled UpdatePDFs binary.')

    # MC-to-Hessian conversion options
    parser.add_argument('--mc2h-neig', type=int, default=50,
                        help='Number of eigenvectors for MC-to-Hessian conversion.')
    parser.add_argument('--mc2h-Q', type=float, default=1.0,
                        help='Energy scale Q (GeV) for MC-to-Hessian conversion.')
    parser.add_argument('--mc2h-epsilon', type=float, default=1000.0,
                        help='Epsilon tolerance for mc2hessian point selection.')
    parser.add_argument('--mc2h-max-nf', type=int, default=3,
                        help='Max number of active flavors for mc2hessian.')
    parser.add_argument('--mc2h-output-dir', type=str, default=None,
                        help='Directory for converted Hessian set (default: auto temp dir).')
    parser.add_argument('--force-mc2h', action='store_true',
                        help='Force MC-to-Hessian conversion even if the set appears to be Hessian.')
    parser.add_argument('--closure-test', action='store_true',
                        help='Run closure test comparing MC and converted Hessian uncertainties.')

    return parser.parse_args()

def parse_measurements(args):
    """
    Parses and merges measurements from CLI flags and/or external files (JSON or text).
    Returns a list of dictionaries, where each measurement dictionary has keys:
    'x', 'Q2', 'value', 'stat', 'uncor_sys', 'cor_sys', 'obs_type', 'flavor', 'xmin', 'xmax', 'nx', 'weight', 'moment'
    """
    measurements = []

    active_config = {
        'obs_type': args.obs_type,
        'flavor': args.flavor,
        'xmin': args.xmin,
        'xmax': args.xmax,
        'nx': args.nx,
        'weight': args.weight,
        'moment': args.moment
    }

    def make_measurement_dict(vals, config):
        m_full = list(vals)
        while len(m_full) < 5:
            m_full.append(0.0)
        
        return {
            'x': m_full[0],
            'Q2': m_full[1],
            'value': m_full[2],
            'stat': m_full[3],
            'uncor_sys': m_full[4],
            'cor_sys': m_full[5:],
            'obs_type': config['obs_type'],
            'flavor': config['flavor'],
            'xmin': config['xmin'],
            'xmax': config['xmax'],
            'nx': config['nx'],
            'weight': config['weight'],
            'moment': config['moment']
        }

    # 1. Parse from repeatable CLI flags (use global CLI configs)
    if args.measurement:
        for m in args.measurement:
            if len(m) < 4:
                print(f"Error: CLI measurement must contain at least [x, Q2, value, stat]. Got: {m}", file=sys.stderr)
                sys.exit(1)
            measurements.append(make_measurement_dict(m, active_config))

    # 2. Parse from text or JSON file
    if args.data_file:
        if not os.path.exists(args.data_file):
            print(f"Error: Data file not found: {args.data_file}", file=sys.stderr)
            sys.exit(1)
            
        if args.data_file.lower().endswith('.json'):
            import json
            try:
                with open(args.data_file, 'r') as f:
                    data_list = json.load(f)
                for item in data_list:
                    item_config = {
                        'obs_type': item.get('obs_type', active_config['obs_type']),
                        'flavor': item.get('flavor', active_config['flavor']),
                        'xmin': float(item.get('xmin', active_config['xmin'])),
                        'xmax': float(item.get('xmax', active_config['xmax'])),
                        'nx': int(item.get('nx', active_config['nx'])),
                        'weight': item.get('weight', active_config['weight']),
                        'moment': int(item.get('moment', active_config['moment']))
                    }
                    vals = [
                        float(item['x']),
                        float(item['Q2']),
                        float(item['value']),
                        float(item['stat']),
                        float(item.get('uncor_sys', 0.0))
                    ]
                    if 'cor_sys' in item:
                        vals.extend([float(v) for v in item['cor_sys']])
                    measurements.append(make_measurement_dict(vals, item_config))
            except Exception as e:
                print(f"Error parsing JSON data file: {e}", file=sys.stderr)
                sys.exit(1)
        else:
            with open(args.data_file, 'r') as f:
                for line in f:
                    line_str = line.strip()
                    if not line_str:
                        continue
                    if line_str.startswith('#') or line_str.startswith('*'):
                        kv_pairs = re.findall(r'(\w+)\s*=\s*(?:"([^"]*)"|\'([^\']*)\'|([^\s]+))', line_str)
                        for key, val1, val2, val3 in kv_pairs:
                            val = val1 or val2 or val3
                            key_lower = key.lower()
                            if key_lower in active_config:
                                if key_lower == 'nx' or key_lower == 'moment':
                                    active_config[key_lower] = int(val)
                                elif key_lower in ['xmin', 'xmax']:
                                    active_config[key_lower] = float(val)
                                else:
                                    active_config[key_lower] = val
                        continue
                        
                    parts = line_str.replace(',', ' ').split()
                    try:
                        vals = [float(p) for p in parts]
                        if len(vals) < 4:
                            print(f"Warning: Ignoring malformed line: {line_str}", file=sys.stderr)
                            continue
                        measurements.append(make_measurement_dict(vals, active_config))
                    except ValueError:
                        print(f"Warning: Ignoring non-numeric line: {line_str}", file=sys.stderr)

    if not measurements:
        print("Error: No measurements supplied. Use --measurement or --data-file.", file=sys.stderr)
        sys.exit(1)

    return measurements

def run_epump(epump_path, base_name):
    """
    Invokes the ePump compiled binary for a given base_name,
    validates the exit status, and parses key profiling metrics from the .out file.
    """
    if not os.path.exists(epump_path):
        print(f"Error: ePump binary not found at: {epump_path}", file=sys.stderr)
        print("Please compile ePump or specify the correct path using --epump-path", file=sys.stderr)
        sys.exit(1)
        
    print(f"Executing ePump profiling: {epump_path} {base_name} ...")
    os.makedirs(base_name, exist_ok=True)
    import subprocess
    result = subprocess.run([epump_path, base_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        print(f"ePump execution failed with exit code {result.returncode}", file=sys.stderr)
        print(f"STDOUT:\n{result.stdout}", file=sys.stderr)
        print(f"STDERR:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)
        
    out_file = f"{base_name}.out"
    if not os.path.exists(out_file):
        print(f"Error: ePump output file {out_file} was not created.", file=sys.stderr)
        sys.exit(1)
        
    print(f"ePump executed successfully. Parsing results from {out_file}...")
    
    chi2_old = None
    chi2_new = None
    
    import re
    with open(out_file, 'r') as f:
        content = f.read()
        
    chi2_old_match = re.search(r'Total\s+Chi2\s*\(z=0\)\s*=\s*([\d\.]+)', content)
    chi2_new_match = re.search(r'Total\s+Chi2\s*\(z=z0\)\s*=\s*([\d\.]+)', content)
    
    if chi2_old_match:
        chi2_old = float(chi2_old_match.group(1))
    if chi2_new_match:
        chi2_new = float(chi2_new_match.group(1))
        
    old_preds = re.findall(r'Old\s+(\d+)\s+([\d\.E+-]+)\s+([\d\.E+-]+)', content)
    new_preds = re.findall(r'New\s+(\d+)\s+([\d\.E+-]+)\s+([\d\.E+-]+)', content)
    
    print("\n============================================================")
    print("                    E-PUMP PROFILING SUMMARY                 ")
    print("============================================================")
    if chi2_old is not None and chi2_new is not None:
        print(f"Initial total Chi2 (before profiling): {chi2_old:.4f}")
        print(f"Profiled total Chi2 (after profiling): {chi2_new:.4f}")
        print(f"Delta Chi2 reduction:                 {chi2_old - chi2_new:.4f}")
    
    if old_preds and new_preds:
        print("\nPredicted Observables & Uncertainties (Before vs After):")
        print(f"  {'Obs ID':<8} | {'Original Value':<15} {'Original Error':<15} | {'Profiled Value':<15} {'Profiled Error':<15}")
        print("-" * 80)
        for old, new in zip(old_preds, new_preds):
            obs_id = old[0]
            old_val, old_err = float(old[1]), float(old[2])
            new_val, new_err = float(new[1]), float(new[2])
            print(f"  {obs_id:<8} | {old_val:<15.6e} {old_err:<15.6e} | {new_val:<15.6e} {new_err:<15.6e}")
    print("============================================================\n")

def report_pdf_comparison(original_set, original_members, updated_set, updated_members, measurements, args):
    """
    Computes original and updated observables for all members,
    calculates their uncertainties using LHAPDF's built-in formulas,
    and prints a detailed comparison report.
    """
    n_members = len(original_members)
    n_obs = len(measurements)
    
    original_vals = np.zeros((n_members, n_obs))
    updated_vals = np.zeros((n_members, n_obs))
    
    # Pre-parse flavor expressions
    parsed_flavors_cache = []
    for m in measurements:
        parsed_terms = parse_flavor_expression(m['flavor'])
        parsed_flavors_cache.append(parsed_terms)
        
    for i in range(n_members):
        orig_member = original_members[i]
        upd_member = updated_members[i]
        
        for j, m in enumerate(measurements):
            x_val = m['x']
            Q2_val = m['Q2']
            parsed_terms = parsed_flavors_cache[j]
            
            # Original observable calculation
            if m['obs_type'] == 'point':
                val_orig = evaluate_pdf_combination(orig_member, parsed_terms, x_val, Q2_val)
            else:
                val_orig = compute_integrated_moment(
                    orig_member, parsed_terms, m['xmin'], m['xmax'], m['nx'], Q2_val,
                    weight_type=m['weight'], moment=m['moment']
                )
            original_vals[i, j] = val_orig
            
            # Updated observable calculation
            if m['obs_type'] == 'point':
                val_upd = evaluate_pdf_combination(upd_member, parsed_terms, x_val, Q2_val)
            else:
                val_upd = compute_integrated_moment(
                    upd_member, parsed_terms, m['xmin'], m['xmax'], m['nx'], Q2_val,
                    weight_type=m['weight'], moment=m['moment']
                )
            updated_vals[i, j] = val_upd
            
    print("\n============================================================")
    print("                PDF COMPARISON & UNCERTAINTY REPORT          ")
    print("============================================================")
    print(f"  {'Obs':<4} | {'Original (LHAPDF)':<30} | {'Profiled (LHAPDF)':<30}")
    print(f"  {'ID':<4} | {'Central':<9} {'Err-':<9} {'Err+':<9} | {'Central':<9} {'Err-':<9} {'Err+':<9}")
    print("-" * 80)
    
    for j, m in enumerate(measurements):
        # Convert numpy array column to clean Python float list for LHAPDF Cython compatibility
        orig_col = [float(x) for x in original_vals[:, j]]
        upd_col = [float(x) for x in updated_vals[:, j]]
        
        orig_unc = original_set.uncertainty(orig_col)
        upd_unc = updated_set.uncertainty(upd_col)
        
        print(f"  {j+1:<4} | {orig_unc.central:<9.5f} -{orig_unc.errminus:<8.5f} +{orig_unc.errplus:<8.5f} | "
              f"{upd_unc.central:<9.5f} -{upd_unc.errminus:<8.5f} +{upd_unc.errplus:<8.5f}")
    print("============================================================\n")

_DEFAULT_EPUMP_PATH = os.path.join(SCRIPT_DIR, 'ePump_kp20221218', 'src', 'UpdatePDFs')

class EProfiler:
    """High-level interface for the ePump PDF profiling pipeline.

    Loads the PDF set on construction so subsequent calls are self-contained.
    For MC replica sets, call convert_to_hessian() before generate_files().
    All file I/O is relative to the working directory at the time of the call.

    Typical usage (Hessian set)::

        ep = EProfiler("CT18NNLO", "my_run")
        ep.add_measurement(x=0.5, Q2=4.0, value=0.151, stat=0.040,
                           obs_type='moment', flavor='u-d',
                           xmin=0.1, xmax=0.7, nx=100)
        ep.generate_files()
        ep.run()
        ep.report()
    """

    def __init__(self, pdf_set_name, run_name, epump_path=None):
        self.pdf_set_name = pdf_set_name
        self.run_name = run_name
        self.epump_path = os.path.abspath(epump_path or _DEFAULT_EPUMP_PATH)
        self.measurements = []
        self.profiled_set = None
        self.profiled_members = None
        setup_lhapdf_path()
        self.pdf_set = lhapdf.getPDFSet(pdf_set_name)
        self.pdf_members = self.pdf_set.mkPDFs()

    def add_measurement(self, x, Q2, value, stat,
                        uncor_sys=0.0, cor_sys=None,
                        obs_type='moment', flavor='u-d',
                        xmin=0.1, xmax=0.9, nx=100,
                        weight='1', moment=0):
        """Append a measurement. All measurements must have equal-length cor_sys lists."""
        self.measurements.append({
            'x': float(x), 'Q2': float(Q2),
            'value': float(value), 'stat': float(stat),
            'uncor_sys': float(uncor_sys),
            'cor_sys': list(cor_sys) if cor_sys is not None else [],
            'obs_type': obs_type, 'flavor': flavor,
            'xmin': float(xmin), 'xmax': float(xmax),
            'nx': int(nx), 'weight': weight, 'moment': int(moment),
        })

    def convert_to_hessian(self, neig=50, Q=1.0, epsilon=1000.0,
                           output_dir=None, max_nf=3):
        """Convert an MC replica PDF set to asymmetric Hessian format in-place.

        If the current set is already Hessian, this is a no-op. When output_dir
        is None, a temporary directory is created and cleaned up at process exit.
        """
        error_type = detect_pdf_error_type(self.pdf_set_name)
        if error_type not in ('replicas', 'mc'):
            print(f"PDF '{self.pdf_set_name}' is already '{error_type}'; skipping conversion.")
            return
        if output_dir is None:
            import tempfile, atexit, shutil
            output_dir = tempfile.mkdtemp(prefix='epump_mc2h_')
            atexit.register(lambda d=output_dir: shutil.rmtree(d, ignore_errors=True))
        else:
            os.makedirs(output_dir, exist_ok=True)
        converted_name, _ = convert_mc_to_hessian(
            self.pdf_set_name, neig=neig, Q=Q,
            epsilon=epsilon, output_dir=output_dir, max_nf=max_nf
        )
        setup_lhapdf_path(custom_path=output_dir)
        self.pdf_set_name = converted_name
        self.pdf_set = lhapdf.getPDFSet(converted_name)
        self.pdf_members = self.pdf_set.mkPDFs()

    def generate_files(self):
        """Write the .in, .data, and .theory files for ePump."""
        if not self.measurements:
            raise ValueError("No measurements added; call add_measurement() first.")
        n_ev_pairs = (len(self.pdf_members) - 1) // 2
        generate_in_file(f"{self.run_name}.in", self.run_name,
                         n_ev_pairs, len(self.measurements), self.pdf_set_name)
        generate_data_file(f"{self.run_name}.data", self.measurements)
        generate_theory_file(f"{self.run_name}.theory",
                             self.pdf_members, self.measurements, self)

    def run(self):
        """Run the ePump binary and load the profiled PDF set."""
        run_epump(self.epump_path, self.run_name)
        self.profiled_set = lhapdf.getPDFSet(self.run_name)
        self.profiled_members = self.profiled_set.mkPDFs()

    def report(self):
        """Print a before/after uncertainty comparison for all measurements."""
        if self.profiled_set is None:
            raise RuntimeError("No profiled set available; call run() first.")
        report_pdf_comparison(
            self.pdf_set, self.pdf_members,
            self.profiled_set, self.profiled_members,
            self.measurements, self
        )

def main():
    args = parse_arguments()
    measurements = parse_measurements(args)

    print("--- Milestone 1 CLI Arguments Parsed Successfully ---")
    print(f"PDF Set: {args.pdf}")
    print(f"Measurements Parsed: {len(measurements)}")

    print("\n--- Milestone 2 Core Computational Functions ---")
    setup_lhapdf_path()

    # Detect PDF error type and convert MC replicas if needed
    original_pdf_name = args.pdf
    temp_dir = None

    error_type = detect_pdf_error_type(args.pdf)
    print(f"Detected ErrorType for '{args.pdf}': {error_type}")

    if error_type in ('replicas', 'mc') or args.force_mc2h:
        if error_type not in ('replicas', 'mc'):
            print(f"Note: --force-mc2h specified but ErrorType is '{error_type}'.")

        if args.mc2h_output_dir:
            output_dir = os.path.abspath(args.mc2h_output_dir)
            os.makedirs(output_dir, exist_ok=True)
        else:
            temp_dir = tempfile.mkdtemp(prefix='epump_mc2h_')
            atexit.register(lambda d=temp_dir: shutil.rmtree(d, ignore_errors=True))
            output_dir = temp_dir

        print(f"Converting MC replica set '{args.pdf}' to Hessian form...")
        converted_name, converted_path = convert_mc_to_hessian(
            pdf_name=args.pdf,
            neig=args.mc2h_neig,
            Q=args.mc2h_Q,
            epsilon=args.mc2h_epsilon,
            output_dir=output_dir,
            max_nf=args.mc2h_max_nf
        )

        # Register converted set directory so LHAPDF can find it
        setup_lhapdf_path(custom_path=output_dir)

        if args.closure_test:
            run_closure_test(original_pdf_name, converted_name, measurements, args)

        args.pdf = converted_name
        print(f"Using converted Hessian set '{converted_name}' for profiling.")

    elif error_type in ('hessian', 'symmhessian'):
        print(f"PDF set '{args.pdf}' is Hessian-type. No conversion needed.")
    else:
        print(f"Warning: Could not determine ErrorType for '{args.pdf}'. Proceeding without conversion.")

    print(f"Loading LHAPDF set '{args.pdf}'...")
    try:
        pdf_set = lhapdf.getPDFSet(args.pdf)
        pdf_members = pdf_set.mkPDFs()
        print(f"Successfully loaded '{args.pdf}' set with {len(pdf_members)} members.")
    except Exception as e:
        print(f"Error loading LHAPDF set '{args.pdf}': {e}", file=sys.stderr)
        sys.exit(1)

    print("\n--- Milestone 3 Dynamic File Generation ---")
    
    in_filepath = f"{args.name}.in"
    data_filepath = f"{args.name}.data"
    theory_filepath = f"{args.name}.theory"
    
    n_ev_pairs = (len(pdf_members) - 1) // 2
    n_obs = len(measurements)
    
    print(f"Generating .in file: {in_filepath}")
    generate_in_file(in_filepath, args.name, n_ev_pairs, n_obs, args.pdf)
    
    print(f"Generating .data file: {data_filepath}")
    generate_data_file(data_filepath, measurements)
    
    print(f"Generating .theory file: {theory_filepath}...")
    try:
        generate_theory_file(theory_filepath, pdf_members, measurements, args)
        print("Dynamic file generation completed successfully.")
    except Exception as e:
        print(f"Error generating theory file: {e}", file=sys.stderr)
        sys.exit(1)

    print("\n--- Milestone 4 ePump Invocation and Output Parsing ---")
    epump_abs_path = os.path.abspath(args.epump_path)
    if not os.path.exists(epump_abs_path):
        epump_abs_path = os.path.abspath(os.path.join(SCRIPT_DIR, args.epump_path))
        
    run_epump(epump_abs_path, args.name)

    print("\n--- Milestone 5 Comparison and Reporting ---")
    print(f"Loading profiled PDF set '{args.name}'...")
    try:
        updated_set = lhapdf.getPDFSet(args.name)
        updated_members = updated_set.mkPDFs()
        print(f"Successfully loaded profiled set '{args.name}' with {len(updated_members)} members.")
    except Exception as e:
        print(f"Error loading profiled PDF set '{args.name}': {e}", file=sys.stderr)
        sys.exit(1)
        
    report_pdf_comparison(pdf_set, pdf_members, updated_set, updated_members, measurements, args)
    
    print("\n============================================================")
    print("                    POSTERIOR PDF STORAGE SUMMARY            ")
    print("============================================================")
    print(f"  Posterior PDFs successfully generated!")
    print(f"  Study Label (Input Label): {args.name}")
    print(f"  Saved LHAPDF Directory:    {os.path.abspath(args.name)}")
    print("============================================================\n")

if __name__ == '__main__':
    main()
