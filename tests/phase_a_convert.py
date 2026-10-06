"""
Phase A: MC-to-Hessian conversion only — timed.

Tests that NNPDF31_nnlo_as_0118 can be converted and the resulting
Hessian set loaded correctly. Run this before phase_b_profile.py.

Usage:
    conda run -n apfelpp python tests/phase_a_convert.py
"""
import sys, time, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from e_profiler import convert_mc_to_hessian, setup_lhapdf_path
import lhapdf

LHAPDF_BASE = '/home/daniel/miniconda3/share/LHAPDF'
OUTPUT_DIR  = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nnpdf31_hessian')
NEIG        = 20

setup_lhapdf_path(LHAPDF_BASE)

print(f"Converting NNPDF31_nnlo_as_0118 → Hessian (neig={NEIG}) ...")
t0 = time.time()
hessian_name, hessian_dir = convert_mc_to_hessian(
    'NNPDF31_nnlo_as_0118',
    neig=NEIG,
    Q=1.0,
    epsilon=1000.0,
    output_dir=OUTPUT_DIR,
    max_nf=3,
)
elapsed = time.time() - t0
print(f"Conversion done in {elapsed:.1f}s  →  '{hessian_name}'")
print(f"Output directory: {hessian_dir}")

lhapdf.setPaths([OUTPUT_DIR] + lhapdf.paths())
s = lhapdf.getPDFSet(hessian_name)
members = s.mkPDFs()
print(f"Loaded: {len(members)} members  (expect {2*NEIG+1} = {2*NEIG+1})")

