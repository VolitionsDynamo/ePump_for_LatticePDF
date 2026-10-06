"""
Phase B: Profile a single observable on the pre-converted Hessian set.

Run phase_a_convert.py first to produce the converted set.

Usage:
    conda run -n apfelpp python tests/phase_b_profile.py
"""
import sys, time, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from e_profiler import EProfiler, setup_lhapdf_path
import lhapdf

LHAPDF_BASE  = '/home/daniel/miniconda3/share/LHAPDF'
HESSIAN_DIR  = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nnpdf31_hessian')
HESSIAN_NAME = 'NNPDF31_nnlo_as_0118_hessian_20'
EPUMP_PATH   = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            '..', 'ePump_kp20221218', 'src', 'UpdatePDFs')
RUN_DIR      = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'phase_b_run')

setup_lhapdf_path(LHAPDF_BASE)
lhapdf.setPaths([HESSIAN_DIR] + lhapdf.paths())

print(f"Loading '{HESSIAN_NAME}' ...")
ep = EProfiler(
    HESSIAN_NAME,
    os.path.join(RUN_DIR, 'phase_b'),
    epump_path=os.path.abspath(EPUMP_PATH),
    lhapdf_path=HESSIAN_DIR,
)
print(f"PDF members: {len(ep.pdf_members)}  (expect 41)")

ep.add_measurement(
    x=0.30, Q2=4.0, value=0.10, stat=0.01,
    obs_type='moment', flavor='u-d',
    xmin=0.25, xmax=0.35, nx=50, weight='gaussian', moment=1,
)

ep.generate_files()

print("Running ePump ...")
t0 = time.time()
ep.run()
elapsed = time.time() - t0
print(f"ePump done in {elapsed:.1f}s")

ep.report()
print("SUCCESS")

