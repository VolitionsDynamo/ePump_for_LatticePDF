"""Phase C: verify EProfiler auto-converts MC replica sets in __init__."""
import sys, time
sys.path.insert(0, '/home/daniel/ePump/ePump_for_LatticePDF')
from e_profiler import EProfiler, setup_lhapdf_path, detect_pdf_error_type
import lhapdf

LHAPDF_BASE = '/home/daniel/miniconda3/share/LHAPDF'
setup_lhapdf_path(LHAPDF_BASE)

MC_PDF = 'NNPDF31_nnlo_as_0118'
RUN_NAME = '/home/daniel/ePump/ePump_for_LatticePDF/tests/phase_c_run/phase_c'

print(f"Original error type: {detect_pdf_error_type(MC_PDF)}")

t0 = time.time()
ep = EProfiler(MC_PDF, RUN_NAME, mc2h_neig=20,
               mc2h_output_dir='/home/daniel/ePump/ePump_for_LatticePDF/tests/nnpdf31_hessian')
t1 = time.time()

print(f"After __init__ ({t1-t0:.1f}s):")
print(f"  pdf_set_name = '{ep.pdf_set_name}'  (expect 'NNPDF31_nnlo_as_0118_hessian_20')")
print(f"  pdf_members  = {len(ep.pdf_members)}  (expect 41)")
assert ep.pdf_set_name == 'NNPDF31_nnlo_as_0118_hessian_20', "Wrong set name after auto-convert"
assert len(ep.pdf_members) == 41, f"Expected 41 members, got {len(ep.pdf_members)}"

ep.add_measurement(
    x=0.30, Q2=4.0, value=0.10, stat=0.01,
    obs_type='moment', flavor='u-d',
    xmin=0.25, xmax=0.35, nx=50, weight='gaussian', moment=1,
)

ep.generate_files()
ep.run()
ep.report()
print(f"\nSUCCESS  (total {time.time()-t0:.1f}s)")

