"""Phase D: verify EProfiler works end-to-end with JAMDiFF23-transversity_lo_nolat."""
import sys, time, warnings
sys.path.insert(0, '/home/daniel/ePump/ePump_for_LatticePDF')
from e_profiler import EProfiler, setup_lhapdf_path, detect_pdf_error_type

JAM_PDF = 'JAMDiFF23-transversity_lo_nolat'
JAM_DIR = '/home/daniel/ePump/ePump_for_LatticePDF'
RUN_NAME = '/home/daniel/ePump/ePump_for_LatticePDF/tests/phase_d_run/phase_d'
OUTPUT_DIR = '/home/daniel/ePump/ePump_for_LatticePDF/tests/jam_hessian'

setup_lhapdf_path(JAM_DIR)
print(f"Original error type: {detect_pdf_error_type(JAM_PDF)}")

# Track RuntimeWarnings during the whole run
t0 = time.time()
caught_warnings = []
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter('always')

    ep = EProfiler(JAM_PDF, RUN_NAME, mc2h_neig=20,
                   mc2h_output_dir=OUTPUT_DIR, mc2h_Q=1.14)

    caught_warnings = [x for x in w if issubclass(x.category, RuntimeWarning)
                       and 'divide' in str(x.message).lower()]

t1 = time.time()

print(f"\nAfter __init__ ({t1-t0:.1f}s):")
print(f"  pdf_set_name = '{ep.pdf_set_name}'")
print(f"  pdf_members  = {len(ep.pdf_members)}  (expect {2*20+1}=41)")
assert len(ep.pdf_members) == 41, f"Expected 41 members, got {len(ep.pdf_members)}"

if caught_warnings:
    print(f"  FAIL: unexpected RuntimeWarning(divide): {caught_warnings[0].message}")
    sys.exit(1)
else:
    print("  No divide RuntimeWarning emitted.")

ep.add_measurement(
    x=0.30, Q2=4.0, value=0.10, stat=0.01,
    obs_type='moment', flavor='u-d',
    xmin=0.25, xmax=0.35, nx=50, weight='gaussian', moment=1,
)

ep.generate_files()

# Verify path length in .in file
with open(f'{RUN_NAME}.in') as f:
    for line in f:
        if 'hessian' in line or './s/' in line:
            parts = line.strip().split()
            pdf_path = parts[0]
            print(f"  PDF path in .in: {repr(pdf_path)} (len={len(pdf_path)})")
            assert len(pdf_path) < 80, f"Path too long for ePump: {len(pdf_path)} chars"
            print("  Path length OK (<80 chars).")

ep.run()
ep.report()
print(f"\nSUCCESS  (total {time.time()-t0:.1f}s)")


