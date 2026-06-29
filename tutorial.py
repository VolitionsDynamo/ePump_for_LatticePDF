# # e_profiler Tutorial and Debugging Guide
#
# This tutorial script steps through every single function in `e_profiler.py` sequentially.
# You can run this directly in python, or easily convert/paste each section into a Jupyter Notebook cell.

# ==========================================
# CELL 1: Imports and Environment Setup
# ==========================================
import sys
import os
import numpy as np
import lhapdf

# Add e_profiler directory to python search path
sys.path.append("/home/daniel/ePump/ePump_for_LatticePDF")

# Import the EProfiler class and any low-level functions used in later cells
import e_profiler
from e_profiler import (
    EProfiler,
    parse_flavor_expression,
    evaluate_pdf_combination,
    compute_integrated_moment,
    run_closure_test,
    detect_pdf_error_type,
)

# ==========================================
# CELL 2: Loading the Starting PDF Set
# ==========================================
print("\n--- Step 2: Loading Starting PDF Set ---")
# EProfiler calls setup_lhapdf_path() and lhapdf.getPDFSet() + mkPDFs() internally,
# so the PDF set is ready as soon as the object is created.
profiler = EProfiler("CT18NNLO", "tutorial_demo")
print(f"Loaded '{profiler.pdf_set_name}' with {len(profiler.pdf_members)} members successfully.")

# ==========================================
# CELL 3: MC-to-Hessian Detection and Conversion
# ==========================================
print("\n--- Step 3: MC-to-Hessian Detection and Conversion ---")

# Check the current set (CT18NNLO is already Hessian)
error_type = detect_pdf_error_type(profiler.pdf_set_name)
print(f"ErrorType for '{profiler.pdf_set_name}': {error_type}")

# Demonstrate conversion with an MC replica set
mc_pdf_name = "NNPDF23_lo_as_0130_qed"
mc_error_type = detect_pdf_error_type(mc_pdf_name)
print(f"ErrorType for '{mc_pdf_name}': {mc_error_type}")

mc_profiler = EProfiler(mc_pdf_name, "tutorial_mc_demo")
mc_profiler.convert_to_hessian(neig=50, Q=1.0, epsilon=1000.0, max_nf=3)
print(f"Converted set: '{mc_profiler.pdf_set_name}' with {len(mc_profiler.pdf_members)} members (expected {2*50+1}=101)")

# ==========================================
# CELL 4: MC-to-Hessian Closure Test
# ==========================================
print("\n--- Step 4: Closure Test ---")

xmin, xmax, nx = 0.1, 0.7, 100
measurements_for_test = [{
    'x': 0.5, 'Q2': 4.0, 'value': 0.151, 'stat': 0.040,
    'uncor_sys': 0.0, 'cor_sys': [],
    'obs_type': 'moment', 'flavor': 'u-d',
    'xmin': xmin, 'xmax': xmax, 'nx': nx,
    'weight': '1', 'moment': 0
}]

# run_closure_test accepts args but never reads from it, so None is safe
run_closure_test(mc_pdf_name, mc_profiler.pdf_set_name, measurements_for_test, None)

# ==========================================
# CELL 5: Parsing Flavor Expressions & Point Evaluation
# ==========================================
print("\n--- Step 5: Parsing Flavor Expressions & Point Evaluation ---")
# Parse flavor expression u - d (flavor name -> PID conversion)
expression = "u - d"
parsed_terms = parse_flavor_expression(expression)
print(f"Parsed Terms for '{expression}':", parsed_terms) # Expected: [(1.0, 2), (-1.0, 1)]

# Evaluate the flavor expression for the central member (member 0) at x=0.5, Q^2=4.0
x_val = 0.5
Q2_val = 4.0
val_point = evaluate_pdf_combination(profiler.pdf_members[0], parsed_terms, x_val, Q2_val)
print(f"Value at x={x_val}, Q2={Q2_val}: {val_point:.6f}")

# ==========================================
# CELL 6: Computing Moments and Specific Weight Functions
# ==========================================
print("\n--- Step 6: Integrated Moments & Specific Weight Functions ---")

# Compute integrated moment with default weight ('1') and moment n=0
val_moment_0 = compute_integrated_moment(
    pdf_member=profiler.pdf_members[0],
    parsed_terms=parsed_terms,
    xmin=xmin,
    xmax=xmax,
    nx=nx,
    Q2=Q2_val,
    weight_type='1',
    moment=0
)
print(f"Integrated Moment (weight='1', n=0) on [{xmin}, {xmax}]: {val_moment_0:.6f}")

# Compute integrated moment with default weight ('1') and moment n=1
val_moment_1 = compute_integrated_moment(
    pdf_member=profiler.pdf_members[0],
    parsed_terms=parsed_terms,
    xmin=xmin,
    xmax=xmax,
    nx=nx,
    Q2=Q2_val,
    weight_type='1',
    moment=1
)
print(f"Integrated Moment (weight='1', n=1) on [{xmin}, {xmax}]: {val_moment_1:.6f}")

# Compute integrated moment with Gaussian weight and moment n=1
val_gaussian = compute_integrated_moment(
    pdf_member=profiler.pdf_members[0],
    parsed_terms=parsed_terms,
    xmin=xmin,
    xmax=xmax,
    nx=nx,
    Q2=Q2_val,
    weight_type='gaussian',
    moment=1
)
print(f"Integrated Moment (weight='gaussian', n=1) on [{xmin}, {xmax}]: {val_gaussian:.6f}")

# ==========================================
# CELL 7: Defining and Adding Measurements
# ==========================================
print("\n--- Step 7: Defining and Adding Measurements ---")
# add_measurement() builds the dict format expected by all downstream generators.
profiler.add_measurement(
    x=0.5, Q2=2.0, value=0.151, stat=0.040,
    obs_type='moment', flavor='u-d',
    xmin=xmin, xmax=xmax, nx=nx, weight='1', moment=0
)
print("Measurement Configuration:", profiler.measurements[0])

# ==========================================
# CELL 8: Dynamic File Generation
# ==========================================
print("\n--- Step 8: Dynamic ePump File Generation ---")
# generate_files() writes the .in, .data, and .theory files in one call.
# No DummyArgs needed — the profiler holds all configuration.
profiler.generate_files()
print("All files generated successfully.")

# ==========================================
# CELL 9: Running ePump Profiling
# ==========================================
print("\n--- Step 9: Executing ePump Profiling ---")
# run() invokes the UpdatePDFs binary and loads the resulting profiled PDF set.
profiler.run()

# ==========================================
# CELL 10: Loading Profiled PDFs & Reporting Comparison
# ==========================================
print("\n--- Step 10: PDF Uncertainty Comparison & Reporting ---")
# report() compares original vs profiled uncertainties for every measurement.
profiler.report()

print("\nTutorial run complete.")
