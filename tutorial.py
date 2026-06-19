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

# Import all core functions from e_profiler
import e_profiler
from e_profiler import (
    setup_lhapdf_path,
    parse_flavor_expression,
    evaluate_pdf_combination,
    compute_integrated_moment,
    generate_in_file,
    generate_data_file,
    generate_theory_file,
    run_epump,
    report_pdf_comparison
)

print("--- Step 1: Environment Setup ---")
setup_lhapdf_path()
print("LHAPDF Search Paths:", lhapdf.paths())

# ==========================================
# CELL 2: Loading the Starting PDF Set
# ==========================================
print("\n--- Step 2: Loading Starting PDF Set ---")
pdf_set_name = "CT18NNLO"
pdf_set = lhapdf.getPDFSet(pdf_set_name)
pdf_members = pdf_set.mkPDFs()
print(f"Loaded '{pdf_set_name}' with {len(pdf_members)} members successfully.")

# ==========================================
# CELL 3: Parsing Flavor Expressions & Point Evaluation
# ==========================================
print("\n--- Step 3: Parsing Flavor Expressions & Point Evaluation ---")
# Parse flavor expression u - d (flavor name -> PID conversion)
expression = "u - d"
parsed_terms = parse_flavor_expression(expression)
print(f"Parsed Terms for '{expression}':", parsed_terms) # Expected: [(1.0, 2), (-1.0, 1)]

# Evaluate the flavor expression for the central member (member 0) at x=0.5, Q^2=4.0
x_val = 0.5
Q2_val = 4.0
val_point = evaluate_pdf_combination(pdf_members[0], parsed_terms, x_val, Q2_val)
print(f"Value at x={x_val}, Q2={Q2_val}: {val_point:.6f}")

# ==========================================
# CELL 4: Computing Moments and Specific Weight Functions
# ==========================================
print("\n--- Step 4: Integrated Moments & Specific Weight Functions ---")
xmin, xmax, nx = 0.1, 0.7, 100

# Compute integrated moment with default weight ('1') and moment n=0
val_moment_0 = compute_integrated_moment(
    pdf_member=pdf_members[0],
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
    pdf_member=pdf_members[0],
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
    pdf_member=pdf_members[0],
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
# CELL 5: Defining and Parsing Measurements
# ==========================================
print("\n--- Step 5: Defining and Parsing Measurements ---")
# Define pseudodata measurements in dict format compatible with the generators.
# Format: [x, Q2, value, statistical_error, uncor_sys, cor_sys...]
raw_measurement = [0.5, 2.0, 0.151, 0.040, 0.0]

measurements = [{
    'x': raw_measurement[0],
    'Q2': raw_measurement[1],
    'value': raw_measurement[2],
    'stat': raw_measurement[3],
    'uncor_sys': raw_measurement[4],
    'cor_sys': [],
    'obs_type': 'moment',
    'flavor': 'u-d',
    'xmin': xmin,
    'xmax': xmax,
    'nx': nx,
    'weight': '1',
    'moment': 0
}]
print("Measurement Configuration:", measurements[0])

# ==========================================
# CELL 6: Dynamic File Generation
# ==========================================
print("\n--- Step 6: Dynamic ePump File Generation ---")
run_name = "tutorial_demo"

in_file = f"{run_name}.in"
data_file = f"{run_name}.data"
theory_file = f"{run_name}.theory"

n_ev_pairs = (len(pdf_members) - 1) // 2
n_obs = len(measurements)

print(f"Generating control file: {in_file}")
generate_in_file(in_file, run_name, n_ev_pairs, n_obs, pdf_set_name)

print(f"Generating data file: {data_file}")
generate_data_file(data_file, measurements)

print(f"Generating theory file: {theory_file}")
# Create dummy namespace class to emulate parsed arguments for theory generator
class DummyArgs:
    obs_type = 'moment'
    xmin = xmin
    xmax = xmax
    nx = nx
    Q2 = 2.0
    flavor = 'u-d'
    weight = '1'
    moment = 0

generate_theory_file(theory_file, pdf_members, measurements, DummyArgs())
print("All files generated successfully.")

# ==========================================
# CELL 7: Running ePump Profiling
# ==========================================
print("\n--- Step 7: Executing ePump Profiling ---")
# Path to compiled UpdatePDFs binary
epump_bin_path = "/home/daniel/ePump/ePump_for_LatticePDF/ePump_kp20221218/src/UpdatePDFs"
run_epump(epump_bin_path, run_name)

# ==========================================
# CELL 8: Loading Profiled PDFs & Reporting Comparison
# ==========================================
print("\n--- Step 8: PDF Uncertainty Comparison & Reporting ---")
# Load profiled PDF set that was generated during step 7
profiled_set = lhapdf.getPDFSet(run_name)
profiled_members = profiled_set.mkPDFs()

# Run full comparative report!
report_pdf_comparison(
    original_set=pdf_set,
    original_members=pdf_members,
    updated_set=profiled_set,
    updated_members=profiled_members,
    measurements=measurements,
    args=DummyArgs()
)

print("\nTutorial run complete. Everything works flawlessly! :)")
