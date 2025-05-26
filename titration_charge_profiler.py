# titration_charge_profiler.py
# 
# CODDM-CPs is a repository of scripts that can be used to study the
# coordinate-origin dependence of dipole moments of charged proteins.
# 
# Copyright (C) 2025  Islam K. Matar
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
titration_charge_profiler.py

This script computes the pH-dependent partial charge profiles of one or more protein structures
(PDB files), using pdb2pqr and PROPKA. It identifies a common or representative formal charge
across proteins, and outputs .pqr files, titration plots, and a CSV summary.

Usage:
    python titration_charge_profiler.py *.pdb --step 0.5 --csv summary.csv
"""

import os
import subprocess
import tempfile
import matplotlib.pyplot as plt
import csv
import shutil
import glob
import argparse

def run_pdb2pqr(input_pdb, pH, output_pqr):
    """
    Runs pdb2pqr to generate a PQR file for a given PDB structure at a specified pH.
    
    Args:
        input_pdb (str): Path to input PDB file.
        pH (float): Target pH for protonation state prediction.
        output_pqr (str): Output filename for the resulting PQR file.
    
    Returns:
        bool: True if pdb2pqr runs successfully, False otherwise.
    """
    try:
        subprocess.run([
            "python", "-m", "pdb2pqr",
            "--ff", "AMBER",
            "--titration-state-method", "propka",
            "--with-ph", str(pH),
            input_pdb, output_pqr
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        print(f"Error at pH {pH} for {input_pdb}")
        return False

def calculate_total_charge(pqr_file):
    """
    Calculates the total partial charge from a PQR file.
    
    Args:
        pqr_file (str): Path to a PQR file.
    
    Returns:
        float: Total charge summed across atoms.
    """
    total_charge = 0.0
    with open(pqr_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    charge = float(line[54:62])
                    total_charge += charge
                except ValueError:
                    continue
    return total_charge

def get_charge_profile(pdb_file, pH_min, pH_max, step):
    """
    Generates a charge profile over a pH range for a single PDB file.
    
    Returns:
        dict: Mapping of pH values to total charges.
    """
    charges = {}
    for i in range(int((pH_max - pH_min) / step) + 1):
        pH = round(pH_min + i * step, 2)
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pqr") as tmpfile:
            pqr_file = tmpfile.name
        if run_pdb2pqr(pdb_file, pH, pqr_file):
            charge = calculate_total_charge(pqr_file)
            charges[pH] = charge
        os.remove(pqr_file)
    return charges

def find_common_formal_charge(charge_profiles):
    """
    Finds the most extreme integer formal charge shared by all charge profiles.
    
    Returns:
        int or None: Most extreme shared formal charge, or None if none exists.
    """
    sets = []
    for charges in charge_profiles.values():
        rounded = set(round(c) for c in charges.values())
        sets.append(rounded)
    common = set.intersection(*sets)
    if not common:
        return None
    return sorted(common, key=lambda x: -abs(x))[0]

def save_best_pqr(pdb_file, pH, output_name):
    """
    Saves the best matching PQR file at a given pH.
    """
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pqr") as tmpfile:
        pqr_file = tmpfile.name
    if run_pdb2pqr(pdb_file, pH, pqr_file):
        shutil.move(pqr_file, output_name)
        print(f"Saved {output_name} at pH {pH}")
    else:
        os.remove(pqr_file)

def plot_charge_curve(pdb_file, charges, target_charge):
    """
    Saves a titration curve plot of total charge vs. pH for a given protein.
    """
    pHs = sorted(charges.keys())
    vals = [charges[pH] for pH in pHs]
    plt.figure()
    plt.plot(pHs, vals, marker='o', label=os.path.basename(pdb_file))
    plt.axhline(target_charge, linestyle='--', color='red', label=f'Target Charge = {target_charge}')
    plt.title(f'Titration Curve: {os.path.basename(pdb_file)}')
    plt.xlabel('pH')
    plt.ylabel('Total Partial Charge')
    plt.legend()
    plt.grid(True)
    out_plot = os.path.splitext(pdb_file)[0] + "_curve.png"
    plt.savefig(out_plot, dpi=300)
    plt.close()

def save_csv_summary(profiles, filename):
    """
    Saves a CSV summary of charge profiles across all proteins.
    """
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        all_pHs = sorted(next(iter(profiles.values())).keys())
        writer.writerow(["pH"] + [os.path.basename(f) for f in profiles])
        for pH in all_pHs:
            row = [pH]
            for f in profiles:
                row.append(round(profiles[f].get(pH, 0.0), 3))
            writer.writerow(row)
    print(f"CSV summary saved to {filename}")

# ---- Entry point ----
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find common extreme formal charge for multiple proteins.")
    parser.add_argument("pdbs", nargs="+", help="List of PDB files (wildcards supported, e.g. *.pdb)")
    parser.add_argument("--step", type=float, default=0.5, help="pH step size (default 0.5)")
    parser.add_argument("--csv", default="multi_charge_summary.csv", help="CSV output")
    args = parser.parse_args()

    # Expand wildcards
    expanded_pdbs = []
    for pattern in args.pdbs:
        expanded_pdbs.extend(glob.glob(pattern))
    args.pdbs = list(set(expanded_pdbs))
    if not args.pdbs:
        print("No PDB files matched the given pattern.\nUsage: python titration_charge_profiler.py *.pdb")
        exit(1)

    # Step 1: Build charge profiles
    charge_profiles = {}
    for pdb in args.pdbs:
        print(f"Scanning {pdb} ...")
        charges = get_charge_profile(pdb, 0.0, 14.0, args.step)
        charge_profiles[pdb] = charges

    # Step 2: Find common or partial matching charge
    common_charge = find_common_formal_charge(charge_profiles)
    matched_pdbs = []

    if common_charge is None:
        print("No single common charge found across all proteins.")
        charge_to_pdbs = {}
        for pdb, profile in charge_profiles.items():
            rounded = set(round(c) for c in profile.values())
            for c in rounded:
                charge_to_pdbs.setdefault(c, set()).add(pdb)
        sorted_charges = sorted(charge_to_pdbs.items(), key=lambda x: (-len(x[1]), -abs(x[0])))
        best_charge, matched_pdbs = sorted_charges[0]
        print(f"Selected charge {best_charge} matches {len(matched_pdbs)}/{len(charge_profiles)} proteins.")
    else:
        best_charge = common_charge
        matched_pdbs = list(charge_profiles.keys())
        print(f"Most extreme common charge found: {common_charge}")

    # Step 3: Generate output
    for pdb in matched_pdbs:
        best_pH = min(charge_profiles[pdb], key=lambda p: abs(round(charge_profiles[pdb][p]) - best_charge))
        output_name = os.path.splitext(pdb)[0] + ".pqr"
        save_best_pqr(pdb, best_pH, output_name)
        plot_charge_curve(pdb, charge_profiles[pdb], best_charge)

    save_csv_summary(charge_profiles, args.csv)
