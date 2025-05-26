# calculate_charge_at_pH7.py
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
calculate_charge_at_pH7.py

This script processes all PDB files in the current directory, computes their total
partial charge using pdb2pqr at pH 7.0, and outputs a CSV summary file.

Usage:
    python calculate_charge_at_pH7.py
"""

import os
import subprocess
import csv
import glob

def run_pdb2pqr(input_pdb, pH, output_pqr):
    """
    Runs pdb2pqr to generate a PQR file for a given PDB structure at the specified pH.

    Args:
        input_pdb (str): Input PDB file path.
        pH (float): Desired pH value.
        output_pqr (str): Output PQR file path.

    Returns:
        bool: True if conversion succeeds, False otherwise.
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
        print(f"Error processing {input_pdb} at pH {pH}")
        return False

def calculate_total_charge(pqr_file):
    """
    Calculates total partial charge from a PQR file.

    Args:
        pqr_file (str): Path to a PQR file.

    Returns:
        float: Total partial charge.
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

def process_all_pdbs(ph=7.0, output_csv="charges_at_pH7.csv"):
    """
    Processes all PDB files in the current directory at the specified pH.

    Args:
        ph (float): Target pH value.
        output_csv (str): Output CSV filename.
    """
    pdb_files = glob.glob("*.pdb")
    results = []

    for pdb_file in pdb_files:
        base_id = os.path.splitext(os.path.basename(pdb_file))[0]
        output_pqr = base_id + f"_pH{int(ph)}.pqr"

        if run_pdb2pqr(pdb_file, ph, output_pqr):
            charge = calculate_total_charge(output_pqr)
            results.append((base_id, round(charge, 3)))
        else:
            results.append((base_id, "Error"))

    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["PDB ID", f"Charge at pH {ph}"])
        writer.writerows(results)

    print(f"Charge report saved to {output_csv}")

# Run as script
if __name__ == "__main__":
    process_all_pdbs()
