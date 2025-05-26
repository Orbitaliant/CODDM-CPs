# dipole_origin_dependence.py
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
dipole_origin_dependence.py

Analyzes coordinate-origin dependence of dipole moments. For each input PQR file:
- Computes the dipole moment from the original structure's center.
- Computes the dipole moment again after aligning to a shifted origin.
- Saves session files and a CSV report of the dipole vectors and magnitudes.

Requires:
- calculate_absolute_dipole.py (PyMOL dipole calculation command)
- centroid.py (utility function for geometric centroid)
"""

import sys
import os
import glob
import csv
import numpy as np
import pymol
from pymol import cmd

# Load local utilities
sys.path.append(r"C:\Users\islam\pymol_scripts")
from centroid import centroid
from calculate_absolute_dipole import calculate_absolute_dipole

def extract_dipole_and_size(selection):
    """
    Computes the dipole vector, magnitude, and charge-weighted radius of gyration for a selection.
    
    Args:
        selection (str): PyMOL object name.
    
    Returns:
        tuple: (dipole_vector_in_debye, dipole_magnitude, charge-weighted Rg)
    """
    atoms = cmd.get_model(selection)
    coords, charges = [], []

    for at in atoms.atom:
        if hasattr(at, "partial_charge"):
            coords.append(at.coord)
            charges.append(at.partial_charge)

    coords = np.array(coords)
    charges = np.array(charges)
    weights = np.abs(charges)

    W_total = np.sum(weights)
    if np.isclose(W_total, 0.0):
        raise ValueError("Sum of absolute charges is zero; Rg undefined.")

    r_cm = np.sum(coords.T * weights, axis=1) / W_total
    dipole = np.sum(charges[:, None] * coords, axis=0)
    squared_distances = np.sum((coords - r_cm) ** 2, axis=1)
    Rg_squared = np.sum(weights * squared_distances) / W_total
    Rg = np.sqrt(Rg_squared)

    dipole_debye = 4.803 * dipole
    dipole_magnitude = np.linalg.norm(dipole_debye)
    return dipole_debye, dipole_magnitude, Rg

# Initialize PyMOL in quiet mode
pymol.finish_launching(['pymol', '-qc'])

pqr_files = glob.glob("*.pqr")
if not pqr_files:
    print("No PQR files found.")
    cmd.quit()
    sys.exit()

summary = []

for pqr in pqr_files:
    base = os.path.splitext(os.path.basename(pqr))[0]
    cmd.load(pqr)
    obj = base

    # ---- First Origin: Center on original centroid ----
    cmd.pseudoatom(obj, name='CRC', resi='1', resn='CRC', chain='ZZ', color='hotpink', pos=centroid(obj))
    cmd.pseudoatom('1st_origin', name='CRC', resi='1', resn='CRC', chain='ZZ', pos='[0,0,0]')
    cmd.align(obj, '1st_origin')
    cmd.remove('resn CRC')
    cmd.delete('1st_origin')

    orig_name = f"{base}_orig"
    cmd.set_name(obj, orig_name)
    calculate_absolute_dipole(orig_name)
    dipole_orig_vec, dipole_orig_mag, size = extract_dipole_and_size(orig_name)
    cmd.save(f"{orig_name}_dipole.pse")

    # ---- Second Origin: Artificial shift to [20, 20, 20] ----
    cmd.delete(f"dipole-{orig_name}")
    cmd.set_name(orig_name, obj)

    cmd.pseudoatom(obj, name='CRC', resi='1', resn='CRC', chain='ZZ', pos=centroid(obj))
    cmd.pseudoatom('2nd_origin', name='CRC', resi='1', resn='CRC', chain='ZZ', pos='[20,20,20]')
    cmd.align(obj, '2nd_origin')
    cmd.remove('resn CRC')
    cmd.delete('2nd_origin')

    mod_name = f"{base}_mod"
    cmd.set_name(obj, mod_name)
    calculate_absolute_dipole(mod_name)
    dipole_mod_vec, dipole_mod_mag, _ = extract_dipole_and_size(mod_name)
    cmd.save(f"{mod_name}_dipole.pse")

    # ---- Summary Row ----
    delta_dipole = abs(dipole_mod_mag - dipole_orig_mag)

    summary.append({
        'Name': base,
        'Dipole_Orig (D)': dipole_orig_mag,
        'Dipole_Mod (D)': dipole_mod_mag,
        'Delta_Dipole (D)': delta_dipole,
        'Charge-weighted Rg (Å)': size,
        'Dipole_Orig_Vector': dipole_orig_vec.tolist(),
        'Dipole_Mod_Vector': dipole_mod_vec.tolist()
    })

# ---- Write CSV Summary ----
with open('dipole_summary.csv', 'w', newline='') as csvfile:
    fieldnames = summary[0].keys()
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in summary:
        writer.writerow(row)

print("Dipole analysis complete. Results saved to dipole_summary.csv")
cmd.quit()
