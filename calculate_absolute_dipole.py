# calculate_absolute_dipole.py
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
calculate_absolute_dipole.py

Defines a PyMOL command to calculate the absolute dipole moment of a selected protein
structure. The dipole vector is calculated using partial charges and atomic coordinates,
and optionally visualized using a CGO arrow (via cgo_arrow.py).

Usage in PyMOL:
    PyMOL> load your_file.pqr
    PyMOL> calculate_absolute_dipole your_object_name
"""

import os
import numpy as np
from pymol import cmd

def calculate_absolute_dipole(selection="all", vector_scale=0.3, cgoarrow_script_path=None):
    """
    Calculates and visualizes the dipole moment of a structure in PyMOL.

    Args:
        selection (str): Name of the object/selection in PyMOL.
        vector_scale (float): Scale factor for the dipole vector visualization.
        cgoarrow_script_path (str): Path to the cgo_arrow.py script (optional).
    """
    if not selection:
        print("No structure selected.")
        return

    cmd.center(selection)
    atoms = cmd.get_model(selection)
    com = np.array(cmd.centerofmass(selection))  # center of mass

    dipole = np.zeros(3)
    total_charge = 0

    for at in atoms.atom:
        coord = np.array(at.coord)
        if hasattr(at, "partial_charge"):
            dipole += at.partial_charge * coord
            total_charge += at.partial_charge
        else:
            print(f"Warning: No charge assigned to atom {at.name} in {selection}, skipping...")

    if total_charge == 0:
        print("Error: No charges found. Dipole moment cannot be computed.")
        return

    ELEMENTARY_CHARGE_C = 1.602176634e-19  # Coulombs
    total_charge_coulombs = total_charge * ELEMENTARY_CHARGE_C
    dipole_magnitude = np.linalg.norm(4.803 * dipole)  # Debye conversion

    print(f"Dipole moment for {selection}: {dipole_magnitude:.3f} Debye")

    # Write output
    output_file = f"{selection}_dipole.txt"
    with open(output_file, "w") as f:
        f.write(f"Selection: {selection}\n")
        f.write(f"Total Charge (e): {total_charge:.6f}\n")
        f.write(f"Total Charge (Coulombs): {total_charge_coulombs:.6e} C\n")
        f.write(f"Dipole Moment: {dipole_magnitude:.3f} Debye\n")
        f.write(f"Dipole Vector: {dipole.tolist()}\n")
        f.write(f"Center of Mass: {com.tolist()}\n")
    print(f"Dipole information saved to {output_file}")

    # Vector for arrow visualization
    v = com + (vector_scale * dipole)

    color_start = "blue"
    color_end = "red"

    # Find cgo_arrow.py if not provided
    if not cgoarrow_script_path:
        home_dir = os.path.expanduser("~")
        cgoarrow_script_path = os.path.join(home_dir, "pymol_scripts", "cgo_arrow.py")

    if not os.path.isfile(cgoarrow_script_path):
        print(f"Warning: Could not find cgo_arrow.py at {cgoarrow_script_path}. Dipole vector visualization skipped.")
        return
    else:
        cmd.do(f"run {cgoarrow_script_path}")

    cgoname = f"dipole-{selection}"
    cmd.delete(cgoname)
    cmd.do(f"cgo_arrow [{com[0]}, {com[1]}, {com[2]}], [{v[0]}, {v[1]}, {v[2]}], gap=1.0, color={color_start} {color_end}, name={cgoname}")
    cmd.enable(cgoname)

# Register command
cmd.extend("calculate_absolute_dipole", calculate_absolute_dipole)
