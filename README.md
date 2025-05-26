# CODDM-CPs: Coordinate-Origin Dependence of Dipole Moments in Charged Proteins

This repository contains Python and PyMOL-based tools for investigating how the choice of coordinate origin affects the calculated dipole moment of charged proteins.

The scripts enable:
- Automated **partial charge profiling** across pH values.
- **Dipole moment calculation** from atomic partial charges.
- Comparison of **dipole vectors across coordinate origins**.
- Consistent, reproducible workflows using PQR structures.

---

## 📂 Repository Structure

| Script | Description |
|--------|-------------|
| `titration_charge_profiler.py` | Scans multiple PDBs across a pH range (0–14), computes formal charge profiles, finds a representative charge state, and saves `.pqr` and titration plots. |
| `calculate_absolute_dipole.py` | PyMOL plugin to compute and visualize dipole moments using charge-weighted atomic positions. |
| `dipole_origin_dependence.py` | Analyzes how dipole magnitude changes when the structure is shifted to a new origin; outputs CSV comparisons and session files. |
| `calculate_charge_at_pH7.py` | Calculates and summarizes partial charges at pH 7 for all `.pdb` files in the current folder. |
| `centroid.py` | Utility script for calculating geometric centroids (required by `dipole_origin_dependence.py`). |
| `cgo_arrow.py` | Custom CGO arrow drawing function for PyMOL (required by `calculate_absolute_dipole.py`). |

---

## 📦 Dependencies

- Python 3.x
- [PyMOL](https://pymol.org/)
- [pdb2pqr](https://github.com/Electrostatics/apbs-pdb2pqr)
- `numpy`, `matplotlib`

Make sure `pdb2pqr` is installed and accessible via `python -m pdb2pqr ...`.

---

## 🧪 Usage Overview

### 1. Charge Profile Scanning

```bash
python titration_charge_profiler.py *.pdb --step 0.5 --csv output.csv
```

- Computes charge vs. pH curves (0–14) for each PDB file.
- Finds a common or best-match formal charge.
- Saves `.pqr` files and titration plots (`.png`).

---

### 2. Dipole Moment Calculation (in PyMOL)

```python
# In PyMOL
load your_file.pqr
calculate_dipole_absolute your_object
```

- Computes the dipole vector based on partial charges.
- Draws a colored CGO arrow (if `cgo_arrow.py` is available).
- Saves dipole information to a `.txt` file.

---

### 3. Dipole Origin Dependence Analysis

```bash
python dipole_origin_dependence.py
```

- Loads all `.pqr` files.
- Calculates dipoles twice: once from original centroid, once from a shifted origin.
- Outputs `dipole_summary.csv` and session files (`.pse`).

---

### 4. Quick Charge Check at pH 7

```bash
python calculate_charge_at_pH7.py
```

- Converts all `.pdb` files to `.pqr` using pH 7.0.
- Computes and records the total partial charge.
- Outputs `charges_at_pH7.csv`.

---

## 📄 License

This repository is licensed under the **GNU General Public License v3.0**.  
See [LICENSE](https://www.gnu.org/licenses/gpl-3.0.html) for details.

---

## ✍️ Author

**Islam K. Matar**  
2025

---

## 🙌 Citation

If you use this code in your research, please cite it appropriately. A citation format (BibTeX or otherwise) will be added here once the work is published.
