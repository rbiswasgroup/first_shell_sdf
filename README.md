# First Shell Spatial Density Calculator

This tool computes the **3D spatial density function (SDF)** of atoms from selected residue types that lie within a **first solvation shell cutoff** around a specified atom from a GROMACS trajectory. It outputs:

* A `.dx` volumetric density map (OpenDX format) centered on the reference atom.
* A `.pdb` file with the reference atom at the origin for visualization.

---

## 📦 Requirements

* Python 3.7+
* [MDAnalysis](https://www.mdanalysis.org/)
* `numpy`
* `tqdm`

Install dependencies:

```bash
pip install MDAnalysis numpy tqdm
```
---
## Installation
You can install it by:
```bash
pip install git+https://github.com/rbiswasgroup/first_shell_sdf.git
```
---
Or after extracting the zip file, inside the folder, just execute
```bash
pip install .
```
---

## 📥 Input Files

* `topol.tpr`: GROMACS topology file
* `traj.xtc`: GROMACS trajectory (see preprocessing section below)
* Atom index (GROMACS-style, 1-based)
* Residue names to include in the shell (e.g., `SOL`, `MLT`, `CHL`)
* Shell cutoff distance (in Å)

---

## ▶️ Usage

```bash
firstshellsdf -t topol.tpr -x traj_fitted.xtc -i 4201 -c 3.5 -r SOL MLT CHL
```

### Optional flags:

* `--spacing`: Grid spacing in Å (default = `0.5`)

### Full command-line help:

```bash
usage: firstshellsdf [-t TPR] [-x XTC] [-i INDEX] [-c CUTOFF] [-r RESNAMES [RESNAMES ...]] [--spacing SPACING]

Compute spatial density of first-shell atoms around a given atom.

optional arguments:
  -t, --tpr       Input .tpr file
  -x, --xtc       Input .xtc trajectory
  -i, --index     GROMACS atom index (1-based)
  -c, --cutoff    First shell cutoff distance (in Å)
  -r, --resnames  Residue names to include (e.g., SOL MLT CHL)
  --spacing       Grid spacing (default: 0.5 Å)
```

---

## 🔄 Preprocessing the Trajectory

The `.xtc` input should be centered and fitted to remove translational motion of the chosen atom. To prepare your trajectory:

1. **Make an index file (`index.ndx`)** with the atom of interest (e.g. atom 4201):

   ```
   [ selected_atom ]
   4201
   ```

2. **Center the trajectory:**

   ```bash
   gmx trjconv -s topol.tpr -f traj.xtc -o traj_centered.xtc \
       -n index.ndx -boxcenter tric -center -ur compact -pbc none -skip 10
   ```

3. **Fit to remove translational motion:**

   ```bash
   gmx trjconv -s topol.tpr -f traj_centered.xtc -o traj_fitted.xtc \
       -fit trans -n index.ndx
   ```

4. **Use `traj_fitted.xtc` as the input to this script.**

---

## 📤 Outputs

* `sdf_<index>_<resnames>.dx`: Spatial density map centered on the reference atom
* `center_atom_<index>.pdb`: Reference atom placed at origin (for alignment)

---

## 🖼 Visualization

1. Load `.dx` into **VMD**, **ChimeraX**, or **PyMOL** as a volumetric map.
2. Load `.pdb` as a VDW sphere.
3. Set isosurfaces or transparency to visualize shell distribution.

---

## 📄 License

MIT License

Copyright (c) [2025] [Rajib Biswas]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

