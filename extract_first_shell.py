import MDAnalysis as mda
import numpy as np
import argparse
from MDAnalysis.lib.distances import distance_array
from tqdm import tqdm
import os

def save_density_as_dx(density, dx_file, grid_spacing, half_box):
    nx, ny, nz = density.shape
    origin = (-half_box, -half_box, -half_box)

    with open(dx_file, 'w') as f:
        f.write(f"object 1 class gridpositions counts {nx} {ny} {nz}\n")
        f.write(f"origin {origin[0]} {origin[1]} {origin[2]}\n")
        f.write(f"delta {grid_spacing} 0 0\n")
        f.write(f"delta 0 {grid_spacing} 0\n")
        f.write(f"delta 0 0 {grid_spacing}\n")
        f.write(f"object 2 class gridconnections counts {nx} {ny} {nz}\n")
        f.write(f"object 3 class array type double rank 0 items {nx*ny*nz} data follows\n")

        flat = density.flatten()
        for i in range(0, len(flat), 3):
            line = " ".join(f"{v:.5e}" for v in flat[i:i+3])
            f.write(line + "\n")

        f.write("object 4 class field\n")
        f.write("component \"positions\" value 1\n")
        f.write("component \"connections\" value 2\n")
        f.write("component \"data\" value 3\n")


def compute_density(tpr, xtc, center_index, cutoff, resnames, grid_spacing=0.5):
    u = mda.Universe(tpr, xtc)
    center_index0 = center_index - 1
    center_atom = u.atoms[center_index0:center_index0+1]

    candidates = u.select_atoms(f"resname {' '.join(resnames)}")
    print(f"Center atom: {center_atom[0].resname} {center_atom[0].name} at index {center_index}")
    print(f"Selected {len(candidates)} atoms from residues: {resnames}")

    # Read box size from the first frame of the trajectory
    u.trajectory[0]
    box_lengths = u.trajectory.ts.dimensions[:3]
    min_box_size = min(box_lengths)
    box_size = min_box_size  # use smallest box dimension to avoid PBC wrap issues
    half_box = box_size / 2.0

    nbins = int(box_size / grid_spacing)
    edges = np.linspace(-half_box, half_box, nbins + 1)
    hist = np.zeros((nbins, nbins, nbins), dtype=np.float32)
    n_frames = 0

    for ts in tqdm(u.trajectory, desc="Processing frames"):
        box = ts.dimensions[:3]
        pos_center = center_atom.positions[0]

        vecs = candidates.positions - pos_center
        vecs -= box * np.round(vecs / box)
        dists = np.linalg.norm(vecs, axis=1)
        shell_atoms = vecs[dists <= cutoff]

        in_box = shell_atoms[(np.abs(shell_atoms) < half_box).all(axis=1)]
        if len(in_box) == 0:
            continue

        h, _ = np.histogramdd(in_box, bins=(edges, edges, edges))
        hist += h
        n_frames += 1

    if n_frames == 0:
        raise RuntimeError("❌ No shell atoms detected. Check cutoff or index.")

    hist /= n_frames

    tag = f"{center_index}_{'_'.join(resnames)}"
    dx_filename = f"sdf_{tag}.dx"
    pdb_filename = f"center_atom_{center_index}.pdb"

    save_density_as_dx(hist, dx_filename, grid_spacing, half_box)
    print(f"✅ Saved DX file: {dx_filename}")

    center_shifted = mda.Merge(center_atom)
    center_shifted.atoms.positions = np.array([[0.0, 0.0, 0.0]])
    center_shifted.atoms.write(pdb_filename)
    print(f"✅ Saved center atom at origin to: {pdb_filename}")


def main():
    parser = argparse.ArgumentParser(description="Compute spatial density of first-shell atoms around a given atom.")
    parser.add_argument("-t", "--tpr", required=True, help="Input .tpr file")
    parser.add_argument("-x", "--xtc", required=True, help="Input .xtc trajectory after centering and fitting for tran+rot motion. See 'gmx spatial' for more info.")
    parser.add_argument("-i", "--index", type=int, required=True, help="GROMACS atom index (1-based)")
    parser.add_argument("-c", "--cutoff", type=float, required=True, help="First shell cutoff distance (in Å)")
    parser.add_argument("-r", "--resnames", nargs='+', required=True, help="Residue names to include (e.g., SOL MLT CHL)")
    parser.add_argument("--spacing", type=float, default=0.5, help="Grid spacing (default: 0.5 Å)")
    args = parser.parse_args()

    compute_density(args.tpr, args.xtc, args.index, args.cutoff, args.resnames,
                    grid_spacing=args.spacing)


if __name__ == "__main__":
    main()

