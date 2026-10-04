# quality checker for molecular dynamics trajectories to ensure no self interaction artifacts

import numpy as np
import MDAnalysis as mda
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

# ---- USER SETTINGS ----
TOPOLOGY = "solvated_system.prmtop"        # your topology file (will have a different ending for diff programs)
TRAJECTORY = "trajectory.dcd"       # again, your trajectory file
SELECTION = "protein"         # just protein is sufficient here unless you're looking at specific residues
CUTOFF = 12.0                 # your nonbonded cutoff, in the SAME units as the box (angstroms by default)
STRIDE = 10                    # analyze every nth frame, 5-10 seems to balance speed and accuracy
# ------------------------

u = mda.Universe(TOPOLOGY, TRAJECTORY)
sel = u.select_atoms(SELECTION)
print(u.trajectory[0].triclinic_dimensions)
print(u.trajectory[-1].triclinic_dimensions)

if len(sel) == 0:
    raise ValueError(f"Selection '{SELECTION}' returned zero atoms - check your selection string.")

times, min_dists = [], []

# looks at the coordinates (ijk) trajectory every nth frame over time
# note, i had originally added the .copy() method to everything with ts because of an issue raised on github but it doesn't seem to change the result
for ts in u.trajectory[::STRIDE]:
    box_vectors = ts.triclinic_dimensions
    positions = sel.unwrap(compound="fragments", reference="com")
    shifts = [
        i * box_vectors[0] + j * box_vectors[1] + k * box_vectors[2]
        for i in (-1, 0, 1)
        for j in (-1, 0, 1)
        for k in (-1, 0, 1)
        if not (i == 0 and j == 0 and k == 0)
    ]

    tree = cKDTree(positions)
    frame_min = np.inf
    for shift in shifts:
        d, _ = tree.query(positions + shift, k=1)
        frame_min = min(frame_min, d.min())

    times.append(ts.time)
    min_dists.append(frame_min)

#sets up the array comparing the minimum distances vs time
times = np.array(times)
min_dists = np.array(min_dists)

# ---- plot ----
plt.figure(figsize=(8, 4.5))
plt.plot(times, min_dists, lw=1)
plt.axhline(CUTOFF, color="red", linestyle="--", label=f"cutoff = {CUTOFF}")
plt.xlabel("Time (ps)")
plt.ylabel("Min. distance to periodic image")
plt.title(f"Self-image distance: '{SELECTION}'")
plt.legend()
plt.tight_layout()
plt.savefig("min_image_distance.png", dpi=150)

#displays the plot, this shows on my other PC running Ubuntu but not through my PC that's using WSL. however it does save the figure and once you close the plot window it does tell you how many violations there were.
plt.show()

n_violations = int((min_dists < CUTOFF).sum())
print(f"{n_violations} / {len(min_dists)} frames have min-image distance below the cutoff ({CUTOFF}).")
if n_violations > 0:
    print("The protein is interacting with its own periodic image in these frames - "
          "consider a bigger box and re-running.")
