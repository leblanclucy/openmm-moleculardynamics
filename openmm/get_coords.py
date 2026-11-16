# purpose: select ligand and determine its coordinates for GIST

import mdtraj as md
import numpy as np

# use the last frame (50,000th frame for me)
with md.open("gist_seed13.dcd") as f:
	n = len(f)
t = md.load_frame("gist_seed13.dcd", n-1, top="solvated_system.prmtop")

# change resnames to whatever corresponds to ligand
resnames = {"4YB", "0YB", "ROH"}
sel = np.array([a.index for a in t.topology.atoms if a.residue.name in resnames])
print(len(sel))
print(sorted({t.topology.atom(i).residue.name for i in sel}))
print(sorted({r.name for r in t.topology.residues
              if not r.is_protein and not r.is_water}))
xyz = t.xyz[0][sel]*10
lo, hi = xyz.min(0), xyz.max(0)

# define grid center and dimensions
center = np.round((lo + hi) / 2, 2)

# provides a 5 angstrom margin
dims = np.ceil((hi - lo + 2 * 5.0) / 0.5).astype(int)

with open("coords.log", "w") as f:
	print("grid center:", center, file=f)
	print("grid dims:", dims, file=f)
