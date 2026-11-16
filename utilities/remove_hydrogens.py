import MDAnalysis as mda

u = mda.Universe('protein.pdb')
heavy_atoms = u.select_atoms('not name H*')

heavy_atoms.write('protein-noH.pdb')
print("Hydrogens removed. Written to protein-noH.pdb")
