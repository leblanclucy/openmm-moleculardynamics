import mdtraj as md

# loads the trajectory and solvated_system file then centers the trajectory so that the complex no longer appears as though it is 'teleporting' when it bumps against the periodic boundary conditions.
# note that although it is possible to unwrap the trajectory.dcd file on the fly in the mdanalysis.py script, there were inconsistencies with physically impossible radial distribution function values (RDF). 
# other readouts like RMSD were unaffected. thus it is highly recommended to unwrap before analysis (and it's required before visualization).
print("Loading trajectory + topology file.")
traj = md.load('trajectory.dcd', top='solvated_system.pdb')

print("Centering protein.")
traj.image_molecules(inplace=True)

print("Saving unwrapped trajectory.")
traj.save_dcd('trajectory_centered.dcd')
print("Done.")
