import mdtraj as md

# loads the trajectory and solvated_system file then centers the trajectory so that the complex no longer appears as though it is 'teleporting' when it bumps against the periodic boundary conditions.
print("Loading trajectory + topology file.")
traj = md.load('trajectory.dcd', top='solvated_system.pdb')

print("Centering protein.")
traj.image_molecules(inplace=True)

print("Saving unwrapped trajectory.")
traj.save_dcd('trajectory_centered.dcd')
print("Done.")
