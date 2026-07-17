import MDAnalysis as mda
from MDAnalysis import transformations
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import mdtraj as md
import numpy as np
import pandas as pd

# load base trajectory. this takes the .prmtop file made during the initial md setup as this helps greatly with accurate hydrogen bond occupancy prediction. using the solvated_system.pdb file can still work for rmsd, rmsf, and other calculations.

u = mda.Universe('solvated_system.prmtop', 'trajectory_centered.dcd')

# optional, used for troubleshooting. verifies that your protein and ligand were properly imported and not subdivided.
print(len(u.select_atoms('protein').fragments))
print(len(u.select_atoms('not (resname HOH WAT NA CL K)').fragments))

# selects only the enzyme and substrate for most downstream analysis rather than the water molecules and ions. probably unnecessary, used for troubleshooting. don't delete, need to look back at older versions and see what this used to be.
protein_and_ligand = u.select_atoms('not (resname HOH WAT NA CL)')

# note the unwrap function. since md is done with periodic boundary conditions, the protein will occasionally cross the boundary of the container and appear on the other side. unwrapping centers/anchors the protein which helps some of the calculations be performed properly. there is a separate unwrapping script if you want to just unwrap the trajectory and view it as is, without going through all this analysis.
# this is only necessary if you did not unwrap the trajectory already with unwrap.py. if you did not do so, change every instance of trajectory_centered.dcd to trajectory.dcd and uncomment the next code line.
# u.trajectory.add_transformations(transformations.unwrap(protein_and_ligand))
protein = u.select_atoms('protein')

# optional, can take out, used this for troubleshooting. verifies that this script can actually read the first 5 frames of your trajectory and times how long that process takes. 
import time
start = time.time()
count = 0
for ts in u.trajectory[:5]:
    count += 1
    print(f"Frame {count}, elapsed: {time.time()-start:.1f}s")

# calculate root mean square deviation (RMSD) which tells you how much your protein moves from frame to frame during the simulation.
R = rms.RMSD(u, select='protein and name CA', ref_frame=0)
R.run()
rmsd_data = R.results.rmsd  # creates an array containing [frame, time, rmsd]

# calculate radius of gyration (Rg) which shows how compact your protein was over the course of the simulation.
times = []
rg_values = []
for ts in u.trajectory:
    times.append(ts.time)
    rg_values.append((protein.radius_of_gyration()))

# load the trajectory into MDTraj for further analysis.
t = md.load('trajectory_centered.dcd', top='solvated_system.prmtop')
protein_indices = t.topology.select('protein')
t_protein = t.atom_slice(protein_indices)

# calculate solvent-accessible surface area (SASA).
sasa_per_atom = md.shrake_rupley(t_protein)
total_sasa_per_frame = sasa_per_atom.sum(axis=1) # Total protein SASA over time

# save time series (frame by frame) data, combining RMSD, Rg, and SASA into one file.
# rmsd_data[:, 2] -> extracts just the RMSD column from MDAnalysis
# rg_values -> a list of Rg numbers per frame
# total_sasa_per_frame -> a 1D numpy array of SASA values

time_series_dict = {
    'Frame': range(len(rg_values)),
    'Time_ps': times,
    'RMSD_A': rmsd_data[:, 2],
    'Radius_of_Gyration_A': rg_values,
    'Total_SASA_nm2': total_sasa_per_frame
}

# convert from a dictionary to a data frame for export.
df_time_series = pd.DataFrame(time_series_dict)

# save to a clean comma-separated text file.
df_time_series.to_csv('trajectory_time_series_summary.csv', index=False)
print("Saved time-series properties to trajectory_time_series_summary.csv")

print("Successfully calculated RMSD, SASA, and radius of gyration.")

# calculate root mean square fluctuation per residue (RMSF) showing how much each individual residue moves over time.
prealigner = align.AlignTraj(u, u, select='protein and name CA', in_memory=True).run()
average = mda.Universe('solvated_system.prmtop', 'trajectory_centered.dcd') # reference
calphas = u.select_atoms('protein and name CA')
rmsfer = rms.RMSF(calphas).run()
rmsf_data = rmsfer.results.rmsf

# export RMSF values by creating an array of residue numbers from N->C terminus, combining residue numbers and RMSF values, and saving as a text file
import numpy as np
residue_numbers = calphas.resids  
rmsf_output = np.column_stack((residue_numbers, rmsf_data))
np.savetxt('protein_rmsf.dat', rmsf_output, 
           header='Residue_ID Fluctuation_A', 
           fmt='%d %.4f')

print("Successfully calculated RMSF.")

# calculate Define Secondary Structure of Proteins (DSSP). the explanation of each code is given at https://biopython.org/docs/1.76/api/Bio.PDB.DSSP.html. note that 'C' is used as a placeholder for 'unable to predict'. need to refine as there have been updates to DSSP (see Hekkelman et al., 2025, doi 10.1002/pro.70208) and some issues with the default naming structure.

dssp_matrix = md.compute_dssp(t_protein, simplified=True)

# exporting DSSP where each row represents another frame of the md simulation
np.savetxt('secondary_structure_dssp.txt', dssp_matrix, fmt='%s', delimiter=' ')

# calculate radial distribution function, g(r) or sometimes noted as rdf, which looks at the spatial distribution of water molecules around the substrate.
ligand = u.select_atoms('resname ROH 4YB 0YB')
is_water_o = u.select_atoms('resname HOH and name O')

rdf = InterRDF(ligand, is_water_o, nbins=75, range=(0.0, 15.0))
rdf.run()

bins = rdf.results.bins
g_r = rdf.results.rdf

# export radial distribution function and DSSP.
rdf_output = np.column_stack((bins, g_r))

np.savetxt('ligand_water_rdf.dat', rdf_output, 
           header='Distance_Angstroms g(r)', 
           fmt='%.3f %.4f')

print("Successfully predicted DSSP and radial distribution function.")

# calculate hydrogen bond occupancy with both ligand and receptor as donors/acceptors given a maximum distance of 3.5 A and a maximum angle of 120 degrees.
protein_sel = 'protein'
ligand_sel = 'resname ROH 4YB 0YB'
hb_protein_donor = HydrogenBondAnalysis(
    universe=u, donors_sel=protein_sel, acceptors_sel=ligand_sel,
    hydrogens_sel='element H', d_a_cutoff=3.5, d_h_a_angle_cutoff=120.0
)
hb_protein_donor.run()
hb_ligand_donor = HydrogenBondAnalysis(
    universe=u, donors_sel=ligand_sel, acceptors_sel=protein_sel,
    hydrogens_sel='element H', d_a_cutoff=3.5, d_h_a_angle_cutoff=120.0
)
hb_ligand_donor.run()

# export hydrogen bond results, sorting by highest occupancy
n_frames = u.trajectory.n_frames
bonds_a = hb_protein_donor.results.hbonds
bonds_b = hb_ligand_donor.results.hbonds
all_bonds = np.vstack((bonds_a, bonds_b))
bond_keys = all_bonds[:, 1:4].astype(int)  # donor_idx, hydrogen_idx, acceptor_idx
unique_bonds, counts = np.unique(bond_keys, axis=0, return_counts=True)
bond_occupancy_percent = (counts / n_frames) * 100

hb_df = pd.DataFrame({
    'Donor_idx': unique_bonds[:, 0],
    'Hydrogen_idx': unique_bonds[:, 1],
    'Acceptor_idx': unique_bonds[:, 2],
    'Occupancy_Percent': bond_occupancy_percent
})
hb_df['Donor_atom'] = [f"{u.atoms[i].resname}{u.atoms[i].resid}-{u.atoms[i].name}" for i in unique_bonds[:, 0]]
hb_df['Acceptor_atom'] = [f"{u.atoms[i].resname}{u.atoms[i].resid}-{u.atoms[i].name}" for i in unique_bonds[:, 2]]
hb_df = hb_df.sort_values(by='Occupancy_Percent', ascending=False)
hb_df.to_csv('hydrogen_bond_occupancy.txt', sep='\t', index=False)

print("Successfully predicted hydrogen bonds.")
