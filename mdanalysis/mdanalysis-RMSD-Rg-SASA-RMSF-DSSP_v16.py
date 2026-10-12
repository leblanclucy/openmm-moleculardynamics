# for 100 ns trajectories (using entire trajectory unless otherwise specified). computes (apo or holo):
# RMSD for the whole protein and for the catalytic domain
# RMSF for the whole protein but using the last 2500 frames only
# radius of gyration (Rg) for the whole protein and for the catalytic domain
# solvent-accessible surface area (SASA) for the whole protein and for the catalytic domain
# secondary structure (DSSP)

# input is solvated_system.prmtop and trajectory_centered.dcd (made by unwrap.py)
# output: trajectory_time_series.csv with RMSD, Rg, SASA, % of each secondary structure; protein_rmsf.csv; secondary_structure.csv

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import mdtraj as md

# set key variables to be used for downstream analysis

TOPOLOGY = 'solvated_system.prmtop'
TRAJECTORY = 'trajectory_centered.dcd'

CATALYTIC_DOMAIN_RESID = "57-300"  # edit for each enzyme corresponding to catalytic domain from MSA

#Dm: 52-294, Ns: 71-313, Cf: 57-300, Da: 59-301. based on homology/MSA from Yoneda et al. 2025. starting at D(V/I)SS

# determine what frames will be analyzed. default all except for DSSP and RMSF.

RMSF_START = 2500
RMSF_STOP = None
RMSF_STRIDE = 1

DSSP_START = 2500
DSSP_STRIDE = 5

# load files and select frames of interest. here, it's the last 2500 frames, given 20 ps/frame, that's the last 50 ns of a trajectory.

u = mda.Universe(TOPOLOGY, TRAJECTORY)
total_frames = u.trajectory.n_frames
 
protein = u.select_atoms('protein')
protein_ca = u.select_atoms('protein and name CA')
catalytic = u.select_atoms(f'protein and resid {CATALYTIC_DOMAIN_RESID}')
catalytic_ca = u.select_atoms(f'protein and name CA and resid {CATALYTIC_DOMAIN_RESID}')
 
print(f"Total frames: {total_frames} | Analyzing all frames, except starting at frame {RMSF_START} for RMSF")
print(f"Protein CA: {len(protein_ca)} | Catalytic domain CA: {len(catalytic_ca)}")
 
assert len(catalytic_ca) > 0, \
    "Catalytic domain selection matched 0 atoms, double check residues assigned"

# calculate RMSD for both the entire protein and catalytic domain.

R_whole = rms.RMSD(u, select='protein and name CA', ref_frame=0)
R_whole.run(backend='multiprocessing', n_workers=8)
 
R_domain = rms.RMSD(
    u, select=f'protein and name CA and resid {CATALYTIC_DOMAIN_RESID}', ref_frame=0
)
R_domain.run(backend='multiprocessing', n_workers=8)
 
assert len(R_whole.results.rmsd) == len(R_domain.results.rmsd), \
    "Whole-protein and catalytic-domain RMSD arrays have different lengths"
    
# calculate Rg for both the entire protein and catalytic domain.

rg_whole = []
rg_domain = []
frame_times = []
for ts in u.trajectory:
    rg_whole.append(protein.radius_of_gyration())
    rg_domain.append(catalytic.radius_of_gyration())
    frame_times.append(ts.time)
 
assert len(rg_whole) == len(R_whole.results.rmsd), \
    "Rg frame count doesn't match RMSD frame count"

# calculate SASA for both the entire protein and catalytic domain.

protein_indices = md.load_topology(TOPOLOGY).select('protein')
t_protein = md.load(TRAJECTORY, top=TOPOLOGY, atom_indices=protein_indices)
 
assert t_protein.n_frames == len(rg_whole), \
    "mdtraj frame count doesn't match MDAnalysis frame count"
 
sasa_per_atom = md.shrake_rupley(t_protein)          # nm^2 per atom
total_sasa_per_frame = sasa_per_atom.sum(axis=1)      # nm^2 per frame
in_domain = np.isin(protein_indices, catalytic.indices)
domain_sasa_per_frame = sasa_per_atom[:, in_domain].sum(axis=1)
 
dssp_matrix = md.compute_dssp(t_protein, simplified=False)
dssp_matrix[dssp_matrix == ' '] = '-'
dssp_codes = ['H', 'B', 'E', 'G', 'I', 'T', 'S', '-']
dssp_pct = {code: (dssp_matrix == code).mean(axis=1) * 100 for code in dssp_codes}
assert np.allclose(sum(dssp_pct.values()), 100), "DSSP percentages don't sum to 100, double check codes"

dssp_frames = np.arange(DSSP_START, total_frames, DSSP_STRIDE)
residue_labels = [f'{res.name}{i + 1}' for i, res in enumerate(t_protein.topology.residues)]
dssp_df = pd.DataFrame(dssp_matrix[dssp_frames], columns=residue_labels)
dssp_df.insert(0, 'Time_ps', np.array(frame_times)[dssp_frames])
dssp_df.insert(0, 'Frame', dssp_frames)
dssp_df.to_csv('secondary_structure_full.csv', index=False)
print(f"Saved secondary_structure_full.csv ({len(dssp_frames)} frames)")

# combine output so far to one file.

time_series_df = pd.DataFrame({
    'Frame': range(len(rg_whole)),
    'Time_ps': frame_times,
    'RMSD_whole_A': R_whole.results.rmsd[:, 2],
    'RMSD_catalytic_domain_A': R_domain.results.rmsd[:, 2],
    'Rg_whole_A': rg_whole,
    'Rg_catalytic_domain_A': rg_domain,
    'SASA_whole_nm2': total_sasa_per_frame,
    'SASA_catalytic_nm2': domain_sasa_per_frame,
    'Pct_H_alpha_helix': dssp_pct['H'],
    'Pct_B_beta_bridge': dssp_pct['B'],
    'Pct_E_strand': dssp_pct['E'],
    'Pct_G_3-10_helix': dssp_pct['G'],
    'Pct_I_pi_helix': dssp_pct['I'],
    'Pct_T_turn': dssp_pct['T'],
    'Pct_S_bend': dssp_pct['S'],
    'Pct_coil': dssp_pct['-'],
})
time_series_df.to_csv('trajectory_time_series.csv', index=False)
print("Saved trajectory_time_series.csv")

# calculate RMSF for the whole protein.

def compute_rmsf(align_select, output_select, out_path):
    u_local = mda.Universe(TOPOLOGY, TRAJECTORY)
    average = align.AverageStructure(
        u_local, u_local, select=align_select, ref_frame=0
    ).run(start=RMSF_START, stop=RMSF_STOP, step=RMSF_STRIDE)
    ref = average.results.universe
    align.AlignTraj(u_local, ref, select=align_select, in_memory=True).run(
        start=RMSF_START, stop=RMSF_STOP, step=RMSF_STRIDE
    )
 
    atoms = u_local.select_atoms(output_select)
    rmsf_result = rms.RMSF(atoms).run(start=RMSF_START, stop=RMSF_STOP, step=RMSF_STRIDE)
 
    output = np.column_stack((atoms.resids, rmsf_result.results.rmsf))
    pd.DataFrame(output, columns=['Residue_ID', 'RMSF_A']).astype({'Residue_ID': int}).to_csv(out_path, index=False)
    return output

compute_rmsf('protein and name CA', 'protein and name CA', 'protein_rmsf.csv')
print("Saved protein_rmsf.csv")
print("\nMD analysis done.")