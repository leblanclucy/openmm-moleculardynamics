# for 100 ns trajectories (using last 2500 frames, or 50 ns). computes the following:
# RMSD for the whole protein and for the catalytic domain
# RMSF for the whole protein and for the catalytic domain
# radius of gyration (Rg)
# solvent-accessible surface area (SASA)
# secondary structure (DSSP)
# substrate-water radial distribution function (RDF)
# hydrogen bond occupancy between protein and substrate (direct)
# hydrogen bond occupancy via water bridges (first order only, less computationally expensive)

# input is solvated_system.prmtop and trajectory_centered.dcd (made by unwrap.py)
# output: trajectory_time_series.csv with RMSD, Rg, SASA, %SS; protein_rmsf.dat, catalytic_domain_rmsf.dat, secondary_structure_full.txt with DSSP, ligand_water_rdf.dat with RDF, hydrogen_bond_occupancy.tsv, and water_bridge_occupancy.tsv

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.analysis.hydrogenbonds import WaterBridgeAnalysis
import mdtraj as md

# set key variables to be used for downstream analysis

TOPOLOGY = 'solvated_system.prmtop'
TRAJECTORY = 'trajectory_centered.dcd'

CATALYTIC_DOMAIN_RESID = "58-300"  # edit for each enzyme corresponding to catalytic domain from MSA

WATER_RESNAMES = "HOH WAT"
ION_RESNAMES = "NA CL"
LIGAND_RESNAMES = "ROH 4YB 0YB"

HBOND_DISTANCE_CUTOFF = 3.5   # Angstroms, donor-acceptor
HBOND_ANGLE_CUTOFF = 150.0    # degrees, D-H...A
ORDER = 1                     # controls whether you look at first order water bridges or beyond

# stride depends on what you're measuring. for example DSSP does not change quickly so a stride of 5 is acceptable, whereas hydrogen bonds are short-lived.

STRIDE_SLOW = 5
STRIDE_HBOND = 1
STRIDE_RDF = 2

# determine what frames will be analyzed.

N_LAST_FRAMES = 2500   # e.g. 2500, analyze only the final 2500 frames
FRAME_START = None     # e.g. 1000, analyze from frame 1000 onward
FRAME_STOP = None      # e.g. 4000, analyze up to (not including) frame 4000

# load files and select frames of interest. here, it's the last 2500 frames, given 20 ps/frame, that's the last 50 ns of a trajectory.

u = mda.Universe(TOPOLOGY, TRAJECTORY)
total_frames = u.trajectory.n_frames

assert not (N_LAST_FRAMES is not None and FRAME_START is not None), \
    "Set either N_LAST_FRAMES or FRAME_START, not both"

if N_LAST_FRAMES is not None:
    assert N_LAST_FRAMES <= total_frames, \
        f"N_LAST_FRAMES ({N_LAST_FRAMES}) exceeds total frame count ({total_frames})"
    FRAME_START = total_frames - N_LAST_FRAMES
    FRAME_STOP = total_frames
else:
    FRAME_START = 0 if FRAME_START is None else FRAME_START
    FRAME_STOP = total_frames if FRAME_STOP is None else FRAME_STOP
 
protein = u.select_atoms('protein')
protein_ca = u.select_atoms('protein and name CA')
catalytic_ca = u.select_atoms(f'protein and name CA and resid {CATALYTIC_DOMAIN_RESID}')
ligand = u.select_atoms(f'resname {LIGAND_RESNAMES}')
water_o = u.select_atoms(f'resname {WATER_RESNAMES} and name O')
 
print(f"Total frames: {total_frames} | Analyzing frames [{FRAME_START}:{FRAME_STOP}] "
      f"({FRAME_STOP - FRAME_START} frames)")
print(f"Protein CA: {len(protein_ca)} | Catalytic domain CA: {len(catalytic_ca)} | "
      f"Ligand atoms: {len(ligand)} | Water O atoms: {len(water_o)}")
 
assert len(catalytic_ca) > 0, \
    "Catalytic domain selection matched 0 atoms -- check CATALYTIC_DOMAIN_RESID"
assert len(ligand) > 0, "Ligand selection matched 0 atoms -- check LIGAND_RESNAMES"
assert len(water_o) > 0, "Water oxygen selection matched 0 atoms -- check WATER_RESNAMES"
assert 0 <= FRAME_START < FRAME_STOP <= total_frames, \
    f"Invalid frame range [{FRAME_START}:{FRAME_STOP}] for a trajectory of {total_frames} frames"

# calculate RMSD for both the entire protein and part of it with a stride of 5.

R_whole = rms.RMSD(u, select='protein and name CA', ref_frame=0)
R_whole.run(start=FRAME_START, stop=FRAME_STOP, step=STRIDE_SLOW, backend='multiprocessing', n_workers=8)
 
R_domain = rms.RMSD(
    u, select=f'protein and name CA and resid {CATALYTIC_DOMAIN_RESID}', ref_frame=0
)
R_domain.run(start=FRAME_START, stop=FRAME_STOP, step=STRIDE_SLOW, backend='multiprocessing', n_workers=8)
 
assert len(R_whole.results.rmsd) == len(R_domain.results.rmsd), \
    "Whole-protein and catalytic-domain RMSD arrays have different lengths"
    
# calculate Rg with a stride of 5.

rg_values = []
frame_times = []
for ts in u.trajectory[FRAME_START:FRAME_STOP:STRIDE_SLOW]:
    rg_values.append(protein.radius_of_gyration())
    frame_times.append(ts.time)
 
assert len(rg_values) == len(R_whole.results.rmsd), \
    "Rg frame count doesn't match RMSD frame count -- stride/range mismatch"

# calculate SASA and DSSP using mdtraj with a stride of 5.

t_full = md.load(TRAJECTORY, top=TOPOLOGY)
t = t_full[FRAME_START:FRAME_STOP:STRIDE_SLOW]
protein_indices = t.topology.select('protein')
t_protein = t.atom_slice(protein_indices)
 
assert t_protein.n_frames == len(rg_values), \
    "mdtraj frame count doesn't match MDAnalysis frame count -- check range/stride"
 
sasa_per_atom = md.shrake_rupley(t_protein)          # nm^2 per atom
total_sasa_per_frame = sasa_per_atom.sum(axis=1)      # nm^2 per frame
 
dssp_matrix = md.compute_dssp(t_protein, simplified=True)  # (n_frames, n_residues)
np.savetxt('secondary_structure_full.txt', dssp_matrix, fmt='%s', delimiter=' ')
 
# Aggregate per-frame % helix / sheet / coil for the main time-series table.
pct_helix = (dssp_matrix == 'H').mean(axis=1) * 100
pct_sheet = (dssp_matrix == 'E').mean(axis=1) * 100
pct_coil = (dssp_matrix == 'C').mean(axis=1) * 100

# combine output so far to one file.

time_series_df = pd.DataFrame({
    'Frame': range(len(rg_values)),
    'Time_ps': frame_times,
    'RMSD_whole_A': R_whole.results.rmsd[:, 2],
    'RMSD_catalytic_domain_A': R_domain.results.rmsd[:, 2],
    'Radius_of_Gyration_A': rg_values,
    'Total_SASA_nm2': total_sasa_per_frame,
    'Pct_Helix': pct_helix,
    'Pct_Sheet': pct_sheet,
    'Pct_Coil': pct_coil,
})
time_series_df.to_csv('trajectory_time_series.csv', index=False)
print("Saved trajectory_time_series.csv")

# calculate RMSF for both the whole protein and the catalytic domain.

def compute_rmsf(align_select, output_select, out_path):
    u_local = mda.Universe(TOPOLOGY, TRAJECTORY)
    average = align.AverageStructure(
        u_local, u_local, select=align_select, ref_frame=0
    ).run(start=FRAME_START, stop=FRAME_STOP, step=STRIDE_SLOW)
    ref = average.results.universe
    align.AlignTraj(u_local, ref, select=align_select, in_memory=True).run(
        start=FRAME_START, stop=FRAME_STOP, step=STRIDE_SLOW
    )
 
    atoms = u_local.select_atoms(output_select)
    rmsf_result = rms.RMSF(atoms).run(start=FRAME_START, stop=FRAME_STOP, step=STRIDE_SLOW)
 
    output = np.column_stack((atoms.resids, rmsf_result.results.rmsf))
    np.savetxt(out_path, output, header='Residue_ID Fluctuation_A', fmt='%d %.4f')
    return output

compute_rmsf('protein and name CA', 'protein and name CA', 'protein_rmsf.dat')
print("Saved protein_rmsf.dat")

compute_rmsf(
    f'protein and name CA and resid {CATALYTIC_DOMAIN_RESID}',
    f'protein and name CA and resid {CATALYTIC_DOMAIN_RESID}',
    'catalytic_domain_rmsf.dat'
)
print("Saved catalytic_domain_rmsf.dat")

# calculate RDF with a stride of 2.

rdf = InterRDF(ligand, water_o, nbins=75, range=(0.0, 15.0))
rdf.run(start=FRAME_START, stop=FRAME_STOP, step=STRIDE_RDF)
 
rdf_output = np.column_stack((rdf.results.bins, rdf.results.rdf))
np.savetxt('ligand_water_rdf.dat', rdf_output,
           header='Distance_Angstroms g(r)', fmt='%.3f %.4f')
print("Saved ligand_water_rdf.dat")

# calculate hydrogen bond occupancy.

protein_sel = 'protein'
ligand_sel = f'resname {LIGAND_RESNAMES}'
 
hb_protein_donor = HydrogenBondAnalysis(
    universe=u, donors_sel=protein_sel, acceptors_sel=ligand_sel,
    hydrogens_sel='element H',
    d_a_cutoff=HBOND_DISTANCE_CUTOFF, d_h_a_angle_cutoff=HBOND_ANGLE_CUTOFF
)
hb_protein_donor.run(start=FRAME_START, stop=FRAME_STOP, step=STRIDE_HBOND)
 
hb_ligand_donor = HydrogenBondAnalysis(
    universe=u, donors_sel=ligand_sel, acceptors_sel=protein_sel,
    hydrogens_sel='element H',
    d_a_cutoff=HBOND_DISTANCE_CUTOFF, d_h_a_angle_cutoff=HBOND_ANGLE_CUTOFF
)
hb_ligand_donor.run(start=FRAME_START, stop=FRAME_STOP, step=STRIDE_HBOND)

n_frames_analyzed = len(range(FRAME_START, FRAME_STOP, STRIDE_HBOND))
 
bonds_a = hb_protein_donor.results.hbonds
bonds_b = hb_ligand_donor.results.hbonds
all_bonds = np.vstack((bonds_a, bonds_b))
bond_keys = all_bonds[:, 1:4].astype(int)  # donor_idx, hydrogen_idx, acceptor_idx
unique_bonds, counts = np.unique(bond_keys, axis=0, return_counts=True)
bond_occupancy_percent = (counts / n_frames_analyzed) * 100
 
hb_df = pd.DataFrame({
    'Donor_idx': unique_bonds[:, 0],
    'Hydrogen_idx': unique_bonds[:, 1],
    'Acceptor_idx': unique_bonds[:, 2],
    'Occupancy_Percent': bond_occupancy_percent,
})
hb_df['Donor_atom'] = [f"{u.atoms[i].resname}{u.atoms[i].resid}-{u.atoms[i].name}"
                        for i in unique_bonds[:, 0]]
hb_df['Acceptor_atom'] = [f"{u.atoms[i].resname}{u.atoms[i].resid}-{u.atoms[i].name}"
                           for i in unique_bonds[:, 2]]
hb_df = hb_df.sort_values(by='Occupancy_Percent', ascending=False)
hb_df.to_csv('hydrogen_bond_occupancy.tsv', sep='\t', index=False)
print("Saved hydrogen_bond_occupancy.tsv")

# calculate hydrogen bonds through water bridges that link the substrate and active site.
# first add additional atom names, beyond just the default, as OpenMM calls water's oxygen atoms just "O" and the carbohydrate substrate's oxygens have special names like O1, O3 etc.

GLYCAM06_DONORS = ('NT', 'OH', 'OW', 'N3', 'N', 'O',
                    'O1', 'O3', 'O4', 'O6', 'N2')       # substrate's actual H bond donors, based on pdb file
GLYCAM06_ACCEPTORS = ('NT', 'OH', 'OW', 'OY', 'OS', 'O', 'SM', 'N', 'O2',
                       'O1', 'O3', 'O4', 'O5', 'O6', 'O2N')  # same but for H bond acceptors.

wb = WaterBridgeAnalysis(
    universe=u,
    selection1=protein_sel,
    selection2=ligand_sel,
    water_selection=f'resname {WATER_RESNAMES}',
    order=ORDER,
    update_water_selection=True,
    distance=HBOND_DISTANCE_CUTOFF,
    angle=HBOND_ANGLE_CUTOFF,
    donors=GLYCAM06_DONORS,
    acceptors=GLYCAM06_ACCEPTORS,
)
wb.run(start=FRAME_START, stop=FRAME_STOP, step=STRIDE_HBOND)

counts = wb.count_by_type()

print(f"\nUnique water bridge types detected: {len(counts)}")

if len(counts) == 0:
    print("\nWARNING: zero water bridges detected. Double check atom names.")

# export water bridge analysis results.

else:
    bridge_df = pd.DataFrame(counts, columns=[
        'Sele1_Index', 'Sele2_Index',
        'Sele1_Resname', 'Sele1_Resid', 'Sele1_Atom',
        'Sele2_Resname', 'Sele2_Resid', 'Sele2_Atom',
        'Occupancy_Fraction',
    ])
    bridge_df['Occupancy_Percent'] = bridge_df['Occupancy_Fraction'] * 100
    bridge_df = bridge_df.sort_values('Occupancy_Percent', ascending=False)
    bridge_df.to_csv('water_bridge_occupancy.tsv', sep='\t', index=False)
    print("\nSaved water_bridge_occupancy.tsv")
    print(bridge_df.head(5))
 
print("\nMD analysis done.")
