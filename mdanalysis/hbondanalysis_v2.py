# initialize and import

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

WATER_RESNAMES = "HOH WAT"
ION_RESNAMES = "NA CL"
LIGAND_RESNAMES = "ROH 4YB 0YB"

HBOND_DISTANCE_CUTOFF = 3.5   # Angstroms, donor-acceptor
HBOND_ANGLE_CUTOFF = 150.0    # degrees, D-H...A

# stride depends on what you're measuring. for example DSSP does not change quickly so a stride of 5 is acceptable, whereas hydrogen bonds are short-lived.

STRIDE_HBOND = 1

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
ligand = u.select_atoms(f'resname {LIGAND_RESNAMES}')
water_o = u.select_atoms(f'resname {WATER_RESNAMES} and name O')
 
print(f"Total frames: {total_frames} | Analyzing frames [{FRAME_START}:{FRAME_STOP}] "
      f"({FRAME_STOP - FRAME_START} frames)")
print(f"Protein CA: {len(protein_ca)} | "
      f"Ligand atoms: {len(ligand)} | Water O atoms: {len(water_o)}")
 
assert len(ligand) > 0, "Ligand selection matched 0 atoms -- check LIGAND_RESNAMES"
assert len(water_o) > 0, "Water oxygen selection matched 0 atoms -- check WATER_RESNAMES"
assert 0 <= FRAME_START < FRAME_STOP <= total_frames, \
    f"Invalid frame range [{FRAME_START}:{FRAME_STOP}] for a trajectory of {total_frames} frames"

# calculate hydrogen bond occupancy between the protein and ligand.

protein_sel = 'protein'
ligand_sel = f'resname {LIGAND_RESNAMES}'

hb = HydrogenBondAnalysis(
    universe=u,
    between=[protein_sel, ligand_sel],
    d_a_cutoff=HBOND_DISTANCE_CUTOFF,
    d_h_a_angle_cutoff=HBOND_ANGLE_CUTOFF,
)
protein_or_ligand_sel = f'({protein_sel}) or ({ligand_sel})'
hb.hydrogens_sel = hb.guess_hydrogens(protein_or_ligand_sel)
hb.acceptors_sel = hb.guess_acceptors(protein_or_ligand_sel)

hb.run(start=FRAME_START, stop=FRAME_STOP, step=STRIDE_HBOND)
n_frames_analyzed = len(range(FRAME_START, FRAME_STOP, STRIDE_HBOND))

id_counts = hb.count_by_ids()

hb_df = pd.DataFrame(id_counts, columns=['Donor_idx', 'Hydrogen_idx', 'Acceptor_idx', 'Count'])
hb_df['Occupancy_Percent'] = hb_df['Count'].astype(float) / n_frames_analyzed * 100

id_to_atom = {atom.id: atom for atom in u.atoms}
hb_df['Donor_atom'] = [f"{id_to_atom[i].resname}{id_to_atom[i].resid}-{id_to_atom[i].name}"
                        for i in hb_df['Donor_idx']]

hb_df['Acceptor_atom'] = [f"{id_to_atom[i].resname}{id_to_atom[i].resid}-{id_to_atom[i].name}"
                           for i in hb_df['Acceptor_idx']]
                           
hb_df = hb_df.sort_values(by='Occupancy_Percent', ascending=False)
hb_df.to_csv('hydrogen_bond_occupancy.tsv', sep='\t', index=False)
print("Saved hydrogen_bond_occupancy.tsv")
