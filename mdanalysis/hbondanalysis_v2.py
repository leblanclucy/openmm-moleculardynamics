# initialize and import packages

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.analysis.hydrogenbonds import WaterBridgeAnalysis
import mdtraj as md

# specify file names
TOPOLOGY = 'solvated_system.prmtop'
TRAJECTORY = 'trajectory_centered.dcd'

# define residue names
WATER_RESNAMES = "HOH"
ION_RESNAMES = "NA CL"
LIGAND_RESNAMES = "ROH 4YB 0YB"

# these are standard cutoffs i see in the literature
HBOND_DISTANCE_CUTOFF = 3.5   # Angstroms, donor-acceptor
HBOND_ANGLE_CUTOFF = 150.0    # degrees, D-H...A

# this is leftover from my other analysis script, i use a stride of 1 for hydrogen bonds since they're short lived and overall a key part of my analysis
# for other parts of mdanalysis like DSSP which does not change as fast, i might use a longer stride
STRIDE_HBOND = 1

# determine what frames will be analyzed.
N_LAST_FRAMES = 2500   # e.g. 2500, analyze only the final 2500 frames
FRAME_START = None     # e.g. 1000, analyze from frame 1000 onward
FRAME_STOP = None      # e.g. 4000, analyze up to (not including) frame 4000

# load files and select frames of interest. here, it's the last 2500 frames, given 20 ps/frame, that's the last 50 ns of a trajectory.
# if you write to the trajectory more frequently (every 10 ps is more standard) you would be covering more frames given the same time

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
