# combines hbond analysis script and fixed water bridge analysis script
# accounts for bond order (shows direct vs. indirect hbonds, which are order 0 vs order != 0)

import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.analysis.hydrogenbonds import WaterBridgeAnalysis

# set key variables to be used for downstream analysis

TOPOLOGY = 'solvated_system.prmtop'
TRAJECTORY = 'trajectory_centered.dcd'

WATER_RESNAMES = "HOH WAT"
ION_RESNAMES = "NA CL"
LIGAND_RESNAMES = "ROH 4YB 0YB"

HBOND_DISTANCE_CUTOFF = 3.5   # Angstroms, donor-acceptor
HBOND_ANGLE_CUTOFF = 150.0    # degrees, D-H...A
ORDER = 1                     # controls whether you look at first order water bridges or beyond

STRIDE_HBOND = 1

N_LAST_FRAMES = 2500   # e.g. 2500, analyze only the final 2500 frames
FRAME_START = None     # e.g. 1000, analyze from frame 1000 onward
FRAME_STOP = None      # e.g. 4000, analyze up to (not including) frame 4000

u = mda.Universe(TOPOLOGY, TRAJECTORY)
total_frames = u.trajectory.n_frames

protein_sel = 'protein'
ligand_sel = f'resname {LIGAND_RESNAMES}'

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

print(f"Protein CA: {len(protein_ca)} | "
      f"Ligand atoms: {len(ligand)} | Water O atoms: {len(water_o)}")

# mdtraj is 0-indexed, this just makes sure that the right frames are being analyzed
print(f'\nStarting at frame {FRAME_START}, ending at frame {FRAME_STOP - 1}, analyzing every {STRIDE_HBOND} frame.')

# hbond analysis using built-in heuristics to determine hbond acceptors and donors
# note this can miss some uncommon hbonds that wbanalysis will find, such as hbonds to O2N on the ligand

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
print(hb_df.head(5))
print("Saved hydrogen_bond_occupancy.tsv")

# explicitly list donors and acceptors based on glycam06 atom IDs
# extends normal list of IDs based on what i see in the pdb

GLYCAM06_DONORS = ('NT', 'OH', 'OW', 'N3', 'N', 'O',
                    'O1', 'O3', 'O4', 'O6', 'N2')       # substrate's actual H bond donors, based on pdb file
GLYCAM06_ACCEPTORS = ('NT', 'OH', 'OW', 'OY', 'OS', 'O', 'SM', 'N', 'O2',
                       'O1', 'O3', 'O4', 'O5', 'O6', 'O2N')  # same but for H bond acceptors.

# distance_type should be set to heavy which is heavy atom to heavy atom, rather than donor hydrogen to heavy atom. the former is how the hbanalysis module works. not sure why that's not consistent with waterbridge analysis by default
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
    distance_type='heavy'
)
wb.run(start=FRAME_START, stop=FRAME_STOP, step=STRIDE_HBOND)

# counts number of hbonds to determine order
# 0 = direct hbond, 1 = first order water bridge

def with_order(current, output, u):
    first, last = current[0], current[-1]
    i1 = first[1] if first[1] is not None else first[0]   # selection 1 end, heavy atom
    i2 = last[3]  if last[3]  is not None else last[2]    # selection 2 end, heavy atom
    s1, s2 = u.atoms[i1], u.atoms[i2]
    output[(s1.resname, s1.resid, s1.name,
            s2.resname, s2.resid, s2.name, len(current) - 1)] += 1

counts = wb.count_by_type(analysis_func=with_order)

print(f'\nUnique water bridge types detected: {len(counts)}')

if len(counts) == 0:
    print('\nWARNING: zero water bridges detected. Double check atom names.')

else:
    bridge_df = pd.DataFrame(
        [(*key, frac) for key, frac in counts],
        columns=['Sele1_Resname', 'Sele1_Resid', 'Sele1_Atom',
                 'Sele2_Resname', 'Sele2_Resid', 'Sele2_Atom',
                 'Order', 'Occupancy_Fraction'],
    )
    bridge_df['Occupancy_Percent'] = bridge_df['Occupancy_Fraction'] * 100
    bridge_df = bridge_df.sort_values('Occupancy_Percent', ascending=False)
    bridge_df.to_csv('water_bridge_occupancy_fixed.tsv', sep='\t', index=False)
    print('\nSaved water_bridge_occupancy_fixed.tsv')
    print(bridge_df.head(5))
    print(bridge_df['Order'].value_counts())
