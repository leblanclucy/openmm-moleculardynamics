import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import mdtraj as md

# please run mm_g-pbsa-prep.py FIRST before you run this. to do MM-GBSA or MM-PBSA, you have to take your system and separate it into three components: enzyme only, substrate only, and whole complex, all without water or ions. 
# this pre-processing occurs to the trajectory file as well.. there are also some formatting incompatibilities that will crash or throw errors if they're not fixed by running that script first (mainly with mm-pbsa which is trickier as that is done through Amber).

# 1. Load the three explicit prmtop files and the trajectory prepared by mm_g-pbsa-prep.py
print("Loading dedicated topologies and trajectory...")
traj = md.load('trajectory_no-solvent.nc', top='complex.prmtop')

top_complex = app.AmberPrmtopFile('complex.prmtop')
top_enzyme  = app.AmberPrmtopFile('enzyme.prmtop')
top_substrate = app.AmberPrmtopFile('substrate.prmtop')

# 2. create the system using the .prmtop files. implicitSolvent has other options, OBC2 is a generalized born model to calculate the electrostatic solvation energy. i also tested OBC1 and had similar results. implicitSolventSaltConc should match your salt molarity set by you in openmmbuild.py for your md simulation.
def build_custom_system(prmtop_obj):
    return prmtop_obj.createSystem(
        nonbondedMethod=app.NoCutoff, 
        implicitSolvent=app.OBC2,          
        implicitSolventSaltConc=0.15*unit.molar
    )

sys_complex   = build_custom_system(top_complex)
sys_enzyme    = build_custom_system(top_enzyme)
sys_substrate = build_custom_system(top_substrate)

# enable GPU usage
platform = mm.Platform.getPlatformByName('CUDA') 
ctx_complex   = mm.Context(sys_complex, mm.VerletIntegrator(1.0*unit.femtoseconds), platform)
ctx_enzyme    = mm.Context(sys_enzyme, mm.VerletIntegrator(1.0*unit.femtoseconds), platform)
ctx_substrate = mm.Context(sys_substrate, mm.VerletIntegrator(1.0*unit.femtoseconds), platform)

# 3. verify atom counts to make sure the calculation will proceed normally, can be helpful for troubleshooting
n_complex = sys_complex.getNumParticles()
n_enzyme = sys_enzyme.getNumParticles()
n_substrate = sys_substrate.getNumParticles()

print(f"Topology configuration: Complex ({n_complex} atoms) = Enzyme ({n_enzyme} atoms) + Substrate ({n_substrate} atoms)")
if n_complex != (n_enzyme + n_substrate):
    raise ValueError("CRITICAL ERROR: Atom counts do not perfectly sum up. Trajectory slices will fail.")

# Set up clean indexing blocks
enzyme_slice = slice(0, n_enzyme)
substrate_slice = slice(n_enzyme, n_complex)

# 4. loop to process the energy of the complex, enzyme alone, and substrate alone, frame by frame
print(f"Processing {traj.n_frames} frames...")
complex_energies = []
enzyme_energies = []
substrate_energies = []
delta_g_binding = []

for frame_idx in range(traj.n_frames):
    positions = traj.openmm_positions(frame_idx)
    
    # A. energy minimize and evaluate full complex. 
    ctx_complex.setPositions(positions)
    mm.LocalEnergyMinimizer.minimize(ctx_complex, maxIterations=50)
    state_c = ctx_complex.getState(getEnergy=True, getPositions=True)
    e_complex = state_c.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    
    # B. extract positions as a raw NumPy array to avoid structural type shifts. may not be necessary
    min_positions = state_c.getPositions(asNumpy=True)
    
    # C. distribute clean slices directly to the target contexts
    ctx_enzyme.setPositions(min_positions[enzyme_slice])
    ctx_substrate.setPositions(min_positions[substrate_slice])
    
    # D. pull exact potential energy fields
    e_enzyme = ctx_enzyme.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    e_substrate = ctx_substrate.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    
    # E. Calculate ΔG
    dg = e_complex - (e_enzyme + e_substrate)
    
    complex_energies.append(e_complex)
    enzyme_energies.append(e_enzyme)
    substrate_energies.append(e_substrate)
    delta_g_binding.append(dg)

# 5. Output results as a file
output_filename = "mm_gbsa_results.txt"
with open(output_filename, "w") as f:
    f.write("="*50 + "\n")
    f.write("             FINAL MM-GBSA RESULTS             \n")
    f.write("="*50 + "\n")
    f.write(f"Avg Total Complex Energy:   {np.mean(complex_energies):12.2f} kcal/mol\n")
    f.write(f"Avg Total Enzyme Energy:    {np.mean(enzyme_energies):12.2f} kcal/mol\n")
    f.write(f"Avg Total Substrate Energy: {np.mean(substrate_energies):12.2f} kcal/mol\n")
    f.write("-"*50 + "\n")
    f.write(f"FINAL DELTA G BINDING:      {np.mean(delta_g_binding):12.2f} +/- {np.std(delta_g_binding):.2f} kcal/mol\n")
    f.write("="*50 + "\n")

print(f"\nMM-GBSA calculations complete. Results saved to {output_filename}")