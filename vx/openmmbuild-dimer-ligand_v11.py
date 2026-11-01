# purpose: set up md simulation with a dimer of the protein of interest as well as three ligands attached (lig 1 on chain a/b, lig 2 on chain a)

from openmm.app import *
from openmm import *
from openmm.unit import *
from rdkit import Chem
from openff.toolkit import Molecule
from openff.units.openmm import to_openmm
import os
import openmmforcefields
from openmmforcefields.generators import GAFFTemplateGenerator
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import sys
import parmed as pmd

# 1. load all ligands as sdf files.
a_lig1 = Molecule.from_file('lig1docked_A_withH_v2.sdf')
b_lig1 = Molecule.from_file('lig1docked_B_withH_v2.sdf')
a_lig2 = Molecule.from_file('ZINC08405079_for-md.sdf')

# 2. prepare enzyme/receptor and protonate it at a given pH using propka. pH may be changed if you're interested in the effects of pH on the system such as histidine protonation. note that this protein .pdb file has no hydrogens, that is because i had to update the coordinates of side chains after flexible docking so the hydrogen naming after exporting from pymol was incompatible with propka.
print("Protonating at pH 7.4...")
os.system("pdb2pqr --ff=AMBER --titration-state-method=propka --with-ph=7.4 dimer-flexABA-noH_v2_md.pdb --ffout=AMBER dimer-flexABA_v2_pH74_md.pdb")

# 2. load protein into topology, check and list all disulfide bonds.
print("Building OpenMM system...")
enzyme_pdb = PDBFile('dimer-flexABA_v2_pH74_md.pdb')
modeller = Modeller(enzyme_pdb.topology, enzyme_pdb.positions)

print("Existing Disulfide Bonds in Topology:")
for bond in modeller.topology.bonds():
    atom1, atom2 = bond
    if atom1.name == 'SG' and atom2.name == 'SG':
        print(f"  Bonded: {atom1.residue.name}{atom1.residue.id} to {atom2.residue.name}{atom2.residue.id}")

# 3. add ligands to protein.
print("Adding lig1 to Chain A...")
top_a = a_lig1.to_topology().to_openmm()
pos_a = to_openmm(a_lig1.conformers[0])
modeller.add(top_a, pos_a)

print("Adding lig1 to Chain B...")
top_b = b_lig1.to_topology().to_openmm()
pos_b = to_openmm(b_lig1.conformers[0])
modeller.add(top_b, pos_b)

print("Adding lig2 inhibitor to Chain A...")
top_lig2 = a_lig2.to_topology().to_openmm()
pos_lig2 = to_openmm(a_lig2.conformers[0])
modeller.add(top_lig2, pos_lig2)

print("All ligands added successfully, applying force fields...")

# 4. apply force fields using GAFF. also, solvate using OPC3 water (a 3-point water model, more recent and accurate than other 3-point models). 
gaff_generator = GAFFTemplateGenerator(molecules=[a_lig1, b_lig1, a_lig2], forcefield="gaff-2.2.20", cache="ligands_lig1-zinc.json")
forcefield = ForceField('amber19-all.xml', 'amber19/opc3.xml')
forcefield.registerTemplateGenerator(gaff_generator.generator)

# add hydrogens. verify enzyme hydrogens.
modeller.addHydrogens(forcefield, pH = 7.4)
print("--- His states of enzyme after addHydrogens ---")
for res in modeller.topology.residues():
    if res.name in ('HIS', 'HID', 'HIE', 'HIP'):
        atom_names = [a.name for a in res.atoms()]
        print(res.name, res.id, atom_names)

# create an octahedral box, add solvent (0.15 M NaCl), and neutralize charges.
# note that model='tip3p' because model in openMM does not take opc as an argument. however, this just tells openMM that you're using a three point water system, and the actual parameters for opc water are set in the forcefield above.
modeller.addSolvent(forcefield, model='tip3p', boxShape='octahedron', padding=1.0*nanometers, positiveIon='Na+', negativeIon='Cl-', ionicStrength=0.15*unit.molar, neutralize=True)

# 4. create topology, add temperature control then barostat, save solvated system. also, create a separate system export to eventually save as a .prmtop file. for some reason, hydrogen bond occupancy didn't work with the solvated_system.pdb file, although other analysis did. constraints and rigidwater have to be removed for that due to incompatibilities between openmm and parmed.
system = forcefield.createSystem(
    modeller.topology, 
    nonbondedMethod=PME, 
    nonbondedCutoff=1.0*nanometers, 
    constraints=HBonds
)

system_for_export = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometers,
    constraints=None,
    rigidWater=False
)

#rename ligands to NAP and lig2 rather than just UNK
lig_names = ['LIG1', 'LIG1', 'LIG2']
k = 0
for res in modeller.topology.residues():
    if res.name == 'UNK':
        res.name = lig_names[k]
        res.id = str(res.index + 1)
        k += 1

from openmm.app import PDBxFile
PDBxFile.writeFile(modeller.topology, modeller.positions, open('solvated_system.cif', 'w'))
structure = pmd.openmm.load_topology(modeller.topology, system_for_export, modeller.positions)
structure.save('solvated_system.prmtop', format='amber')
structure.save('solvated_system.inpcrd', format='rst7')
print(".prmtop file created.")

# 5. restrain heavy atoms by applying a penalty to their harmonic motion, because this allows for the water molecules to settle into the system and can help prevent sudden jumps in energy or unexpected changes in structure.
print("Step 1: Applying positional restraints before equilibration...")
restraint = mm.CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
restraint.addGlobalParameter('k', 100.0 * unit.kilocalories_per_mole/unit.angstrom**2)
restraint.addPerParticleParameter('x0')
restraint.addPerParticleParameter('y0')
restraint.addPerParticleParameter('z0')

restrained_atom_indices = []
for atom in modeller.topology.atoms():
    if atom.residue.name in ('HOH', 'WAT', 'NA', 'CL', 'K'):
        continue
    if atom.element is not None and atom.element.symbol != 'H':
        restraint.addParticle(atom.index, modeller.positions[atom.index])
        restrained_atom_indices.append(atom.index)

system.addForce(restraint)
print(f"Restrained {len(restrained_atom_indices)} heavy atoms.")

# 6. initializes the integrator, which sets initial temperature. change the random seed (set for reproducibility) if desired for doing replicates.
# CUDA ensures that the GPU will be used for computation. mixed precision is sufficient for most cases. higher precision is much more computationally expensive.
integrator = mm.LangevinMiddleIntegrator(0*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
seed=13
integrator.setRandomNumberSeed(seed)
print(f"Random seed: {seed}")
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}
simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setPositions(modeller.positions)

# 7. initialize simulation: perform energy minimization. report energy before and after minimization. 
print("Step 2: Energy minimization and positional restraints...")
state_before = simulation.context.getState(getEnergy=True)
print(f"Potential energy before minimization: {state_before.getPotentialEnergy()}")
simulation.minimizeEnergy(maxIterations=10000)
state_after = simulation.context.getState(getEnergy=True)
print(f"Potential energy after minimization: {state_after.getPotentialEnergy()}")

# 8. gradual heating was chosen for a similar reason as heavy atom restraints. this simulation is now NVT.
# moving the system from its crystallized structure to 300 K abruptly can cause unfolding or sudden unwanted shifts in structure, so this allows the system to settle more gently.
print("Step 3: Gradual heating phase (0 K -> 300 K over 400 ps)...")
heating_steps = 100
steps_per_increment = 2000  # 2000 steps * 2 fs = 4 ps per temperature step, allows the system to slowly heat up
for i in range(heating_steps):
    current_temp = (i / heating_steps) * 300 * unit.kelvin
    integrator.setTemperature(current_temp)
    simulation.step(steps_per_increment)
    if i % 10 == 0:
        print(f"  -> Temperature reached: {current_temp.value_in_unit(unit.kelvin):.1f} K")
integrator.setTemperature(300*unit.kelvin)

# 9. add barostat (which sets the pressure of the system) and release restraints on heavy atoms gradually. 
# do not add the barostat until the temperature is fully set after gradual heating because of pressure/temperature relationship. this simulation is now NPT.
system.addForce(mm.MonteCarloBarostat(1.0*unit.bar, 300*unit.kelvin))
simulation.context.reinitialize(preserveState=True)

# release time of 1 ns, since it's 10 steps, and 100 ps per step
print("Step 4: Releasing positional restraints gradually...")
release_steps = 10
initial_k = 100.0
for i in range(release_steps):
    scale = 1.0 - (i / (release_steps - 1))  # 1.0 -> 0.0
    simulation.context.setParameter('k', initial_k * scale * unit.kilocalories_per_mole/unit.angstrom**2)
    simulation.step(50000)  # 100 ps per stage at 2 fs/step
    print(f"  -> Restraint scale: {scale:.2f}")

print("System equilibrated at 300 K, restraints released. Starting 10 ns production run...")

# 10. set up reporters. downstream analysis such as RMSD and RMSF. will be done on the trajectory.dcd and .prmtop files. run simulation for 5 ns (2,500,000 steps * 0.002 ps). note that this is only a test run, and 100 ns+ (and/or doing replicates of 100 ns+) are strongly recommended for publishable results.
simulation.reporters.append(app.DCDReporter('trajectory.dcd', 10000)) 
simulation.reporters.append(app.CheckpointReporter('checkpoint.chk', 10000))
simulation.reporters.append(app.StateDataReporter(
    sys.stdout, 10000, step=True, potentialEnergy=True, 
    temperature=True, progress=True, totalSteps=2500000, speed=True
))

simulation.step(2500000)
print("Simulation complete.")
