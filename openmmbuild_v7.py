from openmm.app import *
from openmm import *
from openmm.unit import *
from rdkit import Chem
import os
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import sys
import parmed as pmd

# 1. prepare enzyme/receptor and protonate it at a given pH using propka. pH may be changed if you're interested in the effects of pH on the system such as histidine protonation.
print("Protonating enzyme at pH 5.0...")
os.system("pdb2pqr --ff=AMBER --titration-state-method=propka --with-ph=5.0 <enzyme>.pdb --ffout=AMBER <enzyme_pH5>.pdb")

# 2. pick enzyme and substrate, combine into one object, check and list all disulfide bonds, add more disulfide bonds if needed.
print("Building OpenMM system...")
enzyme_pdb = PDBFile('<enzyme_pH5>.pdb')
ligand_pdb = PDBFile('<substrate_topbindingmode>.pdb')

modeller = Modeller(enzyme_pdb.topology, enzyme_pdb.positions)
modeller.add(ligand_pdb.topology, ligand_pdb.positions)

print("Existing Disulfide Bonds in Topology:")
for bond in modeller.topology.bonds():
    atom1, atom2 = bond
    if atom1.name == 'SG' and atom2.name == 'SG':
        print(f"  Bonded: {atom1.residue.name}{atom1.residue.id} to {atom2.residue.name}{atom2.residue.id}")

# optional. in my limited experience the structure already had cysteines close enough to automatically establish them as disulfide bonds, and explicitly declaring them or trying to add them again was not necessary.
n_bonds_before = __builtins__.sum(1 for _ in modeller.topology.bonds())
modeller.topology.createDisulfideBonds(modeller.positions)
n_bonds_after = __builtins__.sum(1 for _ in modeller.topology.bonds())
print(f"Disulfide bonds added: {n_bonds_after - n_bonds_before}")

# 3. apply force fields specific to your enzyme/substrate. also, solvate using OPC water (a 4-point water model). verify the protonation states of histidine residues.
forcefield = ForceField('amber19-all.xml', 'amber14/GLYCAM_06j-1.xml', 'amber19/opc.xml')

ligand_residue_names = {'ROH', '4YB', '0YB'}
ligand_residues = [r for r in modeller.topology.residues() if r.name in ligand_residue_names]

print(f"Found {len(ligand_residues)} ligand residues for hydrogen addition: {[(r.chain.id, r.name) for r in ligand_residues]}")
assert len(ligand_residues) > 0, "No ligand residues found — check chain IDs"

Modeller.loadHydrogenDefinitions('glycam-hydrogens.xml')
modeller.addHydrogens(forcefield)

print("--- His states of enzyme after addHydrogens ---")
for res in modeller.topology.residues():
    if res.name in ('HIS', 'HID', 'HIE', 'HIP'):
        atom_names = [a.name for a in res.atoms()]
        print(res.name, res.id, atom_names)

# note that model='tip4pew' because model in openMM does not take opc as an argument. however, this just tells openMM that you're using a four point water system, and the actual parameters for opc water are set in the forcefield above.
modeller.addSolvent(forcefield, model='tip4pew', boxShape='octahedron', padding=1.1*nanometers, ionicStrength=0.15*unit.molar)

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

# change the random seed if desired. CUDA ensures that the GPU will be used for computation. mixed precision is sufficient for most cases. higher precision is much more computationally expensive.
integrator = mm.LangevinMiddleIntegrator(0*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
seed=13
integrator.setRandomNumberSeed(seed)
print(f"Random seed: {seed}")
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}
simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setPositions(modeller.positions)

from openmm.app import PDBFile
PDBFile.writeFile(modeller.topology, modeller.positions, open('solvated_system.pdb', 'w'))
structure = pmd.openmm.load_topology(modeller.topology, system_for_export, modeller.positions)
structure.save('solvated_system.prmtop', format='amber')
structure.save('solvated_system.inpcrd', format='rst7')
print(".prmtop file created.")

# 5. initialize simulation: perform energy minimization. report energy before and after minimization. restrain heavy atoms by applying a penalty to their harmonic motion, because this allows for the water molecules to settle into the system and can help prevent sudden jumps in energy or unexpected changes in structure. may not be that necessary but doing it just in case for a more defensible method.
print("Step 1: Energy minimization and positional restraints...")
state_before = simulation.context.getState(getEnergy=True)
print(f"Potential energy before minimization: {state_before.getPotentialEnergy()}")
simulation.minimizeEnergy(maxIterations=10000)
state_after = simulation.context.getState(getEnergy=True)
print(f"Potential energy after minimization: {state_after.getPotentialEnergy()}")

print("Applying positional restraints for heating/equilibration...")
restraint = mm.CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
restraint.addGlobalParameter('k', 100.0 * unit.kilocalories_per_mole/unit.angstrom**2)
restraint.addPerParticleParameter('x0')
restraint.addPerParticleParameter('y0')
restraint.addPerParticleParameter('z0')

minimized_positions = simulation.context.getState(getPositions=True).getPositions()

restrained_atom_indices = []
for atom in modeller.topology.atoms():
    # Restrain protein and ligand heavy atoms; skip water/ions and all hydrogens
    if atom.residue.name in ('HOH', 'WAT', 'NA', 'CL', 'K'):
        continue
    if atom.element is not None and atom.element.symbol != 'H':
        restraint.addParticle(atom.index, minimized_positions[atom.index])
        restrained_atom_indices.append(atom.index)

system.addForce(restraint)
simulation.context.reinitialize(preserveState=True)
print(f"Restrained {len(restrained_atom_indices)} heavy atoms.")

# 6. gradual heating was chosen for a similar reason as heavy atom restraints. moving the system from its crystallized structure to 300 K abruptly can cause unfolding or sudden unwanted shifts in structure, so this allows the system to settle more gently.

print("Step 2: Gradual heating phase (0 K -> 300 K over 100 ps)...")
heating_steps = 100
steps_per_increment = 500  # 500 steps * 2 fs = 1 ps per temperature step
for i in range(heating_steps):
    current_temp = (i / heating_steps) * 300 * unit.kelvin
    integrator.setTemperature(current_temp)
    simulation.step(steps_per_increment)
    if i % 10 == 0:
        print(f"  -> Temperature reached: {current_temp.value_in_unit(unit.kelvin):.1f} K")
integrator.setTemperature(300*unit.kelvin)

# 7. add barostat (which sets the pressure of the system) and release restraints on heavy atoms gradually. do not add the barostat until the temperature is fully set after gradual heating because of pressure/temperature relationship. 
system.addForce(mm.MonteCarloBarostat(1.0*unit.bar, 300*unit.kelvin))
simulation.context.reinitialize(preserveState=True)

print("Releasing positional restraints gradually...")
release_steps = 5
initial_k = 100.0
for i in range(release_steps):
    scale = 1.0 - (i / (release_steps - 1))  # 1.0 -> 0.0
    simulation.context.setParameter('k', initial_k * scale * unit.kilocalories_per_mole/unit.angstrom**2)
    simulation.step(50000)  # 100 ps per stage at 2 fs/step
    print(f"  -> Restraint scale: {scale:.2f}")

print("System equilibrated at 300 K, restraints released. Starting 10 ns production run...")

# 8. set up reporters. downstream analysis such as RMSD and RMSF. will be done on the trajectory.dcd and .prmtop files. run simulation for 10 ns (5,000,000 steps * 0.002 ps). note that this is only a test run, and 100 ns+ (or doing replicates of 50 ns+) are strongly recommended for publishable results.
simulation.reporters.append(app.DCDReporter('trajectory.dcd', 10000)) 
simulation.reporters.append(app.CheckpointReporter('checkpoint.chk', 10000))
simulation.reporters.append(app.StateDataReporter(
    sys.stdout, 10000, step=True, potentialEnergy=True, 
    temperature=True, progress=True, totalSteps=5000000, speed=True
))

simulation.step(5000000)
print("Simulation complete.")