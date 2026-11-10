# for doing grid inhomogeneous solvation theory (GIST) to measure solvation entropy
# need to also run a box with just waters since i'm using OPC3
# also do apo runs without ligand, but start from equilibrated holo structure and just delete the ligand
import numpy as np
import mdtraj as md
import openmm as mm
import openmm.app as app
import openmm.unit as unit
# this is for preventing an error that can happen when starting a new trajectory from an existing frame
from openmm.app import AmberPrmtopFile
from openmm.app.internal.unitcell import reducePeriodicBoxVectors


# set random seed and temperature
seed = 13
T = 300 * unit.kelvin

# 1. load last frame of an equilibrated run using the .dcd file
with md.open('trajectory_centered.dcd') as f:
    n_frames = len(f)
frame = md.load_frame('trajectory_centered.dcd', n_frames - 1, top='solvated_system.prmtop')   # 0-indexed

pos = frame.xyz[0]                        # nm, shape (N, 3)
box = frame.unitcell_vectors[0]           # nm, shape (3, 3)

# restrain heavy atoms for apo or holo
heavy = frame.topology.select('(protein or resname 4YB or resname 0YB or resname ROH) and not element H')
print(len(heavy), sorted({frame.topology.atom(i).residue.name for i in heavy}))
# apo: ('protein and not element H')
# holo: '(protein or resname 4YB or resname 0YB or resname ROH) and not element H'

# 2. build the system using the prmtop file as a starting point
prmtop = app.AmberPrmtopFile('solvated_system.prmtop')

# reduce the vectors of the box to prevent an error later
box_vectors = reducePeriodicBoxVectors(prmtop.topology.getPeriodicBoxVectors())
system = prmtop.createSystem(nonbondedMethod=app.PME,
                             nonbondedCutoff=1.0 * unit.nanometer,
                             constraints=app.HBonds)

# add positional restraints on heavy atoms
k_kcal = 10.0                            # kcal/mol/A^2
k = (k_kcal * unit.kilocalories_per_mole / unit.angstrom**2).value_in_unit(
        unit.kilojoules_per_mole / unit.nanometer**2)
restraint = mm.CustomExternalForce('k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)')
restraint.addGlobalParameter('k', k)
for name in ('x0', 'y0', 'z0'):
    restraint.addPerParticleParameter(name)
for i in heavy:
    restraint.addParticle(int(i), pos[i].tolist())
restraint.setForceGroup(1)
system.addForce(restraint)

# add the barostat
barostat = mm.MonteCarloBarostat(1.0 * unit.bar, T)
barostat.setRandomNumberSeed(seed)
system.addForce(barostat)

# add the thermostat
integrator = mm.LangevinMiddleIntegrator(T, 1 / unit.picosecond, 0.002 * unit.picoseconds)
integrator.setRandomNumberSeed(seed)

# initalize GPU platform
platform = mm.Platform.getPlatformByName('HIP')
properties = {'HipPrecision': 'mixed'}
sim = app.Simulation(prmtop.topology, system, integrator, platform, properties)
print(sim.context.getPlatform().getName(), sim.context.getPlatform().getPropertyValue(sim.context, 'HipPrecision'))

# uses reduced box vectors to prevent crash
sim.context.setPeriodicBoxVectors(*box_vectors)
sim.context.setPositions(pos * unit.nanometer)
sim.context.setVelocitiesToTemperature(T, seed)

# 3. allow equilibration for 5 ns, don't save any frames
sim.reporters.append(app.StateDataReporter('gist_equil13.log', 5000, step=True, time=True,
                     potentialEnergy=True, temperature=True, volume=True, density=True))
sim.step(2_500_000)                       # 5 ns at 2 fs per step

r = sim.context.getState(getEnergy=True, groups={1}).getPotentialEnergy()
print('restraint energy after equilibration:', r)

# disable barostat, initialize NVT simulation
# in a previous paper i read they used OpenMM NVT for GIST (Wallerstein et al., 2021)
barostat.setFrequency(0)
sim.context.reinitialize(preserveState=True)

# 4. production MD run for GIST
sim.reporters.clear()
# save every 1000 steps or 2 ps (1000 steps at 2 fs/step)
sim.reporters.append(app.DCDReporter('gist_seed13.dcd', 1000))
sim.reporters.append(app.CheckpointReporter('gist_seed13.chk', 500_000))
sim.reporters.append(app.StateDataReporter('gist_seed13.log', 5000, step=True, time=True,
                     potentialEnergy=True, temperature=True, volume=True, density=True, progress=True, 
                     totalSteps=50000000, speed=True))
# run for 100 ns
sim.step(50_000_000)                      
sim.saveState('gist_final.xml')
