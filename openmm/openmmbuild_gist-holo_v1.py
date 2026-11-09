# for doing grid inhomogeneous solvation theory (GIST) to measure solvation entropy
# need to also run a box with just waters since i'm using OPC3
# also do apo runs without ligand, but start from equilibrated holo structure and just delete the ligand
import numpy as np
import mdtraj as md
import openmm as mm
import openmm.app as app
import openmm.unit as unit

# set random seed and temperature
seed = 13
T = 300 * unit.kelvin

# 1. load last frame of an equilibrated run using the .dcd file
with md.open('prod.dcd') as f:
    n_frames = len(f)
frame = md.load_frame('prod.dcd', n_frames - 1, top='system.prmtop')   # 0-indexed

pos = frame.xyz[0]                        # nm, shape (N, 3)
box = frame.unitcell_vectors[0]           # nm, shape (3, 3)

# restrain heavy atoms for apo or holo
heavy = frame.topology.select('protein and not element H')
# holo: '(protein or resname 4YB or resname 0YB or resname ROH) and not element H'

# 2. build the system using the prmtop file as a starting point
prmtop = app.AmberPrmtopFile('system.prmtop')
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

# NPT: add the barostat
barostat = mm.MonteCarloBarostat(1.0 * unit.bar, T)
barostat.setRandomNumberSeed(seed)
system.addForce(barostat)

# add the thermostat
integrator = mm.LangevinMiddleIntegrator(T, 1 / unit.picosecond, 0.002 * unit.picoseconds)
integrator.setRandomNumberSeed(seed)
sim = app.Simulation(prmtop.topology, system, integrator)

sim.context.setPeriodicBoxVectors(*(v * unit.nanometer for v in box))   # box before positions
sim.context.setPositions(pos * unit.nanometer)
sim.context.setVelocitiesToTemperature(T, seed)

# 3. allow equilibration for 5 ns, don't save any frames
sim.reporters.append(app.StateDataReporter('gist_equil.log', 5000, step=True, time=True,
                     potentialEnergy=True, temperature=True, volume=True, density=True))
sim.step(2_500_000)                       # 5 ns at 2 fs per step

r = sim.context.getState(getEnergy=True, groups={1}).getPotentialEnergy()
print('restraint energy after equilibration:', r)

# 4. production MD run for GIST
sim.reporters.clear()
sim.reporters.append(app.DCDReporter('gist_prod.dcd', 500))   # 1 ps; store every frame
sim.reporters.append(app.StateDataReporter('gist_prod.log', 5000, step=True,
                     potentialEnergy=True, temperature=True, volume=True, density=True))
sim.step(50_000_000)                      # 100 ns
sim.saveState('gist_final.xml')
