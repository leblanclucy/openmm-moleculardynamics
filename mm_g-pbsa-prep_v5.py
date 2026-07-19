import os
import subprocess
import parmed as pmd
from parmed.tools.changeradii import mbondi2

# 1. load the original complex topology (.prmptop file, created during MD run) and coordinates
complex_top = pmd.load_file('solvated_system.prmtop')

# set intrinsic radii to mbondi2. this prevents a crash later with MMPBSA.py.
print("Changing intrinsic radii to mbondi2...")
mbondi2(complex_top)

# define target selection strings so that you can remove the solvent and ligand from the system, isolating just the enzyme.
water_ions = ':HOH,WAT,SOL,NA,CL,K'
ligand_res = ':ROH,4YB,0YB'

# 2. extract the enzyme ONLY, excluding the substrate, water, or ions.
enzyme_top = complex_top[:]
enzyme_top.strip(f'{water_ions},{ligand_res}')
enzyme_top.box = None  # Clear periodic boundary box info
enzyme_top.save('enzyme.prmtop', overwrite=True)
print("Enzyme extracted from topology.")

# 3. extract the substrate ONLY. change "1-294" to reflect the number of amino acids in the enzyme.
substrate_top = complex_top[:]
substrate_top.strip(f':1-294,{water_ions}')  # removes the enzyme and solvent
substrate_top.box = None  # clear periodic boundary box info to prevent a crash later on
substrate_top.save('substrate.prmtop', overwrite=True)
print("Substrate extracted from topology.")

# 4. extract the complex but without any water or ions.
dry_complex_top = complex_top[:]
dry_complex_top.strip(water_ions)
dry_complex_top.box = None  # Clear periodic boundary box info
dry_complex_top.save('complex.prmtop', overwrite=True)
print("Dry complex extracted from topology.")

# 5. use cpptraj to remove the solvent from the trajectory file, which was centered using unwrap.py.
print("Removing water and ions from trajectory...")
cpptraj_input = 'cpptraj_strip.in'
cpptraj_instructions = f"""parm solvated_system.prmtop
trajin trajectory_centered.dcd
strip {water_ions}
trajout trajectory_no-solvent.nc netcdf
run
quit
"""

with open(cpptraj_input, 'w') as f:
    f.write(cpptraj_instructions)

subprocess.run(['cpptraj', '-i', cpptraj_input], stdout=subprocess.DEVNULL)

if os.path.exists(cpptraj_input):
    os.remove(cpptraj_input)
print("Trajectory cleaned successfully (saved as trajectory_no-solvent.nc).")