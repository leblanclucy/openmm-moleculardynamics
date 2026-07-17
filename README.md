OVERALL WORKFLOW
**protein.pdb and ligand.pdb -> openmmbuild.py -> unwrap.py -> mdanalysis.py**

Starting with either a crystal structure or an AlphaFold2-predicted structure of the enzyme (or receptor) as well as the substrate (or ligand) oriented in its top binding mode as predicted by AutoDock Vina or other molecular docking software, this allows you to run molecular dynamics with highly customizable parameters (including protonation states, force fields, water models, temperature) through OpenMM. During the molecular dynamics run, a .prmtop file (which describes the initial state of the solvated system) and a trajectory.dcd file (which describes how the complex moves from frame to frame) are created. 

These are then used for downstream analysis to obtain useful information about the complex, namely the following:
1) Root mean square deviation (RMSD) - how much the protein moves over time
2) Root mean square fluctuation (RMSF) - how much each individual residue moves over time, useful for showing stabilization of certain residues, also shows which regions of the protein are more disordered
3) Define secondary structure of proteins (DSSP) - assigns each residue to a particular secondary structure, such as alpha or pi helices, based on detecting backbone hydrogen bonds. some useful reading as well as possible area of improvement: Hekkelman et al., 2025 (https://doi.org/10.1002/pro.70208)
4) Radius of gyration (Rg) - a measurement of how compact the protein remains over time. a larger Rg indicates that parts of it are moving further and further from its axis of rotation.
5) Solvent-accessible surface area - the surface area of the protein that is accessible to solvent, in this case, water.
6) Hydrogen bond occupancy - given that both the enzyme (or receptor) and ligand (or substrate) has residues that may act as hydrogen bond donors or acceptors, what are the hydrogen bonds formed within the system, and what percent of the time are they (each donor/acceptor pair) close enough to be considered a hydrogen bond?

Future additions:
1) Need to add distance to scissile bond or other more targeted distances worth tracking.
2) Need to add molecular mechanics with Poisson Boltzmann and surface area solvation (MM-PBSA), an end point method that estimates affinity of the receptor for its ligand (old review by Genheden & Ryde in 2015 may be worth a read https://doi.org/10.1517/17460441.2015.1032936). Bonus points if I can add in alchemical perturbation methods as well like free energy perturbation.

Note that my choice of downstream analysis was largely inspired by Karaoli et al., 2026 (https://doi.org/10.1039/d6ra00343e) which did molecular dynamics on thermostable PETases that degrade plastics.
