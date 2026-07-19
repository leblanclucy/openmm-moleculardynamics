OVERALL WORKFLOW
**protein.pdb and ligand.pdb -> openmmbuild.py -> unwrap.py -> mdanalysis.py AND/OR mm_g-pbsa-prep.py**

Starting with either a crystal structure or an AlphaFold2-predicted structure of the enzyme (or receptor) as well as the substrate (or ligand) oriented in its top binding mode as predicted by AutoDock Vina or other molecular docking software, this allows you to run molecular dynamics with highly customizable parameters (including protonation states, force fields, water models, temperature) through OpenMM. During the molecular dynamics run, a .prmtop file (which describes the initial state of the solvated system) and a trajectory.dcd file (which describes how the complex moves from frame to frame) are created. 

These are then used for downstream analysis to obtain useful information about the complex. mdanalysis.py gives the following:
1) Root mean square deviation (RMSD) - how much the protein moves over time
2) Root mean square fluctuation (RMSF) - how much each individual residue moves over time, useful for showing stabilization of certain residues, also shows which regions of the protein are more disordered
3) Define secondary structure of proteins (DSSP) - assigns each residue to a particular secondary structure, such as alpha or pi helices, based on detecting backbone hydrogen bonds. some useful reading as well as possible area of improvement: Hekkelman et al., 2025 (https://doi.org/10.1002/pro.70208)
4) Radius of gyration (Rg) - a measurement of how compact the protein remains over time. a larger Rg indicates that parts of it are moving further and further from its axis of rotation.
5) Solvent-accessible surface area - the surface area of the protein that is accessible to solvent, in this case, water.
6) Hydrogen bond occupancy - given that both the enzyme (or receptor) and ligand (or substrate) has residues that may act as hydrogen bond donors or acceptors, what are the hydrogen bonds formed within the system, and what percent of the time are they (each donor/acceptor pair) close enough to be considered a hydrogen bond?

In addition, you may also run mm_g-pbsa-prep.py. This takes a .prmtop file as input (generated during the molecular dynamics run in addition to the solvated_system.pdb file). It isolates the enzyme alone, substrate alone, and complex, all without water or ions. It also removes water from the trajectory_centered.dcd file. These components are then used for MM-PBSA and/or MM-GBSA, which estimates the binding energy of the complex (more specifically, ΔG). You may read more about these methods at the following references: Genheden & Ryde in 2015 (https://doi.org/10.1517/17460441.2015.1032936), Yau et al. in 2024 (https://doi.org/10.1007/s00894-024-06189-4), or Tuccinardi in 2021 (https://doi.org/10.1080/17460441.2021.1942836).
8) Molecular mechanics/Generalized Born Surface Area (MM-GBSA) is quicker and easier to calculate because it runs through OpenMM using Python script, mm-gbsa.py. Both methods go through the trajectory frame by frame determining energy and ultimately calculating ΔGbinding = ΔGcomplex - ΔGenzyme - ΔGsubstrate.
9) Molecular mechanics/Poisson Boltzmann Surface Area (MM-PBSA) is more computationally expensive and can run into some errors because it's done through MMPBSA.py in Amber, which had some compatibiliy/formatting issues with the workflow I had earlier that I needed to fix, such as removing the periodiic boundary box information. You need to set up a mmpbsa.in file and run it through the command line. However, it is considered to be more accurate by some authors and it gives more detailed information that is helpful for troubleshooting, like showing the van der Waals vs. electrostatic contributions to total energy.
10) The choice of which one is better than the other depends on your system. It's recommended you try both. A similar trend across conditions gives you more confidence in the validity of your system.
11) It's possible to look at the enzyme residue by residue and see which residue contributes most to ΔG but this can take hours whereas taking an overall snapshot takes minutes. Will add in the future.

Future additions:
1) Need to add distance to scissile bond or other more targeted distances worth tracking.
2) Add per residue decomposition for MM-PBSA and MM-GBSA.
3) Alchemical perturbation methods as well like free energy perturbation (FEP).

Note that my choice of downstream analysis was largely inspired by Karaoli et al., 2026 (https://doi.org/10.1039/d6ra00343e) which did molecular dynamics on thermostable PETases that degrade plastics.
