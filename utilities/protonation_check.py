# checks protonation state of your prmtop file compared to the output from propka
import parmed as pmd
parm = pmd.load_file("solvated_system.prmtop")

prot = sum(a.charge for a in parm.atoms if a.residue.name not in {"HOH", "NA", "CL"})
print("protein+ligand charge:", round(prot, 3))

for r in parm.residues:
    if r.name in {"HIS", "HID", "HIE", "HIP"}:
        print(r.idx + 1, r.name, sorted(a.name for a in r.atoms if a.name in {"HD1", "HE2"}))

ss = [(b.atom1.residue.idx + 1, b.atom2.residue.idx + 1) for b in parm.bonds
      if b.atom1.name == "SG" and b.atom2.name == "SG"]
print("disulfides:", ss)

pdb = pmd.load_file("protein_pH7.pdb")

prot2 = sum(a.charge for a in pdb.atoms if a.residue.name not in {"HOH", "NA", "CL"})
print("protein+ligand charge:", round(prot2, 3))

for r in pdb.residues:
    if r.name in {"HIS", "HID", "HIE", "HIP"}:
        print(r.idx + 1, r.name, sorted(a.name for a in r.atoms if a.name in {"HD1", "HE2"}))

# not sure why this showed zero disulfide bonds as they definitely showed in openmm and they also show in the prmtop file. was more concerned about protonation anyway
ss2 = [(b.atom1.residue.idx + 1, b.atom2.residue.idx + 1) for b in pdb.bonds
      if b.atom1.name == "SG" and b.atom2.name == "SG"]
print("disulfides:", ss2)