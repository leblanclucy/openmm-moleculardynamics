from rdkit import Chem
from rdkit.Chem import AllChem

#load the smiles and add hydrogens
smi = "[SMILES HERE]"
molec = Chem.AddHs(Chem.MolFromSmiles(smi))

#generate a conformer using the most recent stable version (v3) couldn't get v4 to work
params = AllChem.ETKDGv3()
params.randomSeed = 13
status = AllChem.EmbedMolecule(molec, params)

#failed if it returns -1, succeeded if returns 0
print("embed status:", status)

#optimize geometry with 5000 iterations
AllChem.MMFFOptimizeMolecule(molec, maxIters=5000)
Chem.MolToMolFile(molec, "6pg_pH74.sdf")
