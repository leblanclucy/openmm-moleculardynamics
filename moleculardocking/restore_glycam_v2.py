# purpose: given the top binding pose, map the heavy atoms back to the reference file from glycamweb. this restores the metadata, namely the IDs of the subunits (ROH/4YB/0YB).

import sys
try:
    from rdkit import Chem
except ImportError:
    print("Please install rdkit first: pip install rdkit")
    sys.exit(1)

# 1. load both files, removing hydrogens to normalize the heavy-atom backbone
print("Reading molecules...")
template = Chem.MolFromPDBFile("GlcNAc4.pdb", removeHs=True)
docked = Chem.MolFromPDBFile("glcnac4_pose1.pdb", removeHs=True)

if template is None or docked is None:
    print("Error: Could not read files.")
    sys.exit(1)

# 2. map the docked heavy atoms to the template using graph topology
print("Mapping atom graphs...")
mapping = template.GetSubstructMatch(docked)

if not mapping or len(mapping) != template.GetNumAtoms():
    print("Error: Topological match failed.")
    sys.exit(1)

# 3. transfer the coordinates from the docked file into the GLYCAM template
template_conf = template.GetConformer()
docked_conf = docked.GetConformer()

for docked_idx, template_idx in enumerate(mapping):
    position = docked_conf.GetAtomPosition(docked_idx)
    template_conf.SetAtomPosition(template_idx, position)

# 4. save the heavy atom structure. may wish to save as sdf to preserve bond order.
output_filename = "docked_glycam_ready.pdb"
Chem.MolToPDBFile(template, output_filename)
print(f"Success! Saved coordinate-mapped structure to {output_filename}")
