# given a log file after executing Autodock GPU on a folder containing small molecules, find the names and affinity of each compound, sort by most negative

import re
import pandas as pd

ligand_pattern = re.compile(r'Ligand file:\s+(\S+)')
score_pattern = re.compile(r'[Bb]est inter \+ intra\s+(-?\d+\.\d+)\s+kcal/mol')

results = []
current_ligand = None

with open('adgpu_results.log') as f:
    for line in f:
        lig_match = ligand_pattern.search(line)
        if lig_match:
            current_ligand = lig_match.group(1)
            continue
        score_match = score_pattern.search(line)
        if score_match and current_ligand is not None:
            results.append((current_ligand, float(score_match.group(1))))
            current_ligand = None

df = pd.DataFrame(results, columns=['Ligand', 'adgpu_Affinity_kcal_mol'])

df['Compound'] = df['Ligand'].str.extract(r'([^/]+)\.pdbqt$')

df = df[['Compound', 'adgpu_Affinity_kcal_mol', 'Ligand']]
df = df.sort_values('adgpu_Affinity_kcal_mol')

df.to_csv('adgpu_results.csv', index=False)
print(f"Parsed {len(df)} results -> adgpu_results_sorted.csv")
print(df.head(10).to_string(index=False))
