import glob
import pandas as pd

rows = []
for f in sorted(glob.glob('results_candidates/*_out.pdbqt')):
    compound_name = f.split('/')[-1].replace('_out.pdbqt', '')
    with open(f) as fh:
        for line in fh:
            if line.startswith('REMARK VINA RESULT'):
                affinity = float(line.split()[3])  # best pose only
                rows.append({'Compound': compound_name, 'affinity_kcal_per_mol': affinity})
                break

df = pd.DataFrame(rows)
df = df.sort_values('affinity_kcal_per_mol')  # ascending: most negative (strongest binder) first
df.to_csv('results_candidates/summary_scores.csv', index=False)
print(df.head(20))