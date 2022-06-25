#!/usr/bin/env python
import pandas as pd

#df = pd.read_csv('smiles_list.csv')
df = pd.read_pickle('smiles_list.pkl')
df2 = df[['SMILES', 'MOL_ID']]
fp = open('rdf_input.txt','w')
num_data = df.shape[0]
for i in range(num_data):
    df0 = df.iloc[i]
    smi0 = df0['SMILES']
    smi = smi0.strip().split()[0]
    mol_id = df0['MOL_ID']
    line_out = '%s %s\n' %(smi, mol_id)
    fp.write(line_out)
fp.close()
