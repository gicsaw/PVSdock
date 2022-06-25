#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
#import cairosvg

def main():

#    list_file = 'select/result_list.csv'
    list_file = 'select/list.txt'
    df = pd.read_csv(list_file,sep='\s+')
#    'MOL_ID,model_id,docking,MCIscore_info,mciscore,mcinfo,SMILES,MCF,MW,LogP,HBD,HBA,TPSA,cyp1a2,cyp2c9,cyp2c19,cyp2d6,cyp3a4,herg'

    size=(500,500)

    fig_dir="fig_mol"
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    num_data = df.shape[0]
    for idx in range(num_data):
        df0 = df.iloc[idx]
        mol_id = df0['MOL_ID']
        smiles = df0['SMILES']
        vmciscore_info = df0['MCIscore_info']
#        MCF = df0['MCF']
#        fig_file=fig_dir+"/M%d_F_%s.svg" %(idx+1,mol_id)
        png_file=fig_dir+"/M%d_F_%s.png" %(idx+1,mol_id)
        mol=Chem.MolFromSmiles(smiles)
        Draw.MolToFile(mol, png_file, size=size)



if __name__ == "__main__":
    main()
