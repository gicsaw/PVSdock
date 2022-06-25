#!/usr/bin/env python

import numpy as np
from openbabel import pybel


def find_hydrophobic_ring(m):
    smarts_pattern_list = ['*1***1', '*1****1', '*1*****1',
                           '*1******1', '*1*******1']
    atoms = m.atoms
    feature_coor_list = list()
    for smarts_pattern in smarts_pattern_list:
        smarts = pybel.Smarts(smarts_pattern)
        atom_idx_group_list = smarts.findall(m)
        print(atom_idx_group_list)
        for atom_idx_group in atom_idx_group_list:
            pseudo_atom = []
            for atom_idx in atom_idx_group:
                coor = np.array(atoms[atom_idx-1].coords)
                pseudo_atom += [coor]
            pseudo_atom = np.array(pseudo_atom)
            pseudo_atom_coor = pseudo_atom.mean(axis=0)
            feature_coor_list += [pseudo_atom_coor]
    return feature_coor_list


def main():

    pdb_chain_id = '3HMMA'
    ligand_id = '855'
    ref_ligand_file = 'fix/%s_%s.pdb' % (pdb_chain_id, ligand_id)
    ms = pybel.readfile('pdb', ref_ligand_file)
    m_ref = list(ms)[0]
    feature_coor_ref_list = find_hydrophobic_ring(m_ref)
    print(feature_coor_ref_list)

    cutoff = 1.5

    pharmacophore_coor_dict = dict()
    for i, feature_coor_ref in enumerate(feature_coor_ref_list):
        pharmacophore_coor_dict[i] = []

    list_file = 'list_final.txt'
    fp = open(list_file)
    lines = fp.readlines()
    fp.close()

    for line in lines:
        lis = line.strip().split()
        pdb_chain_id = lis[0]
        ligand_id = lis[1]
        ligand_file = 'fix/%s_%s.pdb' % (pdb_chain_id, ligand_id)
        ms = pybel.readfile('pdb', ligand_file)
        ms = list(ms)
        for m in ms:
            feature_coor_list = find_hydrophobic_ring(m)
#            print(feature_coor_list)
            for feature_coor in feature_coor_list:
                for i, feature_coor_ref in enumerate(feature_coor_ref_list):
                    diff = np.linalg.norm(feature_coor_ref - feature_coor)
#                    print(i, feature_coor_ref, diff)
                    if diff <= cutoff:
                        pharmacophore_coor_dict[i] += [feature_coor]

    for i in pharmacophore_coor_dict:
        pharmacophore_coor = pharmacophore_coor_dict[i]
#        print(pharmacophore_coor)
        print(np.array(pharmacophore_coor).mean(axis=0))


if __name__ == '__main__':
    main()
