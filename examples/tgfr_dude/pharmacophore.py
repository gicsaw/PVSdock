#!/usr/bin/env python
import os
import numpy as np
from openbabel import pybel
import pandas as pd
from multiprocessing import Manager
from multiprocessing import Process
from multiprocessing import Queue


def find_hydrophobic_ring(m):
    smarts_pattern_list = ['a1:a:a:a:1', 'a1:a:a:a:a1', 'a1:a:a:a:a:a:1',
                           'a1:a:a:a:a:a:a:1', 'a1:a:a:a:a:a:a:a:1']
#    smarts_pattern_list = ['*1***1', '*1****1', '*1*****1',
#                           '*1******1', '*1*******1']

    atoms = m.atoms
    feature_coor_list = list()
    for smarts_pattern in smarts_pattern_list:
        smarts = pybel.Smarts(smarts_pattern)
        atom_idx_group_list = smarts.findall(m)
        for atom_idx_group in atom_idx_group_list:
            pseudo_atom = []
            for atom_idx in atom_idx_group:
                coor = np.array(atoms[atom_idx - 1].coords)
                pseudo_atom += [coor]
            pseudo_atom = np.array(pseudo_atom)
            pseudo_atom_coor = pseudo_atom.mean(axis=0)
            feature_coor_list += [pseudo_atom_coor]
    return feature_coor_list


def creator(q, data, num_sub_proc):
    for d in data:
        idx = d[0]
        q.put((idx, d[1]))

    for i in range(0, num_sub_proc):
        q.put('DONE')


def worker(q, pharmacophore_coor_ref_list, cutoff, return_dict):
    # pid = os.getpid()
    while True:
        qqq = q.get()
        if qqq == 'DONE':
            # print('proc =', os.getpid())
            break
        (idx, mol_id) = qqq

        ligand_file = 'docking/%s/dock_%s.pdb' % (mol_id[0:7], mol_id)
        if not os.path.exists(ligand_file):
            return_dict[idx] = 1.0
            continue
        ms = pybel.readfile('pdb', ligand_file)
        ms = list(ms)
        for m in ms[0:1]:
            feature_coor_list = find_hydrophobic_ring(m)
#            print(feature_coor_list)
            prediction = 0.0
            count = 0
            for i in range(3):
                check = False
                for ref_coor in pharmacophore_coor_ref_list[i]:
                    dif_min = 100
                    for j, feature_coor in enumerate(feature_coor_list):
                        diff = np.linalg.norm(ref_coor - feature_coor)
                        if dif_min > diff:
                            dif_min = diff
                    if dif_min < cutoff:
                        check = True
                if check:
                    count += 1
            if count >= 3:
                prediction = 1.0
        return_dict[idx] = prediction


def main():

    num_sub_proc = 256
    data_file = 'workflow/master/docking.pkl'
    df = pd.read_pickle(data_file)

    cutoff = 2.5
    pharmacophore_coor_ref_list = [[[12.97766667, 64.211, 4.3995]],
                                   [[17.38733333, 68.97966667,  2.449]],
                                   [[15.542, 66.681,  6.613666],
                                       [16.22216667, 67.62283333,  8.757]]]

    num_data = df.shape[0]

    data = list()
    for idx in range(num_data):
        df0 = df.iloc[idx]
        mol_id = df0['MOL_ID']
#        activity = df0['Activity']
        data += [[idx, mol_id]]

    num_data = len(data)
    num_sub_proc = min(num_sub_proc, num_data)

    q1 = Queue()
    manager = Manager()
    return_dict = manager.dict()
    proc_master = Process(target=creator,
                          args=(q1, data, num_sub_proc))
    proc_master.start()

    procs = []
    for sub_id in range(0, num_sub_proc):
        proc = Process(target=worker, args=(q1, pharmacophore_coor_ref_list,
                                            cutoff, return_dict))
        procs.append(proc)
        proc.start()

    q1.close()
    q1.join_thread()
    proc_master.join()
    for proc in procs:
        proc.join()
#    keys = sorted(return_dict.keys())

    prediction_list = list()
    for idx in range(num_data):
        if idx not in return_dict:
            prediction_list += [0.0]
            continue
        prediction = return_dict[idx]
        prediction_list += [prediction]

    df['pharmacophore'] = prediction_list

    df.to_pickle('pharmacophore.pkl')


if __name__ == '__main__':
    main()
