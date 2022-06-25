#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
from sklearn import metrics
from matplotlib import pyplot as plt
import copy


def cal_ef(d, p, rr=0.01):
    num_data = d.shape[0]
    num_act = d.sum()
    num_pos = int(num_data*rr)
    cutoff = -sorted(-p)[num_pos]
    idx_act = d == 1
    idx_pos = p >= cutoff

    idx_tp = idx_act*idx_pos

    num_tp = idx_tp.sum()
    precision = num_tp/num_pos
    active_ratio = num_act/num_data

    ef = precision/active_ratio
    return ef


def main():

    #    data_file = 'workflow/master/docking.pkl'
    data_file = 'pharmacophore.pkl'

    df = pd.read_pickle(data_file)

    num_data = df.shape[0]

    data = list()
    for idx in range(num_data):
        df0 = df.iloc[idx]
        mol_id = df0['MOL_ID']
        activity = df0['Activity']
        docking = df0['Docking1']
        vmciscore = df0['VMCIscore1']
        vmciscore_info = df0['VMCIscore_info1']
        pharmacophore = df0['pharmacophore']

        data += [[activity, docking, vmciscore, vmciscore_info, pharmacophore]]

    num_data = len(data)
    data = np.array(data)
#    fp_out=open("docking.txt","w")
#    for i in range(len(data)):
#        line_out="%.3f %.3f \n" %(data[i,0],data[i,1])
#        fp_out.write(line_out)
#    fp_out.close()
    idx = data[:, 0] == 1
    act = data[:, 0][idx]
    inact = data[:, 1][~idx]
    print(act.shape, inact.shape)
    bins = np.arange(-12, -4, 0.1)
#    plt.hist(inact, bins=bins, density=True, color='blue')
#    plt.hist(act, bins=bins, density=True, color='red')

#    plt.show()

    fpr, tpr, thresholds = metrics.roc_curve(data[:, 0], -data[:, 1])
    roc_auc = metrics.roc_auc_score(data[:, 0], -data[:, 1])
    ef = cal_ef(data[:, 0], -data[:, 1])
    print(roc_auc, ef)

    roc_auc2 = metrics.roc_auc_score(data[:, 0], data[:, 2])
    ef2 = cal_ef(data[:, 0], data[:, 2])
    print(roc_auc, ef2)

    roc_auc3 = metrics.roc_auc_score(data[:, 0], data[:, 3])
    ef3 = cal_ef(data[:, 0], data[:, 3])
    print(roc_auc, ef3)

    pp2 = copy.copy(data[:, 1])
    idx_p = data[:, 4] == 0.0
    pp2[idx_p] += 3.0
    fpr4, tpr4, thresholds4 = metrics.roc_curve(data[:, 0], -pp2)
    roc_auc4 = metrics.roc_auc_score(data[:, 0], -pp2)
    ef4 = cal_ef(data[:, 0], -pp2)
    print(roc_auc4, ef4)

    fp_out = open("roc.txt", "w")
    for i in range(len(fpr)):
        line_out = "%6.2f %6.4f %6.4f\n" % (thresholds[i], fpr[i], tpr[i])
        fp_out.write(line_out)
    fp_out.close()

    plt.rc('font', size=12)
    plt.rc('axes', titlesize=16, labelsize=16)
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('legend', fontsize=12)
    plt.rc('figure', titlesize=24)

    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.title('ROC')
    plt.plot(fpr, tpr)
    plt.plot(fpr4, tpr4)

    plt.plot([0, 1], [0, 1], '-')

    plt.legend(['vina_%.2f_%.2f' % (roc_auc, ef),
               'vina+a_%.2f_%.2f' % (roc_auc4, ef4)])
#    plt.show()
    plt.savefig('TGFR_DUDE_ROC.png', dpi=300)


#    activity = df['Activity']
#    prediction = df['fragment']

#    idx_act = activity == 1
#    idx_pos = prediction == 1
#    idx_tp = idx_act*idx_pos
#    num_act = idx_act.sum()
#    num_pos = idx_pos.sum()
#    num_tp = idx_tp.sum()
#    print(num_tp, num_pos, num_act, num_data)
#    print(num_tp/num_pos)
#    print((num_tp/num_pos)/(num_act/num_data))

if __name__ == '__main__':
    main()
