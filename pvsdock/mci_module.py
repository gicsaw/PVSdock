import numpy as np
from mci import mcifinder


class mci_class(object):
    def __init__(self, docking_params):
        self.use_mciscore = docking_params['use_mciscore']
        if self.use_mciscore:
            self.pharma = docking_params['pharma']
            self.use_mcinfo = self.pharma.use_mcinfo
            self.vmciscore_w1 = docking_params['vmciscore_w1']
            self.vmciscore_w2 = docking_params['vmciscore_w2']

    def simulation_process(self, idx, mol_id, smi, pid, smi_p,
                           out_dock_dir1, docking_pdb_file, result_dict):
        if self.use_mciscore:
            if self.output_save:
                pf_ligand_file = '%s/pf_%s.txt' % (out_dock_dir1, mol_id)
                interaction_file = '%s/mcinter_%s.txt' % (
                    out_dock_dir1, mol_id)
            else:
                pf_ligand_file = None
                interaction_file = None
            try:
                result_mciscore = self.pharma.cal_mciscore(
                    docking_pdb_file, pf_ligand_file, interaction_file)
                total_count_dict, total_count_info_dict = result_mciscore
                keys = list(total_count_dict.keys())
                mciscore = list()
                if self.use_mcinfo:
                    mcinfo = list()
                for key in keys:
                    count_dict = total_count_dict[key]
                    mcis = count_dict['Score']
                    mciscore += [mcis]
                    if self.use_mcinfo:
                        count_info_dict = total_count_info_dict[key]
                        pin = count_info_dict['Score']
                        mcinfo += [pin]
                mciscore = np.array(mciscore, dtype=np.float32)
                if self.use_mcinfo:
                    mcinfo = np.array(mcinfo, dtype=np.float32)

            except Exception as e:
                print(e, 'mciscore', idx, mol_id, smi_p, flush=True)
                mciscore = np.array([0], dtype=np.float32)
                if self.use_mcinfo:
                    mcinfo = np.array([0], dtype=np.float32)

            result_dict['mciscore'] = mciscore
            if self.rescoring:
                docking_rescore = result_dict['docking_re']
                vmciscore = (self.vmciscore_w1 * (-docking_rescore)
                             + (1.0 - self.vmciscore_w1) * mciscore)
            else:
                docking_score = result_dict['docking']
                vmciscore = (self.vmciscore_w1 * (-docking_score)
                             + (1.0 - self.vmciscore_w1) * mciscore)

            result_dict['vmciscore'] = vmciscore

            if self.use_mcinfo:
                result_dict['mcinfo'] = mcinfo
                if self.rescoring:
                    docking_rescore = result_dict['docking_re']
                    vmciscore_info = (self.vmciscore_w1 * (-docking_rescore)
                                      + (1.0 - self.vmciscore_w1)
                                      * (self.vmciscore_w2 * mciscore
                                         + (1.0 - self.vmciscore_w1) * mcinfo))
                else:
                    docking_score = result_dict['docking']
                    vmciscore_info = (self.vmciscore_w1 * (-docking_score)
                                      + (1.0 - self.vmciscore_w1)
                                      * (self.vmciscore_w2 * mciscore
                                         + (1.0 - self.vmciscore_w1) * mcinfo))

                result_dict['vmciscore_info'] = vmciscore_info

    def predict(self, smiles_list, result_dict, return_dict):
        data = list(enumerate(smiles_list))
        num_data = len(data)

        keys = sorted(return_dict.keys())
        if self.use_mciscore:
            mciscore_list = list()
            vmciscore_list = list()
            if self.use_mcinfo:
                mcinfo_list = list()
                vmciscore_info_list = list()

        for key in range(num_data):
            if key in keys:
                result_dict0 = return_dict[key]
                if self.use_mciscore:
                    if 'mciscore' in result_dict0:
                        mciscore = result_dict0['mciscore']
                    else:
                        mciscore = np.array([0], dtype=np.float32)

                    if 'vmciscore' in result_dict0:
                        vmciscore = result_dict0['vmciscore']
                    else:
                        vmciscore = np.array([0], dtype=np.float32)
                    if self.use_mcinfo:
                        if 'mcinfo' in result_dict0:
                            mcinfo = result_dict0['mcinfo']
                        else:
                            mcinfo = np.array([0], dtype=np.float32)
                        if 'vmciscore_info' in result_dict0:
                            vmciscore_info = result_dict0['vmciscore_info']
                        else:
                            vmciscore_info = np.array([0], dtype=np.float32)
            else:
                if self.use_mciscore:
                    mciscore = np.array([0], dtype=np.float32)
                    vmciscore = np.array([0], dtype=np.float32)
                    if self.use_mcinfo:
                        mcinfo = np.array([0], dtype=np.float32)
                        vmciscore_info = np.array([0], dtype=np.float32)
            if self.use_mciscore:
                mciscore_list += [mciscore]
                vmciscore_list += [vmciscore]
                if self.use_mcinfo:
                    mcinfo_list += [mcinfo]
                    vmciscore_info_list += [vmciscore_info]
        if self.use_mciscore:
            result_dict['mciscore_list'] = mciscore_list
            result_dict['vmciscore_list'] = vmciscore_list
            if self.use_mcinfo:
                result_dict['mcinfo_list'] = mcinfo_list
                result_dict['vmciscore_info_list'] = vmciscore_info_list
        return result_dict


def parser_arg(parser):
    parser.add_argument('--mciscore_receptor', required=False,
                        default=None,
                        help='input receptor pdb file for mciscore')
    parser.add_argument('--pf_receptor', required=False,
                        default=None,
                        help='output for pharmacophoric feature of receptor')
    parser.add_argument('--mcinfo_ligand', required=False,
                        default=None, help='template ligand file for mcinfo')
    parser.add_argument('--pf_receptor_info', required=False, default=None,
                        help='output for template feature of receptor')
    parser.add_argument('--mci_cutoff', type=float, default=6.5, required=False,
                        help='pharmacophoric interaction cutoff distance')
    parser.add_argument('--vmciscore_w1', type=float, required=False,
                        default=0.2,
                        help='weight_1 of VMCIscore\n' +
                        'VMCIscore = w1*vina_score + (1-w1)*mciscore')
    parser.add_argument('--vmciscore_w2', type=float, required=False,
                        default=0.8,
                        help='weight_2 of VMCIscore_info\n' +
                        'VMCIscore_info = w1*vina_score + (1-w1)*(w2*mciscore'
                        + '(1-w2)*mcinfo)')
    parser.add_argument('--include_hydrophobic', action='store_true',
                        required=False,
                        help='include hydrophobic feature for template')


def arg_to_params(parser, docking_params):

    args = parser.parse_args()
    use_mciscore = False
    if args.mcinfo_ligand is not None or args.pf_receptor_info is not None:
        use_mciscore = True
        pharma = mcifinder.set_mcifinder(args)
        docking_params['use_mciscore'] = use_mciscore
        docking_params['pharma'] = pharma
        docking_params['use_mcinfo'] = pharma.use_mcinfo
        docking_params['vmciscore_w1'] = args.vmciscore_w1
        docking_params['vmciscore_w2'] = args.vmciscore_w2
    return docking_params


def mci_score_to_df(df, docking_params, result_dict):
    if docking_params['use_mciscore']:
        mciscore = result_dict['mciscore_list']
        vmciscore = result_dict['vmciscore_list']
        vmciscore_1 = [x.max() for x in vmciscore]
        df['VMCIscore1'] = vmciscore_1

        if docking_params['use_mcinfo']:
            mcinfo = result_dict['mcinfo_list']
            vmciscore_info = result_dict['vmciscore_info_list']
            vmciscore_info_1 = [x.max() for x in vmciscore_info]

            df['VMCIscore_info1'] = vmciscore_info_1

        df['VMCIscore'] = vmciscore
        if docking_params['use_mcinfo']:
            df['VMCIscore_info'] = vmciscore_info
        df['MCIscore'] = mciscore
        if docking_params['use_mcinfo']:
            df['MCinfo'] = mcinfo
