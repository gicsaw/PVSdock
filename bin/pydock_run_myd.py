#!/usr/bin/env python
import sys
import os
import argparse
from pvsdock import pydock_myd as pydock
import pandas as pd
# import modin.pandas as pd


class LoadFromConfig(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        with values as f:
            parser.parse_args(f.read().split(), namespace)


class ExtendAction(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest) or []
        items.extend(values)
        setattr(namespace, self.dest, items)


def parser_arg(parser):
    # docking parameter

    parser.register('action', 'extend', ExtendAction)
    parser.add_argument('--arg_file', type=open, required=False, default=None,
                        action=LoadFromConfig, help='argment file')
    parser.add_argument('--dock_config', type=str, required=False,
                        default=None, help='docking config file ')
    parser.add_argument('-v', '--docking_program', type=str, required=False,
                        default='rbdock',
                        help='select rdock, rbdock')
    parser.add_argument('--my_module', type=str, required=False,
                        default=None,
                        help='set user python module path (for pifinder)')
    parser.add_argument('--neutralize', action='store_true',
                        required=False, help='neutralize smiles ')
    parser.add_argument('--pH', type=float, default=None,
                        required=False, help='protonate state for pH 7.4 ')
    parser.add_argument('--output_save', action='store_true', required=False,
                        help='default output pdbqt is temp file ')
    parser.add_argument('--gen3d_dir', type=str, required=False, default='tmp',
                        help='3d initial conformation directory')
    parser.add_argument('--dock_dir', type=str, required=False,
                        default='tmp', help='binding conformation directory')
    parser.add_argument('--num_sub_proc', type=int, required=False,
                        default=10, help=' --num_sub_proc 10')
    parser.add_argument('--timeout_gen3d', type=int, required=False,
                        default=1, help=' --timeout_gen3d 1')
    parser.add_argument('--timeout_dock', type=int, required=False,
                        default=120, help=' --timeout_dock 120')
    parser.add_argument('--tlen', type=int, default='7', required=False,
                        help='lenth of sub directory name, default: 7')
    parser.add_argument('--pout', type=int, default='0', required=False,
                        help='print processing out: 0 or number, default: 0')
    parser.add_argument('--rescoring_program', type=str, required=False,
                        default='smina', help='smina path')
    parser.add_argument('--rescoring_config', type=str, required=False,
                        default=None, help='docking config file for rescoring')

    parser.add_argument('--tether_ref_lig', type=str, required=False,
                        default=None,
                        help='reference ligand for tether docking')
    parser.add_argument('--tether_SMARTS', type=str, required=False,
                        default=None, help='SMARTS pattern for tether docking')
    parser.add_argument('--tether_ref_coor_file', type=str, required=False,
                        default=None,
                        help='reference coordinate file for tether docking')

    parser.add_argument('--exhaustiveness', type=int, required=False,
                        default=10,
                        help='exhaustiveness for rdock')

    return


def arg_to_params(parser):

    use_my_module = False
    for i, m in enumerate(sys.argv):
        if m == '--my_module':
            my_module_path = sys.argv[i+1]
            use_my_module = True
            my_module_dir = os.path.dirname(my_module_path)
            sys.path.append(my_module_dir)
            import my_module
            my_module.parser_arg(parser)

    args = parser.parse_args()

    docking_program = args.docking_program
    num_sub_proc = args.num_sub_proc
    timeout_gen3d = args.timeout_gen3d
    timeout_dock = args.timeout_dock
    output_save = args.output_save
    gen3d_dir = args.gen3d_dir
    dock_dir = args.dock_dir
    dock_config_file = args.dock_config

    tlen = args.tlen
    pout = args.pout
    neutralize = args.neutralize
    pH = args.pH

    tether_ref_lig = args.tether_ref_lig
    tether_SMARTS = args.tether_SMARTS
    tether_ref_coor_file = args.tether_ref_coor_file

    exhaustiveness = args.exhaustiveness

    rescoring = False
    rescoring_config_file = args.rescoring_config
    rescoring_program = args.rescoring_program
    if rescoring_config_file is not None:
        rescoring = True

    docking_params = dict()
    docking_params['docking_program'] = docking_program
    docking_params['gen3d_dir'] = gen3d_dir
    docking_params['dock_dir'] = dock_dir
    docking_params['num_sub_proc'] = num_sub_proc
    docking_params['timeout_gen3d'] = timeout_gen3d
    docking_params['timeout_dock'] = timeout_dock
    docking_params['output_save'] = output_save
    docking_params['tlen'] = tlen
    docking_params['pout'] = pout
    docking_params['neutralize'] = neutralize
    docking_params['pH'] = pH
    docking_params['dock_config_file'] = dock_config_file
    docking_params['rescoring'] = rescoring
    docking_params['rescoring_program'] = rescoring_program
    docking_params['rescoring_config_file'] = rescoring_config_file
    docking_params['tether_ref_lig'] = tether_ref_lig
    docking_params['tether_SMARTS'] = tether_SMARTS
    docking_params['tether_ref_coor_file'] = tether_ref_coor_file
    docking_params['exhaustiveness'] = exhaustiveness

    my_module_path = args.my_module
    docking_params['use_my_module'] = use_my_module
    docking_params['my_module_path'] = my_module_path

    if use_my_module:
        docking_params = my_module.arg_to_params(parser, docking_params)

    return args, docking_params


def main():

    parser = argparse.ArgumentParser(description='docking with multi process')
    parser.add_argument('-l', '--ligand_list_file', type=str, required=False,
                        default='smiles.txt',
                        help=' --ligand_list_file smiles.txt')
    parser.add_argument('-o', '--output_file', type=str, required=False,
                        default='docking.txt',
                        help=' --output_file docking.txt')

    parser_arg(parser)

    args, docking_params = arg_to_params(parser)
    ligand_list_file = args.ligand_list_file
    output_file = args.output_file

    if len(sys.argv) < 2:
        parser.print_usage()
        sys.exit()
    if args.dock_config is None:
        parser.print_usage()
        print('dock_config is missing')
        sys.exit()

    df = pd.read_pickle(ligand_list_file)

#    num_data = df.shape[0]
    fkey = df.keys()[0]
    if fkey.startswith('#'):
        df.rename(columns={fkey: fkey[1:]}, inplace=True)
    smiles_list = df[['MOL_ID', 'SMILES']].values.tolist()

#    smiles_list = smiles_list[0:10]
    docking = pydock.Docking(docking_params)

    result_dict = docking.predict(smiles_list)
    docking_score_list = result_dict['docking']
    docking_min = [x[0] for x in docking_score_list]
    df['Docking1'] = docking_min
    df['Docking'] = docking_score_list
    if docking_params['rescoring']:
        rescoring = result_dict['docking_re']
        df['Docking_re'] = rescoring

    use_my_module = docking_params['use_my_module']
    my_module_path = docking_params['my_module_path']

    if use_my_module:
        my_module_dir = os.path.dirname(my_module_path)
        sys.path.append(my_module_dir)
        import my_module
        my_module.my_score_to_df(df, docking_params, result_dict)

    df.to_pickle(output_file)


if __name__ == "__main__":
    main()
