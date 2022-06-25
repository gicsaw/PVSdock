#!/usr/bin/env python
import sys
import os
import numpy as np
import argparse
import time
from filelock import FileLock
import pandas as pd
# import modin.pandas as pd
import subprocess


def get_job_list(list_dir):
    job_idx_list = list()
    list_file = list_dir + '/list.txt'
    if not os.path.exists(list_file):
        return job_idx_list
    freeze_lock = FileLock('{}.lock'.format(list_file))
    with freeze_lock.acquire(timeout=30):
        with open(list_file, 'r') as fp:
            lines = fp.readlines()
        for line in lines:
            job_idx_list += [int(line.strip())]
    job_idx_list2 = sorted(job_idx_list)

    return job_idx_list2


def set_job_list(job_idx_list, list_dir):
    list_file = list_dir + '/list.txt'
    freeze_lock = FileLock('{}.lock'.format(list_file))
    with freeze_lock.acquire(timeout=30):
        if os.path.exists(list_file):
            with open(list_file, 'r') as fp:
                lines = fp.readlines()
        else:
            lines = list()
        with open(list_file, 'w') as fp:
            for line in lines:
                fp.write(line)
            for job_idx in job_idx_list:
                line = job_idx + '\n'
                fp.write(line)

    return


def data_split(params_dict):

    file_size = params_dict['file_size']
    master_dir = params_dict['master_dir']
    todo_dir = params_dict['todo_dir']

    remain_smi_file = master_dir + '/' + 'remain.pkl'

    try:
        df_remain = pd.read_pickle(remain_smi_file)
    except pd.errors.EmptyDataError:
        return False

#    0: 'MOL_IDX', 1: 'MOL_ID', 2: 'SMILES' 3:, 'Pred', 4: 'uncertainty'
    fkey = df_remain.keys()[0]
    if fkey.startswith('#'):
        df_remain.rename(columns={fkey: fkey[1:]}, inplace=True)

    num_remain = df_remain.shape[0]
    if num_remain == 0:
        return False

    sample_smi_file = '%s/sample.pkl' % (master_dir)
    os.replace(remain_smi_file, sample_smi_file)

    todo_list = list()

    num_data = df_remain.shape[0]
    num_file = int(np.ceil(num_data/file_size))
    for idx in range(0, num_file):
        job_idx = '%d' % (idx)
        ini = idx * file_size
        fin = (idx+1) * file_size
        df_todo = df_remain[ini:fin]
        job_todo_file = todo_dir + '/' + job_idx + '.pkl'
        df_todo.to_pickle(job_todo_file)
        todo_list += [job_idx]

    set_job_list(todo_list, master_dir)
    set_job_list(todo_list, todo_dir)

    return True


def write_stage(stage_list, stage_file):
    fp_stage = open(stage_file, 'w')
    for stage in stage_list:
        line_out = '%s\n' % (stage)
        fp_stage.write(line_out)
    fp_stage.close()


def check_done(params_dict):

    check = False
    master_dir = params_dict['master_dir']
    done_dir = params_dict['done_dir']

    master_job_idx_list = get_job_list(master_dir)
    done_job_idx_list = get_job_list(done_dir)

    master_job = set(master_job_idx_list)
    done_job = set(done_job_idx_list)
    running_job = master_job - done_job
    if len(running_job) == 0:
        check = True
    return check


def gather_result(params_dict):

    master_dir = params_dict['master_dir']
    done_dir = params_dict['done_dir']

#    head_dict = {0: 'MOL_IDX', 1: 'MOL_ID',
#                 2: 'SMILES', 3: 'Docking1', 4: 'Docking'}

    docking_smi_file = master_dir + '/' + 'docking.pkl'

    all_data = list()
    if os.path.exists(docking_smi_file):
        df_docking = pd.read_pickle(docking_smi_file)
        all_data.append(df_docking)
    done_job_idx_list = get_job_list(done_dir)
    for idx in done_job_idx_list:
        job_idx = '%d' % (idx)
        job_done_file = done_dir + '/' + job_idx + '.pkl'
        df_done = pd.read_pickle(job_done_file)

#        df_done.rename(columns=head_dict, inplace=True)
        all_data.append(df_done)
    df = pd.concat(all_data, axis=0, ignore_index=True)
    df.to_pickle(docking_smi_file)
    return


def gapjil(params_dict):
    '''
        master gapjil
    '''
    sleep_cycle = params_dict['sleep_cycle']
    master_dir = params_dict['master_dir']
    auto_sub = params_dict['auto_sub']
    sub_dock_script = params_dict['sub_dock_script']
    num_sub = params_dict['num_sub']

    stage_file = master_dir + '/' + 'stage.txt'

    stage_list = []
    check_stage = None
    if os.path.exists(stage_file):
        stage_list = [x.strip() for x in open(stage_file)]
        if len(stage_list) != 0:
            if stage_list[-1] != 'done':
                check_stage = 'current'
            else:
                check_stage = 'done'
        line_out = 'restart simulation'
        print(line_out, flush=True)

    if check_stage is None:
        check_split = data_split(params_dict)
        if not check_split:
            line_out = 'strange: data is null'
            print(line_out, flush=True)
            sys.exit()
        check_stage = 'current'
        stage_list += ['current']
        write_stage(stage_list, stage_file)
        line_out = 'Start simulation'
        print(line_out, flush=True)

        if auto_sub == 'shell':
            run_line = 'bash ' + sub_dock_script
            for i in range(num_sub):
                subprocess.Popen(run_line.split())
            line_out = '%d sub_dock is generation' % num_sub
            print(line_out, flush=True)

        elif auto_sub == 'batch':
            run_line = 'sbatch ' + sub_dock_script
            for i in range(num_sub):
                subprocess.Popen(run_line.split())
            line_out = '%d sub_dock_batch is submitted' % num_sub
            print(line_out, flush=True)

    while check_stage != 'done':
        check_docking = check_done(params_dict)
        if check_docking:
            check_stage = 'done'
            stage_list += ['done']
            write_stage(stage_list, stage_file)
            break

        if sleep_cycle is None:
            line_out = 'stop master: sleep_cycle is None '
            print(line_out, flush=True)
            break
        time.sleep(sleep_cycle)

    gather_result(params_dict)
    line_out = 'End simulation'
    print(line_out, flush=True)

    return


class LoadFromConfig(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        with values as f:
            parser.parse_args(f.read().split(), namespace)


class ExtendAction(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest) or []
        items.extend(values)
        setattr(namespace, self.dest, items)


def main():

    parser = argparse.ArgumentParser(description='Master of puppets')
    parser.register('action', 'extend', ExtendAction)
    parser.add_argument('--arg_file', type=open, required=False, default=None,
                        action=LoadFromConfig, help='argment file')
    parser.add_argument('--work_dir', type=str, required=False,
                        default='workflow', help='workflow directory')
    parser.add_argument('--sleep_cycle', type=int, required=False, default=None,
                        help='sleep master for xx seconds per loop\n' +
                        'if None (default), escape without repeating the loop')
    parser.add_argument('--file_size', type=int, required=False,
                        default='1000', help='smiles per file')
    parser.add_argument('--auto_sub', type=str, required=False,
                        default=None,
                        help='shell (bash), batch (slurm) default is None')
    parser.add_argument('--sub_dock_script', type=str, required=False,
                        default='run_sub_dock.sh',
                        help='run_sub_dock.sh or slurm_sub_dock.sh')
    parser.add_argument('--num_sub', type=int, required=False, default=1,
                        help='number of docking subprocess')

    if len(sys.argv) < 2:
        parser.print_usage()
        sys.exit()
    args = parser.parse_args()

    work_dir = args.work_dir
    master_dir = work_dir + '/master'
    todo_dir = work_dir + '/todo'
    current_dir = work_dir + '/current'
    done_dir = work_dir + '/done'

    file_size = args.file_size
    sleep_cycle = args.sleep_cycle

    auto_sub = args.auto_sub
    sub_dock_script = args.sub_dock_script
    num_sub = args.num_sub

    params_dict = dict()
    params_dict['work_dir'] = work_dir
    params_dict['master_dir'] = master_dir
    params_dict['todo_dir'] = todo_dir
    params_dict['current_dir'] = current_dir
    params_dict['done_dir'] = done_dir
    params_dict['file_size'] = file_size
    params_dict['sleep_cycle'] = sleep_cycle
    params_dict['auto_sub'] = auto_sub
    params_dict['sub_dock_script'] = sub_dock_script
    params_dict['num_sub'] = num_sub

    gapjil(params_dict)


if __name__ == "__main__":
    main()
