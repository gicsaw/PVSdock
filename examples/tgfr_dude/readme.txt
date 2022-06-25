################################################
# VHTS for batch system
# prepare
# set pdb dir 
# prepare receptor_pdbqt file

# if usepharmacophore model: 
# find_pf_receptor.py -r $receptor_file.pdb -p $pf_receptor.txt -t $ref_ligand.pdb -u $pf_receptor_info.txt -v config.txt(docking box)
# check $pf_receptor_info.txt file and select feature

# transform input list file to pkl
csv2pkl.py -i smiles_list.csv -o smiles_list.pkl

### option1 vhts version #######################
# VHTS batch system
mkdir workflow
cd workflow
mkdir current done master todo
cd ..

cp smiles_list.pkl workflow/master/remain.pkl
# master
# check master_config.txt
nohup master_dock.py --arg_file master_config.txt > mm.txt 2>&1 &

# sub_dock
# check subdock_config.txt and config.txt

# for batch system

sbatch slurm_sub_dock.sh 
# repeat 4 (number of cpu node) 
# check slurm: squeue, sinfo

# or background (no batch ) 

# nohup sub_dock.py --arg_file subdock_config.txt --dock_config config.txt --use_mci --num_sub_proc 256  > a.txt 2>&1 &
### option1 vhts version end ####################


### option2 simple version ######################
# (simple version)
# check pydock_config.txt
nohup pydock_run.py -l smiles_list.pkl -o docking.pkl --arg_file pydock_config.txt --dock_config config.txt --use_mci --num_sub_proc=256 > a.txt 2>&1 &
### option2 simple version end #################


################################################
#####       post process                   #####

### Medicinal Chemistry Filter #################
# rd_filter: 
python rdf_input.py
rd_filters filter --in rdf_input.txt --prefix rdf --rules rules.json --np 32



#### option1: reference fragment coordinate filter  ######
# reference fragment check 

# find reference fragment coordinate 
# example: find_fragment_pdb.py
# Check if fragment matches reference coordinates 
# example: fragment.py 

# Modify these two files according to the problem 
sbatch slurm_frag.sh

# Rescoring by reflecting fragment reference matching 
python rescoring_fragment.py

# select molecule and pose
python split_model_fragment.py

#### option1: reference fragment coordinate filter end #####


#### option2: pharmacophore rescoring #############
#### (Vina + MCIScore + info)
python rescoring_vmcisinfo.py
# select molecule and pose
python split_model_vmcisinfo.py

#### option2 pharmacophore rescoring end ##########

# see select dir


# draw mol figure
python fig.py
# draw pharmacophore
./draw.sh

# draw building block (only for Enamine REAL)
python synthesis.py







