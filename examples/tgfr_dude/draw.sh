#!/bin/bash
list_file=select/list.txt
#list_file=select/a.txt
mkdir pymol
idx=-1
while read line
do
    idx=$((idx+1))
    if [ $idx -eq 0 ]
    then
        continue
    fi
    
    array=($line)
    mol_id=${array[0]}
    model_id=${array[1]}
    m_id=$mol_id"_"$model_id
    receptor_pdb=pdb/2ZDXA_model_missing.pdb
    receptor_pf=pdb/2ZDXA_missing_pf_receptor.txt
    pi_draw.py -r $receptor_pdb -p $receptor_pf -l select/dock_$m_id.pdb -q select/pf_$m_id.txt -i select/pinter_$m_id.txt -i select/pinter_$m_id.txt -o pymol/M$idx"_"$m_id.pse -f pymol/M$idx"_"$m_id.png --show_cartoon
done < $list_file


