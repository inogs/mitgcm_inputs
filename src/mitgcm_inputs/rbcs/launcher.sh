#! /bin/bash

## clean
#rm *static*

URL=https://medeaf.ogs.it/internal-validation/gbolzon/MER/Domain_static_data
domains="GOT GSN ION LIG NAD SAD SAR SIC TYR"

for dom in ${domains}
do
        mkdir ${dom}
        wget ${URL}/${dom}/MIT_static.nc
        wget ${URL}/${dom}/rivers_positions.json
        mv MIT_static.nc ${dom}/MIT_static_${dom}.nc
        mv rivers_positions.json ${dom}/rivers_positions_${dom}.json
done

## jsons
base_path=/g100_scratch/userexternal/fgiordan/MER/
script_path=/g100_scratch/userexternal/fgiordan/MER/MITgcm_BC_IC/mitgcm/bc_ic/static-data/masks/MER

#python ${script_path}read_excellistaFIUMI_do_json_for_7domains.py ${base_path} ${script_path}
python ${script_path}/scarichi_json_gen.py ${base_path} ${script_path}
echo "Excel read!"

## download and make binaries 


for dom in ${domains}; do
        python -u ${script_path}/rbcs_gen.py ${base_path} ${dom}
        echo ${dom}" done!"
done

