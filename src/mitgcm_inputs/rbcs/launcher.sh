#! /bin/bash

## clean
#rm *static*

JSON_DIR=$PWD/JSON_DIR
mkdir -p $JSON_DIR

for dom in ${domains}
do
        mkdir -p ${dom}
        wget ${URL}/${dom}/MIT_static.nc ###
        wget ${URL}/${dom}/rivers_positions.json ###
        mv MIT_static.nc ${dom}/MIT_static_${dom}.nc ###
        mv rivers_positions.json ${dom}/rivers_positions_${dom}.json ###
done

python scarichi_json_gen.py -i $PWD/Allegato_1_Capitolato_Tecnico_B32_B35_scarichi.xlsx -o $JSON_DIR ###
python fiumi_json_gen.py -i $PWD/listaFIUMI_MER.xlsx -o $JSON_DIR ###

URL=https://medeaf.ogs.it/internal-validation/gbolzon/MER/Domain_static_data
domains="GOT GSN ION LIG NAD SAD SAR SIC TYR"


for dom in ${domains}
do
        python -u rbcs_gen.py -i $JSON_DIR -d ${dom} --domdir $PWD/${dom}
        echo ${dom}" done!"
done

