#! /bin/bash

## clean
#rm *static*

. ./profile.inc

JSON_DIR=$PWD/JSON_DIR
mkdir -p $JSON_DIR

domains="ION LIG NAD SAD SAR SIC TYR GOT GSN" 
URL=https://medeaf.ogs.it/internal-validation/gbolzon/MER/Domain_static_data

for dom in ${domains}
do
        mkdir -p ${dom}
        wget ${URL}/${dom}/MIT_static.nc ###
        wget ${URL}/${dom}/rivers_positions.json ###
        wget https://raw.githubusercontent.com/inogs/bathytools/refs/heads/main/mer_domains/rivers/${dom}.json
        mv MIT_static.nc ${dom}/MIT_static_${dom}.nc ###
        mv rivers_positions.json ${dom}/rivers_positions_${dom}.json ###
        mv ${dom}.json ${JSON_DIR}/${dom}.json
done

my_prex_or_die "python scarichi_json_gen.py -i $PWD/Allegato_1_Capitolato_Tecnico_B32_B35_scarichi.xlsx -o $JSON_DIR"

for dom in ${domains} ; do
        my_prex_or_die "python -u rbcs_gen.py -i $JSON_DIR -d ${dom} --domdir $PWD/${dom}"
done

