#! /bin/bash

## clean
#rm *static*

. ./profile.inc

domains="ION LIG NAD SAD SAR SIC TYR GOT GSN" 
domains="GSN"
URL=https://medeaf.ogs.it/internal-validation/gbolzon/MER/Domain_static_data

for dom in ${domains} ; do
        domdir=${PWD}/${dom}
        mkdir -p $domdir
        #wget ${URL}/${dom}/MIT_static.nc -O ${domdir}/MIT_static.nc ###
        #wget ${URL}/${dom}/meshmask.nc -O ${domdir}/meshmask.nc ###
        #wget ${URL}/${dom}/rivers_positions.json -O ${domdir}/rivers_positions.json ###
        curl "https://raw.githubusercontent.com/inogs/bathytools/refs/heads/main/mer_domains/rivers/${dom}.json" > ${domdir}/rivers.json


        rivers_json=${domdir}/rivers.json
        maskfile=${domdir}/meshmask.nc
        sewers_json=${domdir}/PointSource.json
        my_prex_or_die "python3 scarichi_json_gen.py -i $PWD/Allegato_1_Capitolato_Tecnico_B32_B35_scarichi.xlsx -d ${domdir} -m $maskfile -o $sewers_json"
        my_prex_or_die "python3 -u rbcs_gen.py -s $sewers_json -r $rivers_json --domdir ${domdir}"
done

