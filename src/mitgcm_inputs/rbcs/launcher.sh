#! /bin/bash

## clean
#rm *static*

. ./profile.inc

domains="ION LIG NAD SAD SAR SIC TYR GOT GSN"
domains="GSN"
URL=https://medeaf.ogs.it/internal-validation/gbolzon/MER/Domain_static_data
URL=/leonardo_work/OGS_test2528_0/MER/bathytools/

for dom in ${domains} ; do
        domdir=${PWD}/${dom}
        mkdir -p $domdir
	cp $URL/${dom}/MIT_static.nc ${domdir}
	cp ${URL}/${dom}/rivers_positions.json ${domdir}
	cp ${URL}/${dom}/meshmask.nc ${domdir}
        #wget ${URL}/${dom}/MIT_static.nc -O ${domdir}/MIT_static.nc ###
        #wget ${URL}/${dom}/meshmask.nc -O ${domdir}/meshmask.nc ###
        #wget ${URL}/${dom}/rivers_positions.json -O ${domdir}/rivers_positions.json ###
        curl "https://raw.githubusercontent.com/inogs/bathytools/refs/heads/main/mer_domains/rivers/${dom}.json" > ${domdir}/rivers.json


        rivers_json=${domdir}/rivers.json
        maskfile=${domdir}/meshmask.nc
        sewers_json=${domdir}/PointSource.json
        my_prex_or_die "python scarichi_json_gen.py -i $PWD/Allegato_1_Capitolato_Tecnico_B32_B35_scarichi.xlsx -d ${domdir} -m $maskfile -o $sewers_json"
        my_prex_or_die "python -u rbcs_gen.py -s $sewers_json -r $rivers_json --domdir ${domdir}"
done
