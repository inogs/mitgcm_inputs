#! /bin/bash

## clean
#rm -rf ???/

. ./profile.inc

domains="ION LIG NAD SAD SAR SIC TYR GOT GSN ISO CAG"

URL=https://medeaf.ogs.it/internal-validation/gbolzon/MER/Domain_static_data
URL=/leonardo_work/OGS_test2528_0/pub/gbolzon/MER/Domain_static_data/

curl "https://raw.githubusercontent.com/inogs/bathytools/refs/heads/main/mer_domains/rivers/main.json" > rivers_main.json
for dom in ${domains} ; do
        domdir=${PWD}/${dom}
        mkdir -p $domdir
	if [[ ! $URL =~ ^https ]]; then
		my_prex_or_die "cp $URL/${dom}/MIT_static.nc ${domdir}/MIT_static.nc"
		my_prex_or_die "cp ${URL}/${dom}/rivers_positions.json ${domdir}"
		my_prex_or_die "gzip -dc ${URL}/${dom}/meshmask.nc.gz > ${domdir}/meshmask.nc"
	else
            wget ${URL}/${dom}/MIT_static.nc -O ${domdir}/MIT_static.nc ###
            wget ${URL}/${dom}/meshmask.nc.gz -O ${domdir}/meshmask.nc.gz
            gzip -dc ${domdir}/meshmask.nc.gz > ${domdir}/meshmask.nc && rm -f ${domdir}/meshmask.nc.gz
            wget ${URL}/${dom}/rivers_positions.json -O ${domdir}/rivers_positions.json ###
        fi


        maskfile=${domdir}/meshmask.nc
        sewers_json=${domdir}/PointSource.json
        my_prex_or_die "python scarichi_json_gen.py -i $PWD/Allegato_1_Capitolato_Tecnico_B32_B35_scarichi.xlsx -d ${domdir} -m $maskfile -o $sewers_json"
        my_prex_or_die "python -u rbcs_gen.py -s $sewers_json -r rivers_main.json --domdir ${domdir}"
done
