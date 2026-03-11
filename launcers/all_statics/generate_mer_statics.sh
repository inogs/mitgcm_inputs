#!/bin/bash

set -e

OUTPUTDIR="/dev/shm/miotest"

BATHYTOOLSDIR="/home/spiani/projects/bathytools"
MITGCMINPUTSDIR="/home/spiani/projects/mitgcm_inputs"

DOMAINS_REG="NAD SAD ION SIC TYR LIG SAR GOT GSN"
DOMAINS_HR="HGT GOR CON PES LAM PAN NAP GAE FOL CAG ISO"

DOMAINS="$DOMAINS_REG $DOMAINS_HR"

mkdir -p ${OUTPUTDIR}

for domain in ${DOMAINS}
  do
    domaindir=${OUTPUTDIR}/${domain}
    mkdir ${domaindir}

    echo
    echo "----------^^^-----------  DOMAIN ${domain}  ----------^^^-----------"
    echo

    # Setta il valore dello SPONGE
    if [[ " $DOMAINS_HR " =~ " $domain " ]]; then
        SPONGE=10
    else
        SPONGE=16
    fi

    echo SPONGE $SPONGE

    cd $BATHYTOOLSDIR
    # Check where is the correct configuration file for the domain (in the mer_domains
    # directory or in the mer_domains/high_reso dir)
    domain_descriptor="$BATHYTOOLSDIR/mer_domains/${domain,,}.yaml"
    if [ ! -f "${domain_descriptor}" ]; then
        domain_descriptor="$BATHYTOOLSDIR/mer_domains/high_reso/${domain,,}.yaml"
    fi

    # Run bathytools
    poetry run bathytools -c ${domain_descriptor} -o ${domaindir} --mer

    # Generate fluxes.tar.gz
    cd ${MITGCMINPUTSDIR}
    poetry run mitgcm_inputs FLUXES -m ${domaindir}/meshmask.nc -o ${domaindir}/fluxes.tar.gz

    # Check if this domain has rivers to write the rbcs and ob_indices command
    # line options
    rpositions=${domaindir}/rivers_positions.json
    rcustom=${BATHYTOOLSDIR}/mer_domains/rivers/${domain}.json
    rbcs_options=""
    ob_indices_options=""

    if [ -f "$rpositions" ]; then
        rbcs_options="-p $rpositions -r $BATHYTOOLSDIR/mer_domains/rivers/main.json"
        ob_indices_options="-p $rpositions"
        if [ -f "${rcustom}" ]; then
            rbcs_options="${rbcs_options} -d ${rcustom}"
        fi
    fi

    poetry run mitgcm_inputs ob_indices -m ${domaindir}/meshmask.nc -r ${domaindir}/additional_variables.nc -s $SPONGE -o ${domaindir}/rivers_open_boundaries.txt -n ${domaindir}/nudging_indices.txt ${ob_indices_options}

    poetry run mitgcm_inputs rbcs -m ${domaindir}/meshmask.nc -o ${domaindir} ${rbcs_options}
  done
