#!/bin/bash

set -e

OUTPUTDIR="/dev/shm/miotest"

BATHYTOOLSDIR="/home/spiani/projects/bathytools"
MITGCMINPUTSDIR="/home/spiani/projects/mitgcm_inputs"

DOMAINS_REG="NAD SAD ION SIC TYR LIG SAR GOT GSN"
DOMAINS_HR="HGT GOR CON PES LAM PAN NAP GAE FOL CAG ISO"

DOMAINS="$DOMAINS_REG $DOMAINS_HR"

MAXJOBS=5   # numero massimo di domini processati in parallelo

mkdir -p "${OUTPUTDIR}"

process_domain () {
    set -e
    domain=$1

    domaindir=${OUTPUTDIR}/${domain}
    rm -rf "${domaindir}"
    mkdir -p "${domaindir}"

    echo
    echo "----------^^^-----------  DOMAIN ${domain}  ----------^^^-----------"
    echo

    # Setta SPONGE
    if [[ " $DOMAINS_HR " =~ " $domain " ]]; then
        SPONGE=10
    else
        SPONGE=16
    fi

    cd "$BATHYTOOLSDIR"

    domain_descriptor="$BATHYTOOLSDIR/mer_domains/${domain,,}.yaml"
    if [ ! -f "${domain_descriptor}" ]; then
        domain_descriptor="$BATHYTOOLSDIR/mer_domains/high_reso/${domain,,}.yaml"
    fi

    poetry run bathytools -c "${domain_descriptor}" -o "${domaindir}" --mer

    cd "$MITGCMINPUTSDIR"

    poetry run mitgcm_inputs FLUXES \
        -m "${domaindir}/meshmask.nc" \
        -o "${domaindir}/fluxes.tar.gz"

    rpositions="${domaindir}/rivers_positions.json"
    rcustom="${BATHYTOOLSDIR}/mer_domains/rivers/${domain}.json"

    rbcs_options=""
    ob_indices_options=""

    if [ -f "$rpositions" ]; then
        rbcs_options="-p $rpositions -r $BATHYTOOLSDIR/mer_domains/rivers/main.json"
        ob_indices_options="-p $rpositions -r ${domaindir}/additional_variables.nc"

        if [ -f "${rcustom}" ]; then
            rbcs_options="${rbcs_options} -d ${rcustom}"
        fi
    fi

    poetry run mitgcm_inputs ob_indices \
        -m "${domaindir}/meshmask.nc" \
        -s "$SPONGE" \
        -o "${domaindir}/rivers_open_boundaries.txt" \
        -n "${domaindir}/nudging_indices.txt" \
        ${ob_indices_options}

    poetry run mitgcm_inputs rbcs \
        -m "${domaindir}/meshmask.nc" \
        -o "${domaindir}" \
        ${rbcs_options}

}

export -f process_domain
export OUTPUTDIR BATHYTOOLSDIR MITGCMINPUTSDIR DOMAINS_HR

if command -v parallel >/dev/null 2>&1; then
    parallel --lb --tag --halt soon,fail=1 -j ${MAXJOBS} process_domain ::: ${DOMAINS}
else
    echo "GNU parallel not found, falling back to serial execution." >&2
    for domain in $DOMAINS
    do
        process_domain "$domain"
    done
fi
