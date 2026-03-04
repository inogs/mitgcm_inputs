OUTPUTDIR="/dev/shm/miotest"

BATHYTOOLSDIR="/home/spiani/projects/bathytools"
MITGCMINPUTSDIR="/home/spiani/projects/mitgcm_inputs"


DOMAINS="NAD SAD ION SIC TYR LIG SAR GOT"
DOMAINS="$DOMAINS HGT GOT CON PES LAM PAN NAP GAE FOL CAG ISO"

mkdir -p ${OUTPUTDIR}

for domain in ${DOMAINS}
  do
    domaindir=${OUTPUTDIR}/${domain}
    mkdir ${domaindir}

    echo
    echo "----------^^^-----------  DOMAIN ${domain}  ----------^^^-----------"
    echo

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

    # Check if this domain has rivers to write the rbcs command line options
    rpositions=${domaindir}/rivers_positions.json
    rcustom=${BATHYTOOLSDIR}/mer_domains/rivers/${domain}.json
    rbcs_options=""

    if [ -f "$rpositions" ]; then
        rbcs_options="-p $rpositions -r $BATHYTOOLSDIR/mer_domains/rivers/main.json"
        if [ -f "${rcustom}" ]; then
            rbcs_options="${rbcs_options} -d ${rcustom}"
        fi
    fi

    echo poetry run mitgcm_inputs rbcs -m ${domaindir}/meshmask.nc -o ${domaindir} ${rbcs_options}
    poetry run mitgcm_inputs rbcs -m ${domaindir}/meshmask.nc -o ${domaindir} ${rbcs_options}
  done
