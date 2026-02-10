#! /bin/bash


. ./profile.inc
domains="ION LIG NAD SAD SAR SIC TYR GOT GSN ISO CAG"

PUBDIR=/leonardo_work/OGS_test2528_0/pub/
# generation of files for https://medeaf.ogs.it/internal-validation/gbolzon/MER/Domain_static_data
HERE=$PWD
for dom in ${domains} ; do
        cd ${dom}
        echo ${dom}
        mkdir -p conc relax
        my_prex_or_die "mv conc*bin conc/"

        my_prex_or_die "mv bottom*bin relax/"
        my_prex_or_die "tar -cf conc.tar conc/"
        my_prex_or_die "tar -cf relax.tar relax/"
        my_prex_or_die "gzip -c conc.tar > conc.tar.gz && rm -f conc.tar"
        my_prex_or_die "gzip -c relax.tar > relax.tar.gz && rm -f relax.tar"
        my_prex_or_die "mv conc.tar.gz  relax.tar.gz $PUBDIR/gbolzon/MER/Domain_static_data/${dom}/"
        cd $HERE
done
