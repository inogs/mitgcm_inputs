#! /bin/bash

# valid only for g100

. ./profile.inc
domains="ION LIG NAD SAD SAR SIC TYR GOT GSN" 
# generation of files for https://medeaf.ogs.it/internal-validation/gbolzon/MER/Domain_static_data
HERE=$PWD
for dom in ${domains} ; do
        cd ${dom}
        echo ${dom}
        mkdir -p conc relax
        my_prex "mv conc*bin conc/"
        my_prex "mv bottom*bin relax/"
        my_prex_or_die "tar -cf conc.tar conc/"
        my_prex_or_die "tar -cf relax.tar relax/"
        my_prex_or_die "gzip -c conc.tar > conc.tar.gz && rm -f conc.tar"
        my_prex_or_die "gzip -c relax.tar > relax.tar.gz && rm -f relax.tar"
        cd $HERE
done

for dom in ${domains} ; do
   my_prex_or_die "mv {dom}/conc.tar.gz /g100_work/OGS_devC/Benchmark/pub/gbolzon/MER/Domain_static_data/${dom}"
   my_prex_or_die "mv {dom}/relax.tar.gz /g100_work/OGS_devC/Benchmark/pub/gbolzon/MER/Domain_static_data/${dom}"
done