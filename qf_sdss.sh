#!/bin/usr/
setup idlutils; setup -r ~/primus/; setup -r ~/pro/
idl -e ".r build_target_sample.pro"
idl -e ".r build_env_sample.pro"
idl -e ".r build_environment_cylinder.pro"
idl -e ".r build_vmax_avail.pro"

#echo "run="; read run
#echo "N_ransack="; read Nransack
#echo "N_random="; read Nrandom
run=$1
Nransack=$2
Nrandom=$3
##########################################################################################
# Ransack, Random, and V_max,avail: 
########################################################################################
#sh build_ransack.sh sdss $Nransack $Nrandom

#idl -e "build_vmax_avail,"$run","$Nransack","$Nrandom",/sdss"
##########################################################################################
# Build target and environment defining population:  
##########################################################################################
#idl -e "build_target_sample,"$run","$Nransack","$Nrandom",/literature,/sdss"
#
#idl -e "build_env_sample"
#
idl -e "build_environment_cylinder,"$run","$Nransack",/sdss,/literature"
#idl -e "build_smf,"$run",/sdss,/literature,/prank"
#idl -e "build_smf,"$run",/sdss,/literature"
