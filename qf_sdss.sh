#!/bin/usr/
#setup idlutils; setup -r ~/primus/; setup -r ~/pro/
idl -e ".r build_target_sample.pro"
idl -e ".r build_env_sample.pro"
idl -e ".r build_environment_cylinder.pro"
idl -e ".r build_vmax_avail.pro"

run=$1
Nransack=$2
Nrandom=$3
echo "run= "$run", Nransack= "$Nransack", Nrandom= "$Nrandom
##########################################################################################
# Ransack, Random, and V_max,avail: 
########################################################################################
#sh build_ransack.sh sdss $Nransack $Nrandom
idl -e "build_vmax_avail,"$run","$Nransack","$Nrandom",/sdss"

##########################################################################################
# Build target and environment defining population:  
##########################################################################################
idl -e "build_target_sample,"$run","$Nransack","$Nrandom",/sdss"

idl -e "build_environment_cylinder,"$run","$Nransack",/sdss"

#idl -e "build_smf_im_mf_vmax,"$run
