#!/bin/bash
set -euo pipefail

readonly Es=(34500 35000 35500 36000 36223 37000 37600 37688 38000 38688 38814 41010)
readonly R=$(dirname "$(readlink -f "$0")")

readonly TIME=90

readonly E="$1"
readonly K="$2"

echo "Launching $K, E = $E"
sleep 2s;
cd $R/$E
qsub -A QMC-AE -n 512 --mode c32 -t $TIME --env STRATT_ROOT=$STRATT_ROOT --jobname "$E-$K:1"\
    $(readlink -f $R/../HOCH_Geodesics_MPI) -e $E -t 1-16000 -$K
qsub -A QMC-AE -n 512 --mode c32 -t $TIME --env STRATT_ROOT=$STRATT_ROOT --jobname "$E-$K:2"\
    $(readlink -f $R/../HOCH_Geodesics_MPI) -e $E -t 16001-32000 -$K
qsub -A QMC-AE -n 512 --mode c32 -t $TIME --env STRATT_ROOT=$STRATT_ROOT --jobname "$E-$K:3"\
    $(readlink -f $R/../HOCH_Geodesics_MPI) -e $E -t 32001-48000 -$K
	 
cd $R
