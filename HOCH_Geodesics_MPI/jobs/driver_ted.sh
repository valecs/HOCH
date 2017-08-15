#!/bin/bash
set -euo pipefail

readonly R=$(dirname "$(readlink -f "$0")")

readonly E="$1"
readonly K="$2"

readonly OUT="$R/$E/sge"

echo "Launching $K, E = $E"
sleep 2s;
cd "$R/$E"

# qsub -A QMC-AE -n 512 --mode c32 -t $TIME --env STRATT_ROOT=$STRATT_ROOT --jobname "$E-$K:1"\
#     $(readlink -f $R/../../shared/bin/HOCH_Geodesics_MPI) -$K -e $E

qsub -b yes -shell no -j yes -o $OUT -wd $OUT -N "HOCH$E$K" -p -1 -q stratt,chemistry -tc 20\
     -v LD_LIBRARY_PATH="$HOME/toor/lib:$HOME/toor/gcc62/lib64:/share/apps/lib",STRATT_ROOT=$STRATT_ROOT\
     -t 1-5000\
     $(readlink -f $R/../../shared/bin/HOCH_Geodesics_MPI) -e $E -G -$K

qsub -b yes -shell no -j yes -o $OUT -wd $OUT -N "HOCH$E$K" -p -1 -q stratt,chemistry -tc 20\
     -v LD_LIBRARY_PATH="$HOME/toor/lib:$HOME/toor/gcc62/lib64:/share/apps/lib",STRATT_ROOT=$STRATT_ROOT\
     -t 25001-30000\
     $(readlink -f $R/../../shared/bin/HOCH_Geodesics_MPI) -e $E -G -$K

cd -

echo "Don't forget to rebalance! e.g.:"
echo "qstat | REBALANCE_MAX=110 rebalance.awk"
