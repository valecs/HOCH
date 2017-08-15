#!/bin/bash
set -euo pipefail

readonly R=$(dirname "$(readlink -f "$0")")

readonly E="$1"
readonly K="$2"

readonly I="$R/incompletes/$E.$K.incomplete"
readonly P="$I.part"

if [ -f $I ]
then
    echo "Reading incompletes from $I"
else
    echo "Cannot read $I"
    exit
fi

awk '$1<=5000 || ($1>=25000 && $1<=30000)' "$I" > "$P"
readonly N=$(wc -l < "$P")

echo "Have $N in $P"

readonly OUT="$R/$E/sge"

echo "Launching $K, E = $E, incompletes"
sleep 2s;
cd "$R/$E"

qsub -b yes -shell no -j yes -o "$OUT" -wd "$OUT" -N "HOCH$E$K"I -p -1 -q stratt,chemistry -tc 20\
     -v LD_LIBRARY_PATH="$HOME/toor/lib:$HOME/toor/gcc62/lib64:/share/apps/lib",STRATT_ROOT="$STRATT_ROOT"\
     -t "1-$N"\
     $(readlink -f $R/../../shared/bin/HOCH_Geodesics_MPI) -e $E -G -$K -F $P

cd -

echo "Don't forget to rebalance! e.g.:"
echo "qstat | REBALANCE_MAX=110 rebalance.awk"
