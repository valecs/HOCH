#!/bin/bash

set -euo pipefail

readonly f_name=$1

readonly exe=$(readlink -f ../../../shared/bin/HOCH_Geodesics_MPI)

if [ ! -f "$f_name" ]
then
    echo "$f_name is not a file; exiting!"
    exit 1
fi

shift 1

readonly F=$(readlink -f "$f_name")

readonly out_name="./${f_name}_o"
mkdir "$out_name" || true

readonly out=$(readlink -f "$out_name")

readonly N=$(wc -l < "$F")

echo "submitting $F"
echo "writing to $out"
echo "$N tasks"

sleep 3s

qsub -b yes -shell no -j yes -o "$out" -wd "$out" -N "j$f_name" -p -1 -q stratt,chemistry,gaussian,doll -tc 140 \
     -v LD_LIBRARY_PATH="$HOME/toor/lib:$HOME/toor/gcc62/lib64:/share/apps/lib",STRATT_ROOT="$STRATT_ROOT" \
     -t "1-$N" "$@" "$exe" -G -F "$F"
