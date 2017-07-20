#!/bin/bash
set -euo pipefail

readonly N=40000
readonly Es=(34500 35000 35500 36000 36223 37688 38688 38814 41010)

for E in "${Es[@]}"
do
    echo "----- ENERGY = $E -----"
    ./HOCH_MicroCanonical $E $N
done
