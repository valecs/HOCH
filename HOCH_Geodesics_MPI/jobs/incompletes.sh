#!/bin/bash
set -euo pipefail

#sort <(seq -w 1 16000) <(awk '/SUCCESS/ { print $1 }' < 1132999.error) | uniq -u | wc

readonly E=$1
readonly R=$(dirname "$(readlink -f "$0")")
readonly P="$PWD"

readonly X="incomplete"
readonly COBALT="cobaltlog"
readonly OUTPUT="error"

(# Subshell
    cd "$R/$E"

    for K in D A R
    do
	fname="$P/$E.$K.$X"
	rm "$fname" 2>/dev/null && echo "deleting $fname" || echo "will create $fname"
	for F in *.$COBALT
	do
	    if grep "$K:" "$F" > /dev/null
	    then
		N=$(awk '/Command:/ { R=$(NF-1); sub("-", " ", R); gsub("'\''", "", R); print R }' "$F")
		echo "$F matches for $K: $N"

		S=$(basename "$F" .$COBALT).$OUTPUT

		# intentionally allowing $N to split
		sort <(seq -w $N) <(awk '/SUCCESS/ { print $1 }' "$S")\
		    | uniq -u >> "$fname"
	    fi
	done
	sort -o "$fname" "$fname"
    done
)
