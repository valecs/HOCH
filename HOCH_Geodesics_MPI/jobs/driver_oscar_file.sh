#!/bin/bash
set -euo pipefail

readonly R=$(dirname "$(readlink -f "$0")")

if [ ! -f "$1" ]
then
    echo "$1 is not a file; exiting!"
    exit 1
fi

readonly F=$(readlink -f "$1")

readonly suffix=$(basename $F)

shift 1

readonly OUT="$R/mixed/slurm"

readonly N=$(wc -l < "$F")

echo "Have $N in $F"

echo "Launching mixed job: $F, have $N"
sleep 2s;

cd "$OUT"

for start in $(seq 1 1000 $N)
do
    export MY_ARRAY_START=$((start -1))
    echo "launching $MY_ARRAY_START"
    #sleep 4s

    if [ "$(((N/1000)*1000))" -eq "$MY_ARRAY_START" ]
    then
	M=$((N-MY_ARRAY_START))
    else
	M=1000
    fi
    
    sbatch --array=1-$M --workdir="$OUT" --hint=compute_bound\
	   --job-name "HOCH-$suffix-$MY_ARRAY_START" --time 2-0 --mem-per-cpu=256 "$@" <<EOF
#!/bin/bash
SLURM_ARRAY_TASK_ID="\$((SLURM_ARRAY_TASK_ID + MY_ARRAY_START))"
export SLURM_ARRAY_TASK_ID
srun $(readlink -f $R/../../shared/bin/HOCH_Geodesics_MPI) -G -F $F
EOF
done

cd -

## sbatch flags of note

## look into checkpointing...

# --array=<indexes> Submit a job array, multiple jobs to be executed
#     with identical parameters.  The indexes specification identifies
#     what array index values should be used. Multiple values may be
#     specified using a comma separated list and/or a range of values
#     with a "-" separator. For example, "--array=0-15" or
#     "--array=0,6,16-32".  A step function can also be specified with a
#     suffix containing a colon and number. For example,
#     "--array=0-15:4" is equivalent to "--array=0,4,8,12".  A maximum
#     number of simultaneously running tasks from the job array may be
#     specified using a "%" separator.  For example "--array=0-15%4"
#     will limit the number of simultaneously running tasks from this
#     job array to 4.  The minimum index value is 0.  the maximum value
#     is one less than the configuration parameter MaxArraySize.
