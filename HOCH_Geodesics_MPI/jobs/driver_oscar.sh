#!/bin/bash
set -euo pipefail

readonly R=$(dirname "$(readlink -f "$0")")

readonly E="$1"
readonly K="$2"
shift 2

readonly OUT="$R/$E/slurm"

echo "Launching $K, E = $E"
sleep 2s;
cd "$OUT"

sbatch --array=1-2%2 --workdir="$OUT" --export=STRATT_ROOT --hint=compute_bound\
       --job-name "HOCH$E$K" --time 8-0 $@ <<EOF
#!/bin/sh
srun $(readlink -f $R/../../shared/bin/HOCH_Geodesics_MPI) -e $E -G -$K
EOF

## R.B. Ted jobs did 1-5000 and 25001-30000

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
