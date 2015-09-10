#!/bin/sh

#SBATCH --job-name=test
#SBATCH --ntasks=21

#SBATCH --mem-per-cpu=2048

#SBATCH --exclude=synch
#SBATCH --ntasks-per-node=4
#SBATCH --partition=normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=thiagos@lptl.jussieu.fr

# for large numbers of tasks the controlling node will have a large number
# of processes, so it will be necessary to change the user process limit
#ulimit -u 10000

# --delay .2 prevents overloading the controlling node
# -j is the number of tasks parallel runs so we set it to $SLURM_NTASKS
# --joblog makes parallel create a log of tasks that it has already run
# --resume makes parallel use the joblog to resume from where it has left off
# the combination of --joblog and --resume allow jobs to be resubmitted if
# necessary and continue from where they left off
parallel="/users/lptl/thiagos/softs/parallel/bin/parallel --delay 0.2 -j $SLURM_NTASKS"
