#!/bin/bash
export DISPLAY=:0 # Or localhost:10.0
#SBATCH --account=NN9895K
#SBATCH --job-name AlphaTherapy
#SBATCH --mail-type=ALL
#SBATCH --mail-user=idapro@fys.uio.no
#SBATCH --time=0-00:10:00
#SBATCH --mem-per-cpu=2000
#SBATCH --ntasks=1 --cpus-per-task=10
#SBATCH --array=5
##SBATCH --array=1,3,5,10,25,50,75,100,150

EXECUTABLE=main

ACTIVITY=$SLURM_ARRAY_TASK_ID
OUTPUTFILE="TerminalOutput_C4_2_${ACTIVITY}kBq.txt"

CELLLINE=C4_2
NUMBERITERATIONS=10


./$EXECUTABLE $CELLLINE $ACTIVITY $NUMBERITERATIONS  > $OUTPUTFILE