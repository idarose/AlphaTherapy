#!/bin/bash
#SBATCH --account=NN9895K
#SBATCH --job-name AlphaTherapy
#SBATCH --mail-type=ALL
#SBATCH --mail-user=idapro@fys.uio.no
#SBATCH --time=1-10:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1 --cpus-per-task=40
#SBATCH --array=25,50,75,100,150

EXECUTABLE=main

ACTIVITY=$SLURM_ARRAY_TASK_ID
OUTPUTFILE="TerminalOutput_PC3_PIP_${ACTIVITY}kBq.txt"

CELLLINE=PC3_PIP
NUMBERITERATIONS=40


./$EXECUTABLE $CELLLINE $ACTIVITY $NUMBERITERATIONS  > $OUTPUTFILE