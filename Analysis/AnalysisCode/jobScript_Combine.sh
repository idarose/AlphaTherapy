#!/bin/bash
#SBATCH --account=NN9895K
#SBATCH --job-name AlphaTherapy
#SBATCH --mail-type=ALL
#SBATCH --mail-user=idapro@fys.uio.no
#SBATCH --time=1-10:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --ntasks=1 --cpus-per-task=10
#SBATCH --array=1,3,5,10,25,50,75,100,150

EXECUTABLE=combine

ACTIVITY=$SLURM_ARRAY_TASK_ID

CELLGEOMETRY=D12CP

./$EXECUTABLE $CELLGEOMETRY $ACTIVITY