#!/bin/bash
#SBATCH --account=NN9895K
#SBATCH --job-name VF_KCWLI
#SBATCH --mail-type=ALL
#SBATCH --mail-user=k.c.w.li@fys.uio.no
#SBATCH --time=0-12:00:00
#SBATCH --mem-per-cpu=1500
#SBATCH --ntasks=1 --cpus-per-task=20
#SBATCH --array=245,461

EXECUTABLE=VF_scenario${SLURM_ARRAY_TASK_ID}_4plusH1
OUTPUTFILE=/cluster/projects/nn9895k/kcwli/VeridicalFit/VeridicalFitResults/12C_2021/scenario${SLURM_ARRAY_TASK_ID}_4plusH1
OUTPUTFILE="${OUTPUTFILE}_TerminalOutput.txt"

./$EXECUTABLE > $OUTPUTFILE
Skriv til Kevin Li

