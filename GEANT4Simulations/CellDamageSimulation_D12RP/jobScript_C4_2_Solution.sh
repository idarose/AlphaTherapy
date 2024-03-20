#!/bin/bash
#SBATCH --account=NN9895K
#SBATCH --job-name AlphaTherapy
#SBATCH --mail-type=ALL
#SBATCH --mail-user=idapro@fys.uio.no
#SBATCH --time=0-00:10:00
#SBATCH --mem-per-cpu=2000
#SBATCH --ntasks=1 --cpus-per-task=40
#SBATCH --array=0-359

# export G4RADIOACTIVEDATA="/Users/kcwli/Academic/Codes/Installations/GEANT4/geant4.10.07.p03-data/RadioactiveDecay5.6"
# export G4RADIOACTIVEDATA="/Users/kcwli/Academic/Codes/Installations/GEANT4/geant4.10.07.p03-data/RadioactiveDecay5.6_212Pb_212Bi"
# export G4RADIOACTIVEDATA="/Users/kcwli/Academic/Codes/Installations/GEANT4/geant4.10.07.p03-data/RadioactiveDecay5.6_212Bi_208Tl"
# export G4RADIOACTIVEDATA="/Users/kcwli/Academic/Codes/Installations/GEANT4/geant4.10.07.p03-data/RadioactiveDecay5.6_212Bi_212Po"
# export G4RADIOACTIVEDATA="/Users/kcwli/Academic/Codes/Installations/GEANT4/geant4.10.07.p03-data/RadioactiveDecay5.6_212Pb_212Bi"
# export G4RADIOACTIVEDATA="/Users/kcwli/Academic/Codes/Installations/GEANT4/geant4.10.07.p03-data/RadioactiveDecay5.6_212Po_208Pb"



EXECUTABLE=exampleB4a

OUTPUTFILE=$SLURM_ARRAY_TASK_ID
#OUTPUTFILE=Sim$SLURM_ARRAY_TASK_ID
OUTPUTFILE="TerminalOutput_${OUTPUTFILE}.txt"

MACRO=Run_C4_2_Solution_Case_$SLURM_ARRAY_TASK_ID
MACRO="${MACRO}.mac"

./$EXECUTABLE -m $MACRO > $OUTPUTFILE