#!/bin/bash
#SBATCH --job-name=pyGenomeTracks				 # Job name
#SBATCH --mail-type=ALL				 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tajo5912@colorado.edu		 # Where to send mail
#SBATCH --nodes=1					 # Number of cores job will run on
#SBATCH --ntasks=1					 # Number of CPU (processers, tasks)
#SBATCH --time=00:10:00					 # Time limit hrs:min:sec
#SBATCH --partition short				 # Job queue
#SBATCH --mem=4gb					 # Memory limit
#SBATCH --output=/scratch/Users/tajo5912/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/tajo5912/eofiles/%x_%j.err

################### SET VARIABLES ######################################

OUT='/scratch/Users/tajo5912/peptide/igv_combined'
INI='/scratch/Users/tajo5912/peptide/igv_combined/ini'
ROOTNAME='chip_peptide_r1Nut_btg2'

########################################################################
################### LOAD NECESSARY MODULES #############################

module load bedtools ##this is needed to add the genome track
source ~/venv_pygenometracks/bin/activate

########################################################################
################## PRINT JOB INFO ######################################

printf "Output Dir: $OUT"
printf "\nINI File: $INI"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

#######################################################################

pyGenomeTracks \
 --tracks ${INI}/${ROOTNAME}.ini \
 --region chr1:203,305,393-203,313,476 \
 --outFileName $OUT/${ROOTNAME}.svg

printf "complete"
