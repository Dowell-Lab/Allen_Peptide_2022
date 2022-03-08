#!/bin/bash
#SBATCH --job-name=bedtoolsmulticov                     # Job name
#SBATCH --mail-type=ALL                           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tajo5912@colorado.edu                       # Where to send mail
#SBATCH --nodes=1                               # Number of nodes requested
#SBATCH --ntasks=1                              # Number of CPUs (processor cores/tasks)
#SBATCH --mem=8gb                               # Memory limit
#SBATCH --time=05:00:00                         # Time limit hrs:min:sec
#SBATCH --partition=short                       # Partition/queue requested on server
#SBATCH --output=/scratch/Users/tajo5912/peptide/ChIP-seq/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/tajo5912/peptide/ChIP-seq/eofiles/%x_%j.out 

### Displays the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

pwd; hostname; date

HOME=/scratch/Users/tajo5912/peptide/ChIP-seq

### Loads modules
module load bedtools/2.25.0

### <SOFTWARE SPECIFICS>
bedtools multicov -bams \
${HOME}/mapped/bams/DMSO_R1_NEW_S3_R1_001.sorted.bam \
${HOME}/mapped/bams/DMSO_R2_NEW_S4_R1_001.sorted.bam \
${HOME}/mapped/bams/DMSO_pep_R1_S1_R1_001.sorted.bam \
${HOME}/mapped/bams/DMSO_pep_R2_S2_R1_001.sorted.bam \
${HOME}/mapped/bams/Nutlin_R1_NEW_S15_R1_001.sorted.bam \
${HOME}/mapped/bams/Nutlin_R2_NEW_S17_R1_001.sorted.bam \
${HOME}/mapped/bams/Nutlin_R2_S16_R1_001.sorted.bam \
${HOME}/mapped/bams/Nutlin_pep_R1_S13_R1_001.sorted.bam \
${HOME}/mapped/bams/Nutlin_pep_R2_S14_R1_001.sorted.bam \
${HOME}/mapped/bams/in_DMSO_R1_S7_R1_001.sorted.bam \
${HOME}/mapped/bams/in_DMSO_R2_S8_R1_001.sorted.bam \
${HOME}/mapped/bams/in_DMSO_pep_R1_S5_R1_001.sorted.bam \
${HOME}/mapped/bams/in_DMSO_pep_R2_S6_R1_001.sorted.bam \
${HOME}/mapped/bams/in_Nutlin_R1_S11_R1_001.sorted.bam \
${HOME}/mapped/bams/in_Nutlin_R2_S12_R1_001.sorted.bam \
${HOME}/mapped/bams/in_Nutlin_pep_R1_S9_R1_001.sorted.bam \
${HOME}/mapped/bams/in_Nutlin_pep_R2_S10_R1_001.sorted.bam \
-bed ${HOME}/annotations/hg38_longest_isoform_shuffle.bed \
> ${HOME}/counts/hg38_longest_isoform_shuffle_counts.txt


echo Job finished at `date +"%T %a %d %b %Y"`
