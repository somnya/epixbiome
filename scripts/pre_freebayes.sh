#!/bin/bash
# Name: pre_freebayes
# Author: Eddy Mendoza
# Date: 02/02/24

# === SLURM Config ===
#SBATCH -J pre_bayes
#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --array=1-169%13
#SBATCH -o scripts/exegfiles/pre_freebayes.out
#SBATCH -e scripts/exegfiles/pre_freebayes.err

# === SET UP ENVIRONMENT ===

module load slurm_setup
source activate env_g
cd taller/ra42hux/

# === ARRAY INSTRUCTIONS ===

KEY="er_data/tes/key.txt"

line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY} )

# === JOB ===

# sort and index every file

samtools sort -@ 16 er_data/snps_rev_bay/masked/${line}.bam > er_data/snps_rev_bay/masked/sorted_${line}.bam

samtools index -@ 16 er_data/snps_rev_bay/masked/sorted_${line}.bam