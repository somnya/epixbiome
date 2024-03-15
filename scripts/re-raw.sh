#!/bin/bash
# Name: re-raw BAM file
# Author: Eddy Mendoza
# Date: 29/01/23

# === SLURM Config ===
#SBATCH -J re-raw
#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --array=1-180%20
#SBATCH -o scripts/exegfiles/re-raw.out
#SBATCH -e scripts/exegfiles/re-raw.err

# === SET UP ENVIRONMENT ===

module load slurm_setup
source activate env_g
cd taller/ra42hux/

# === ARRAY INSTRUCTIONS ===

KEY="er_data/snps/key.txt"

line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY} ) 

# === WORK COMMANDS ===

name=$(echo $line | sed -E 's/.+am.(\w+)R.+/\1/g')
echo "WORKING WITH ${name}"

# split bam into 2 fastq files, removing non-paired and non-mapped reads

samtools fastq -@ 16 er_data/snps/bam/eR${name}.bam \
    -1 er_data/fastq/eR${name}_1.fastq.gz \
    -2 er_data/fastq/eR${name}_2.fastq.gz \
    -s er_data/fastq/eR${name}_1.fastq.gz \
    -0 er_data/fastq/eR${name}_1.fastq.gz \
    -n