#!/bin/bash
# Name: revelio
# Author: Eddy Mendoza
# Date: 02/02/24

# === SLURM Config ===
#SBATCH -J revelio
#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --array=1-169%13
#SBATCH -o scripts/exegfiles/revelio.out
#SBATCH -e scripts/exegfiles/revelio.err

# === SET UP ENVIRONMENT ===

module load slurm_setup
source activate env_g
cd taller/ra42hux/

# === ARRAY INSTRUCTIONS ===

KEY="er_data/tes/key.txt"

line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY} )


# === WORK COMMANDS ===

python er_data/snps_rev_bay/revelio.py \
    -f ref/genome_CM.fa \
    -t work/ \
    --threads 16 \
    er_data/snps/bam/${line}.bam \
    er_data/snps_rev_bay/masked/${line}.bam