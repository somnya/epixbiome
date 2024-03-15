#!/bin/bash
# Name: freebayes
# Author: Eddy Mendoza
# Date: 02/02/24

# === SLURM Config ===
#SBATCH -J snp_bayes
#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --array=1-169%13
#SBATCH -o scripts/exegfiles/freebayes.out
#SBATCH -e scripts/exegfiles/freebayes.err

# === SET UP ENVIRONMENT ===

module load slurm_setup
source activate env_g
cd taller/ra42hux/

# === ARRAY INSTRUCTIONS ===

KEY="er_data/tes/key.txt"

line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY} )

# === JOB ===

#freebayes \
#    -f ref/genome_CM.fa \
#    -q 20 -m 20 -C 4 --min-coverage 4 \
#    -L er_data/snps_rev_bay/bam.list > er_data/snps_rev_bay/vcf/freebayes_raw.vcf
# freebayes alone runs in one thread, which would take forever !!!!

# run multithread parallel jobs for each sample

freebayes-parallel <(fasta_generate_regions.py ref/genome_CM.fa 10000) 32 \
    -f ref/genome_CM.fa \
    -q 20 -m 20 -C 4 --min-coverage 4 \
    er_data/snps_rev_bay/masked/sorted_${line}.bam > er_data/snps_rev_bay/vcf/${line}.vcf

# compress and index
# bgzip --threads 32 er_data/snps_rev_bay/vcf/${line}.vcf
# bcftools index --threads 32 er_data/snps_rev_bay/vcf/${line}.vcf.gz