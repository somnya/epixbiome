#!/bin/bash
# Author: Eddy Mendoza
# Date: 06/02/2024

#########################
#------MethylRille------#
#########################

##### AIM: Process MethylScore output to get the raw methylation calls per sample per context per DMR
##### This can be done for different matrices and the same bed file and viceversa

# === SLURM Config ===
#SBATCH -J rille
#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -o scripts/exegfiles/rille_%j.out
#SBATCH -e scripts/exegfiles/rille_%j.err

# === SET UP ENVIRONMENT ===

module load slurm_setup
source activate env_g
cd taller/ra42hux/

MAT=er_data/all/03matrix # MethylScore 03matrix output directory, no slash
DMRS=er_data/exp/05DMRs # MethylScore 05DMRs output directory, no slash
OUT=er_data/all # Output directory, no slash

# === WORK COMMANDS ===

# Clean and split the raw mC calling matrix

rm -R $OUT/rille/
mkdir $OUT/rille

gunzip -c $MAT/genome_matrix.tsv.gz > $OUT/rille/int_genome_matrix.tsv
cat $OUT/rille/int_genome_matrix.tsv | sed -E 's/\w+\/(\w+\.\w+)\/\w+\/\w+/\1/g' > $OUT/rille/mc_rate.tsv

#### for each context:
#### Split the matrix and convert it to bed
#### Duplicate the second column and rename header 
#### the header includes the sample info, always keep!

for CNTXT in $(echo "CG" "CHG" "CHH")
do

echo "Working with ${CNTXT} information"

paste \
        <(awk -v CNTXT="$CNTXT" '$3 == CNTXT || $3 == "class" ' $OUT/rille/mc_rate.tsv | \
                awk '{print $1, $2, $2, $3 }' | sed 's/#//g' | sed 's/pos class/pos_rep context/g') \
        <(awk -v CNTXT="$CNTXT" '$3 == CNTXT || $3 == "class" ' $OUT/rille/mc_rate.tsv | cut -f 1-4 --complement) \
         | sed -E 's/ /\t/g' > $OUT/rille/int_mc_rate_${CNTXT}.bed 
         

#### Clean the DMR bed files
ND=$(wc -l $DMRS/all/DMRs.${CNTXT}.bed | tr ' ' '\n' | head -1)
# awk '{ printf $0 "\t" "dmr_%d_CG\n", NR }'
cat $DMRS/all/DMRs.${CNTXT}.bed | cut -f 1-3 | awk -v cntxt="$CNTXT" '{ printf $0 "\t" "dmr_%d_" cntxt "\n", NR }' > $OUT/rille/int_dmrs_${CNTXT}.bed


#### Get the intersections and include the header
bedtools intersect -a <(tail -n +2 $OUT/rille/int_mc_rate_${CNTXT}.bed) -b $OUT/rille/int_dmrs_${CNTXT}.bed -wb | \
    cat <( paste <(head -1 $OUT/rille/int_mc_rate_${CNTXT}.bed) <(echo -e "chr_dmr\tstart_dmr\tend_dmr\tid_dmr") ) - \
    > $OUT/rille/int_for_surco_${CNTXT}.tsv

echo "Done working with ${CNTXT}"
echo "--------------------------"
    
done

# Merge BED files
cat $OUT/rille/int_dmrs_CG.bed $OUT/rille/int_dmrs_CHG.bed $OUT/rille/int_dmrs_CHH.bed > $OUT/rille/all_DMRs.bed   

# Merge matrix bed files
cat $OUT/rille/int_for_surco_CG.tsv \
    <(tail -n +2 $OUT/rille/int_for_surco_CHG.tsv) \
    <(tail -n +2 $OUT/rille/int_for_surco_CHH.tsv) > $OUT/rille/for_surco.tsv
# All three files have a header so I remove them from the last 2

# clean 
rm $OUT/rille/int*

# Run the Rscript and get the mean tables

cd $OUT/rille/
export TMPDIR=./ # set tmp dir or R won't work
surco.R

# Then make the tree
perl -S tab_to_phy.pl "phylo.tsv" 
iqtree -nt 16 -s phylo.tsv.fa -B 1000