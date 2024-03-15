#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(tibble)
setDTthreads(16)

# columns resulting from bedtools intersect: chr, pos, pos_rep, context, sample_i, sample_j, chr_dmr, start_dmr, end_dmr, id_dmr
# where the number of samples can be infinite (sample_i ... sample _n)

# function to process the data
surco_dmr = function(matrix) {
  
          # move and remove columns
          matrix = matrix %>% dplyr::relocate(c(chr_dmr, start_dmr, end_dmr, id_dmr), .before = context)
          matrix = matrix %>% dplyr::select(-c(chr, pos, pos_rep)) 
          
          # identify the columns with sample info
          # this is important so it doesn't matter how many samples we have
          sample_names = colnames(matrix[,-c(1:5)])

          # calculate mead methylation rate per DMR, iterate over samples
          matrix = as.data.table(matrix)
          matrix[, (sample_names) := lapply(.SD, mean), .SDcols = sample_names, by=.(id_dmr)]

          # get one row per dmr
          matrix = unique(matrix, by = "id_dmr")
          
          # make binary file, using the threshold of <0.15 to set unmethylated DMRs
          # this threshold was based on a visual inspection of the distribution of %
          bin = matrix[,-c(1:5)]
          bin[bin <= 0.15] = 0 
          bin[bin > 0.15] = 1
          bin = as.data.frame(bin)
          rownames(bin) = matrix$id_dmr
          
          # save one table for phylogenetic tree and one for binary matrix
          mat_bin = data.frame(matrix[,1:5], bin)
          phy_bin = as.data.frame(t(bin)) %>% mutate(sample = sample_names) %>% relocate(sample) %>% remove_rownames()
          gwas_mat = as.data.frame(t(matrix[,-c(1:5)])) %>% mutate(sample = sample_names) %>% relocate(sample) %>% remove_rownames()
          colnames(gwas_mat) = c("sample", matrix$id_dmr)
          
          # print results
          list(matrix, mat_bin, phy_bin, gwas_mat)
}

# READ INPUT FILE AND RUN THE FUNCTIONS

input = read.table("for_surco.tsv", header=T) # READ THE HEADER!!

tables = surco_dmr(input)

# Save output
write.table(tables[[1]], file = "mean_matrix.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(tables[[2]], file = "binary_matrix.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(tables[[3]], file = "gwas_binary_matrix.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(tables[[4]], file = "gwas_matrix.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

# File for the tree must not have a header
write.table(tables[[3]], file = "phylo.tsv", sep = "\t", col.names = F, row.names = F, quote = F)