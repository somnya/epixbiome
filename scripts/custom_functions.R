# function to calculate the mean based on a BETA distribution

mean_beta = function(vector) {
  
  # In case the DMR has no methylation
  if (sum(vector) == 0) {
    print(0)
    
    # If there is methylation     
  } else {
    vector = as.numeric(vector)
    vector = round(vector/100, digits = 2)
    vector = vector[!is.na(vector)]
    #using the formula {alpha / (alpha + beta)} for the mean of a beta distribution
    as.numeric(ebeta(vector)$parameters[1] / (ebeta(vector)$parameters[1] + ebeta(vector)$parameters[2]) )
  }
}

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
  # for(i in sample_names) matrix[, (i) := mean(get(i)), by=.(id_dmr)] 
  # the last achieves the same results, but is probably slower when many samples
  
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


# To calculate pi

pi_epi = function(binary_table, b_size){
  
  # this function requires the binary table to be genome coordinate sorted
  binary_table = as.data.frame(binary_table)
  
  # assing blocks based on a value (how many DMRs per block)
  n_blocks = ceiling(nrow(binary_table) / b_size) # round to the nearest integer
  binary_table = binary_table %>% mutate(block = rep(1:n_blocks, each=b_size)[1:nrow(binary_table)])
  
  # create empty list to save the values
  genome_epi_D = list()
  
  # iterate the D calculation over each block
  for(b in 1:n_blocks) {
    
    # select data and remove categorical columns
    df = filter(binary_table, block==b)[,-c(1:5)] %>% dplyr::select(-c("block"))
    
    # set variables
    n = ncol(df)
    S = nrow(df)
    
    # create empty matrix
    pairwise_differences = matrix(0, nrow = n, ncol = n)
    
    # compute pairwise differences
    for (i in 1:n ) {
      for (j in 1:n ) {
        pairwise_differences[i, j] = sum(df[,i] != df[,j])
      }
    }
    
    # Remove lower triangle
    pairwise_differences[lower.tri(pairwise_differences, diag = T)] = 0
    
    # calcule pi and w
    epi_pi = sum(pairwise_differences)/((n * (n - 1)) / 2)
    
    # get block info
    info = filter(binary_table, block==b) %>% group_by(block) %>% summarise(chr = min(chr_dmr),
                                                                            start = min(start_dmr), # this basically takes the bigger chr number in case the block is between two
                                                                            chr_end = max(chr_dmr),
                                                                            end = max(end_dmr),
                                                                            mean_length = mean(end_dmr-start_dmr),
                                                                            context) %>% 
      ungroup() %>% unique()
    
    # save the values and info
    genome_epi_D[[b]] = data.frame(info,
                                   epi_pi = epi_pi)
  }
  
  # print all values
  genome = do.call(rbind, genome_epi_D)
  genome
  
}


#### function to calculate the effect of each factor in beta and alpha diversity values

diversity_factor_lm = function(diversity_data) {
    # this requires each measurement in one column and the columns of the experimental design
  data.frame(p_val = c(
    round(glance(lm(phenotype ~ sample, data = diversity_data))$p.value, digits = 3),
    round(glance(lm(phenotype ~ block, data = diversity_data))$p.value, digits = 3),
    round(glance(lm(phenotype ~ position, data = diversity_data))$p.value, digits = 3),
    round(glance(lm(phenotype ~ harv_date, data = diversity_data))$p.value, digits = 3),
    round(glance(lm(phenotype ~ harv_pers, data = diversity_data))$p.value, digits = 3)
  ),
  factor = c("sample", "block", "position", "harv_date", "harv_pers")
  )
}

#### Functions to calculate heritability and other contributions from a Random Effects model

# to get the contributions, only for epigenotypes

heritator = function(model) {
  # get the variance from the model
  int = as.data.frame(VarCorr(model))[, c("grp", "vcov")] 
  
  # get variance contributions
  heritability = round(int[int$grp == "sample", 2]*100 / sum(int$vcov), digits = 0)
  residual = round(int[int$grp == "Residual", 2]*100 / sum(int$vcov), digits = 0)
  experimental = round(sum(int[int$grp != c("sample", "Residual"), 2])*100 / sum(int$vcov), digits = 0)
  
  # extract the name of the independent variable tested
  name_y = as.character(model@call$formula[[2]])
  
  # merge data and print to output
  data.frame(variable = c("heritability", "experimental", "error"),
             per_variance = c(heritability, experimental, residual)) %>% 
    mutate(name = name_y)
}

# genome scan and testing

DMR_H = function(model_data) {
  
  # heritability analysis
  model = lmer(phenotype ~ (1|mC) + (1|pool) + (1|block) + (1|position) + (1|harv_date) + (1|harv_pers), data = model_data, REML = T)
  heritability = DMR_heritator(model, model_data) # heritator reuires the model_data object
  heritability = heritability %>% mutate(id_dmr = model_data$id_dmr[1], phenotype = model_data$pheno_name[1])
  heritability
}

DMR_heritator = function(model, model_data) {
  # this assumes there is a "model_data" object in the temp environment
  
  # get the variance from the model
  int = as.data.frame(VarCorr(model))[, c("grp", "vcov")] 
  
  # get variance contributions
  heritability = round(int[int$grp == "mC", 2]*100 / sum(int$vcov), digits = 0)
  residual = her = round(int[int$grp == "Residual", 2]*100 / sum(int$vcov), digits = 0)
  experimental = 100 - heritability - residual
  
  # perform a loglikelihood test against a reduced model
  reduced_model = lmer(phenotype ~ + (1|pool) + (1|block) + (1|position) + (1|harv_date) + (1|harv_pers), data = model_data, REML = T)
  ratio = 2*abs(as.numeric(logLik(reduced_model)) - as.numeric(logLik(model)))
  ratio_test_p = pchisq(ratio, df = 1, lower.tail = FALSE)
  
  # evaluate residuals
  residuals_norm_p = shapiro.test(residuals(model))$p.val
  
  # evaluate heteroscedasticity
  het_test_p = as.double(check_heteroscedasticity(model))
  
  # merge data and print to output
  data.frame(heritability = heritability,
             experimental = experimental,
             residuals = residual,
             ratio_test_p = ratio_test_p,
             residuals_norm_p = residuals_norm_p,
             het_test_p = het_test_p)
}


# genome scan, also for data.table

DMR_GWAS = function(model_data) {

  # assuming the data is split into phenotype x epigenotype combinations
  model = lmer(phenotype ~ mC + pool + (1|block) + (1|position) + (1|harv_date) + (1|harv_pers), data = model_data, REML = T)
  
  # test model and save value, using a Chi-Squared test
  test_p = car::Anova(model)[1,3]
  
  # evaluate residuals
  residuals_norm_p = shapiro.test(residuals(model))$p.val
  
  # evaluate heteroscedasticity
  het_test_p = as.double(check_heteroscedasticity(model))

  # merge data and print to output
  data.frame(pheno_name = model_data$pheno_name[[1]],
             test_p = test_p,
             residuals_norm_p = residuals_norm_p,
             het_test_p = het_test_p
             )
}


epi_H = function(model_data, formula) {
  
  # heritability analysis, including all factors
  # this one uses the variable "sample" as the epigenotype
  
  model = lmer(formula = paste("phenotype", formula, sep = " "), data = model_data, REML = T)
  
  # get the variance from the model
  int = as.data.frame(VarCorr(model))[, c("grp", "vcov")] 
  
  # get variance contributions
  heritability = round(int[int$grp == "sample", 2]*100 / sum(int$vcov), digits = 0)
  residual = her = round(int[int$grp == "Residual", 2]*100 / sum(int$vcov), digits = 0)
  experimental = 100 - heritability - residual
  
  # perform a loglikelihood test against a reduced model
  reduced_formula = str_replace(formula, "\\(1\\|sample\\) [+] ", "")
  reduced_model = lmer(formula = paste("phenotype", reduced_formula, sep = " "), data = model_data, REML = T)
  ratio = 2*abs(as.numeric(logLik(reduced_model)) - as.numeric(logLik(model)))
  ratio_test_p = pchisq(ratio, df = 1, lower.tail = FALSE)
  
  # evaluate residuals
  residuals_norm_p = shapiro.test(residuals(model))$p.val
  
  # evaluate heteroscedasticity
  het_test_p = as.double(check_heteroscedasticity(model))
  
  # merge data and print to output
  data.frame(heritability = heritability,
             experimental = experimental,
             residuals = residual,
             ratio_test_p = ratio_test_p,
             residuals_norm_p = residuals_norm_p,
             het_test_p = het_test_p)

}


epi_AS = function(model_data, formula) {
  
  # assuming the data is split into phenotype x epigenotype combinations
  model = lmer(paste("phenotype", formula, sep = " "), data = model_data, REML = T)
  
  # test model and save value, using a Chi-Squared test
  test_p = car::Anova(model)[1,3]
  
  # evaluate residuals
  residuals_norm_p = shapiro.test(residuals(model))$p.val
  
  # evaluate heteroscedasticity
  het_test_p = as.double(check_heteroscedasticity(model))
  
  # merge data and print to output
  # not necessary to add the pheno_name as data.frame will keep it
  data.frame(test_p = test_p,
             residuals_norm_p = residuals_norm_p,
             het_test_p = het_test_p
  )
}

#### 
# compute alpha and diversity, using relative abundance and counts

alpha_beta_div = function(phyloseq_object) {
  # get relative abundances
  tss = otu_table(as.data.frame(lapply(as.data.frame(phyloseq_object@otu_table), function(x) x*100/sum(x)) ), taxa_are_rows =F)
  rownames(tss) = rownames(phyloseq_object@otu_table)
  
  # get alpha div estimates, Chao1 and Fisher don't accept digits
  tss_alpha = estimate_richness(tss, 
                          measures=c("Shannon", "Simpson")) %>% 
    mutate(ID = rownames(as.data.frame(phyloseq_object@otu_table))) %>% 
    relocate(ID) %>%
    left_join(phyloseq_object@sam_data, by="ID") %>%
    dplyr::rename("TSS_Shannon" = "Shannon",
                  "TSS_Simpson" = "Simpson")
  
  # now without the transformation
  alpha = estimate_richness(phyloseq_object@otu_table, 
                            measures=c("Shannon", "Simpson", "Chao1", "Fisher")) %>% 
    mutate(ID = rownames(as.data.frame(phyloseq_object@otu_table))) %>%
    dplyr::select(-c(se.chao1))
  
  # calculate beta diversity, two methods, two components each and both relative and raw counts
  # Bray Curtis 
  beta = as.data.frame(ordinate(phyloseq_object, method = "MDS", distance = "bray")$vectors[,1:2]) %>% 
    rownames_to_column("ID") %>% 
    dplyr::rename("Bray-Curtis_MDS1" = "Axis.1",
                  "Bray-Curtis_MDS2" = "Axis.2") %>% left_join(
                    
                    # DCA with phylogenetic-aware distances
                    as.data.frame(ordinate(phyloseq_object, distance = "wunifraq")[["rproj"]])[,1:2] %>% 
                      rownames_to_column("ID") %>% 
                      dplyr::rename("W-UniFrac_DCA1"="DCA1",
                                    "W-UniFrac_DCA2"="DCA2"), by = "ID") %>% left_join(
                                      
                                      # Now with relative abundances
                                      as.data.frame(ordinate(tss, method = "MDS", distance = "bray")$vectors[,1:2]) %>% 
                                        rownames_to_column("ID") %>% 
                                        dplyr::rename("TSS_Bray-Curtis_MDS1" = "Axis.1",
                                                      "TSS_Bray-Curtis_MDS2" = "Axis.2"), by = "ID") %>% left_join(
                                                        
                                                        as.data.frame(ordinate(tss, distance = "wunifraq")[["rproj"]])[,1:2] %>% 
                                                          rownames_to_column("ID") %>% 
                                                          dplyr::rename("TSS_W-UniFrac_DCA1"="DCA1",
                                                                        "TSS_W-UniFrac_DCA2"="DCA2"), by = "ID")
  # merge all
  tss_alpha %>% left_join(alpha, by="ID") %>% left_join(beta, by="ID")
}

alpha_beta_div_its = function(phyloseq_object) {
  # get relative abundances
  tss = otu_table(as.data.frame(lapply(as.data.frame(phyloseq_object@otu_table), function(x) x*100/sum(x)) ), taxa_are_rows =F)
  rownames(tss) = rownames(phyloseq_object@otu_table)
  
  # get alpha div estimates, Chao1 and Fisher don't accept digits
  tss_alpha = estimate_richness(tss, 
                                measures=c("Shannon", "Simpson")) %>% 
    mutate(ID = rownames(as.data.frame(phyloseq_object@otu_table))) %>% 
    relocate(ID) %>%
    left_join(phyloseq_object@sam_data, by="ID") %>%
    dplyr::rename("TSS_Shannon" = "Shannon",
                  "TSS_Simpson" = "Simpson")
  
  # now without the transformation
  alpha = estimate_richness(phyloseq_object@otu_table, 
                            measures=c("Shannon", "Simpson")) %>% 
    mutate(ID = rownames(as.data.frame(phyloseq_object@otu_table))) 
  
  # calculate beta diversity, two methods, two components each and both relative and raw counts
  # Bray Curtis 
  beta = as.data.frame(ordinate(phyloseq_object, method = "MDS", distance = "bray")$vectors[,1:2]) %>% 
    rownames_to_column("ID") %>% 
    dplyr::rename("Bray-Curtis_MDS1" = "Axis.1",
                  "Bray-Curtis_MDS2" = "Axis.2") %>% left_join(
                    
                    # DCA with phylogenetic-aware distances
                    as.data.frame(ordinate(phyloseq_object, distance = "wunifraq")[["rproj"]])[,1:2] %>% 
                      rownames_to_column("ID") %>% 
                      dplyr::rename("W-UniFrac_DCA1"="DCA1",
                                    "W-UniFrac_DCA2"="DCA2"), by = "ID") %>% left_join(
                                      
                                      # Now with relative abundances
                                      as.data.frame(ordinate(tss, method = "MDS", distance = "bray")$vectors[,1:2]) %>% 
                                        rownames_to_column("ID") %>% 
                                        dplyr::rename("TSS_Bray-Curtis_MDS1" = "Axis.1",
                                                      "TSS_Bray-Curtis_MDS2" = "Axis.2"), by = "ID") %>% left_join(
                                                        
                                                        as.data.frame(ordinate(tss, distance = "wunifraq")[["rproj"]])[,1:2] %>% 
                                                          rownames_to_column("ID") %>% 
                                                          dplyr::rename("TSS_W-UniFrac_DCA1"="DCA1",
                                                                        "TSS_W-UniFrac_DCA2"="DCA2"), by = "ID")
  # merge all
  tss_alpha %>% left_join(alpha, by="ID") %>% left_join(beta, by="ID")
}


#####
# plot results from heritability

post_H = function(alpha_table, threshold, colin){
  # colin is the palette for valid/non-valid model distinction
  
  # make a vector to arrange the order later
  top_pheno = dplyr::arrange(alpha_table, heritability)$pheno_name
  # transform
  table = alpha_table %>% mutate(ratio_test_p = round(ratio_test_p, digits = 3),
                                 p_val = ifelse(ratio_test_p == 0, "P < 0.001", paste("P = ", ratio_test_p, sep = "") ),
                                 performance = ifelse(residuals_norm_p > threshold & het_test_p > threshold, "valid", "biased")
  ) %>%
    
    melt(id.vars = c("pheno_name", "p_val", "performance"), 
         measure.vars = c("heritability", "residuals", "experimental") ) %>% 
    mutate(value = ifelse(value==0, NA, value))
  # apply order
  table$variable = factor(table$variable, levels = c("residuals", "experimental", "heritability"))
  table$pheno_name = factor(table$pheno_name, levels = top_pheno)
  
  # plot
  ggplot(table) +
    geom_bar(aes(y=pheno_name, x=value, fill=variable), stat = "identity", color="black", width = 0.6) +
    geom_text(aes(y=pheno_name, x=value, label = value), size = 4, position = position_stack(vjust = 0.5), color="white") +
    geom_text(aes(y=pheno_name, x=110, label=p_val, color=performance), vjust=0) +
    scale_fill_manual(values = palette_H,
                      limits = c("heritability", "experimental", "residuals")) +
    scale_color_manual(values = colin) +
    scale_x_continuous(breaks = c(0, 50, 100)) +
    labs(x= "Variance explained (%)", y= NULL) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text = element_text(size=20),
          axis.title = element_text(size=20),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.x.bottom = element_line(color = 'black')
    ) 
}

# evaluate and plot Association analysis

post_AS = function(input_table, threshold, colin, names=TRUE){
  
  # colin is the palette for valid/non-valid model distinction
  
  # compute bias
  table = input_table %>% mutate(performance = ifelse(residuals_norm_p > threshold & het_test_p > threshold, "valid", "non-valid"))
  
  # plot
  
  if (names == TRUE) {
    ggplot(table) +
      geom_point(aes(y=pheno_name, x=-log10(test_p), color=performance), size=8) +
      geom_segment(aes(y=pheno_name, yend=pheno_name, x=0, xend=-log10(test_p), color=performance), size=4) +
      geom_vline(xintercept = -log10(threshold), size=1, color="gray40", linetype = "dashed")  +
      scale_color_manual(values = colin) +
      labs(x= expression(-log[10](P)), y= NULL) +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(size=20),
            axis.title = element_text(size=20),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line.y = element_line(color = 'black'),
            axis.line.x.bottom = element_line(color = 'black')
      ) 
    
  } else {
    ggplot(table) +
      geom_point(aes(y=pheno_name, x=-log10(test_p), color=performance), size=8) +
      geom_segment(aes(y=pheno_name, yend=pheno_name, x=0, xend=-log10(test_p), color=performance), size=4) +
      geom_vline(xintercept = -log10(threshold), size=1, color="gray40", linetype = "dashed")  +
      scale_color_manual(values = colin) +  
      labs(x= expression(-log[10](P)), y= NULL) +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(size=20),
            axis.title = element_text(size=20),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line.y = element_line(color = 'black'),
            axis.line.x.bottom = element_line(color = 'black')
      ) 
  }
  
}

### Variance Stabilization to OTU counts
VSTPhy = function(phyloseq_object) {
  table = t(DESeq2::varianceStabilizingTransformation( t(as.matrix(as.data.frame(phyloseq_object@otu_table)))+1))
  as.data.frame(table) %>% rownames_to_column("ID")
}

#### For the selection analysis, calculate mean phenotype values and SE from a LMM
LS_estimate = function(table) {
  as.data.frame(emmeans(lmer(phenotype ~ sample + (1|position) + (1|block) + (1|harv_date) + (1|harv_pers), data = table, REML = TRUE), specs = "sample"))[,1:2]
}
# when including the sequencing run
LS_estimate_pool = function(table) {
  as.data.frame(emmeans(lmer(phenotype ~ sample * pool+ (1|position) + (1|block) + (1|harv_date) + (1|harv_pers), data = table, REML = TRUE), specs = "sample"))[,1:2]
}



##### Selection Entropy

s_entropy = function(p, q, n) {
  # calculate shannon entropy from the selection contrast
  if (any(p == 0 | q == 0)) {
    0
  } else {
    (p/n)*log2(1/(p/n))+(q/n)*log2(1/(q/n))
  }
}


permutation_entropy = function(in_entropy, n) {
  # perform a permutation test to see if the observed entropy is higher than what is expected by chance
  
  if (any(in_entropy == 0)) {
    1
  } else {
    comparisons = vector(length = 1000)
    for (i in 1:1000) {
      random = sample(c("positive", "negative", NA), n, replace = T)
      p = length(which(random == "positive"))
      q = length(which(random == "negative"))
      comparisons[i] = in_entropy > s_entropy(p,q,n)
    }
    length(which(comparisons == TRUE))/1000
  }
}