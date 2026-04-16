library(TwoSampleMR)

input_dir <- "./data/"
output_dir <- "./results/Individual_MR/"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

exposure_file <- paste0(input_dir, "exposure.F.csv")
outcome_file <- paste0(input_dir, "outcome.csv")
gene_file <- paste0(output_dir, "../Consensus_Genes.txt") 
outcome_name <- "AAA"

rt <- read.csv(exposure_file, header = TRUE, sep = ",", check.names = FALSE)
sig_genes <- read.table(gene_file, header = FALSE, sep = "\t", check.names = FALSE)

for(gene_name in as.vector(sig_genes[, 1])) {
  
  safe_name <- gsub("\\%|\\/|\\:", "_", gene_name)
  single_exposure_file <- paste0(output_dir, safe_name, "_exposure.csv")
  
  exposure_set <- rt[rt$exposure == gene_name, ]
  write.csv(exposure_set, file = single_exposure_file, row.names = FALSE)
  
  exposure_dat <- read_exposure_data(
    filename = single_exposure_file,
    sep = ",",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    pval_col = "pval.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure",
    phenotype_col = "exposure",
    id_col = "id.exposure",
    samplesize_col = "samplesize.exposure",
    chr_col = "chr.exposure", 
    pos_col = "pos.exposure",
    clump = FALSE
  )
  
  outcome_data <- read_outcome_data(
    snps = exposure_dat$SNP,
    filename = outcome_file, 
    sep = ",",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome"
  )
  
  outcome_data$outcome <- outcome_name
  dat <- harmonise_data(exposure_dat, outcome_data)
  
  out_tab <- dat[dat$mr_keep == TRUE, ]
  write.csv(out_tab, file = paste0(output_dir, safe_name, "_table_SNP.csv"), row.names = FALSE)
  
  mr_result <- mr(dat)
  mr_tab <- generate_odds_ratios(mr_result)
  write.csv(mr_tab, file = paste0(output_dir, safe_name, "_table_MRresult.csv"), row.names = FALSE)
  
  heter_tab <- mr_heterogeneity(dat)
  write.csv(heter_tab, file = paste0(output_dir, safe_name, "_table_heterogeneity.csv"), row.names = FALSE)
  
  pleio_tab <- mr_pleiotropy_test(dat)
  write.csv(pleio_tab, file = paste0(output_dir, safe_name, "_table_pleiotropy.csv"), row.names = FALSE)
  
  pdf(file = paste0(output_dir, safe_name, "_scatter_plot.pdf"), width = 7, height = 6.5)
  print(mr_scatter_plot(mr_result, dat))
  dev.off()
  
  res_single <- mr_singlesnp(dat)      
  pdf(file = paste0(output_dir, safe_name, "_forest.pdf"), width = 6.5, height = 5)
  print(mr_forest_plot(res_single))
  dev.off()
  
  pdf(file = paste0(output_dir, safe_name, "_funnel_plot.pdf"), width = 6.5, height = 6)
  print(mr_funnel_plot(singlesnp_results = res_single))
  dev.off()
  
  pdf(file = paste0(output_dir, safe_name, "_leaveoneout.pdf"), width = 6.5, height = 5)
  print(mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat)))
  dev.off()
  
}

writeLines(capture.output(sessionInfo()), paste0(output_dir, "SessionInfo.txt"))
