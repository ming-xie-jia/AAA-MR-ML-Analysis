library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)

exposure_file <- "eQTL_exposure.csv"
outcome_file <- "finngen_R12_I9_ABAORTANEUR"
outcome_name <- "AAA_eQTL"

exposure_dat <- read_exposure_data(
  filename = exposure_file,
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

outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outcome_file, 
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval",
  eaf_col = "af_alt"
)
write.csv(outcome_dat, file = "outcome_eQTL.csv", row.names = FALSE)

outcome_dat$outcome <- outcome_name
dat <- harmonise_data(exposure_dat, outcome_dat)
dat <- dat[dat$pval.outcome > 5e-06, ]

out_tab <- dat[dat$mr_keep == TRUE, ]
write.csv(out_tab, file = "table.SNP_eQTL.csv", row.names = FALSE)

mr_result <- mr(dat)
mr_tab <- generate_odds_ratios(mr_result)
write.csv(mr_tab, file = "table.MRresult_eQTL.csv", row.names = FALSE)

heter_tab <- mr_heterogeneity(dat)
write.csv(heter_tab, file = "table.heterogeneity_eQTL.csv", row.names = FALSE)

pleio_tab <- mr_pleiotropy_test(dat)
write.csv(pleio_tab, file = "table.pleiotropy_eQTL.csv", row.names = FALSE)

pdf(file = "pic.scatter_plot_eQTL.pdf", width = 7.5, height = 7)
mr_scatter_plot(mr_result, dat)
dev.off()

res_single <- mr_singlesnp(dat)      

pdf(file = "pic.forest_eQTL.pdf", width = 7, height = 5.5)
mr_forest_plot(res_single)
dev.off()

pdf(file = "pic.funnel_plot_eQTL.pdf", width = 7, height = 6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

pdf(file = "pic.leaveoneout_eQTL.pdf", width = 7, height = 5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()
