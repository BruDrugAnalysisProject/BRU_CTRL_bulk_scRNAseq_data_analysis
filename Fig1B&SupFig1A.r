rm(list = ls())

setwd("~/bulk/1.result_figs/bulk_TCGALIHC_survival_analysis")

library(tidyverse)
library(survival)
library(survminer)
library(tidyestimate)

source("./functions.r")

# load data

exp <- readRDS("~/public_data/TCGADownload/TCGAGDC_UCSCXena_preprocessed/TCGA-LIHC_explog2_add1_processed.rds") %>%
  as.data.frame()
sur_phen <- readRDS("~/public_data/TCGADownload/TCGAGDC_UCSCXena_preprocessed/TCGA-LIHC_sur_phen_processed.rds") %>%
  as.data.frame()
estimate <- readRDS("~/public_data/TCGADownload/TCGAGDC_UCSCXena_preprocessed/TCGA-LIHC_estimate_deconv.rds") %>%
  t() %>% 
  as.data.frame()
cd8_timer <- readRDS("~/public_data/TCGADownload/TCGAGDC_UCSCXena_preprocessed/TCGA-LIHC_timer_deconv.rds") %>%
  as.data.frame() %>%
  column_to_rownames(var = "cell_type") %>%
  t() %>%
  as.data.frame()
cd8_cibersort <- readRDS("~/public_data/TCGADownload/TCGAGDC_UCSCXena_preprocessed/TCGA-LIHC_cibersort_deconv.rds")%>%
  as.data.frame()
cd8_xcell <- readRDS("~/public_data/TCGADownload/TCGAGDC_UCSCXena_preprocessed/TCGA-LIHC_xcell_deconv.rds")%>%
  as.data.frame() %>%
  column_to_rownames(var = "cell_type") %>%
  t() %>%
  as.data.frame()

tumor_exp <- exp[,rownames(sur_phen)[sur_phen$sample_type.samples == "Primary Tumor"]]
dim(tumor_exp)
normal_exp <- exp[,rownames(sur_phen)[sur_phen$sample_type.samples == "Solid Tissue Normal"]]
dim(normal_exp)

# adjust NFE2L2(NRF2) expression level in tumor of bulk RNA-seq data -- TCGA-LIHC

gene_NFE2L2_corrected <- correct_gene_expression(gene_name = "NFE2L2", expr_matrix = tumor_exp, purity_df = estimate, normal_expr = normal_exp)
gene_NFE2L2_corrected <- gene_NFE2L2_corrected %>%
  column_to_rownames(var = "sample_id")

# prepare survival analysis data

cp_ge_sur_matrix <- gene_NFE2L2_corrected
cp_ge_sur_matrix$CD8A <- as.numeric(tumor_exp["CD8A",rownames(cp_ge_sur_matrix)])
cp_ge_sur_matrix$CD8Tc_TIMER <- as.numeric(cd8_timer[rownames(cp_ge_sur_matrix),"T cell CD8+"])
cp_ge_sur_matrix$CD8Tc_CIBERSORT <- as.numeric(cd8_cibersort[rownames(cp_ge_sur_matrix),"T cells CD8"])
cp_ge_sur_matrix$CD8Tc_xCell <- as.numeric(cd8_xcell[rownames(cp_ge_sur_matrix),"T cell CD8+"])
colnames(cp_ge_sur_matrix)[2] <- "NFE2L2"
cp_ge_sur_matrix$OS.time <- sur_phen[rownames(cp_ge_sur_matrix),]$OS.time
cp_ge_sur_matrix$OS <- sur_phen[rownames(cp_ge_sur_matrix),]$OS
cp_ge_sur_matrix_bak <- cp_ge_sur_matrix
cp_ge_sur_matrix <- cp_ge_sur_matrix[!str_detect(rownames(cp_ge_sur_matrix),"-01B"),] # filter 01B sample
cp_ge_sur_matrix <- cp_ge_sur_matrix[cp_ge_sur_matrix$OS.time>1,] # filter survival time 1 day

# survival analysis
color_list <- c("#DC2819","#4d4e75")
survival_result_single <- km_survival_single(gene = "NFE2L2",cutoff = 50,dat = cp_ge_sur_matrix,auto_cutoff = FALSE, col1 = color_list[1], col2 = color_list[2])
pdf(paste0("sur_TCGALIHC_NFE2L2_adjust_TP.pdf"),width = 2.8,height = 3.5,onefile = FALSE)
print(survival_result_single)
dev.off()

survival_result_single_na <- km_survival_single(gene = "origial_expr",cutoff = 50,dat = cp_ge_sur_matrix,auto_cutoff = FALSE, col1 = color_list[1], col2 = color_list[2])
pdf(paste0("sur_TCGALIHC_NFE2L2_notadjust_TP.pdf"),width = 2.8,height = 3.5,onefile = FALSE)
print(survival_result_single_na)
dev.off()

for(i in c("CD8A","CD8Tc_TIMER","CD8Tc_CIBERSORT","CD8Tc_xCell")){
  survival_results <- km_survival(genea = "NFE2L2",geneb = i,cutoff = 50,dat = cp_ge_sur_matrix, auto_cutoff = FALSE,col1 = color_list[1],col2 = color_list[2])
  pdf(paste0("sur_TCGALIHC_",i,"_high_","NFE2L2_adjust_TP.pdf"),width = 2.8,height = 3.5,onefile = FALSE)
  print(survival_results[[1]])
  dev.off()
  pdf(paste0("sur_TCGALIHC_",i,"_low_","NFE2L2_adjust_TP.pdf"),width = 2.8,height = 3.5,onefile = FALSE)
  print(survival_results[[2]])
  dev.off()
  print(i)
}

