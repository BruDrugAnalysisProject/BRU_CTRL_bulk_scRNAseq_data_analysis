rm(list = ls())

library(colocr)
library(Seurat)
library(tidyverse)
library(EBImage)
library(spdep)
library(tidyverse)
library(raster)
library(spatialEco)

source("~functions.r")

load_packages()

setwd("~/mIF/stain_images")

img_file_list <- list.files("./",pattern = "-4.jpg|-5.jpg|-6.jpg|-4-PD-L1.jpg|-5-Nrf2.jpg|-6-CD8.jpg")
img_file_list
names(img_file_list) <- c(rep("BRV-2mg-3-1",3),
                          rep("BRV-2mg-3-2",3),
                          rep("BRV-2mg-3-3",3),
                          rep("BRV-2mg-3-4",3),
                          rep("BRV-3-1",3),
                          rep("BRV-5-1",3),
                          rep("BRV-5-2",3),
                          rep("BRV-5-3",3),
                          rep("ctrl-2-1",3),
                          rep("ctrl-2-2",3),
                          rep("ctrl-5-1",3),
                          rep("ctrl-5-2",3),
                          rep("ctrl-5-3",3),
                          rep("ctrl-5-4",3),
                          rep("ctrl-5-5",3),
                          rep("ctrl-5-6",3),
                          rep("ctrl-5-7",3),
                          rep("ctrl-5-8",3),
                          rep("ctrl-5-9",3),
                          rep("ctrl-8-1",3))

for(i in unique(names(img_file_list))){
  file_paths <- img_file_list[names(img_file_list)==i]
  protein_names <- c("PD-L1","NRF2","CD8a")
  channels <- read_protein_channels(
    file_paths = file_paths,
    protein_names = protein_names,
    remove_scale = TRUE,
    scale_height = 100,  # 例如：50像素
    scale_position = "top"
  )
  df <- create_protein_expression_df(
    channels = channels,
    sample_fraction = 1
  )
  processed_data <- preprocess_expression_data(
    df = df,
    remove_background = TRUE,
    normalize = TRUE
  )
  df_processed <- processed_data$processed_data
  pdl1_nrf2_cor <- cor.test(df_processed$`PD-L1`,df_processed$`NRF2`)
  pdl1_cd8a_cor <- cor.test(df_processed$`PD-L1`,df_processed$`CD8a`)
  cd8a_nrf2_cor <- cor.test(df_processed$`CD8a`,df_processed$`NRF2`)
  if(i == unique(names(img_file_list))[1]){
    results_cor_df <- data.frame("sample_id" = i,
                                 "pdl1_nrf2_cor" = pdl1_nrf2_cor$estimate,
                                 "pdl1_cd8a_cor" = pdl1_cd8a_cor$estimate,
                                 "cd8a_nrf2_cor" = cd8a_nrf2_cor$estimate)
  }else{
    results_cor_df <- rbind(results_cor_df,
                            data.frame("sample_id" = i,
                                      "pdl1_nrf2_cor" = pdl1_nrf2_cor$estimate,
                                      "pdl1_cd8a_cor" = pdl1_cd8a_cor$estimate,
                                      "cd8a_nrf2_cor" = cd8a_nrf2_cor$estimate))
  }
  print(i)
}

View(results_cor_df)

### plot results

results_cor_df$group <- ifelse(str_detect(results_cor_df$sample_id,"BRV"),"BRU-treated","Control")
colnames(results_cor_df)
results_cor_df$group <- factor(results_cor_df$group,levels = c("Control","BRU-treated"))

for(cor_i in colnames(results_cor_df)[2:4]){
  print(cor_i)
  results_cor_df$tmp_col <- results_cor_df[,colnames(results_cor_df)%in%cor_i]
  if(cor_i == "pdl1_nrf2_cor"){
    p_cor <- ggboxplot(results_cor_df, x = "group", y = "tmp_col",
                       color = "group",add = "jitter",add.params = list(size = 0.5),
                       width = 0.5,ylab = "Pearson correlation coefficient of PD-L1 and NRF2",
                       palette = c("#4393C3","#D6604D"))+
      stat_compare_means(comparisons = list(c("Control","BRU-treated")), label = "p.signif")
    pdf(paste0("~/mIF/1.results/img_",cor_i,"_with_psignif.pdf"),width = 3, height = 4.2)
    print(p_cor)
    dev.off()
    p_cor <- ggboxplot(results_cor_df, x = "group", y = "tmp_col",
                       color = "group",add = "jitter",add.params = list(size = 0.5),
                       width = 0.5,ylab = "Pearson correlation coefficient of PD-L1 and NRF2",
                       palette = c("#4393C3","#D6604D"))+
      stat_compare_means(comparisons = list(c("Control","BRU-treated")))
    pdf(paste0("~/mIF/1.results/img_",cor_i,".pdf"),width = 3, height = 4.2)
    print(p_cor)
    dev.off()
  }else if(cor_i == "pdl1_cd8a_cor"){
    p_cor <- ggboxplot(results_cor_df, x = "group", y = "tmp_col",
                       color = "group",add = "jitter",add.params = list(size = 0.5),
                       width = 0.5,ylab = "Pearson correlation coefficient of PD-L1 and CD8a",
                       palette = c("#4393C3","#D6604D"))+
      stat_compare_means(comparisons = list(c("Control","BRU-treated")), label = "p.signif")
    pdf(paste0("~/mIF/1.results/img_",cor_i,"_with_psignif.pdf"),width = 3, height = 4.2)
    print(p_cor)
    dev.off()
    p_cor <- ggboxplot(results_cor_df, x = "group", y = "tmp_col",
                       color = "group",add = "jitter",add.params = list(size = 0.5),
                       width = 0.5,ylab = "Pearson correlation coefficient of PD-L1 and CD8a",
                       palette = c("#4393C3","#D6604D"))+
      stat_compare_means(comparisons = list(c("Control","BRU-treated")))
    pdf(paste0("~/mIF/1.results/img_",cor_i,".pdf"),width = 3, height = 4.2)
    print(p_cor)
    dev.off()
  }else if(cor_i == "cd8a_nrf2_cor"){
    p_cor <- ggboxplot(results_cor_df, x = "group", y = "tmp_col",
                       color = "group",add = "jitter",add.params = list(size = 0.5),
                       width = 0.5,ylab = "Pearson correlation coefficient of CD8a and NRF2",
                       palette = c("#4393C3","#D6604D"))+
      stat_compare_means(comparisons = list(c("Control","BRU-treated")), label = "p.signif")
    pdf(paste0("~/mIF/1.results/img_",cor_i,"_with_psignif.pdf"),width = 3, height = 4.2)
    print(p_cor)
    dev.off()
    p_cor <- ggboxplot(results_cor_df, x = "group", y = "tmp_col",
                       color = "group",add = "jitter",add.params = list(size = 0.5),
                       width = 0.5,ylab = "Pearson correlation coefficient of CD8a and NRF2",
                       palette = c("#4393C3","#D6604D"))+
      stat_compare_means(comparisons = list(c("Control","BRU-treated")))
    pdf(paste0("~/mIF/1.results/img_",cor_i,".pdf"),width = 3, height = 4.2)
    print(p_cor)
    dev.off()
  }

}

results_cor_df_turn <- results_cor_df[,2:5] %>%
  pivot_longer(
    cols = c("pdl1_nrf2_cor","pdl1_cd8a_cor","cd8a_nrf2_cor"),
    names_to = "Protein-Protein",
    values_to = "Pearson correlation coefficient" 
  )
results_cor_df_turn$`Protein-Protein` <- ifelse(results_cor_df_turn$`Protein-Protein` == "pdl1_nrf2_cor","PD-L1 - NRF2",results_cor_df_turn$`Protein-Protein`)
results_cor_df_turn$`Protein-Protein` <- ifelse(results_cor_df_turn$`Protein-Protein` == "pdl1_cd8a_cor","PD-L1 - CD8a",results_cor_df_turn$`Protein-Protein`)
results_cor_df_turn$`Protein-Protein` <- ifelse(results_cor_df_turn$`Protein-Protein` == "cd8a_nrf2_cor","CD8a - NRF2",results_cor_df_turn$`Protein-Protein`)

p_cor_all <- ggboxplot(results_cor_df_turn, 
                       x = "group", y = "Pearson correlation coefficient",
                       facet.by = "Protein-Protein",
                   color = "group",add = "jitter",add.params = list(size = 0.5),
                   width = 0.5,ylab = "Pearson correlation coefficient",
                   palette = c("#4393C3","#D6604D"))+
  stat_compare_means(comparisons = list(c("Control","BRU-treated")))

pdf(paste0("~/mIF/1.results/img_",cor_i,"all.pdf"),width = 6, height = 4.8)
print(p_cor_all)
dev.off()

p_cor_all <- ggboxplot(results_cor_df_turn, 
                       x = "group", y = "Pearson correlation coefficient",
                       facet.by = "Protein-Protein",
                       color = "group",add = "jitter",add.params = list(size = 0.5),
                       width = 0.5,ylab = "Pearson correlation coefficient",
                       palette = c("#4393C3","#D6604D"))+
  stat_compare_means(comparisons = list(c("Control","BRU-treated")), label = "p.signif")

pdf(paste0("~/mIF/1.results/img_",cor_i,"all_psignf.pdf"),width = 6, height = 4.8)
print(p_cor_all)
dev.off()

