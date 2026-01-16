rm(list = ls())

setwd("~/bulk/")
source("./functions.r")

library(tidyverse)
library(Seurat)
library(AUCell)
library(org.Mm.eg.db)
library(clusterProfiler)
library(presto)
library(DESeq2)
library(factoextra)
library(ImmuCellAImouse)
library(ggrepel)
library(RColorBrewer)

# read batch1 and batch2 gene count dat & preprocess

dat1 <- data.table::fread("~/bulk/0.source_data/bulk_bru_batch1.txt") %>%
  as.data.frame() %>%
  column_to_rownames(var = "SYMBOL")
colnames(dat1) <- c("BRU_1","BRU_2","CTRL_1","CTRL_2")

dat2 <- data.table::fread("~/bulk/0.source_data/bulk_bru_batch2.txt") %>%
  as.data.frame()

dat2 <- dat2 %>% 
  select("Transcriptid","Geneid","GeneName","BRU_1","BRU_2","Ctrl_1","Ctrl_2") %>%
  group_by(GeneName) %>% 
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames(var = "GeneName") %>%
  mutate(row_sum = rowSums(.)) %>%
  filter(row_sum!=0) %>%
  select(-row_sum)

colnames(dat2) <- c("BRU_3","BRU_4","CTRL_3","CTRL_4")

intersect_genes <- intersect(rownames(dat1),rownames(dat2))

dat <- cbind(dat1[intersect_genes,],dat2[intersect_genes,])

dat <- dat[,c("BRU_1","BRU_2","BRU_3","BRU_4","CTRL_1","CTRL_2","CTRL_3","CTRL_4")]

saveRDS(dat,file = "./1.results_data/0.merged_batch1_2_bru_ctrl_bulk_dat.rds")
write.table(dat,file = "./1.results_data/0.merged_batch1_2_bru_ctrl_bulk_dat.txt",quote = F,col.names = T,row.names = T,sep = "\t")

# remove batch
dat_matrix <- as.matrix(dat)
dat_matrix <- sva::ComBat_seq(dat_matrix, batch=c(1,1,2,2,1,1,2,2), group=c(rep(1,4),rep(2,4)))
dat_rb <- as.data.frame(dat_matrix)
dat_rb <- dat_rb[,c("CTRL_1","CTRL_2","CTRL_3","CTRL_4","BRU_1","BRU_2","BRU_3","BRU_4")]

saveRDS(dat_rb,file = "./1.results_data/0.merged_batch1_2_bru_ctrl_bulk_dat_remove_batcheffect.rds")
write.table(dat_rb,file = "./1.results_data/0.merged_batch1_2_bru_ctrl_bulk_dat_remove_batcheffect.txt",quote = F,col.names = T,row.names = T,sep = "\t")

# DESeq2 deg

sample_df <- as.data.frame(matrix(nrow = ncol(dat_rb),c(1:(2*ncol(dat_rb))))) %>%
      dplyr::rename("Sample" ="V1","Condition" = "V2") %>%
      dplyr::mutate("Condition" = factor(c(rep("CTRL",4),rep("BRU",4)),levels = c("CTRL","BRU"))) %>%
      dplyr::mutate("Sample" = colnames(dat_rb)) %>%
      column_to_rownames(var = "Sample")

diff_dat<- DESeqDataSetFromMatrix(countData = dat_rb, colData = sample_df, design = ~Condition) %>%
  DESeq() %>%
  results() %>%
  as.data.frame() %>%
  filter(!is.na(log2FoldChange)) %>%
  filter(!is.na(padj)) %>%
  filter(!is.na(pvalue))

diff_dat$label <- ifelse(diff_dat$log2FoldChange>=1&diff_dat$padj<=0.05,"Up","Stable")
diff_dat$label <- ifelse(diff_dat$log2FoldChange<=-1&diff_dat$padj<=0.05,"Down",diff_dat$label)

saveRDS(diff_dat,file = "./1.results_data/1.merged_batch1_2_bru_ctrl_bulk_dat_deg.rds")
write.table(diff_dat,file = "./1.results_data/1.merged_batch1_2_bru_ctrl_bulk_dat_deg.txt",quote = F,col.names = T,row.names = T,sep = "\t")

diff_dat <- diff_dat[order(diff_dat$log2FoldChange,decreasing = TRUE),]

p_deg_bru_ctrl_vol <- vol_plot(dat = diff_dat,selected_genes = c("Chil3","Cxcl12","Cxcl2","Nos2","Tnn","Saa3","Gzma","Gzme","Gzmb","Gzmc","Akr1c6","Rfx6","Ttr","Gcg","Colec11","Zfhx2os","Gsta3","Agxt2","Akr1c18"))
pdf("./1.result_figs/deg_bru_ctrl_vol.pdf",width = 6,height = 6)
print(p_deg_bru_ctrl_vol)
dev.off()

# deg GO analysis
# up
up_entrezid_genes=bitr(rownames(diff_dat)[diff_dat$label=="Up"], fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
go_up_gene <- enrichGO(gene = unique(up_entrezid_genes$ENTREZID),
                      OrgDb = org.Mm.eg.db,
                      keyType = 'ENTREZID',
                      ont = "BP",
                      minGSSize = 1,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,
                      readable = TRUE)
go_up_gene <- clusterProfiler::simplify(go_up_gene, cutoff=0.5, by="p.adjust", select_fun=min)
go_up_gene=go_up_gene@result

select_go_id <- c("GO:0010717","GO:0050678","GO:0045765","GO:0062012","GO:0140507","GO:0002709","GO:0002286","GO:0001906","GO:0050863")

go_up_gene_selected <- go_up_gene[select_go_id,]
order1 <- sort(go_up_gene_selected$pvalue,index.return = TRUE,decreasing = T)
go_up_gene_selected$Description <- factor(go_up_gene_selected$Description,levels = go_up_gene_selected$Description[order1$ix])
p <- ggplot(data = go_up_gene_selected,aes(Description,-log10(pvalue),fill = -log10(pvalue))) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low= "#FED976",high="#BD0026") +
  theme_classic() + coord_flip() +
  theme(legend.position="none",axis.text=element_text(size=12)) +
  ggtitle("Biological processes of up-regulation genes in BRU")
p
pdf("./1.result_figs/go_bulk_bru_vs_ctrl_up_genes.pdf",width = 8.5,height = 3.5)
print(p)
dev.off()

write.table(rownames(diff_dat)[diff_dat$label=="Up"],file = "./1.results_data/up_gene_list.txt",quote = F,col.names = F,row.names = F,sep = "\t")

######### GSEA analysis ##########
##################################

library(enrichplot)
library(msigdbr)  # Molecular Signatures Database

# rank genes
logFC <- diff_dat %>%
  mutate(SYMBOL = rownames(.)) %>%
  as.data.frame() %>%
  dplyr::select(SYMBOL,log2FoldChange)

logFC.sort <- arrange(logFC, desc(log2FoldChange))
geneList <- unlist(logFC.sort$log2FoldChange)
names(geneList)<- unlist(logFC.sort$SYMBOL)

# extract mm df
mm_msigdbr <- msigdbr(species="Mus musculus")
mm_df <- as.data.frame(mm_msigdbr)

library(GseaVis)

for(i in select_go_id){
  
  print(i)
  
  select_ID <- i
  
  GO_term_id <- select_ID
  
  GO_term_term <- names(select_go_id)[select_go_id==select_ID]
  
  outab <- as.data.frame(matrix(c(1:3),nrow = 1))
  
  colnames(outab) <- c("Term_id","pvalue","NES")
  
  outab <- outab[-1,]
  
  GO_DATA <- get_GO_data("org.Mm.eg.db", "ALL", "SYMBOL")
  
  go_1_1<-unlist(GO_DATA$PATHID2EXTID[GO_term_id]) %>%
    as.vector()
  
  go_1_1<-cbind(rep(GO_term_term,length(go_1_1)),go_1_1) %>%
    as.data.frame()
  
  colnames(go_1_1)<-c("ont","gene")
  
  if(i == select_go_id[1]){
    go_df_rbind <- go_1_1
  }else{
    go_df_rbind <- rbind(go_df_rbind, go_1_1)
  }
  
  GSEA_result <- GSEA(geneList, TERM2GENE = go_1_1, verbose=FALSE,minGSSize = 1, pvalueCutoff = 1)
  pvalue <- GSEA_result@result[,"p.adjust"]
  NES <- GSEA_result@result[,"NES"]
  
  # fgsea plot
  gsea_plot <- gseaplot2(GSEA_result,1,title = GO_term_term,pvalue_table = TRUE)
  pdf(paste0("~/bulk/1.result_figs/GSEA/GSEA_",GO_term_term,"_fgsea.pdf"),width = 7,height = 6)
  print(gsea_plot)
  dev.off()

  # gsea junjun plot
  gsea_plot <- gseaNb(object = GSEA_result,
                      geneSetID = GO_term_term,
                      subPlot = 3,
                      #addGene = mygene,
                      rmSegment = TRUE,
                      addPval = T,
                      pvalX = 0.6,pvalY = 0.8,
                      pCol = 'black',
                      pHjust = 0)

  pdf(paste0("~/bulk/1.result_figs/GSEA/GSEA_",GO_term_term,".pdf"),width = 7,height = 6)
  print(gsea_plot)
  dev.off()
  
}

colnames(go_df_rbind) <- c("gs_name","gene")

GSEA_all_result <- GSEA(geneList, TERM2GENE = go_df_rbind, verbose=TRUE,minGSSize = 1, pvalueCutoff = 300)

geneSetID = GSEA_all_result@result[["Description"]]

# all plot
gsea_plot <- gseaNb(object = GSEA_all_result,
       geneSetID = geneSetID,
       curveCol = scPalette(17)
       )

pdf(paste0("~/bulk/1.result_figs/GSEA/GSEA_all_change_color",".pdf"),width = 9,height = 6)
print(gsea_plot)
dev.off()
