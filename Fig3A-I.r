rm(list = ls())

setwd("~/scRNAseq/1.result_figs/20250814/")

# load functions
source("~/functions.r")

library(tidyverse)
library(Seurat)
library(AUCell)
library(org.Mm.eg.db)
library(clusterProfiler)
library(presto)
library(SCpubr)
library(scRNAtoolVis)
library(CellChat)
library(ggpubr)
library(ggrepel)
library(AUCell)
library(scDblFinder)
library(harmony)

# read data
sc_dat_bru <- Read10X_h5("~/scRNAseq/0.source_data/raw_count/BRU/output_bru_raw_filtered_seurat.h5")
sc_dat_ctrl <- Read10X_h5("~/scRNAseq/0.source_data/raw_count/Ctrl/output_ctrl_raw_filtered_seurat.h5")

# bru
sc_dat_bru_sce <- SingleCellExperiment(list(counts=sc_dat_bru))
sc_dat_bru_sce <- scDblFinder(sc_dat_bru_sce)
db <- data.frame("Barcode" = rownames(colData(sc_dat_bru_sce)), "scDblFinder_DropletType" = sc_dat_bru_sce$scDblFinder.class, "scDblFinder_Score" = sc_dat_bru_sce$scDblFinder.score)
rownames(db) <- db$Barcode
sc_dat_bru <- CreateSeuratObject(sc_dat_bru,min.cells = 0, min.features = 0)
sc_dat_bru@meta.data <- cbind(sc_dat_bru@meta.data,db[rownames(sc_dat_bru@meta.data),])
# ctrl
sc_dat_ctrl_sce <- SingleCellExperiment(list(counts=sc_dat_ctrl))
sc_dat_ctrl_sce <- scDblFinder(sc_dat_ctrl_sce)
db <- data.frame("Barcode" = rownames(colData(sc_dat_ctrl_sce)), "scDblFinder_DropletType" = sc_dat_ctrl_sce$scDblFinder.class, "scDblFinder_Score" = sc_dat_ctrl_sce$scDblFinder.score)
rownames(db) <- db$Barcode
sc_dat_ctrl <- CreateSeuratObject(sc_dat_ctrl,min.cells = 0, min.features = 0)
sc_dat_ctrl@meta.data <- cbind(sc_dat_ctrl@meta.data,db[rownames(sc_dat_ctrl@meta.data),])

sc_dat_bru@meta.data$group <- "BRU"
sc_dat_ctrl@meta.data$group <- "CTRL"
sc_dat <- merge(sc_dat_bru,sc_dat_ctrl)

sc_dat[["percent.mt"]] <- PercentageFeatureSet(sc_dat,pattern = "^mt-")

sc_dat_filtered <- subset(sc_dat, subset = nFeature_RNA >= 700 & nCount_RNA >= 700 &
                            nCount_RNA <= 30000 & percent.mt <= 20 &scDblFinder_DropletType == "singlet") 

sc_dat_filtered@meta.data$cellbarcode <- rownames(sc_dat_filtered@meta.data)

sc_dat_filtered <- SCTransform(sc_dat_filtered,vars.to.regress = "percent.mt",return.only.var.genes = TRUE)

sc_dat_filtered <- sc_dat_filtered %>%
  RunPCA(., verbose = TRUE) %>%
  RunHarmony(.,reduction.use = "pca",group.by.vars = "group") %>%
  RunUMAP(.,reduction = "harmony", dims = 1:30, verbose = TRUE) %>%
  FindNeighbors(.,reduction = "harmony", dims = 1:30, verbose = TRUE) 

sc_dat_filtered <- FindClusters(sc_dat_filtered, resolution = 1.2)
sc_dat_filtered <- FindClusters(sc_dat_filtered, resolution = 0.7)
sc_dat_filtered <- FindClusters(sc_dat_filtered, resolution = 0.1)

sc_dat_filtered_markers <- FindAllMarkers(sc_dat_filtered,only.pos = TRUE)

major_celltype <- c("CD8+ T cell",
                    "M1 macrophage",
                    "Neutrophils", 
                    "M1 macrophage", 
                    "M2 macrophage", 
                    "CD8+ T cell", 
                    "Mitochondrial hyperactive macrophage", 
                    "CD8+ T cell", 
                    "Treg",
                    "NK cell",
                    "B cell", 
                    "epc", 
                    "pDC")
Idents(sc_dat_filtered) <- sc_dat_filtered@meta.data$seurat_clusters
names(major_celltype) <- levels(sc_dat_filtered)
sc_dat_filtered <- RenameIdents(sc_dat_filtered,major_celltype)
sc_dat_filtered@meta.data$major_celltype <- Idents(sc_dat_filtered)

DimPlot(sc_dat_filtered,label = TRUE)

sc_dat_filtered <- subset(sc_dat_filtered,major_celltype%in%c("CD8+ T cell",
                                            "M1 macrophage",
                                            "Neutrophils", 
                                            "M1 macrophage", 
                                            "M2 macrophage", 
                                            "CD8+ T cell", 
                                            "Mitochondrial hyperactive macrophage", 
                                            "CD8+ T cell", 
                                            "Treg",
                                            "NK cell",
                                            "B cell", 
                                            "pDC"))
sc_dat_bak <- sc_dat_filtered
sc_dat <- sc_dat_filtered

# 1. sc_dat dimplot -- vis all cells
# col list
col_list <- scPalette(17)[1:length(unique(sc_dat@meta.data$major_celltype))]
names(col_list) <- unique(sc_dat@meta.data$major_celltype)
  
p_dimplot_group <- clusterCornerAxes(object = sc_dat,
                       reduction = 'umap',
                       noSplit = F,
                       groupFacet = 'orig.ident',
                       clusterCol = "major_celltype",
                       aspect.ratio = 1,
                       pSize = 0.01,
                       relLength = 0.5)+scale_color_manual(values = col_list)
p_dimplot_group
pdf("./dimplot_celltype_split_group.pdf",width = 10,height = 6)
print(p_dimplot_group)
dev.off()

p_dimplot <- clusterCornerAxes(object = sc_dat,
                               reduction = 'umap',
                               #noSplit = F,
                               #groupFacet = 'orig.ident',
                               clusterCol = "major_celltype",
                               aspect.ratio = 1,
                               pSize = 0.01,
                               relLength = 0.5)+scale_color_manual(values = scPalette(17))
p_dimplot
pdf("./dimplot_celltype.pdf",width = 7,height = 6)
print(p_dimplot)
dev.off()

# 2. sc_dat marker genes -- vis all cells dotplot
Idents(sc_dat) <- "major_celltype"

p <- DotPlot(sc_dat,features = c("Cd8a", "Tox", "Grap2", "Cxcr6","Ifng","Gzmb",
                                 "S100a8","S100a9",
                                 "Lcn2","Cxcl9","Ccl24","H2-DMb1","Arg1",
                                 "Ccl2","mt-Co1","mt-Co2","mt-Co3",
                                 "Ncr1","Klrb1c",
                                 "Foxp3","Cd4","Ms4a1",
                                 "Siglech","Tcf4","Ccr9","Bcl11a","Rnase6"),
             #cols = c(brewer.pal(11,name = "YlGnBu")[c(1,8)]),
             scale = TRUE) + 
  theme_test()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black"),
        axis.text.y = element_text(color = "black"))+
  ylab("Cell type")+
  scale_y_discrete(limits = c("CD8+ T cell",
                              "Neutrophils",
                              "M1 macrophage",
                              "M2 macrophage",
                              "Mitochondrial hyperactive macrophage",
                              "NK cell",
                              "Treg",
                              "B cell",
                              "pDC"))
p
pdf("./dotplot_genes_major_celltypes.pdf",width = 10,height = 4)
print(p)
dev.off()

# calculate cell type propration of BRU and CTRL

metadata <- sc_dat@meta.data

celltype_group <- metadata %>%
  group_by(orig.ident)%>%
  dplyr::count(major_celltype)
celltype_group$proportion <- c(celltype_group$n[1:9]/sum(celltype_group$n[1:9]),celltype_group$n[10:18]/sum(celltype_group$n[10:18]))
celltype_group$major_celltype <- factor(celltype_group$major_celltype, levels = names(col_list))

View(celltype_group)

colnames(celltype_group) <- c("sample","major celltype","cell number","proportion")

write.table(celltype_group,file = "summary_celltype_persample.txt",sep = "\t", quote = F,col.names = T,row.names = F)

library(ggpubr)

p <- ggbarplot(celltype_group, "major_celltype", "proportion",
               fill = "orig.ident", 
               xlab = "Celltype",ylab = "Proportion",
               label = FALSE,
               position = position_dodge(0.9))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_manual(values = c("#D6604D","#4393C3"))+
  scale_x_discrete(limits = c("CD8+ T cell",
                              "Neutrophils",
                              "M1 macrophage",
                              "M2 macrophage",
                              "Mitochondrial hyperactive macrophage",
                              "NK cell",
                              "Treg",
                              "B cell",
                              "pDC"))
p
pdf("./barplot_sc_dat_celltype_group.pdf",width = 8,height = 6)
print(p)
dev.off()

p_celltype_percent_barplot <- celltype_percent_barplot(dat = celltype_group,
                                     x_name = "orig.ident",
                                     y_name = "proportion",
                                     fill_name = "major_celltype",
                                     col_list = col_list)
p_celltype_percent_barplot

pdf("./barplot_sc_dat_celltype_ggalluvial_group.pdf",width = 6,height = 5)
print(p_celltype_percent_barplot)
dev.off()

### analysis Cd8+ t cells 

sc_dat_cd8 <- subset(sc_dat,subset = major_celltype=="CD8+ effector T cell")

DimPlot(sc_dat_cd8)

Idents(sc_dat_cd8) <- "orig.ident"

DimPlot(sc_dat_cd8)

View(sc_dat_cd8@meta.data)

diff_markers <- FindAllMarkers(sc_dat_cd8,min.pct = 0,logfc.threshold = 0,return.thresh = 1)

cd8_diff_bru_ctrl <- diff_markers

saveRDS(diff_markers,file = "../../1.results_data/cd8_deg_20250814.rds")

View(cd8_diff_bru_ctrl)

cd8_diff_bru <- cd8_diff_bru_ctrl[cd8_diff_bru_ctrl$cluster=="BRU",]
cd8_diff_bru <- cd8_diff_bru[order(cd8_diff_bru$avg_log2FC,decreasing = TRUE),]

p_rankplot <- function_rankplot(dat = cd8_diff_bru,plot_title="The highly expressed genes ofCD8+ T cells in BRU and control groups",input_type="FindAllMarkers",top_n = 15, down_n = 15,top_log2fc = 0.25,down_log2fc=-0.25)
p_rankplot
pdf("./deg_bru_ctrl_rank_plot.pdf",width = 8,height = 7)
print(p_rankplot)
dev.off()

genes.use <- cd8_diff_bru$gene[cd8_diff_bru$avg_log2FC>=0.25]

library(org.Mm.eg.db)
library(clusterProfiler)

# clusterProfiler GO enrichment

entrezid_genes=bitr(genes.use, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
go_cd8_df <- enrichGO(gene = unique(entrezid_genes$ENTREZID),
                      OrgDb = org.Mm.eg.db,
                      keyType = 'ENTREZID',
                      ont = "BP",
                      minGSSize = 1,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.5,
                      readable = TRUE)
go_cd8_df=go_cd8_df@result
View(go_cd8_df)

saveRDS(go_cd8_df, file = "../../1.results_data/go_cd8_deg_df_20250814.rds")

selected_go_id <- c("GO:0002429","GO:0031341","GO:0050851","GO:0032609","GO:0002699","GO:0002697","GO:0002424","GO:0002418","GO:0002302","GO:0002703","GO:0036037","GO:0140507")
go_cd8_df_selected <- go_cd8_df[selected_go_id,]
names(selected_go_id) <- go_cd8_df_selected$Description

selected_go_id <- selected_go_id[c("regulation of cell killing",
                                   "immune response-activating cell surface receptor signaling pathway",
                                   "immune response to tumor cell",
                                   "positive regulation of immune effector process",
                                   "type II interferon production",
                                   "granzyme-mediated programmed cell death signaling pathway")]

go_cd8_df_selected <- go_cd8_df[selected_go_id,]

order1 <- sort(go_cd8_df_selected$pvalue,index.return = TRUE,decreasing = T)
go_cd8_df_selected$Description <- factor(go_cd8_df_selected$Description,levels = go_cd8_df_selected$Description[order1$ix])
p <- ggplot(data = go_cd8_df_selected,aes(Description,-log10(pvalue),fill = -log10(pvalue))) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low="#FED976",high="#FD8D3C") +
  theme_classic() + coord_flip() +
  theme(legend.position="none",axis.text=element_text(size=12)) +
  ggtitle("Effector CD8+ T cells BRU vs CTRL")
p

pdf("./go_barplot_cd8_bru_vs_ctrl.pdf",width = 8,height = 2.5)
print(p)
dev.off()

# AUC analysis
dir.create("./AUCell_analysis_pathways")

# granzyme_mediated_cell_death_genes
granzyme_mediated_cell_death_genes <- function_get_pathway_gene_list(species = "mouse", GO_term_id = "GO:0140507", GO_term_term = "granzyme-mediated programmed cell death signaling pathway")

expMatrix <- sc_dat_cd8@assays$RNA@counts
cells_rankings <- AUCell_buildRankings(expMatrix, nCores=1, plotStats=TRUE) 

geneSets <- list(c(granzyme_mediated_cell_death_genes$gene))
names(geneSets) <- "geneSet"
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,aucMaxRank = ceiling(0.05*nrow(cells_rankings)))
aucMatrix <- as.data.frame(t(getAUC(cells_AUC)))
sc_dat_cd8$granzyme_mediated_cell_death_genes <- aucMatrix$geneSet

sc_dat_cd8@meta.data$Group <- sc_dat_cd8@meta.data$orig.ident

p_gzm_score <- 
  ggboxplot(sc_dat_cd8@meta.data, x = "Group", y="granzyme_mediated_cell_death_genes",width = 0.65,
          fill = "Group",xlab = "Group", ylab = "Granzyme-mediated programmed cell death\n signaling pathway enrichment score")+
  stat_compare_means()+
  scale_fill_manual(values = c("#D6604D","#4393C3"))
p_gzm_score
pdf(paste0("./AUCell_analysis_pathways/granzyme_mediated_cell_death_genes_enrichment_score.pdf"),width = 2.8,height = 4)
print(p_gzm_score)
dev.off()

sc_dat_cd8@meta.data <- sc_dat_cd8@meta.data[,-c(30:ncol(sc_dat_cd8@meta.data))]

# other pathways

other_pathways <- selected_go_id[-6]

for(i in other_pathways){
  pathway_genes <- function_get_pathway_gene_list(species = "mouse", GO_term_id = i, GO_term_term = names(selected_go_id)[selected_go_id==i])
  
  expMatrix <- sc_dat_cd8@assays$RNA@counts
  cells_rankings <- AUCell_buildRankings(expMatrix, nCores=1, plotStats=TRUE) 
  
  geneSets <- list(c(pathway_genes$gene))
  
  names(geneSets) <- "geneSet"
  
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,aucMaxRank = ceiling(0.05*nrow(cells_rankings)))
  
  aucMatrix <- as.data.frame(t(getAUC(cells_AUC)))
  
  sc_dat_cd8$pathway_genes <- aucMatrix$geneSet
  
  sc_dat_cd8$tmp_pathway_genes <- sc_dat_cd8$pathway_genes
  
  colnames(sc_dat_cd8@meta.data)[ncol(sc_dat_cd8@meta.data)] <- names(selected_go_id)[selected_go_id==i]
  
  sc_dat_cd8@meta.data$Group <- sc_dat_cd8@meta.data$orig.ident
  
  p_gzm_score <- 
    ggboxplot(sc_dat_cd8@meta.data, x = "Group", y="pathway_genes",width = 0.65,
              fill = "Group",xlab = "Group", ylab = paste0(names(selected_go_id)[selected_go_id==i]," enrichment score"))+
    stat_compare_means()+
    scale_fill_manual(values = c("#D6604D","#4393C3"))
  
  pdf(paste0("./AUCell_analysis_pathways/",names(selected_go_id)[selected_go_id==i],"_genes_enrichment_score.pdf"),width = 2.8,height = 4)
  print(p_gzm_score)
  dev.off()
  
  print(i)
  
}

# analysis CD8+ T cells

# subset CD8+ T cells

cd8_t <- subset(sc_dat, subset = rename_celltype %in% "CD8+ T cell")
cd8_t <- cd8_t %>%
  SCTransform(., vars.to.regress = "percent.mt",verbose = FALSE) %>%
  RunPCA(., npcs = 30, verbose = TRUE) %>%
  RunUMAP(.,reduction = "pca", dims = 1:30, verbose = TRUE) %>%
  FindNeighbors(.,reduction = "pca", dims = 1:30, verbose = TRUE) 
cd8_t <- FindClusters(cd8_t, resolution = 1.5)
cd8_t <- FindClusters(cd8_t, resolution = 0.7)
cd8_t <- FindClusters(cd8_t, resolution = 0.3)
DimPlot(cd8_t,label = TRUE)
DimPlot(cd8_t,label = TRUE,group = "orig.ident")
DimPlot(cd8_t,label = TRUE,split.by = "orig.ident")

Idents(cd8_t) <- "SCT_snn_res.0.3"

cd8_markers <- FindAllMarkers(cd8_t)

# annotation

minor_cd8t_optimized <- c(
  "Effector CD8+ T cell (Ccl5-hi)",       
  "Effector CD8+ T cell (Ccl5-hi)", 
  "Proliferating CD8+ T cell",            
  "Activated/Stressed CD8+ T cell",       
  "Tcf7+ Progenitor-like CD8+ T cell",    
  "Effector CD8+ T cell (Ccl5-hi)"
)

Idents(cd8_t) <- cd8_t@meta.data$SCT_snn_res.0.3
names(minor_cd8t_optimized) <- levels(cd8_t)
cd8_t <- RenameIdents(cd8_t,minor_cd8t_optimized)
cd8_t@meta.data$minor_cd8t_celltype_optimized <- Idents(cd8_t)
DimPlot(cd8_t,label = TRUE)

cd8_t@meta.data$minor_cd8t_celltype_optimized <- factor(cd8_t@meta.data$minor_cd8t_celltype_optimized,
                                                        levels = c("Effector CD8+ T cell (Ccl5-hi)",
                                                                   "Proliferating CD8+ T cell",
                                                                   "Activated/Stressed CD8+ T cell",
                                                                   "Tcf7+ Progenitor-like CD8+ T cell"))

celltype_cd8t_cols <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
                 "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
                 "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91",
                 "#EE4C97", "#C6DBEF", "#08519C", "#FED976", "#D94801")[1:length(unique(cd8_t@meta.data$minor_cd8t_celltype_optimized))]

names(celltype_cd8t_cols) <- unique(minor_cd8t_optimized)

p_cd8_dimplot <- clusterCornerAxes(object = cd8_t,
                               reduction = 'umap',
                               cellLabelSize  = 3,
                               noSplit = T,
                               cellLabel = TRUE,
                               #groupFacet = 'orig.ident',
                               clusterCol = "minor_cd8t_celltype_optimized",
                               aspect.ratio = 1,
                               pSize = 0.001,
                               relLength = 0.5)+
  scale_color_manual(values = celltype_cd8t_cols)
p_cd8_dimplot
pdf("CD8_minor_celltypes_dimplot_20251213.pdf",width = 8,height = 6)
print(p_cd8_dimplot)
dev.off()

p <- SCpubr::do_ExpressionHeatmap(sample = cd8_t, 
                                  features = c("Tcf7", "Tnfsf8",
                                               "Pclaf", "Stmn1", "Mki67",
                                               "Hspa1a", "Hspa1b", "Hsph1", "Ifng","Ccl4","Gzmb","Ccl5"),
                                  features.order = c("Tcf7", "Tnfsf8",
                                                     "Pclaf", "Stmn1", "Mki67",
                                                     "Hspa1a", "Hspa1b", "Hsph1", "Ifng","Ccl4","Gzmb","Ccl5"),
                                  groups.order = unique(minor_cd8t_optimized),
                                  use_viridis = TRUE,
                                  font.size = 7,
                                  legend.width = 1,
                                  legend.length = 8,
                                  #max.cutoff = 4.45,
                                  viridis.palette = "magma")

p
pdf("CD8_minor_celltypes_marker_genes_ExpressionHeatmap_20251213.pdf",width = 4,height = 3)
print(p)
dev.off()

# calculate percentage
metadata <- cd8_t@meta.data
celltype_group <- metadata %>%
  group_by(orig.ident)%>%
  dplyr::count(minor_cd8t_celltype_optimized)
celltype_group$proportion <- c(celltype_group$n[1:sum(celltype_group$orig.ident=="BRU")]/sum(sc_dat@meta.data$orig.ident=="BRU"),
                               celltype_group$n[(sum(celltype_group$orig.ident=="BRU")+1):nrow(celltype_group)]/sum(sc_dat@meta.data$orig.ident=="CTRL"))

celltype_group$minor_cd8t_celltype_optimized <- factor(celltype_group$minor_cd8t_celltype_optimized, 
                                                       levels = c("Effector CD8+ T cell (Ccl5-hi)",
                                                                  "Proliferating CD8+ T cell",
                                                                  "Activated/Stressed CD8+ T cell",
                                                                  "Tcf7+ Progenitor-like CD8+ T cell"))

p <- ggbarplot(celltype_group, "minor_cd8t_celltype_optimized", "proportion",
               fill = "orig.ident", #color = "orig.ident", 
               xlab = "Celltype",ylab = "Proportion",
               label = FALSE,
               position = position_dodge(0.9))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_manual(values = c("#D6604D","#4393C3"))+
  scale_x_discrete(limits = c("Effector CD8+ T cell (Ccl5-hi)",
                              "Proliferating CD8+ T cell",
                              "Activated/Stressed CD8+ T cell",
                              "Tcf7+ Progenitor-like CD8+ T cell"))
p
pdf("./barplot_cd8t_sc_dat_celltype_group20251213.pdf",width = 4,height = 5)
print(p)
dev.off()

# calculate pathways

load("~/scRNAseq/1.results_data/immune_genesets.RData")
for(g_s in names(immune_geneset)){
  if(g_s == names(immune_geneset)[1]){
    immune_geneset_df <- data.frame("source" = g_s,"target" = immune_geneset[[g_s]])
  }else{
    immune_geneset_df <- rbind(immune_geneset_df,data.frame("source" = g_s,"target" = immune_geneset[[g_s]]))
  }
  print(g_s)
}

library(babelgene)
mouse_orthologs <- orthologs(genes = immune_geneset_df$target, species = "mouse")
mouse_orthologs <- merge(mouse_orthologs,immune_geneset_df,by.x = "human_symbol", by.y = "target")
for(g_s in unique(mouse_orthologs$source)){
  if(g_s == unique(mouse_orthologs$source)[1]){
    immune_geneset_df_unique <- data.frame("source" = g_s,"target" = unique(mouse_orthologs$symbol[mouse_orthologs$source%in%g_s]))
  }else{
    immune_geneset_df_unique <- rbind(immune_geneset_df_unique,data.frame("source" = g_s,"target" = unique(mouse_orthologs$symbol[mouse_orthologs$source%in%g_s])))
  }
  print(g_s)
}

for(g_s in unique(immune_geneset_df_unique$source)){
  print(g_s)
  if(g_s == unique(immune_geneset_df_unique$source)[1]){
    immune_geneset_list_unique <- list(unique(immune_geneset_df_unique$target[immune_geneset_df_unique$source%in%g_s]))
    names(immune_geneset_list_unique) <- g_s
  }else{
    immune_geneset_list_unique <- c(immune_geneset_list_unique,list(unique(immune_geneset_df_unique$target[immune_geneset_df_unique$source%in%g_s])))
    names(immune_geneset_list_unique)[length(immune_geneset_list_unique)] <- g_s 
  }
}

cd8t_scored <- AddModuleScore(cd8_t,features = immune_geneset_list_unique)
colnames(cd8t_scored@meta.data)[22:ncol(cd8t_scored@meta.data)] <- names(immune_geneset_list_unique)

expMatrix <- cd8_t@assays$RNA@counts
cells_rankings <- AUCell_buildRankings(expMatrix, nCores=1, plotStats=TRUE) 
cells_AUC <- AUCell_calcAUC(immune_geneset_list_unique_selected, cells_rankings,aucMaxRank = ceiling(0.05*nrow(cells_rankings)))
aucMatrix <- as.data.frame(t(getAUC(cells_AUC)))
cd8_t@meta.data <- cd8_t@meta.data[,c(1:21)]
cd8_t@meta.data <- cbind(cd8_t@meta.data,aucMatrix[rownames(cd8_t@meta.data),])

for(i in c("Pan-interferon response",
           "Antigen processing and presentation",
           "Antigen presentation by MHC-I",
           "Cytotoxic activity")){
  input_df <- cd8_t@meta.data[cd8_t@meta.data$minor_cd8t_celltype_optimized=="Effector CD8+ T cell (Ccl5-hi)",]
  wilcox_results <- pairwise.wilcox.test(
    input_df[,colnames(input_df)%in%i],
    input_df$orig.ident,
    p.adjust.method = "BH"  # 可选：BH, bonferroni, holm等
  )
  p_value <- as.numeric(wilcox_results$p.value)
  fc <- mean(input_df[,colnames(input_df)%in%i][input_df$orig.ident=="BRU"])/mean(input_df[,colnames(input_df)%in%i][input_df$orig.ident=="CTRL"])
  if(i == "Pan-interferon response"){
    summary_df <- data.frame("pathway_name" = i,"p_value" = p_value, "fc" = fc)
  }else{
    summary_df_tmp <- data.frame("pathway_name" = i,"p_value" = p_value, "fc" = fc)
    summary_df <- rbind(summary_df,summary_df_tmp)
  }
  print(i)
}
summary_df <- summary_df %>%
  arrange(fc) %>%
  mutate(
    pathway_name = factor(pathway_name, levels = unique(pathway_name))
  ) 
p_wrap <- ggplot(summary_df, aes(x = pathway_name, y = fc, color = fc)) +
  geom_point(aes(size = -log10(p_value))) +
  scale_color_viridis()+
  #scale_color_gradient(low="#C7E9B4",high="#1D91C0") +
  coord_flip() +
  theme_classic() +
  labs(
    y = "Fold change (Mean score)"
  )
p_wrap
pdf("~/scRNAseq/1.result_figs/20251213/fold_change_BRU_Effect_cd8_t_20260103.pdf",width = 5.5,height = 2.5)
print(p_wrap)
dev.off()
