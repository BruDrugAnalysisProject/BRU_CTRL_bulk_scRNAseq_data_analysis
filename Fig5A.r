######################## MHC-I ###################################
library(GseaVis)
GO_DATA <- get_GO_data("org.Mm.eg.db", "ALL", "SYMBOL")
go_1_1<-unlist(GO_DATA$PATHID2EXTID["GO:0050851"])
go_1_1<-cbind(rep("MHC-I-dependent antigen processing and presentation",length(go_1_1)),as.data.frame(go_1_1))
colnames(go_1_1)<-c("ont","gene")
GSEA_result <- GSEA(geneList, TERM2GENE = go_1_1, verbose=FALSE, pvalueCutoff = 300)
NES <- GSEA_result@result[,"NES"]
gsea_plot <- gseaNb(object = GSEA_result,
       geneSetID = "MHC-I-dependent antigen processing and presentation",
       curveCol = "#5EABD6",
       htHeight = 1,
       htCol = c("#1450A3","#8C1007"),
       base_size = 8,
       rankCol = "#FFC436",
       addPval = T,
       pvalX = 0.9,pvalY = 0.8
       )
gsea_plot
pdf(paste0("~/1.BRU_YanZhu_WorkSpace_20250807/bulk/1.result_figs/GSEA/GSEA_MHC-I-dependent antigen processing and presentation.pdf"),width = 4,height = 4.5)
print(gsea_plot)
dev.off()
