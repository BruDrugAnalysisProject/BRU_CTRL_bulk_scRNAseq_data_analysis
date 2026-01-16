km_survival <- function(genea,geneb,cutoff,dat,auto_cutoff,col1,col2){
  surdat <- dat[,c("OS","OS.time",as.character(genea),as.character(geneb))]
  colnames(surdat) <- c("OS","OS.time","genea","geneb")
  surdat[,"genea"] <- as.numeric(surdat[,"genea"])
  surdat[,"geneb"] <- as.numeric(surdat[,"geneb"])
  surdat[,"OS.time"] <- ((as.numeric(surdat[,"OS.time"]))/365)*12
  if (auto_cutoff==FALSE) {
    surdat <- surdat[order(surdat$geneb,decreasing = T),]
    surdat$geneblabel <- c(rep("geneb_high",floor(nrow(surdat)*0.01*cutoff)),rep("geneb_low",(nrow(surdat) - (floor(nrow(surdat)*0.01*cutoff)))))
    genebh <- surdat %>%
      dplyr::filter(geneblabel == "geneb_high")
    genebl <- surdat %>%
      dplyr::filter(geneblabel == "geneb_low")
    
    ### geneb high expression KM survival analysis of genea
    genebh$genealabel <- ifelse(genebh$genea>median(genebh$genea),"high","low")
    fit=survfit(Surv(OS.time, OS) ~ genealabel , data = genebh)
    diff <- survdiff(Surv(OS.time, OS) ~ genealabel , data = genebh)
    pValue = 1-pchisq(diff$chisq,df = 1)
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("P=",sprintf("%.03f",pValue))
    }
    HR = (diff$obs[1]/diff$exp[1])/(diff$obs[2]/diff$exp[2])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
    ### Plot Survival Curve
    surPlot_genebh_genea = ggsurvplot(fit,
                                      data=genebh,
                                      #pval=pValue,
                                      #pval.size=4.5,
                                      font.pval = c(size = 3, face = "italic", color = "black"),
                                      pval.coord = c(3,0.15),
                                      legend.labs=paste(levels(factor(genebh$genealabel)),"expression"),
                                      #                       legend.title=paste("The survival analysis of",genea,"when the",geneb,"highly expressed"),
                                      palette =  c(col1,col2),
                                      #font.legend=4.5,
                                      legend = c(0.3,0.09),
                                      size = 0.25,
                                      censor.size = 1,
                                      font.x = 7,
                                      font.y = 7,
                                      fontsize = 7,
                                      conf.int = FALSE,
                                      font.legend = list(size = 7, color = "black"),
                                      #legend.title = none,
                                      xlab="Time(Months)",
                                      ylab="Surivial Probability",
                                      font.tickslab = c(7,"plain"),
                                      break.time.by = 12,
                                      risk.table.title="",
                                      risk.table=T,     
                                      risk.table.height=.25,
                                      risk.table.y.text = F, 
                                      ggtheme = theme_test())
    surPlot_genebh_genea$plot <- surPlot_genebh_genea$plot +
      ggplot2::annotate(
        "text",
        x = Inf, y = Inf,
        vjust = 1, hjust = 1,
        label = paste0(pValue,"\n","HR = ",round(HR,digits = 2)," (CI ",round(low95,digits = 2)," - ",round(up95,digits = 2),")"),
        size = 3
      ) +
      theme(legend.title=element_blank())
    #surPlot_genebh_genea <- ggpubr::ggpar(surPlot_genebh_genea,font.legend = list(size = 14, color = "black"))
    
    ### geneb low expression KM survival analysis of genea
    genebl$genealabel <- ifelse(genebl$genea>median(genebl$genea),"high","low")
    fit=survfit(Surv(OS.time, OS) ~ genealabel , data = genebl)
    diff <- survdiff(Surv(OS.time, OS) ~ genealabel , data = genebl)
    pValue = 1-pchisq(diff$chisq,df = 1)
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("p=",sprintf("%.03f",pValue))
    }
    HR = (diff$obs[1]/diff$exp[1])/(diff$obs[2]/diff$exp[2])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
    ### Plot Survival Curve
    surPlot_genebl_genea = ggsurvplot(fit,
                                      data=genebl,
                                      #pval=pValue,
                                      font.pval = c(size = 3, face = "italic", color = "black"),
                                      pval.coord = c(3,0.15),
                                      legend.labs=paste(levels(factor(genebl$genealabel)),"expression"),
                                      #                                    legend.title=paste("The survival analysis of",genea,"when the",geneb,"lowly expressed"),
                                      palette =  c(col1,col2),
                                      size = 0.25,
                                      font.x = 7,
                                      font.y = 7,
                                      fontsize = 7,
                                      conf.int = FALSE,
                                      font.legend = list(size = 7, color = "black"),
                                      legend = c(0.3,0.09),
                                      censor.size = 1,
                                      #legend.title = "",
                                      font.tickslab = 7,
                                      xlab="Time(Months)",
                                      ylab="Surivial Probability",
                                      break.time.by = 12,
                                      risk.table.title="",
                                      risk.table=T,     
                                      risk.table.height=.25,
                                      risk.table.y.text = F, 
                                      ggtheme = theme_test())
    surPlot_genebl_genea$plot <- surPlot_genebl_genea$plot+
      ggplot2::annotate(
        "text",
        x = Inf, y = Inf,
        vjust = 1, hjust = 1,
        label = paste0(pValue,"\n","HR = ",round(HR,digits = 2)," (CI ",round(low95,digits = 2)," - ",round(up95,digits = 2),")"),
        size = 3
      ) +
      theme(legend.title=element_blank())
    #surPlot_genebl_genea <- ggpubr::ggpar(surPlot_genebl_genea,font.legend = list(size = 14, color = "black"))
    survival_plotlist <- list(surPlot_genebh_genea,surPlot_genebl_genea)
    return(survival_plotlist)
  }else{
    dat.cut=surv_cutpoint(surdat, time = "OS.time", event = "OS",variables ="geneb")
    dat.cat=surv_categorize(dat.cut)
    genebh <- surdat[c(rownames(dat.cat[dat.cat$geneb == "high",])),]
    genebl <- surdat[c(rownames(dat.cat[dat.cat$geneb == "low",])),]
    
    ### geneb high expression KM survival analysis of genea
    dat.cut_genebh=surv_cutpoint(genebh, time = "OS.time", event = "OS",variables ="genea")
    dat.cat_genebh=surv_categorize(dat.cut_genebh)
    fit=survfit(Surv(OS.time, OS) ~ genea , data = dat.cat_genebh)
    diff <- survdiff(Surv(OS.time, OS) ~ genea , data = dat.cat_genebh)
    pValue = 1-pchisq(diff$chisq,df = 1)
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("P=",sprintf("%.03f",pValue))
    }
    HR = (diff$obs[1]/diff$exp[1])/(diff$obs[2]/diff$exp[2])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
    ### Plot Survival Curve
    surPlot_genebh_genea = ggsurvplot(fit,
                                      data=dat.cat_genebh,
                                      #pval=pValue,
                                      #pval.size=4.5,
                                      font.pval = c(size = 4.5, face = "italic", color = "black"),
                                      pval.coord = c(3,0.15),
                                      legend.labs=paste(levels(factor(dat.cat_genebh$genea)),"expression"),
                                      #                       legend.title=paste("The survival analysis of",genea,"when the",geneb,"highly expressed"),
                                      palette =  c(col1,col2),
                                      conf.int = FALSE,
                                      font.legend = list(size = 14, color = "black"),
                                      font.x = 16,
                                      font.y = 16,
                                      legend = c(0.3,0.09),
                                      size = 1.5,
                                      #legend.title = "",
                                      font.tickslab = 12,
                                      xlab="Time(Months)",
                                      ylab="Surivial Probability",
                                      break.time.by = 24,
                                      risk.table.title="",
                                      risk.table=T,     
                                      risk.table.height=.25,
                                      risk.table.y.text = F, 
                                      ggtheme = theme_test())
    surPlot_genebh_genea$plot <- surPlot_genebh_genea$plot +
      ggplot2::annotate(
        "text",
        x = Inf, y = Inf,
        vjust = 1, hjust = 1,
        label = paste0(pValue,"\n","HR = ",round(HR,digits = 2)," (CI ",round(low95,digits = 2)," - ",round(up95,digits = 2),")"),
        size = 5
      ) +
      theme(legend.title=element_blank()) 
    #surPlot_genebh_genea <- ggpubr::ggpar(surPlot_genebh_genea,font.legend = list(size = 14, color = "black"))
    
    
    ### geneb low expression KM survival analysis of genea
    dat.cut_genebl=surv_cutpoint(genebl, time = "OS.time", event = "OS",variables ="genea")
    dat.cat_genebl=surv_categorize(dat.cut_genebl)
    fit=survfit(Surv(OS.time, OS) ~ genea , data = dat.cat_genebl)
    diff <- survdiff(Surv(OS.time, OS) ~ genea , data = dat.cat_genebl)
    pValue = 1-pchisq(diff$chisq,df = 1)
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("P=",sprintf("%.03f",pValue))
    }
    HR = (diff$obs[1]/diff$exp[1])/(diff$obs[2]/diff$exp[2])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
    ### Plot Survival Curve
    surPlot_genebl_genea = ggsurvplot(fit,
                                      data=dat.cat_genebl,
                                      conf.int = FALSE,
                                      #pval=pValue,
                                      #pval.size=4.5,
                                      font.pval = c(size = 4.5, face = "italic", color = "black"),
                                      pval.coord = c(3,0.15),
                                      legend.labs=paste(levels(factor(dat.cat_genebl$genea)),"expression"),
                                      #                                    legend.title=paste("The survival analysis of",genea,"when the",geneb,"lowly expressed"),
                                      palette =  c(col1,col2),
                                      font.legend = list(size = 14, color = "black"),
                                      legend = c(0.3,0.09),
                                      size = 1.5,
                                      font.x = 16,
                                      font.y = 16,
                                      #legend.title = "",
                                      font.tickslab = c(7,"plain"),
                                      xlab="Time(Months)",
                                      ylab="Surivial Probability",
                                      break.time.by = 24,
                                      risk.table.title="",
                                      risk.table=T,     
                                      risk.table.height=.25,
                                      risk.table.y.text = F, 
                                      ggtheme = theme_test()) 
    surPlot_genebl_genea$plot <- surPlot_genebl_genea$plot +
      ggplot2::annotate(
        "text",
        x = Inf, y = Inf,
        vjust = 1, hjust = 1,
        label = paste0(pValue,"\n","HR = ",round(HR,digits = 2)," (CI ",round(low95,digits = 2)," - ",round(up95,digits = 2),")"),
        size = 3
      ) +
      theme(legend.title=element_blank()) 
    #surPlot_genebl_genea <- ggpubr::ggpar(surPlot_genebl_genea,font.legend = list(size = 14, color = "black"))
    survival_plotlist <- list(surPlot_genebh_genea,surPlot_genebl_genea)
    return(survival_plotlist)
  }
}

km_survival_single <- function(gene,cutoff,dat,auto_cutoff,col1,col2){
  surdat <- dat[,c("OS.time","OS",gene)]
  colnames(surdat) <- c("OS.time","OS","gene")
  surdat$gene <- as.numeric(surdat$gene)
  surdat$OS.time <- ((as.numeric(surdat$OS.time))/365)*12
  surdat <- surdat[order(surdat$gene,decreasing = T),]
  if(auto_cutoff==FALSE){
    surdat$genelabel <- c(rep("high",floor(nrow(surdat)*0.01*cutoff)),rep("low",(nrow(surdat) - (floor(nrow(surdat)*0.01*cutoff)))))
    fit=survfit(Surv(OS.time, OS) ~ genelabel , data = surdat)
    diff <- survdiff(Surv(OS.time, OS) ~ genelabel , data = surdat)
    pValue = 1-pchisq(diff$chisq,df = 1)
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("P=",sprintf("%.03f",pValue))
    }
    HR = (diff$obs[1]/diff$exp[1])/(diff$obs[2]/diff$exp[2])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
    ### Plot Survival Curve
    surPlot = ggsurvplot(fit,
                         data=surdat,
                         #pval=pValue,
                         #pval.size=4.5,
                         font.pval = c(size = 3, face = "italic", color = "black"),
                         pval.coord = c(3,0.15),
                         legend.labs=paste(levels(factor(surdat$genelabel)),"expression"),
                         #                       legend.title=paste("The survival analysis of",genea,"when the",geneb,"highly expressed"),
                         palette =  c(col1,col2),
                         censor.size = 1,
                         font.legend = list(size = 7, color = "black"),
                         legend = c(0.3,0.09),
                         conf.int = FALSE,
                         size = 0.25,
                         font.x = 7,
                        font.y = 7,
                        fontsize = 1,
                         #legend.title = "",
                         font.tickslab = c(7,"plain"),
                         xlab="Time(Months)",
                         ylab="Surivial Probability",
                         break.time.by = 36,
                         risk.table.title="",
                         risk.table=T,     
                         risk.table.height=.25,
                         risk.table.y.text = F, 
                         ggtheme = theme_test())
    surPlot$plot <- surPlot$plot +
      ggplot2::annotate(
        "text",
        x = Inf, y = Inf,
        vjust = 1, hjust = 1,
        label = paste0(pValue,"\n","HR = ",round(HR,digits = 2)," (CI ",round(low95,digits = 2)," - ",round(up95,digits = 2),")"),
        size = 3
      ) +
      theme(legend.title=element_blank()) 
    #surPlot <- ggpubr::ggpar(surPlot,font.legend = list(size = 14, color = "black"))
  }else{
    dat.cut=surv_cutpoint(surdat, time = "OS.time", event = "OS",variables ="gene")
    dat.cat=surv_categorize(dat.cut)
    fit=survfit(Surv(OS.time, OS) ~ gene , data = dat.cat)
    diff <- survdiff(Surv(OS.time, OS) ~ gene , data = dat.cat)
    pValue = 1-pchisq(diff$chisq,df = 1)
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("P=",sprintf("%.03f",pValue))
    }
    HR = (diff$obs[1]/diff$exp[1])/(diff$obs[2]/diff$exp[2])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
    ### Plot Survival Curve
    surPlot = ggsurvplot(fit,
                         data=dat.cat,
                         #pval=pValue,
                         #pval.size=4.5,
                         font.pval = c(size = 3, face = "italic", color = "black"),
                         pval.coord = c(3,0.15),
                         legend.labs=paste(levels(factor(dat.cat$gene)),"expression"),
                         #                       legend.title=paste("The survival analysis of",genea,"when the",geneb,"highly expressed"),
                         palette =  c(col1,col2),
                         font.legend = list(size = 7, color = "black"),
                         legend = c(0.3,0.09),
                         conf.int = FALSE,
                         size = 0.25,
                         font.x = 7,
                        font.y = 7,
                        fontsize = 1,
                        censor.size = 1,
                         #legend.title = "",
                         font.tickslab = c(7,"plain"),
                         xlab="Time(Months)",
                         ylab="Surivial Probability",
                         break.time.by = 36,
                         risk.table.title="",
                         risk.table=T,     
                         risk.table.height=.25,
                         risk.table.y.text = F, 
                         ggtheme = theme_test())
    surPlot$plot <- surPlot$plot +
      ggplot2::annotate(
        "text",
        x = Inf, y = Inf,
        vjust = 1, hjust = 1,
        label = paste0(pValue,"\n","HR = ",round(HR,digits = 2)," (CI ",round(low95,digits = 2)," - ",round(up95,digits = 2),")"),
        size = 3
      ) +
      theme(legend.title=element_blank()) 
    #surPlot <- ggpubr::ggpar(surPlot,font.legend = list(size = 14, color = "black"))
  }
  return(surPlot)
}

correct_gene_expression <- function(gene_name, expr_matrix, purity_df, normal_expr = NULL) {
  gene_expr <- expr_matrix[gene_name, ]
  common_samples <- intersect(colnames(gene_expr), rownames(purity_df))
  gene_expr <- gene_expr[common_samples]
  purity_vec <- purity_df$TumorPurity[match(common_samples, rownames(purity_df))]
  if (!is.null(normal_expr)) {
    c_val <- median(as.numeric(normal_expr[gene_name, ]), na.rm = TRUE)
  } else {
    c_val <- min(gene_expr, na.rm = TRUE)
  }
  corrected_expr <- (gene_expr - c_val * (1 - purity_vec)) / purity_vec
  corrected_expr[corrected_expr < 0] <- 0
  return(data.frame(
    sample_id = common_samples,
    origial_expr = as.numeric(gene_expr),
    corrected_expr = as.numeric(corrected_expr),
    purity = purity_vec
  ))
}

vol_plot <- function(dat,selected_genes){
  diff_dat <- dat
  if (length(selected_genes) != 0 & !is.null(selected_genes)) {
    diff_dat$selected_genes <- ""
    gene_list <- selected_genes
    for (i in 1:length(gene_list)) {
      diff_dat[gene_list[i],"selected_genes"] <- gene_list[i]
    }
    
    up <- diff_dat %>%
      dplyr::filter(selected_genes!="") %>%
      dplyr::filter(label == "Up")
    
    down <- diff_dat %>%
      dplyr::filter(selected_genes!="") %>%
      dplyr::filter(label == "Down")
    
    abs_xlim <- ifelse(max(abs(diff_dat$log2FoldChange)) == max(diff_dat$log2FoldChange),max(diff_dat$log2FoldChange)+floor(max(diff_dat$log2FoldChange)/2)+1,max(abs(diff_dat$log2FoldChange))+floor(max(diff_dat$log2FoldChange)/2)+1)
    
    p <- ggplot(diff_dat,aes(x = log2FoldChange, y = -log10(padj), color = label)) +
      geom_point(size = 1) + 
      theme_test() +
      xlim(-abs_xlim,abs_xlim) + 
      scale_color_manual(values = c("#4d4e75","grey60","#DC2819")) + 
      geom_hline(yintercept =-log10(0.05),linetype = 3,size = 0.7,color = "black",lty = "dashed") +
      geom_vline(xintercept =c(-1,1),linetype = 3,size = 0.7,color = "black",lty = "dashed") +
      theme(legend.position = "none") +
      labs(x = "Log2(FC)", y = "-Log10(adjusted P-value)") +
      ggtitle(label = "The Volcano Plot Of DEGs") +
      theme_test() +
      theme(legend.background=element_blank(), legend.key=element_blank(),
            plot.title = element_text(hjust = 0.5), # Center-justify plot title
            axis.text=element_text(size=7),
            axis.title = element_text(size = 7),
            # legend.title = element_blank(),
            legend.text = element_text(size = 7),
            plot.margin = margin(15,5.5,5.5,5.5),
            panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
      # ggrepel::geom_text_repel(aes(label = diff_dat$selected_genes),
      #                                 size = 3.5,
      #                                 box.padding = 0.5,
      #                                 min.segment.length = 0.5,
      #                                 max.overlaps = 100000,
      #                                 segment.color = "black",show.legend = T) + coord_fixed(ratio = 1.1)
    p <- p + 
    # geom_point(data = up,
    #                     aes(x = log2FoldChange, y = -log10(padj)),
    #                     color = '#EB4232', size = 2.5, alpha = 0.2) +
             geom_text_repel(data = up,
                              aes(x = log2FoldChange, y = -log10(padj), label = selected_genes),
                              seed = 233,
                              size = 3.5,
                              color = 'black',
                              min.segment.length = 0,
                              force = 3,
                              force_pull = 2,
                              box.padding = 0.1,
                              max.overlaps = Inf,
                              segment.linetype = 3, #线段类型,1为实线,2-6为不同类型虚线
                              segment.color = 'black', #线段颜色
                              segment.alpha = 0.5, #线段不透明度
                              nudge_x = 3 + up$log2FoldChange, #标签x轴起始位置调整
                              direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
                              hjust = 0) + #对齐标签：0右对齐，1左对齐，0.5居中
            #  geom_point(data = down,
            #             aes(x = log2FoldChange, y = -log10(padj)),
            #             color = '#2DB2EB', size = 2.5, alpha = 0.2) +
             geom_text_repel(data = down,
                              aes(x = log2FoldChange, y = -log10(padj), label = selected_genes),
                              seed = 233,
                              size = 3.5,
                              color = 'black',
                              min.segment.length = 0,
                              force = 3,
                              force_pull = 2,
                              box.padding = 0.1,
                              max.overlaps = Inf,
                              segment.linetype = 3, #线段类型,1为实线,2-6为不同类型虚线
                              segment.color = 'black', #线段颜色
                              segment.alpha = 0.5, #线段不透明度
                              nudge_x = -3 + down$log2FoldChange, #标签x轴起始位置调整
                              direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
                              hjust = 0)
  }else if(length(selected_genes) == 0 | is.null(selected_genes)){
    abs_xlim <- ifelse(max(abs(diff_dat$log2FoldChange)) == max(diff_dat$log2FoldChange),max(diff_dat$log2FoldChange)+1,max(abs(diff_dat$log2FoldChange))+1)
    p <- ggplot(diff_dat,aes(x = log2FoldChange, y = -log10(padj), color = label)) +
      geom_point(size = 1) + 
      theme_test() +
      xlim(-abs_xlim,abs_xlim) + 
      scale_color_manual(values = c("#4d4e75","grey60","#DC2819")) + 
      geom_hline(yintercept =-log10(0.05),linetype = 3,size = 0.7,color = "black",lty = "dashed") +
      geom_vline(xintercept =c(-1,1),linetype = 3,size = 0.7,color = "black",lty = "dashed") +
      theme(legend.position = "none") +
      labs(x = "Log2(FC)", y = "-Log10(adjusted P-value)") +
      ggtitle(label = "The Volcano Plot Of DEGs") +
      theme_test() +
      theme(legend.background=element_blank(), legend.key=element_blank(),
            plot.title = element_text(hjust = 0.5), # Center-justify plot title
            axis.text=element_text(size=7),
            axis.title = element_text(size = 7),
            # legend.title = element_blank(),
            legend.text = element_text(size = 7),
            plot.margin = margin(15,5.5,5.5,5.5),
            panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  }
  return(p)
}

library(DOSE)
library(GOSemSim)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(dplyr)
library(GO.db)
#
get_GO_data <- function(OrgDb, ont, keytype) {
  GO_Env <- get_GO_Env()
  use_cached <- FALSE
  
  if (exists("organism", envir=GO_Env, inherits=FALSE) &&
      exists("keytype", envir=GO_Env, inherits=FALSE)) {
    
    org <- get("organism", envir=GO_Env)
    kt <- get("keytype", envir=GO_Env)
    
    if (org == DOSE:::get_organism(OrgDb) &&
        keytype == kt &&
        exists("goAnno", envir=GO_Env, inherits=FALSE)) {
      ## https://github.com/GuangchuangYu/clusterProfiler/issues/182
      ## && exists("GO2TERM", envir=GO_Env, inherits=FALSE)){
      
      use_cached <- TRUE
    }
  }
  
  if (use_cached) {
    goAnno <- get("goAnno", envir=GO_Env)
  } else {
    OrgDb <- GOSemSim:::load_OrgDb(OrgDb)
    kt <- keytypes(OrgDb)
    if (! keytype %in% kt) {
      stop("keytype is not supported...")
    }
    
    kk <- keys(OrgDb, keytype=keytype)
    goAnno <- suppressMessages(
      AnnotationDbi::select(OrgDb, keys=kk, keytype=keytype,
                            columns=c("GOALL", "ONTOLOGYALL")))
    
    goAnno <- unique(goAnno[!is.na(goAnno$GOALL), ])
    
    assign("goAnno", goAnno, envir=GO_Env)
    assign("keytype", keytype, envir=GO_Env)
    assign("organism", DOSE:::get_organism(OrgDb), envir=GO_Env)
  }
  
  if (ont == "ALL") {
    GO2GENE <- unique(goAnno[, c(2,1)])
  } else {
    GO2GENE <- unique(goAnno[goAnno$ONTOLOGYALL == ont, c(2,1)])
  }
  
  GO_DATA <- DOSE:::build_Anno(GO2GENE, get_GO2TERM_table())
  
  goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
  goOnt <- goOnt.df[,2]
  names(goOnt) <- goOnt.df[,1]
  assign("GO2ONT", goOnt, envir=GO_DATA)
  return(GO_DATA)
}

get_GO_Env <- function () {
  if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".GO_clusterProfiler_Env", new.env(), envir=envir)
  }
  get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
}

get_GO2TERM_table <- function() {
  GOTERM.df <- get_GOTERM()
  GOTERM.df[, c("go_id", "Term")] %>% unique
}

get_GOTERM <- function() {
  pos <- 1
  envir <- as.environment(pos)
  if (!exists(".GOTERM_Env", envir=envir)) {
    assign(".GOTERM_Env", new.env(), envir)
  }
  GOTERM_Env <- get(".GOTERM_Env", envir = envir)
  if (exists("GOTERM.df", envir = GOTERM_Env)) {
    GOTERM.df <- get("GOTERM.df", envir=GOTERM_Env)
  } else {
    GOTERM.df <- toTable(GOTERM)
    assign("GOTERM.df", GOTERM.df, envir = GOTERM_Env)
  }
  return(GOTERM.df)
}




load_packages <- function() {
  packages <- c("tidyverse", "EBImage", "tiff", "corrplot", 
                "patchwork", "spdep", "spatialEco", "FNN", "magick")
  
  cran_packages <- packages[!packages %in% c("EBImage")]
  for(pkg in cran_packages) {
    if(!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
  
  if(!require("EBImage", character.only = TRUE)) {
    if(!require("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("EBImage")
    library(EBImage)
  }
}

read_protein_channels <- function(file_paths, protein_names, 
                                  remove_scale = TRUE, 
                                  scale_height = NULL,
                                  scale_position = "bottom",
                                  crop_region = NULL) {
  if(length(file_paths) != 3 || length(protein_names) != 3) {
    stop("three files")
  }
  
  channels <- map2(file_paths, protein_names, ~{
    cat(sprintf("read %s channel: %s\n", .y, basename(.x)))
    
    img <- if(grepl("\\.tif(f)?$", .x, ignore.case = TRUE)) {
      readTIFF(.x, all = FALSE)
    } else {
      readImage(.x)
    }
    
    if(length(dim(img)) == 3) {
      dims <- dim(img)
      if(dims[3] >= 3) {
        img <- 0.299 * img[,,1] + 0.587 * img[,,2] + 0.114 * img[,,3]
      }
    }
    
    if(remove_scale) {
      img <- remove_scale_from_image(img, scale_height, scale_position, crop_region)
    }
    
    img <- img / max(img, na.rm = TRUE)
    
    return(img)
  }) %>% set_names(protein_names)
  
  dims <- map(channels, dim)
  min_rows <- min(map_dbl(dims, 1))
  min_cols <- min(map_dbl(dims, 2))
  
  if(length(unique(dims)) > 2) {
    warning("size different")
    channels <- map(channels, ~.[1:min_rows, 1:min_cols])
  }
  
  cat("working...\n")
  return(channels)
}

remove_scale_from_image <- function(img, scale_height = NULL, 
                                    scale_position = "bottom",
                                    crop_region = NULL) {
  dims <- dim(img)
  rows <- dims[1]
  cols <- dims[2]
  
  # 方法1：如果指定了裁剪区域
  if(!is.null(crop_region) && length(crop_region) == 4) {
    x1 <- crop_region[1]
    y1 <- crop_region[2]
    x2 <- crop_region[3]
    y2 <- crop_region[4]
    
    if(x1 >= 1 && x2 <= cols && y1 >= 1 && y2 <= rows) {
      return(img[y1:y2, x1:x2])
    } else {
      warning("skip")
    }
  }
  
  if(is.null(scale_height)) {
    scale_height <- auto_detect_scale_height(img, scale_position)
  }
  
  if(scale_height > 0) {
    cat(sprintf("  check: %d pixel\n", scale_height))
    
    if(scale_position == "bottom") {
      if(scale_height < rows) {
        img <- img[1:(rows - scale_height), ]
      }
    } else if(scale_position == "top") {
      if(scale_height < rows) {
        img <- img[(scale_height + 1):rows, ]
      }
    } else if(scale_position == "left") {
      if(scale_height < cols) {
        img <- img[, (scale_height + 1):cols]
      }
    } else if(scale_position == "right") {
      if(scale_height < cols) {
        img <- img[, 1:(cols - scale_height)]
      }
    }
  }
  
  return(img)
}

auto_detect_scale_height <- function(img, scale_position = "bottom", 
                                     threshold_ratio = 0.1) {
  dims <- dim(img)
  rows <- dims[1]
  cols <- dims[2]
  
  if(scale_position == "bottom") {

    bottom_rows <- max(1, floor(rows * 0.9)):rows
    bottom_section <- img[bottom_rows, ]
    

    row_means <- apply(bottom_section, 1, mean, na.rm = TRUE)
    

    if(length(row_means) > 10) {

      gradients <- diff(row_means)
      significant_gradients <- which(abs(gradients) > sd(gradients, na.rm = TRUE) * 2)
      
      if(length(significant_gradients) > 0) {

        scale_start <- max(significant_gradients)
        scale_height <- length(row_means) - scale_start
        return(max(5, min(scale_height, rows * 0.2))) 
      }
    }
  }
  
  if(scale_position == "top") {

    bottom_rows <- max(1, floor(rows * 0.1)):rows
    bottom_section <- img[bottom_rows, ]
    

    row_means <- apply(bottom_section, 1, mean, na.rm = TRUE)
    

    if(length(row_means) > 10) {

      gradients <- diff(row_means)
      significant_gradients <- which(abs(gradients) > sd(gradients, na.rm = TRUE) * 2)
      
      if(length(significant_gradients) > 0) {

        scale_start <- max(significant_gradients)
        scale_height <- scale_start
        return(max(5, min(scale_height, rows * 0.2))) 
      }
    }
  }
  

  cat("  error\n")
  return(0)
}

create_protein_expression_df <- function(channels, sample_fraction = 1, seed = 123) {
  dims <- dim(channels[[1]])
  total_pixels <- prod(dims)
  if(sample_fraction < 1) {
    set.seed(seed)
    n_samples <- round(total_pixels * sample_fraction)
    sampled_indices <- sample(1:total_pixels, n_samples)
  } else {
    sampled_indices <- 1:total_pixels
  }
  y_coords <- rep(1:dims[1], times = dims[2])
  x_coords <- rep(1:dims[2], each = dims[1])
  df <- data.frame(
    x = x_coords[sampled_indices],
    y = y_coords[sampled_indices]
  )
  
  for(protein in names(channels)) {
    df[[protein]] <- as.vector(channels[[protein]])[sampled_indices]
  }
  
  return(df)
}

preprocess_expression_data <- function(df, remove_background = TRUE, normalize = TRUE) {
  protein_names <- setdiff(colnames(df), c("x", "y"))
  processed_df <- df
  
  if(remove_background) {
    
    otsu_threshold <- function(x) {
      x <- x[is.finite(x)]
      if(length(x) < 2) return(0)
      
      h <- hist(x, breaks = 100, plot = FALSE)
      counts <- h$counts
      breaks <- h$breaks
      
      total <- sum(counts)
      if(total == 0) return(0)
      
      sumB <- wB <- maximum <- 0
      threshold <- 0
      
      for(i in seq_along(counts)) {
        wB <- wB + counts[i]
        if(wB == 0) next
        
        wF <- total - wB
        if(wF == 0) break
        
        sumB <- sumB + i * counts[i]
        mB <- sumB / wB
        mF <- (sum(seq_along(counts) * counts) - sumB) / wF
        
        between <- wB * wF * (mB - mF)^2
        
        if(between > maximum) {
          threshold <- breaks[i]
          maximum <- between
        }
      }
      
      return(threshold)
    }
    
    for(protein in protein_names) {
      threshold <- otsu_threshold(df[[protein]])
      processed_df[[protein]] <- ifelse(df[[protein]] < threshold, 0, df[[protein]])
    }
  }
  
  # 标准化到0-1范围
  if(normalize) {
    processed_df[protein_names] <- map_df(processed_df[protein_names], ~{
      min_val <- min(.x, na.rm = TRUE)
      max_val <- max(.x, na.rm = TRUE)
      if(max_val - min_val == 0) return(rep(0, length(.x)))  
      (.x - min_val) / (max_val - min_val)
    })
  }
  
  stats <- processed_df[protein_names] %>%
    map_df(~data.frame(
      Mean = mean(.x, na.rm = TRUE),
      SD = sd(.x, na.rm = TRUE),
      Median = median(.x, na.rm = TRUE),
      Min = min(.x, na.rm = TRUE),
      Max = max(.x, na.rm = TRUE),
      Zero_Percent = 100 * sum(.x == 0, na.rm = TRUE) / length(.x)
    ), .id = "Protein")
  
  print(stats)
  
  return(list(
    processed_data = processed_df,
    protein_names = protein_names,
    stats = stats
  ))
}
