if (T) {
  dir.create('scripts')
  dir.create('files')
  dir.create('PDFs')
  dir.create('results')
  dir.create('raw_dat/tcga', recursive = T)
}

plot_RS_Event_heatmap <- function(dat, genes, position = c(0, 1) ,col=c("#00468BFF", "#ED0000FF")) {
  # dat：数据框，行为样本，列为特征（生存时间、生存状态、风险得分）
  library(ggsci)
  library(pheatmap)
  library(dplyr)
  library(ggplotify)
  library(ggplot2)
  dat <- as.data.frame(dat)
  df <- dat[, c("OS.time", "OS", 'RiskScore')]
  df$Type <- ifelse(df$RiskScore > median(df$RiskScore), 'High Risk', 'Low Risk')
  df$Event <- as.factor(ifelse(df$OS == 0, 'Alive', 'Dead'))
  # print(head(df))
  df <- arrange(df, RiskScore)
  df$Order <- seq(1, length(df$RiskScore))
  
  # 生存时间转为年
  if (max(df$OS.time) > 365) {
    df$OS.time <- df$OS.time / 365
  } else if (max(df$OS.time) > 24) {
    df$OS.time <- df$OS.time / 12
  } else {
    df$OS.time <- df$OS.time
  }
  
  # 绘制 RiskScore 的分布
  p1 <- ggplot(df,aes(x = Order,y = RiskScore,color = Type)) +
    geom_point(lwd=1.1) + 
    geom_hline(aes(yintercept=median(RiskScore),), 
               colour="black", linetype="dashed") +
    geom_vline(aes(xintercept=median(Order),), 
               colour="black", linetype="dashed") +
    theme_few() + 
    ylab("RiskScore") + xlab('Samples') +
    scale_color_lancet() +
    theme(legend.position=c(0, 1),
          legend.justification=c(0, 1),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.title = element_blank())
  # 绘制生存状态的分布
  p2 <- ggplot(df, aes(x=Order, y=OS.time, color=Event)) + 
    geom_point(lwd=1.1) + theme_few() + 
    geom_vline(aes(xintercept=median(Order),), 
               colour="black", linetype="dashed") +
    ylab("Time(years)") + xlab('Samples') +
    scale_color_manual(values = col) +
    # scale_color_npg() +
    theme(legend.position=position,
          legend.justification=position,
          legend.background = element_rect(fill = NA, colour = NA),
          legend.title = element_blank())
  
  # 绘制模型基因的热图 #####
  df1 <- dat[, c(genes, 'RiskScore')]
  # print(head(df1))
  df1 <- arrange(df1, RiskScore)
  annotation_col <- data.frame(RiskScore = df1$RiskScore)
  rownames(annotation_col) <- rownames(df1)
  ann_colors = list(RiskScore = c('#EFC000', "white", "#0073C2"))
  bk=unique(c(seq(-1.2, 1.2, length=100)))
  p3 <- pheatmap(t(df1[, genes]), 
                 scale = 'row', 
                 breaks = bk,
                 annotation_col = annotation_col,
                 annotation_legend = FALSE,
                 cluster_cols = F, cluster_rows = F,
                 show_rownames = T, show_colnames = F,
                 # color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                 annotation_colors =ann_colors,
                 border_color = NA,
                 silent = T
  )
  p3 <- as.ggplot(p3)
  
  # 合并文件
  p <- cowplot::plot_grid(p1, p2, p3,
                          ncol=1, labels=LETTERS[1:3], 
                          align = 'hv')
  return(p)
}
ggsurvplotKM <- function(dat, title = 'Groups', 
                         lables = c(), col = c('lancet'),
                         risk.table = TRUE,
                         tables.height = 0.25) {
  # dat：数据框，行为样本，列为时间、状态以及分组
  # the color palette to be used. Allowed values include "hue" for the default hue color scale; "grey" for grey color palettes; brewer palettes e.g. "RdBu", "Blues", ...; or custom color palette e.g. c("blue", "red"); and scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty". See details section for more information. Can be also a numeric vector of length(groups); in this case a basic color palette is created using the function palette.
  library(ggplot2)
  library(survminer)
  library(survival)
  library(ggthemes)
  library(ggplotify)
  library(ggsci)
  
  colnames(dat) <- c("OS.time", "OS", "Groups")
  # 将时间转换成年份
  if (max(dat[, 1]) > 365) {
    dat[, 1] <- dat[, 1] / 365
  } else if (max(dat[, 1]) > 24) {
    dat[, 1] <- dat[, 1] / 12
  } else {
    dat <- dat
  }
  fit <- survfit(Surv(OS.time, OS) ~ Groups,data=dat)
  surv.fit <- ggsurvplot(fit, data = dat, palette = col,
                         pval = TRUE, 
                         pval.method = T,
                         pval.method.size = 4,
                         pval.method.coord = c(0, 0.15),
                         surv.median.line='hv',
                         linetype = 1, 
                         pval.coord=c(0, 0.05), 
                         pval.size = 4,
                         risk.table = risk.table,
                         legend.title = title,
                         legend.labs = lables,
                         xlab = 'Time(years)',
                         ggtheme=theme_bw(),
                         risk.table.y.text = FALSE,
                         tables.height = tables.height)
  # 将图形转换为 ggplot 对象
  if (risk.table) {
    surv.fit1 <- surv.fit$plot + 
      theme(legend.position=c(1,1), 
            legend.justification=c(1,1),
            plot.margin=unit(c(0.1, 0.15, 0, 0.15), "inches"),
            legend.background = element_rect(fill = NA, colour = NA)
            # ,
            # axis.text.x=element_blank(),
            # axis.title.x=element_blank()
      )
    
    surv.fit2 <- surv.fit$table + 
      theme(plot.title=element_blank(),
            plot.margin=unit(c(0, 0.15, 0, 0.15), "inches")) +
      ylab('')
    surv.fit <- ggpubr::ggarrange(surv.fit1,
                                  surv.fit2, 
                                  ncol = 1, 
                                  nrow = 2,
                                  heights = c(1 - tables.height, 
                                              tables.height),
                                  align = "hv")
    # surv.fit <- cowplot::plot_grid(surv.fit1, 
    #                                surv.fit2,
    #                                align = 'hv',
    #                                ncol = 1,
    #                                rel_heights = c(1 - tables.height, tables.height))
  } else {
    surv.fit <- surv.fit$plot + 
      theme(legend.position=c(1,1), 
            legend.justification=c(1,1),
            # plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"),
            legend.background = element_rect(fill = NA, colour = NA))
  }
  return(surv.fit)
}
ggplotROC <- function(dat, breaks=c(1,3,5), color = color8) {
  # time：生存时间
  # status：生存状态
  # RS：风险得分或其他特征
  library(timeROC)
  library(ggplot2)
  library(survival)
  # 将时间转换成年份
  if (max(dat$OS.time) > 365) {
    dat$OS.time <- dat$OS.time / 365
  } else if (max(dat$OS.time) > 24) {
    dat$OS.time <- dat$OS.time / 12
  } else {
    dat$OS.time <- dat$OS.time
  }
  
  test <- timeROC(T=dat$OS.time,
                  delta=dat$OS
                  ,marker=dat$RiskScore,
                  cause=1,weighting="marginal",
                  times=breaks,
                  iid=TRUE)
  test_TP <- test$TP
  test_FP <- test$FP
  test_AUC <- test$AUC
  test_CI <- confint(test)$CI_AUC
  ROC_dat <- rbind()
  for (i in 1:length(test_AUC)) {
    # test_dat1 <- data.frame(TP = test_TP[, i],
    #                         FP = test_FP[, i],
    #                         Times =paste0(breaks[i], 
    #                                       ' Years AUC=', round(test_AUC[i], 2), ',95%CI(',
    #                                       round(test_CI[i,1] / 100, 2), '-',
    #                                       round(test_CI[i,2] / 100, 2), ')')
    # )
    test_dat1 <- data.frame(TP = test_TP[, i],
                            FP = test_FP[, i],
                            Times =paste0(breaks[i], 
                                          ' Years AUC: ', round(test_AUC[i], 2))
    )
    ROC_dat <- rbind(ROC_dat, test_dat1)
  }
  ggplot(ROC_dat, aes(x=FP,y=TP, fill=Times))+
    geom_line(aes(colour=Times),lwd=0.75)+
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    xlab('False positive rate')+ylab('True positive rate')+ 
    scale_colour_manual(values = color) +
    theme_bw() +
    theme(legend.position=c(1,0),legend.justification=c(1,0),
          # plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches"),
          legend.background = element_rect(fill = NA, colour = NA)) 
}
ggplot_bar <- function(data, width = 0.5, cols = color9,
                       xlab = '', ylab = '', name = 'Groups') {
  # data：数据框，行和列为不同的分类特征，数据框的数据是特征的统计量
  for (i in 1:ncol(data)) {
    data[, i] = data[, i] / sum(data[, i])
  }
  data=reshape2::melt(data)
  colnames(data)=c('Type','variable','Perc')
  data[, 1] = as.character(data[, 1])
  data[, 2] = as.character(data[, 2])
  bar = ggplot(data, aes(x=variable, y=Perc, fill=Type))+geom_bar(stat = "identity", width=width) + theme_get()
  bar = bar+geom_text(data=data,aes(label=sprintf("%0.2f", round(Perc, digits = 2))),position=position_stack(vjust=0.5)) 
  bar = bar+labs(x=xlab, y=ylab)+scale_fill_manual(name = name, values = cols)+theme(legend.position = "top")# + scale_fill_discrete(name = name)
  return(bar)
}
bar_plot <- function(data, width = 0.6, cols = color9,
                     xlab = '', ylab = '', name = 'Groups',
                     pvalue = T) {
  oridat <- data
  for (i in 1:ncol(data)) {
    data[, i] = data[, i] / sum(data[, i])
  }
  data=reshape2::melt(data)
  colnames(data)=c('Type','variable','Perc')
  data[, 1] = as.character(data[, 1])
  data[, 2] = as.character(data[, 2])
  bar = ggplot(data, aes(x=variable, y=Perc, fill=Type))+geom_bar(stat = "identity", width=width) + theme_get()
  bar = bar+geom_text(data=data,aes(label=sprintf("%0.2f", round(Perc, digits = 2))),position=position_stack(vjust=0.5))
  if (pvalue) {
    pvalue <- round(chisq.p.value(oridat), 4)
    bar = bar+labs(x=paste0('Chi-squared test P=', pvalue), y=ylab)+scale_fill_manual(name = name, values = cols)+theme(legend.position = "top")# + scale_fill_discrete(name = name)
  } else {
    bar = bar+labs(x=xlab, y=ylab)+scale_fill_manual(name = name, values = cols)+theme(legend.position = "top")# + scale_fill_discrete(name = name)
  }
  
  return(bar)
}
ggplotViolin_group <- function(dat, groups, cols = color9,
                               title = 'Groups',
                               xlab = '',
                               ylab = '',
                               angle = 30) {
  # dat: 行为样本，列为特征
  # groups：样本的分类,分类可以是两类也可以是多类
  # 这个是可以将相同分组的不同特征数据画在同一张图上
  library(ggplot2)
  library(ggpubr)
  features <- c()
  values <- c()
  Groups <- c()
  if (nrow(dat) == length(groups)) {
    for (fea in colnames(dat)) {
      features <- c(features, rep(fea, length(groups)))
      values <- c(values, as.numeric(dat[, fea]))
      Groups <- c(Groups, groups)
    }
    res <- data.frame(Groups = Groups,
                      Features = features,
                      Values = values)
    p <- ggplot(res, aes(x=Features, y=Values,fill=Groups)) +
      stat_compare_means(aes(group=Groups), label = "p.signif") + 
      # geom_violin(trim=F, color="white") + 
      geom_boxplot(position=position_dodge(0.95))+ 
      scale_fill_manual(values = cols)+
      theme_bw()+ 
      theme(axis.text.x=element_text(angle = angle, vjust = 1,  hjust = 1,colour="black"),
            # panel.border = element_blank(),axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),   #不显示网格线
            panel.grid.minor = element_blank(),
            legend.position = 'top') +  #不显示网格线
      labs(fill = title) + xlab(xlab) + ylab(ylab)
    return(p)
  } else {
    print('dat error')
    return(NULL)
  }
}

library(scales)
library(ggsci)
mypal <- pal_lancet(palette = c("lanonc"), alpha =1)(9)




genetag <- read.delim('e:/public/GeneTag.genecode.v32.txt', 
                      header = T, stringsAsFactors = F)
genetag <- genetag[!duplicated(genetag$ENSGID), ]
rownames(genetag) <- genetag$ENSGID
table(genetag$TYPE)
protein_genes <- genetag[genetag$TYPE == 'protein_coding', ]$SYMBOL
protein_genes

library(data.table)
# 1、数据整理 #############
# 1.0、缺氧基因的整理 ##############
HYPOXIA_genes_dat <- read.delim('raw_dat/MSigDB7.3_HALLMARK_HYPOXIA.txt')
HYPOXIA_genes <- HYPOXIA_genes_dat$HALLMARK_HYPOXIA
intersect(protein_genes, HYPOXIA_genes)

# 1.1、TCGA 数据预处理 ##########
luad_cli <- read.table('raw_dat/Clinical BCR XML.merge.txt', header = T, sep = '\t')
luad_cli <- luad_cli[, c("A0_Samples","A1_OS", "A2_Event", "age_at_initial_pathologic_diagnosis","A3_T", "A4_N", "A5_M", "A6_Stage", "tobacco_smoking_history", "A18_Sex")]
luad_cli <- reshape::rename(luad_cli, c(age_at_initial_pathologic_diagnosis = 'Age',A18_Sex = 'Gender', tobacco_smoking_history = 'Smoking', A1_OS = 'OS.time', A2_Event = 'OS'))
luad_cli$A3_T <- gsub('[ab]', '', luad_cli$A3_T)
luad_cli$A4_N[luad_cli$A4_N == ''] <- 'NX'
luad_cli$A5_M <- gsub('[ab]', '', luad_cli$A5_M)
luad_cli$A5_M[luad_cli$A5_M == ''] <- 'MX'
luad_cli$A6_Stage <- gsub('[AB]', '', luad_cli$A6_Stage)
luad_cli$A6_Stage <- gsub('Stage ', '', luad_cli$A6_Stage)
luad_cli$A6_Stage[luad_cli$A6_Stage == ''] <- 'X'
luad_cli$Age[luad_cli$Age == 'Not Available'] <- NA
luad_cli$Smoking[luad_cli$Smoking == 'Not Available'] <- c(7)
luad_cli$Smoking[luad_cli$Smoking == 'Unknown'] <- c(7)
luad_cli <- luad_cli[!is.na(luad_cli$OS.time) & luad_cli$OS.time > 0, ]
luad_cli$OS <- ifelse(luad_cli$OS == 'Alive', 0, 1)
luad_cli$Gender <- ifelse(luad_cli$Gender == 'FEMALE', 'Female', 'Male')
luad_cli$A0_Samples <- paste0(luad_cli$A0_Samples, '-01')
rownames(luad_cli) <- luad_cli$A0_Samples

library(data.table)
luad_tpm <- fread('raw_dat/TCGA_LUAD_TPM_genesymbol.txt', header = T, data.table = F)
rownames(luad_tpm) <- luad_tpm$ENSG_IDs
luad_tpm <- luad_tpm[, -1]
luad_tmr_samples <- intersect(luad_cli$A0_Samples, colnames(luad_tpm))
luad_cli <- luad_cli[luad_tmr_samples, ]
luad_tpm <- luad_tpm[, luad_tmr_samples]
luad_tpm_filtered <- luad_tpm[which(apply(luad_tpm,1,function(x){return(sum(x>=1))})>=0.5*ncol(luad_tpm)),]
luad_tpm_filtered_log <- log2(luad_tpm_filtered + 1)
dim(luad_tpm_filtered_log)
boxplot(luad_tpm_filtered_log[, 1:10], las = 2)

HYPOXIA_genes_filtered <- intersect(rownames(luad_tpm_filtered_log), HYPOXIA_genes)

# 2、TCGA 数据的一致性聚类 ###########
# 2.1、缺氧基因的单因素分析 ######
HYPOXIA_genes_cox <- Sig_cox(luad_tpm_filtered_log[HYPOXIA_genes_filtered, luad_cli$A0_Samples],
                             luad_cli$OS.time, 
                             luad_cli$OS)
write.csv(HYPOXIA_genes_cox,
          file = 'results/HYPOXIA_genes_cox.csv',
          row.names = F, quote = F)
table(HYPOXIA_genes_cox$P.Value < 0.05)
HYPOXIA_cox_genes <- HYPOXIA_genes_cox[HYPOXIA_genes_cox$P.Value < 0.05, ]$Genes

write.csv(HYPOXIA_cox_genes,
          file = 'results/S1.csv',
          row.names = F, col.names = F, quote = F)

# 2.2、一致性聚类分析 ######
library(ConsensusClusterPlus)
?ConsensusClusterPlus
tcga_results <- ConsensusClusterPlus(as.matrix(luad_tpm_filtered_log[HYPOXIA_cox_genes,]),
                                     maxK=10, 
                                     reps=100,
                                     pItem=0.8, 
                                     pFeature=1, 
                                     title = "ConsensusClusterPlus",
                                     clusterAlg="km", 
                                     distance="euclidean",
                                     corUse = "everything", 
                                     innerLinkage='ward.D2',
                                     seed=1984, 
                                     plot="pdf", 
                                     writeTable = T)
tcga_cluster <- data.frame(sample = names(tcga_results[[2]]$consensusClass),
                           Cluster = tcga_results[[2]]$consensusClass)
tcga_cluster$Cluster <- paste0('C', tcga_cluster$Cluster)
table(tcga_cluster$Cluster)

tcga_cluster_cli <- cbind(luad_cli, tcga_cluster[luad_cli$A0_Samples, ])



Cluster_km <- ggsurvplotKM(tcga_cluster_cli[, c("OS.time", "OS", "Cluster")],
                           title = 'Cluster',
                           lables = c('C1', 'C2'),
                           risk.table = T)
Cluster_km
ggsave(plot = Cluster_km,
       filename = 'PDFs/tcga_HYPOXIA_cluster_km.pdf',
       width = 5, height = 5)

# 2.3、t-SNE analysis ###############
library(Rtsne)
library(survminer)
citation('survminer')
set.seed(1984)
tsne_out<-Rtsne(t(luad_tpm_filtered_log[HYPOXIA_cox_genes, tcga_cluster_cli$A0_Samples]))

tsne_data<-data.frame(tsne_out$Y, tcga_cluster_cli$Cluster)
colnames(tsne_data)<-c("Y1","Y2","Cluster")

tsne_plot <- ggplot(tsne_data, aes(Y1, Y2, fill=Cluster)) + 
  geom_point(size=3, colour="black", alpha=1, shape=21) + 
  scale_fill_manual(values=c("#00468BFF", "#ED0000FF")) +
  # scale_fill_npg() + 
  xlab('Coordinate 1') + ylab('Coordinate 2') +
  theme_bw() + 
  theme(legend.direction = 'horizontal', legend.position = 'top')
tsne_plot
ggsave(plot = tsne_plot,
       filename = 'PDFs/tcga_cluster_tsne_plot.pdf',
       width = 5, height = 5)


# 2.4、缺氧基因在亚型中的热图 ###############
library(pheatmap)
library(dplyr)
tcga_cluster <- arrange(tcga_cluster, Cluster)
annotation_col  <- data.frame(Cluster = factor(tcga_cluster$Cluster))
rownames(annotation_col) <- tcga_cluster$sample
tcga_siggene_heatmap <- log2(luad_tpm[HYPOXIA_cox_genes, tcga_cluster$sample] + 1)
bk=unique(c(seq(-1.2, 1.2, length=100)))
ann_colors = list(
  Cluster = c(C1 = '#00468BFF', C2 = '#ED0000FF')
)
table(tcga_cluster$Cluster)
pdf('tcga_HYPOXIA_cox_genes_heatmap.pdf',width = 6,height = 6)
pheatmap(tcga_siggene_heatmap, 
         scale = 'row', 
         breaks = bk,
         annotation_col = annotation_col,
         cluster_cols = F, cluster_rows = T,
         show_rownames = F, show_colnames = F,
         gaps_col = 257, 
         annotation_colors =ann_colors, clustering_method = 'ward.D'
)
dev.off()


# 2.5、亚型缺氧通路的ssGSEA得分比较 ###########
library(GSVA)
library(GSEABase)
hpathway <- getGmt('e:/public/msigdb_v7.3_GMTs/h.all.v7.3.symbols.gmt')

tcga_h_ssgsea <- gsva(as.matrix(log2(luad_tpm[, tcga_cluster_cli$A0_Samples] + 1)), 
                   hpathway,
                   method='ssgsea',
                   kcdf='Gaussian',
                   abs.ranking=TRUE)
tcga_h_ssgsea <- as.data.frame(t(tcga_h_ssgsea))
tcga_h_ssgsea$Cluster <- tcga_cluster_cli[rownames(tcga_h_ssgsea), ]$Cluster

ggplotViolin_muti_group <- function(dat, cols = color8,
                                    title = 'Groups',
                                    xlab = '',
                                    ylab = '',
                                    angle = 0) {
  # dat: 行为样本，列为特征：列特征可以是两类也可以是多类
  # 这个是一个图片只能画单个特征，但是可以分成多组
  library(ggplot2)
  library(ggpubr)
  res <- dat
  colnames(res) <- c('Groups', 'Values')
  p <- ggplot(res, aes(x=Groups, y=Values)) +
    stat_compare_means(aes(group=Groups)) + 
    geom_violin(aes(fill = Groups)) + 
    geom_boxplot(width=0.15, position=position_dodge(0.95)) +
    scale_fill_manual(values = cols)+
    theme_bw()+ 
    theme(axis.text.x=element_text(angle = angle, vjust = 0.5,  hjust = 0.5,colour="black"),
          # panel.border = element_blank(),axis.line = element_line(colour = "black"),
          # panel.grid.major = element_blank(),   #不显示网格线
          # panel.grid.minor = element_blank(),
          legend.position = 'top') +  #不显示网格线
    labs(fill = title) + xlab(xlab) + ylab(ylab)
  return(p)
}

tcga_h_ssgsea_Violin <- ggplotViolin_muti_group(tcga_h_ssgsea[, c("Cluster", "HALLMARK_HYPOXIA")], 
                        title = 'Cluster', ylab = 'HALLMARK_HYPOXIA', 
                        cols = c('#00468BFF', '#ED0000FF'))
tcga_h_ssgsea_Violin
ggsave(plot = tcga_h_ssgsea_Violin,
       filename = 'PDFs/tcga_h_ssgsea_Violin.pdf', 
       width = 5, height = 5)


# 2.6、亚型之间差异基因的鉴定 ###############
tcga_cluster_deg <- DEGs_limma(luad_tpm_filtered_log[, tcga_cluster_cli$A0_Samples],
                               tcga_cluster_cli[, c("A0_Samples", "Cluster")],
                               'C2', 'C1')
tcga_cluster_deg_filtered <- tcga_cluster_deg[tcga_cluster_deg$adj.P.Val < 0.05 & abs(tcga_cluster_deg$logFC) > 1, ]
dim(tcga_cluster_deg_filtered)
table(tcga_cluster_deg_filtered$logFC > 0)
write.csv(tcga_cluster_deg_filtered,
          file = 'results/S3.csv')
tcga_cluster_deg_up <- rownames(tcga_cluster_deg_filtered[tcga_cluster_deg_filtered$logFC > 0, ])
tcga_cluster_deg_dn <- rownames(tcga_cluster_deg_filtered[tcga_cluster_deg_filtered$logFC < 0, ])


# 2.6.1、亚型差异基因的火山图与热图绘制 ############
wb_volcano <- function(gene, p.value, logFC,
                       cut_p = 0.05, cut_logFC = 1,
                       xlab = "log2FC", ylab = "-log10(FDR)",
                       title = 'State') {
  library(ggplot2)
  dat <- data.frame(Gene = gene, p.value = p.value, logFC = logFC)
  dat$type = ifelse(dat$p.value < cut_p & abs(dat$logFC) >= cut_logFC, ifelse(dat$logFC > cut_logFC, 'Up', 'Down'), 'Stable')

  p <- ggplot(dat, aes(x = logFC,
                       y = -log10(p.value))) +
    geom_point(aes(color = type),alpha=0.4, size=2) +
    scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
    # 辅助线
    geom_vline(xintercept=c(-cut_logFC,cut_logFC),
               lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(cut_p),lty=4,col="black",lwd=0.8) +
    # 坐标轴
    labs(x=xlab, y=ylab) + xlim(-3, 3) +
    theme_bw()+
    # 图例
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="top")
  return(p)
}

tcga_cluster_diffgene_volcano <- wb_volcano(rownames(tcga_cluster_deg), 
           tcga_cluster_deg$adj.P.Val, 
           tcga_cluster_deg$logFC, 
           cut_p = 0.05, 
           cut_logFC = 1)
tcga_cluster_diffgene_volcano
ggsave(plot = tcga_cluster_diffgene_volcano,
       filename = 'PDFs/tcga_cluster_diffgene_volcano.pdf',
       width = 5, height = 5)



library(pheatmap)
annotation_col  <- data.frame(Cluster = factor(tcga_cluster$Cluster))
rownames(annotation_col) <- tcga_cluster$sample
tcga_cluster_diff_exp_gene_heatmap <- luad_tpm_filtered_log[rownames(tcga_cluster_deg_filtered),  
                                                            tcga_cluster$sample]
bk=unique(c(seq(-1.2, 1.2, length=100)))
ann_colors = list(
  Cluster = c(C1 = '#00468BFF', C2 = '#ED0000FF')
)
pdf('tcga_subtype_diff_exp_gene_heatmap.pdf',width = 6, height = 6)
pheatmap(tcga_cluster_diff_exp_gene_heatmap, 
         scale = 'row', 
         breaks = bk,
         annotation_col = annotation_col,
         cluster_cols = F, cluster_rows = T,
         show_rownames = F, show_colnames = F,
         gaps_col = 257, #gaps_row = 50,
         # cellwidth = 0.8, cellheight = 0.35,
         # color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
         annotation_colors =ann_colors, clustering_method = 'ward.D')
dev.off()

# 2.7、亚型之间的 GSEA 分析 ###############
library(enrichplot)
library(ggplot2)
library(clusterProfiler)

plot_gsea1 <- function(dat, gene) {
  tmp <- as.data.frame(dat)
  for (i in 1:nrow(tmp)) {
    print(i)
    if (dat[i, "enrichmentScore"] > 0) {
      col <- 'red'
    } else {
      col <- 'green'
    }
    p <- gseaplot2(dat, geneSetID = i, 
                   title = dat$Description[i],
                   pvalue_table = F,
                   base_size = 10,
                   color = col,
                   rel_heights = c(1.5, 0.3, 1))
    anno <- dat[i, c("enrichmentScore", "pvalue", "NES", "p.adjust")]
    colnames(anno) <- c("ES", "pvalue", "NES", "FDR")
    class(anno)
    lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
    # print(lab)
    p <- p + annotate("text", 0.1, 0.5, label = lab, hjust = 0, vjust = 0)
    ggsave(p, filename = paste0(gene, '/', dat[i, "Description"], '.png'), width = 8, height = 6)
    ggsave(p, filename = paste0(gene, '/', dat[i, "Description"], '.pdf'), width = 8, height = 6)
  }
}

Pathways <- read.gmt("e:/public/msigdb_v7.3_GMTs/h.all.v7.3.symbols.gmt")
limma_gene <- arrange(tcga_cluster_deg, desc(logFC))
Ranks_gene <- limma_gene$logFC
names(Ranks_gene) <- rownames(limma_gene)
egmt2 <- GSEA(Ranks_gene,
              TERM2GENE=Pathways,
              verbose=FALSE,
              pvalueCutoff = 1)
clusterProfiler_res <- as.data.frame(egmt2)
dir.create('Cluster.GSEA')
write.csv(clusterProfiler_res,
          file = paste0('Cluster.GSEA', '/', 'Cluster.GSEA', '_GSEA_result.csv'))
plot_gsea1(egmt2, 'Cluster.GSEA')






# 3、免疫评分分组 ##################
tcga_estimate <- estimate_score(luad_tpm[, tcga_cluster_cli$A0_Samples])
class(tcga_estimate)
rownames(tcga_estimate) <- gsub('\\.', '-', rownames(tcga_estimate))

tcga_cluster_cli <- cbind(tcga_cluster_cli, tcga_estimate[tcga_cluster_cli$A0_Samples, ])



library(survival)
library(cutoff)
library(ggpubr)
library(survminer)
res.cut <-surv_cutpoint(tcga_cluster_cli, time = "OS.time", event = "OS",variables = c("ImmuneScore"))
summary(res.cut)
pdf('PDFs/ImmuneScore_cutoff.pdf', width = 8, height = 6)
plot(res.cut, "ImmuneScore", palette = "npg")
dev.off()


res.cat <- surv_categorize(res.cut)
head(res.cat)

tcga_cluster_cli$ImmuneGroups <- res.cat[tcga_cluster_cli$A0_Samples, ]$ImmuneScore

ggsurvplotKM(res.cat[, c("OS.time", "OS", "ImmuneScore")],
             # title = 'Cluster',
             # lables = c('C1', 'C2'),
             risk.table = T)



ImmuneGroups_km <- ggsurvplotKM(tcga_cluster_cli[, c("OS.time", "OS", "ImmuneGroups")],
                                title = 'ImmuneGroups',
                                lables = c('ImmuneHigh', 'ImmuneLow'),
                                risk.table = T)
ImmuneGroups_km
ggsave(plot = ImmuneGroups_km,
       filename = 'tcga_ImmuneGroups_km.pdf',
       width = 5, height = 5)


# 3.1、免疫分组的差异基因鉴定
tcga_cluster_cli$ImmuneGroups
tcga_ImmuneGroups_deg <- DEGs_limma(luad_tpm_filtered_log[, tcga_cluster_cli$A0_Samples],
                               tcga_cluster_cli[, c("A0_Samples", "ImmuneGroups")],
                               'low', 'high')
tcga_ImmuneGroups_deg_filtered <- tcga_ImmuneGroups_deg[tcga_ImmuneGroups_deg$adj.P.Val < 0.05 & abs(tcga_ImmuneGroups_deg$logFC) > 1, ]
table(tcga_ImmuneGroups_deg_filtered$logFC > 0)
dim(tcga_ImmuneGroups_deg_filtered)
write.csv(tcga_ImmuneGroups_deg_filtered,
          file = 'results/S4.csv')
tcga_ImmuneGroups_deg_up <- rownames(tcga_ImmuneGroups_deg_filtered[tcga_ImmuneGroups_deg_filtered$logFC > 0, ])
tcga_ImmuneGroups_deg_dn <- rownames(tcga_ImmuneGroups_deg_filtered[tcga_ImmuneGroups_deg_filtered$logFC < 0, ])

# 3.1.1、亚型差异基因的火山图与热图绘制 ############
tcga_ImmuneGroups_diffgene_volcano <- wb_volcano(rownames(tcga_ImmuneGroups_deg), 
                                            tcga_ImmuneGroups_deg$adj.P.Val, 
                                            tcga_ImmuneGroups_deg$logFC, 
                                            cut_p = 0.05, 
                                            cut_logFC = 1)
tcga_ImmuneGroups_diffgene_volcano
ggsave(plot = tcga_ImmuneGroups_diffgene_volcano,
       filename = 'PDFs/tcga_ImmuneGroups_diffgene_volcano.pdf',
       width = 5, height = 5)



library(pheatmap)
tcga_cluster_cli <- arrange(tcga_cluster_cli, ImmuneGroups)
table(tcga_cluster_cli$ImmuneGroups)
annotation_col  <- data.frame(ImmuneGroups = factor(tcga_cluster_cli$ImmuneGroups))
rownames(annotation_col) <- tcga_cluster_cli$A0_Samples
tcga_ImmuneGroups_diff_exp_gene_heatmap <- luad_tpm_filtered_log[rownames(tcga_ImmuneGroups_deg_filtered),  
                                                            tcga_cluster_cli$A0_Samples]
bk=unique(c(seq(-1.2, 1.2, length=100)))
ann_colors = list(
  ImmuneGroups = c(high = '#00468BFF', low = '#ED0000FF')
)
pdf('tcga_ImmuneGroups_diff_exp_gene_heatmap.pdf',width = 6, height = 6)
pheatmap(tcga_ImmuneGroups_diff_exp_gene_heatmap, 
         scale = 'row', 
         breaks = bk,
         annotation_col = annotation_col,
         cluster_cols = F, cluster_rows = T,
         show_rownames = F, show_colnames = F,
         gaps_col = 181,
         annotation_colors =ann_colors, clustering_method = 'ward.D2')
dev.off()




# 4、缺氧免疫分组的鉴定以及分析 ###########
tcga_cluster_cli$HypoxiaImmuneGroups <- ''
tcga_cluster_cli$HypoxiaImmuneGroups[tcga_cluster_cli$Cluster == 'C1' & tcga_cluster_cli$ImmuneGroups == 'high'] <- 'Hypoxia-low & Immune-high'
tcga_cluster_cli$HypoxiaImmuneGroups[tcga_cluster_cli$Cluster == 'C2' & tcga_cluster_cli$ImmuneGroups == 'low'] <- 'Hypoxia-high & Immune-low'
tcga_cluster_cli$HypoxiaImmuneGroups[tcga_cluster_cli$HypoxiaImmuneGroups == ''] <- 'Mixed'


HypoxiaImmuneGroups_km <- ggsurvplotKM(tcga_cluster_cli[, c("OS.time", "OS", "HypoxiaImmuneGroups")],
                                title = 'HypoxiaImmuneGroups',
                                lables = c('Hypoxia-high & Immune-low', 'Hypoxia-low & Immune-high', 'Mixed'),
                                risk.table = T)
HypoxiaImmuneGroups_km
ggsave(plot = HypoxiaImmuneGroups_km,
       filename = 'tcga_HypoxiaImmuneGroups_km.pdf',
       width = 5, height = 5)

# 4.1、缺氧免疫分组的差异基因鉴定
tcga_hypoImmune_cli <- tcga_cluster_cli[tcga_cluster_cli$HypoxiaImmuneGroups != 'Mixed', ]
table(tcga_hypoImmune_cli$HypoxiaImmuneGroups)
tcga_HypoxiaImmuneGroups_deg <- DEGs_limma(luad_tpm_filtered_log[, tcga_hypoImmune_cli$A0_Samples],
                                    tcga_hypoImmune_cli[, c("A0_Samples", "HypoxiaImmuneGroups")],
                                    'Hypoxia-high & Immune-low', 'Hypoxia-low & Immune-high')
tcga_HypoxiaImmuneGroups_deg_filtered <- tcga_HypoxiaImmuneGroups_deg[tcga_HypoxiaImmuneGroups_deg$adj.P.Val < 0.05 & abs(tcga_HypoxiaImmuneGroups_deg$logFC) > 1, ]
table(tcga_HypoxiaImmuneGroups_deg_filtered$logFC > 0)
dim(tcga_HypoxiaImmuneGroups_deg_filtered)

write.csv(tcga_HypoxiaImmuneGroups_deg_filtered,
          file = 'results/S5.csv')

tcga_HypoxiaImmuneGroups_deg_up <- rownames(tcga_HypoxiaImmuneGroups_deg_filtered[tcga_HypoxiaImmuneGroups_deg_filtered$logFC > 0, ])
tcga_HypoxiaImmuneGroups_deg_dn <- rownames(tcga_HypoxiaImmuneGroups_deg_filtered[tcga_HypoxiaImmuneGroups_deg_filtered$logFC < 0, ])

# 4.1.1、亚型差异基因的火山图与热图绘制 ############
tcga_HypoxiaImmuneGroups_diffgene_volcano <- wb_volcano(rownames(tcga_HypoxiaImmuneGroups_deg), 
                                                 tcga_HypoxiaImmuneGroups_deg$adj.P.Val, 
                                                 tcga_HypoxiaImmuneGroups_deg$logFC, 
                                                 cut_p = 0.05, 
                                                 cut_logFC = 1)
tcga_HypoxiaImmuneGroups_diffgene_volcano
ggsave(plot = tcga_HypoxiaImmuneGroups_diffgene_volcano,
       filename = 'PDFs/tcga_HypoxiaImmuneGroups_diffgene_volcano.pdf',
       width = 5, height = 5)



library(pheatmap)
tcga_hypoImmune_cli <- arrange(tcga_hypoImmune_cli, HypoxiaImmuneGroups)
table(tcga_hypoImmune_cli$HypoxiaImmuneGroups)
annotation_col  <- data.frame(HypoxiaImmuneGroups = factor(tcga_hypoImmune_cli$HypoxiaImmuneGroups))
rownames(annotation_col) <- tcga_hypoImmune_cli$A0_Samples
tcga_HypoxiaImmuneGroups_diff_exp_gene_heatmap <- luad_tpm_filtered_log[rownames(tcga_HypoxiaImmuneGroups_deg_filtered),  
                                                                 tcga_hypoImmune_cli$A0_Samples]
bk=unique(c(seq(-1.2, 1.2, length=100)))
ann_colors = list(
  HypoxiaImmuneGroups = c(`Hypoxia-high & Immune-low` = '#00468BFF', `Hypoxia-low & Immune-high` = '#ED0000FF')
)
pdf('tcga_HypoxiaImmuneGroups_diff_exp_gene_heatmap.pdf',width = 6, height = 6)
pheatmap(tcga_HypoxiaImmuneGroups_diff_exp_gene_heatmap, 
         scale = 'row', 
         breaks = bk,
         annotation_col = annotation_col,
         cluster_cols = F, cluster_rows = T,
         show_rownames = F, show_colnames = F,
         gaps_col = 147, 
         annotation_colors =ann_colors, clustering_method = 'ward.D2')
dev.off()



# 5、保护风险基因以及危险风险基因的鉴定以及功能富集分析
# 5.1、保护风险基因以及危险风险基因的鉴定以及功能富集分析 ##############
ProgDEGs <- unique(c(intersect(tcga_HypoxiaImmuneGroups_deg_dn, tcga_ImmuneGroups_deg_dn), 
                     intersect(tcga_HypoxiaImmuneGroups_deg_dn, tcga_cluster_deg_dn)))

library(venn)
# venn 绘图方法
pdf('PDFs/ProgDEGs-Venn.pdf', width = 5, height = 5)
venn(list(`Hypoxia-Immune` = tcga_HypoxiaImmuneGroups_deg_dn,
          Immune = tcga_ImmuneGroups_deg_dn,
          Hypoxia = tcga_cluster_deg_dn),
     zcolor = color8[1:3], 
     box = FALSE)
dev.off()
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
enrich_ProgDEGs <- bitr(ProgDEGs, fromType = "SYMBOL",
                          toType = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                          OrgDb = org.Hs.eg.db)

enrich_ProgDEGs_BP <- enrichGO(gene          = enrich_ProgDEGs$ENTREZID,
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 1,
                                 qvalueCutoff  = 1,
                                 readable      = TRUE)
head(enrich_ProgDEGs_BP)
ProgDEGsBP <- dotplot(enrich_ProgDEGs_BP,
                    showCategory = 10,
                    # color = "pvalue",
                    title = "ProgDEGs BP")
ProgDEGsBP

enrich_ProgDEGs_BP_dat <- as.data.frame(enrich_ProgDEGs_BP)


enrich_ProgDEGs_MF <- enrichGO(gene          = enrich_ProgDEGs$ENTREZID,
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "MF",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 1,
                                 qvalueCutoff  = 1,
                                 readable      = TRUE)
head(enrich_ProgDEGs_MF)
ProgDEGsMF <- dotplot(enrich_ProgDEGs_MF,
                    showCategory = 10,
                    # color = "pvalue",
                    title = "ProgDEGs MF")
ProgDEGsMF

enrich_ProgDEGs_MF_dat <- as.data.frame(enrich_ProgDEGs_MF)


enrich_ProgDEGs_CC <- enrichGO(gene          = enrich_ProgDEGs$ENTREZID,
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "CC",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 1,
                                 qvalueCutoff  = 1,
                                 readable      = TRUE)
head(enrich_ProgDEGs_CC)
ProgDEGsCC <- dotplot(enrich_ProgDEGs_CC,
                    showCategory = 10,
                    # color = "pvalue",
                    title = "ProgDEGs CC")
ProgDEGsCC

enrich_ProgDEGs_CC_dat <- as.data.frame(enrich_ProgDEGs_CC)


enrich_ProgDEGs_KEGG <- enrichKEGG(gene         = enrich_ProgDEGs$ENTREZID,
                                     organism     = 'hsa',
                                     pvalueCutoff = 1,
                                     pAdjustMethod = "BH")
head(enrich_ProgDEGs_KEGG)
ProgDEGskegg <- dotplot(enrich_ProgDEGs_KEGG,
                      showCategory = 10,
                      # color = "pvalue",
                      title = "ProgDEGs KEGG")
ProgDEGskegg
enrich_ProgDEGs_KEGG_dat <- setReadable(enrich_ProgDEGs_KEGG, 
                                          OrgDb = org.Hs.eg.db, 
                                          keyType="ENTREZID")
enrich_ProgDEGs_KEGG_dat <- as.data.frame(enrich_ProgDEGs_KEGG_dat)



colnames(enrich_ProgDEGs_KEGG_dat)
colnames(enrich_ProgDEGs_BP_dat)
colnames(enrich_ProgDEGs_MF_dat)
colnames(enrich_ProgDEGs_CC_dat)

enrich_ProgDEGs_KEGG_dat$TYPE <- 'KEGG'
enrich_ProgDEGs_BP_dat$TYPE <- 'BP'
enrich_ProgDEGs_MF_dat$TYPE <- 'MF'
enrich_ProgDEGs_CC_dat$TYPE <- 'CC'

enrich_ProgDEGs_GO_KEGG <- rbind(enrich_ProgDEGs_KEGG_dat,
                                   enrich_ProgDEGs_BP_dat,
                                   enrich_ProgDEGs_MF_dat,
                                   enrich_ProgDEGs_CC_dat)
table(enrich_ProgDEGs_GO_KEGG$TYPE)

write.csv(enrich_ProgDEGs_GO_KEGG,
          file = 'files/enrich_ProgDEGs_GO_KEGG.csv',
          row.names = F)

write.csv(enrich_ProgDEGs_GO_KEGG,
          file = 'results/enrich_ProgDEGs_GO_KEGG.csv',
          row.names = F)

ProgDEGsgo_kegg <- plot_grid(ProgDEGsBP, 
                          ProgDEGsMF,
                          ProgDEGsCC,
                          ProgDEGskegg, 
                          ncol=2, 
                          labels = LETTERS[1:4],
                          align = 'hv')

ProgDEGsgo_kegg
library(ggplot2)
ggsave(plot = ProgDEGsgo_kegg,
       filename = 'PDFs/ProgDEGsgo_kegg.pdf',
       width = 18,height = 10)


# 5.2、危险风险基因的鉴定以及功能富集分析 ##############
RiskDEGs <- unique(c(intersect(tcga_HypoxiaImmuneGroups_deg_up, tcga_ImmuneGroups_deg_up), 
                     intersect(tcga_HypoxiaImmuneGroups_deg_up, tcga_cluster_deg_up)))
pdf('PDFs/RiskDEGs-Venn.pdf', width = 5, height = 5)
venn(list(`Hypoxia-Immune` = tcga_HypoxiaImmuneGroups_deg_up,
          Immune = tcga_ImmuneGroups_deg_up,
          Hypoxia = tcga_cluster_deg_up),
     zcolor = color8[1:3], 
     box = FALSE)
dev.off()

enrich_RiskDEGs <- bitr(RiskDEGs, fromType = "SYMBOL",
                        toType = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                        OrgDb = org.Hs.eg.db)

enrich_RiskDEGs_BP <- enrichGO(gene          = enrich_RiskDEGs$ENTREZID,
                               OrgDb         = org.Hs.eg.db,
                               ont           = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 1,
                               qvalueCutoff  = 1,
                               readable      = TRUE)
head(enrich_RiskDEGs_BP)
RiskDEGsBP <- dotplot(enrich_RiskDEGs_BP,
                      showCategory = 10,
                      # color = "pvalue",
                      title = "RiskDEGs BP")
RiskDEGsBP

enrich_RiskDEGs_BP_dat <- as.data.frame(enrich_RiskDEGs_BP)


enrich_RiskDEGs_MF <- enrichGO(gene          = enrich_RiskDEGs$ENTREZID,
                               OrgDb         = org.Hs.eg.db,
                               ont           = "MF",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 1,
                               qvalueCutoff  = 1,
                               readable      = TRUE)
head(enrich_RiskDEGs_MF)
RiskDEGsMF <- dotplot(enrich_RiskDEGs_MF,
                      showCategory = 10,
                      # color = "pvalue",
                      title = "RiskDEGs MF")
RiskDEGsMF

enrich_RiskDEGs_MF_dat <- as.data.frame(enrich_RiskDEGs_MF)


enrich_RiskDEGs_CC <- enrichGO(gene          = enrich_RiskDEGs$ENTREZID,
                               OrgDb         = org.Hs.eg.db,
                               ont           = "CC",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 1,
                               qvalueCutoff  = 1,
                               readable      = TRUE)
head(enrich_RiskDEGs_CC)
RiskDEGsCC <- dotplot(enrich_RiskDEGs_CC,
                      showCategory = 10,
                      # color = "pvalue",
                      title = "RiskDEGs CC")
RiskDEGsCC

enrich_RiskDEGs_CC_dat <- as.data.frame(enrich_RiskDEGs_CC)


enrich_RiskDEGs_KEGG <- enrichKEGG(gene         = enrich_RiskDEGs$ENTREZID,
                                   organism     = 'hsa',
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "BH")
head(enrich_RiskDEGs_KEGG)
RiskDEGskegg <- dotplot(enrich_RiskDEGs_KEGG,
                        showCategory = 10,
                        # color = "pvalue",
                        title = "RiskDEGs KEGG")
RiskDEGskegg
enrich_RiskDEGs_KEGG_dat <- setReadable(enrich_RiskDEGs_KEGG, 
                                        OrgDb = org.Hs.eg.db, 
                                        keyType="ENTREZID")
enrich_RiskDEGs_KEGG_dat <- as.data.frame(enrich_RiskDEGs_KEGG_dat)



colnames(enrich_RiskDEGs_KEGG_dat)
colnames(enrich_RiskDEGs_BP_dat)
colnames(enrich_RiskDEGs_MF_dat)
colnames(enrich_RiskDEGs_CC_dat)

enrich_RiskDEGs_KEGG_dat$TYPE <- 'KEGG'
enrich_RiskDEGs_BP_dat$TYPE <- 'BP'
enrich_RiskDEGs_MF_dat$TYPE <- 'MF'
enrich_RiskDEGs_CC_dat$TYPE <- 'CC'

enrich_RiskDEGs_GO_KEGG <- rbind(enrich_RiskDEGs_KEGG_dat,
                                 enrich_RiskDEGs_BP_dat,
                                 enrich_RiskDEGs_MF_dat,
                                 enrich_RiskDEGs_CC_dat)
table(enrich_RiskDEGs_GO_KEGG$TYPE)

write.csv(enrich_RiskDEGs_GO_KEGG,
          file = 'files/enrich_RiskDEGs_GO_KEGG.csv',
          row.names = F)

write.csv(enrich_RiskDEGs_GO_KEGG,
          file = 'results/enrich_RiskDEGs_GO_KEGG.csv',
          row.names = F)

RiskDEGsgo_kegg <- plot_grid(RiskDEGsBP, 
                             RiskDEGsMF,
                             RiskDEGsCC,
                             RiskDEGskegg, 
                             ncol=2, 
                             labels = LETTERS[1:4],
                             align = 'hv')

RiskDEGsgo_kegg
library(ggplot2)
ggsave(plot = RiskDEGsgo_kegg,
       filename = 'PDFs/RiskDEGsgo_kegg.pdf',
       width = 17,height = 10)


# 6、模型数据准备 #################
# 6.1、tcga_model_dat #############
boxplot(luad_tpm_filtered_log[, 1:10])
rownames(tcga_cluster_cli) <- tcga_cluster_cli$A0_Samples
tcga_model_dat <- cbind(tcga_cluster_cli[, c("OS.time", "OS")],
                        t(luad_tpm_filtered_log[c(ProgDEGs, RiskDEGs), 
                                        tcga_cluster_cli$A0_Samples]))
colnames(tcga_model_dat) <- gsub('-', '__', colnames(tcga_model_dat))
class(tcga_model_dat)
tcga_model_dat <- as.data.frame(tcga_model_dat, stringsAsFactors = F)
str(tcga_model_dat)

# 6.2、GSE30219_model_dat ##################
GSE30219_model_dat <- read.delim('raw_dat/GSE30219/GSE30219_model_dat.txt', row.names = 1)

# 6.3、GSE31210_model_dat ##################
GSE31210_model_dat <- read.delim('raw_dat/GSE31210/GSE31210_model_dat.txt', row.names = 1)

# 6.4、GSE50081_model_dat ##################
GSE50081_model_dat <- read.delim('raw_dat/GSE50081/GSE50081_model_dat.txt', row.names = 1)





model_genes <- intersect(intersect(colnames(tcga_model_dat), 
                                   colnames(GSE30219_model_dat)),
                         intersect(colnames(GSE31210_model_dat), 
                                   colnames(GSE50081_model_dat)))
model_genes
tcga_model_dat <- tcga_model_dat[, model_genes]








# 7、模型筛选 #################
num <- 31
# 7.1、开始验证 #######
tra.samples <- rownames(read.delim(paste0('model_select/traindat_median_',num,'.txt'), 
                                   header = T, row.names = 1, stringsAsFactors = F))
test.samples <- rownames(read.delim(paste0('model_select/testdat_median_',num,'.txt'),
                                    header = T, row.names = 1, stringsAsFactors = F))

tra.data <- tcga_model_dat[tra.samples, ]
test.data <- tcga_model_dat[test.samples, ]


trainSingleCox <- Sig_cox(dat = t(tra.data[, -c(1, 2)]), 
                          time = tra.data[, 1], 
                          status = tra.data[, 2])
trainSingleCox <- na.omit(trainSingleCox)
dim(trainSingleCox)
cutCox <- trainSingleCox[which(trainSingleCox[, 2] < 0.05), ]
dim(cutCox)
write.csv(cutCox,
          file = 'results/S6.csv',
          quote = F, row.names = F)
SingleGene <- rownames(cutCox)



library(glmnet)
set.seed(num)
fit1=glmnet(as.matrix(tra.data[,SingleGene]),
            cbind(time=tra.data$OS.time,
                  status=tra.data$OS),
            family="cox",
            nlambda=100, 
            alpha=1) 

cv.fit<-cv.glmnet(as.matrix(tra.data[,SingleGene]),
                  cbind(time=tra.data$OS.time,
                        status=tra.data$OS),
                  family="cox",
                  nlambda=100, 
                  alpha=1)
sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
cv.fit$lambda.min
filter_genes <- names(sig.coef)

pdf('PDFs/lasso.pdf',width = 10,height = 5)
par(mfrow=c(1, 2))
plot(fit1, xvar="lambda")
plot(cv.fit)
dev.off()

library(survival)
# 基因的多因素
tcga_dat1 <- cbind(time=tra.data$OS.time,
                   status=tra.data$OS,
                   tra.data[,filter_genes])

fmla <- as.formula(paste0("Surv(time, status) ~"
                          ,paste0(filter_genes,collapse = '+')))


cox <- coxph(fmla, data =as.data.frame(tcga_dat1))



cox1 <- step(cox, trace = 0)
lan <- coef(cox1)
round(lan, 3)
genes <- names(cox1$coefficients)



summary(cox1)
pdf('PDFs/model_genes_multi_cox.pdf', width = 8, height = 4, onefile = FALSE)
survminer::ggforest(cox1,data=tcga_dat1)
dev.off()


Siggene_KM <- function(dat, gene) {
  p <- ggsurvplotKM(data.frame(dat$OS.time,
                               dat$OS,
                               ifelse(dat[,gene] > median(dat[,gene]),
                                      'High', 'Low')),
                    title = gene,
                    lables = c('High Exp', 'Low Exp'),
                    risk.table = F)
  return(p)
}
KM_MS4A1 <- Siggene_KM(tra.data, "MS4A1")
KM_MS4A1
KM_CPA3 <- Siggene_KM(tra.data, "CPA3")
KM_CPA3
KM_FSCN1 <- Siggene_KM(tra.data, "FSCN1")
KM_FSCN1
KM_PTPRH <- Siggene_KM(tra.data, "PTPRH")
KM_PTPRH
KM_DKK1 <- Siggene_KM(tra.data, "DKK1")
KM_DKK1



gene_KM <- cowplot::plot_grid(KM_MS4A1,
                              KM_CPA3,
                              KM_FSCN1,
                              KM_PTPRH,
                              KM_DKK1,
                              ncol=3, nrow = 2, labels=LETTERS[1:5],
                              align = 'hv')
gene_KM
library(ggplot2)
ggsave(plot = gene_KM,
       filename = 'PDFs/gene_KM.pdf',
       width = 12, height = 8)

library(ggthemes)
# 5.2、数据集的ROC和KM曲线 #####
risk.tra <- as.numeric(lan%*%as.matrix(t(tra.data[,genes])))
tra.data.plot <- tra.data[, c("OS.time", "OS", genes)]
tra.data.plot$RiskScore <- risk.tra
tra.data.plot$RiskType <- ifelse(tra.data.plot$RiskScore > median(tra.data.plot$RiskScore), 'High', 'Low')
tra_plot <- plot_RS_model(tra.data.plot, genes, position = c(0, 1), breaks = c(1,3,5))
tra_plot
ggsave(plot = tra_plot,
       filename = 'PDFs/tra_plot.pdf',
       width = 8, height = 10)

# 测试集
risk.test <- as.numeric(lan%*%as.matrix(t(test.data[,genes])))

test.data.plot <- test.data[, c("OS.time", "OS", genes)]
test.data.plot$RiskScore <- risk.test
test.data.plot$RiskType <- ifelse(test.data.plot$RiskScore > median(test.data.plot$RiskScore), 'High', 'Low')
test_plot <- plot_RS_model(test.data.plot, genes, position = c(0, 1), breaks = c(1,3,5))
test_plot
ggsave(plot = test_plot,
       filename = 'PDFs/test_plot.pdf',
       width = 8, height = 10)

# 全部数据集
risk.all <- as.numeric(lan%*%as.matrix(t(tcga_model_dat[,genes])))

all.data.plot <- tcga_model_dat[, c("OS.time", "OS", genes)]
all.data.plot$RiskScore <- risk.all
all.data.plot$RiskType <- ifelse(all.data.plot$RiskScore > median(all.data.plot$RiskScore), 'High', 'Low')
dim(all.data.plot)
all_plot <- plot_RS_model(all.data.plot, genes, position = c(0, 1), breaks = c(1,3,5))
all_plot
ggsave(plot = all_plot,
       filename = 'PDFs/all_plot.pdf',
       width = 8, height = 10)

library(ggthemes)
# GSE30219 数据集
risk.GSE30219 <- as.numeric(lan%*%as.matrix(t(GSE30219_model_dat[,genes])))


GSE30219.data.plot <- GSE30219_model_dat[, c("OS.time", "OS", genes)]
GSE30219.data.plot$RiskScore <- risk.GSE30219
GSE30219.data.plot$RiskType <- ifelse(GSE30219.data.plot$RiskScore > median(GSE30219.data.plot$RiskScore), 'High', 'Low')
dim(GSE30219.data.plot)
GSE30219_plot <- plot_RS_model(GSE30219.data.plot, genes, position = c(0, 1), breaks = c(1,3,5))
GSE30219_plot
ggsave(plot = GSE30219_plot,
       filename = 'PDFs/GSE30219_plot.pdf',
       width = 8, height = 10)


# GSE31210 数据集
cox.GSE31210 <- coxph(as.formula(paste0("Surv(OS.time, OS) ~"
                                        ,paste0(genes,collapse = '+'))),
                      data =as.data.frame(GSE31210_model_dat))
lan.GSE31210 <- coef(cox.GSE31210)
risk.GSE31210 <- as.numeric(lan.GSE31210%*%as.matrix(t(GSE31210_model_dat[,genes])))


GSE31210.data.plot <- GSE31210_model_dat[, c("OS.time", "OS", genes)]
GSE31210.data.plot$RiskScore <- risk.GSE31210
GSE31210.data.plot$RiskType <- ifelse(GSE31210.data.plot$RiskScore > median(GSE31210.data.plot$RiskScore), 'High', 'Low')
dim(GSE31210.data.plot)
GSE31210_plot <- plot_RS_model(GSE31210.data.plot, genes, position = c(0, 1), breaks = c(1,3,5))
GSE31210_plot
ggsave(plot = GSE31210_plot,
       filename = 'PDFs/GSE31210_plot.pdf',
       width = 8, height = 10)




# GSE50081 数据集
cox.GSE50081 <- coxph(as.formula(paste0("Surv(OS.time, OS) ~"
                                        ,paste0(genes,collapse = '+'))),
                      data =as.data.frame(GSE50081_model_dat))
lan.GSE50081 <- coef(cox.GSE50081)
risk.GSE50081 <- as.numeric(lan.GSE50081%*%as.matrix(t(GSE50081_model_dat[,genes])))

# risk.GSE50081 <- as.numeric(lan%*%as.matrix(t(GSE50081_model_dat[,genes])))


GSE50081.data.plot <- GSE50081_model_dat[, c("OS.time", "OS", genes)]
GSE50081.data.plot$RiskScore <- risk.GSE50081
GSE50081.data.plot$RiskType <- ifelse(GSE50081.data.plot$RiskScore > median(GSE50081.data.plot$RiskScore), 'High', 'Low')
dim(GSE50081.data.plot)
GSE50081_plot <- plot_RS_model(GSE50081.data.plot, genes, position = c(0, 1), breaks = c(1,3,5))
GSE50081_plot
ggsave(plot = GSE50081_plot,
       filename = 'PDFs/GSE50081_plot.pdf',
       width = 8, height = 10)



# TCGA 临床特征的风险得分比较 ###########
tcga_cli_RS <- tcga_cluster_cli[rownames(tcga_model_dat), ]
tcga_cli_RS$RiskScore <- risk.all
tcga_cli_RS$RiskType <- ifelse(tcga_cli_RS$RiskScore > median(tcga_cli_RS$RiskScore), 'High', 'Low')
fivenum(as.numeric(tcga_cli_RS$Age))
tcga_cli_RS$Age1 <- ifelse(as.numeric(tcga_cli_RS$Age) > 65, '>65', '<=65')

A3_T.dat <- tcga_cli_RS[, c("A3_T", "RiskScore")]
A3_T.dat <- A3_T.dat[A3_T.dat$A3_T != 'TX', ]
A3_T.violin <- ggplotViolin_muti_group(A3_T.dat, 
                                  title = 'T.Stage', 
                                  ylab = 'RiskScore')
A3_T.violin

A4_N.dat <- tcga_cli_RS[, c("A4_N", "RiskScore")]
A4_N.dat <- A4_N.dat[A4_N.dat$A4_N != 'NX' & A4_N.dat$A4_N != 'N3', ]
A4_N.violin <- ggplotViolin_muti_group(A4_N.dat, 
                                  title = 'N.Stage', 
                                  ylab = 'RiskScore')
A4_N.violin

A5_M.dat <- tcga_cli_RS[, c("A5_M", "RiskScore")]
A5_M.dat <- A5_M.dat[A5_M.dat$A5_M != 'MX', ]
A5_M.violin <- ggplotViolin_muti_group(A5_M.dat, 
                                  title = 'M.Stage', 
                                  ylab = 'RiskScore')
A5_M.violin

A6_Stage.dat <- tcga_cli_RS[, c("A6_Stage", "RiskScore")]
A6_Stage.dat <- A6_Stage.dat[A6_Stage.dat$A6_Stage != 'X', ]
A6_Stage.violin <- ggplotViolin_muti_group(A6_Stage.dat, 
                                      title = 'Stage', 
                                      ylab = 'RiskScore')
A6_Stage.violin

Gender.dat <- tcga_cli_RS[, c("Gender", "RiskScore")]
Gender.violin <- ggplotViolin_muti_group(Gender.dat, 
                                    title = 'Gender', 
                                    ylab = 'RiskScore')
Gender.violin

Smoking.dat <- tcga_cli_RS[, c("Smoking", "RiskScore")]
Smoking.dat <- Smoking.dat[Smoking.dat$Smoking != 5 & Smoking.dat$Smoking != 7, ]
Smoking.violin <- ggplotViolin_muti_group(Smoking.dat, 
                                     title = 'Smoking', 
                                     ylab = 'RiskScore')
Smoking.violin

Age1.dat <- tcga_cli_RS[, c("Age1", "RiskScore")]
Age1.dat <- na.omit(Age1.dat)
Age1.violin <- ggplotViolin_muti_group(Age1.dat, 
                                  title = 'Age', 
                                  ylab = 'RiskScore')
Age1.violin

Cluster.dat <- tcga_cli_RS[, c("Cluster", "RiskScore")]
Cluster.violin <- ggplotViolin_muti_group(Cluster.dat, 
                                     title = 'Cluster', 
                                     ylab = 'RiskScore')
Cluster.violin

ImmuneGroups.dat <- tcga_cli_RS[, c("ImmuneGroups", "RiskScore")]
ImmuneGroups.violin <- ggplotViolin_muti_group(ImmuneGroups.dat, 
                                          title = 'ImmuneGroups', 
                                          ylab = 'RiskScore')
ImmuneGroups.violin

hypoxiaImmuneGroups.dat <- tcga_cli_RS[, c("HypoxiaImmuneGroups", "RiskScore")]
HypoxiaImmuneGroups.violin <- ggplotViolin_muti_group(hypoxiaImmuneGroups.dat, 
                                               title = 'HypoxiaImmuneGroups', 
                                               ylab = 'RiskScore')
HypoxiaImmuneGroups.violin



tcga_cli_RS_compare <- cowplot::plot_grid(A3_T.violin,
                                          A4_N.violin,
                                          A5_M.violin,
                                          A6_Stage.violin, 
                                          Cluster.violin,
                                          ImmuneGroups.violin,
                                          HypoxiaImmuneGroups.violin,
                                          Smoking.violin, 
                                          Gender.violin,
                                          Age1.violin,
                                          ncol=4, labels=LETTERS[1:10], 
                                          align = 'hv')
tcga_cli_RS_compare
ggsave("PDFs/tcga_cli_RS_compare.pdf", 
       tcga_cli_RS_compare, 
       width = 16, 
       height = 12)


# TCGA 风险分组在临床特征上的表现 ###########
age_km <- tcga_cli_RS[, c("OS.time", "OS", "Age1", "RiskType")]
table(age_km$Age1)
age_km1 <- age_km[age_km$Age1 == '>65', ]
age_km1 <- na.omit(age_km1)
age_km1 <- ggsurvplotKM(data.frame(age_km1$OS.time,
                                   age_km1$OS,
                                   age_km1$RiskType),
                        title = 'Age > 65',
                        lables = c('High', 'Low'))
age_km1

age_km2 <- age_km[age_km$Age1 == '<=65', ]
age_km2 <- na.omit(age_km2)
age_km2 <- ggsurvplotKM(data.frame(age_km2$OS.time,
                                   age_km2$OS,
                                   age_km2$RiskType),
                        title = 'Age ≤ 65',
                        lables = c('High', 'Low'))
age_km2


Gender_km <- tcga_cli_RS[, c("OS.time", "OS", "Gender", "RiskType")]
table(Gender_km$Gender)
Gender_km1 <- Gender_km[Gender_km$Gender == 'Female', ]
Gender_km1 <- na.omit(Gender_km1)
Gender_km1 <- ggsurvplotKM(data.frame(Gender_km1$OS.time,
                                      Gender_km1$OS,
                                      Gender_km1$RiskType),
                           title = 'Female',
                           lables = c('High', 'Low'))
Gender_km1

Gender_km2 <- Gender_km[Gender_km$Gender == 'Male', ]
Gender_km2 <- na.omit(Gender_km2)
Gender_km2 <- ggsurvplotKM(data.frame(Gender_km2$OS.time,
                                      Gender_km2$OS,
                                      Gender_km2$RiskType),
                           title = 'Male',
                           lables = c('High', 'Low'))
Gender_km2



A3_T_km <- tcga_cli_RS[, c("OS.time", "OS", "A3_T", "RiskType")]
table(A3_T_km$A3_T)
A3_T_km1 <- A3_T_km[A3_T_km$A3_T == 'T1' | A3_T_km$A3_T == 'T2', ]
A3_T_km1 <- na.omit(A3_T_km1)
A3_T_km1 <- ggsurvplotKM(data.frame(A3_T_km1$OS.time,
                                    A3_T_km1$OS,
                                    A3_T_km1$RiskType),
                         title = 'T1-T2',
                         lables = c('High', 'Low'))
A3_T_km1

A3_T_km2 <- A3_T_km[A3_T_km$A3_T == 'T3' | A3_T_km$A3_T == 'T4', ]
A3_T_km2 <- na.omit(A3_T_km2)
A3_T_km2 <- ggsurvplotKM(data.frame(A3_T_km2$OS.time,
                                    A3_T_km2$OS,
                                    A3_T_km2$RiskType),
                         title = 'T3-T4',
                         lables = c('High', 'Low'))
A3_T_km2


A4_N_km <- tcga_cli_RS[, c("OS.time", "OS", "A4_N", "RiskType")]
table(A4_N_km$A4_N)
A4_N_km1 <- A4_N_km[A4_N_km$A4_N == 'N0', ]
A4_N_km1 <- na.omit(A4_N_km1)
A4_N_km1 <- ggsurvplotKM(data.frame(A4_N_km1$OS.time,
                                    A4_N_km1$OS,
                                    A4_N_km1$RiskType),
                         title = 'N0',
                         lables = c('High', 'Low'))
A4_N_km1

A4_N_km2 <- A4_N_km[A4_N_km$A4_N == 'N1' | A4_N_km$A4_N == 'N2' | A4_N_km$A4_N == 'N3', ]
A4_N_km2 <- na.omit(A4_N_km2)
A4_N_km2 <- ggsurvplotKM(data.frame(A4_N_km2$OS.time,
                                    A4_N_km2$OS,
                                    A4_N_km2$RiskType),
                         title = 'N1-N3',
                         lables = c('High', 'Low'))
A4_N_km2


A5_M_km <- tcga_cli_RS[, c("OS.time", "OS", "A5_M", "RiskType")]
table(A5_M_km$A5_M)
A5_M_km1 <- A5_M_km[A5_M_km$A5_M == 'M0', ]
A5_M_km1 <- na.omit(A5_M_km1)
A5_M_km1 <- ggsurvplotKM(data.frame(A5_M_km1$OS.time,
                                    A5_M_km1$OS,
                                    A5_M_km1$RiskType),
                         title = 'M0',
                         lables = c('High', 'Low'))
A5_M_km1

A5_M_km2 <- A5_M_km[A5_M_km$A5_M == 'M1', ]
A5_M_km2 <- na.omit(A5_M_km2)
A5_M_km2 <- ggsurvplotKM(data.frame(A5_M_km2$OS.time,
                                    A5_M_km2$OS,
                                    A5_M_km2$RiskType),
                         title = 'M1',
                         lables = c('High', 'Low'))
A5_M_km2


A6_Stage_km <- tcga_cli_RS[, c("OS.time", "OS", "A6_Stage", "RiskType")]
table(A6_Stage_km$A6_Stage)
A6_Stage_km1 <- A6_Stage_km[A6_Stage_km$A6_Stage == 'I' | A6_Stage_km$A6_Stage == 'II', ]
A6_Stage_km1 <- na.omit(A6_Stage_km1)
A6_Stage_km1 <- ggsurvplotKM(data.frame(A6_Stage_km1$OS.time,
                                        A6_Stage_km1$OS,
                                        A6_Stage_km1$RiskType),
                             title = 'I-II',
                             lables = c('High', 'Low'))
A6_Stage_km1

A6_Stage_km2 <- A6_Stage_km[A6_Stage_km$A6_Stage == 'III' | A6_Stage_km$A6_Stage == 'IV', ]
A6_Stage_km2 <- na.omit(A6_Stage_km2)
A6_Stage_km2 <- ggsurvplotKM(data.frame(A6_Stage_km2$OS.time,
                                        A6_Stage_km2$OS,
                                        A6_Stage_km2$RiskType),
                             title = 'III-IV',
                             lables = c('High', 'Low'))
A6_Stage_km2



Smoking_km <- tcga_cli_RS[, c("OS.time", "OS", "Smoking", "RiskType")]
table(Smoking_km$Smoking)
Smoking_km1 <- Smoking_km[Smoking_km$Smoking == 1, ]
Smoking_km1 <- na.omit(Smoking_km1)
Smoking_km1 <- ggsurvplotKM(data.frame(Smoking_km1$OS.time,
                                       Smoking_km1$OS,
                                       Smoking_km1$RiskType),
                            title = 'Smoking 1',
                            lables = c('High', 'Low'))
Smoking_km1

Smoking_km2 <- Smoking_km[Smoking_km$Smoking == 2 | Smoking_km$Smoking == 3 | Smoking_km$Smoking == 4, ]
Smoking_km2 <- na.omit(Smoking_km2)
Smoking_km2 <- ggsurvplotKM(data.frame(Smoking_km2$OS.time,
                                       Smoking_km2$OS,
                                       Smoking_km2$RiskType),
                            title = 'Smoking 2-4',
                            lables = c('High', 'Low'))
Smoking_km2

tcga_cli_RS_KM <- cowplot::plot_grid(A3_T_km1,
                                     A3_T_km2,
                                     A4_N_km1,
                                     A4_N_km2,
                                     A5_M_km1,
                                     A5_M_km2, 
                                     A6_Stage_km1, 
                                     A6_Stage_km2,
                                     age_km1,
                                     age_km2,
                                     Gender_km1,
                                     Gender_km2,
                                     Smoking_km1,
                                     Smoking_km2, 
                                     nrow = 4, ncol = 4,
                                     align = 'hv',
                                     labels = LETTERS[1:14])
tcga_cli_RS_KM
ggsave("PDFs/tcga_cli_RS_KM.pdf", tcga_cli_RS_KM, width = 20, height = 20)



# TCGA 单因素与多因素分析 ##############
tcga_cox_dat <- tcga_cli_RS

library(forestplot)
library(survcomp)
tcga_cox_dat[tcga_cox_dat == ''] <- NA
table(tcga_cox_dat$RiskType)
tcga_cox_dat$RiskType <- ifelse(tcga_cox_dat$RiskType == 'High', '1High', '0Low')


table(tcga_cox_dat$Age1)
tcga_cox_dat$Age1 <- ifelse(tcga_cox_dat$Age1 == '>65', '1>65', '0≤65')


table(tcga_cox_dat$A3_T)
tcga_cox_dat$A3_T[tcga_cox_dat$A3_T == 'T1' | tcga_cox_dat$A3_T == 'T2'] <- 'T1+T2'
tcga_cox_dat$A3_T[tcga_cox_dat$A3_T == 'T3' | tcga_cox_dat$A3_T == 'T4'] <- 'T3+T4'
tcga_cox_dat$A3_T[tcga_cox_dat$A3_T == 'TX'] <- NA

table(tcga_cox_dat$A4_N)
tcga_cox_dat$A4_N[tcga_cox_dat$A4_N == 'N1' | tcga_cox_dat$A4_N == 'N2' | tcga_cox_dat$A4_N == 'N3'] <- 'N1+N2+N3'
tcga_cox_dat$A4_N[tcga_cox_dat$A4_N == 'NX'] <- NA

table(tcga_cox_dat$A5_M)
tcga_cox_dat$A5_M[tcga_cox_dat$A5_M == 'MX'] <- NA

table(tcga_cox_dat$A6_Stage)
tcga_cox_dat$A6_Stage[tcga_cox_dat$A6_Stage == 'I' | tcga_cox_dat$A6_Stage == 'II'] <- 'I+II'
tcga_cox_dat$A6_Stage[tcga_cox_dat$A6_Stage == 'III' | tcga_cox_dat$A6_Stage == 'IV'] <- 'III+IV'
tcga_cox_dat$A6_Stage[tcga_cox_dat$A6_Stage == 'X'] <- NA

table(tcga_cox_dat$Gender)


table(tcga_cox_dat$Smoking)
tcga_cox_dat$Smoking[tcga_cox_dat$Smoking == 2 | tcga_cox_dat$Smoking == 3 | tcga_cox_dat$Smoking == 4] <- '2+3+4'
tcga_cox_dat$Smoking[tcga_cox_dat$Smoking == 5 | tcga_cox_dat$Smoking == 7] <- NA



library(survival)
age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age1,
                             data=tcga_cox_dat))
age_sig_cox_dat <- data.frame(Names=rownames(age_sig_cox[[8]]),
                              HR = round(age_sig_cox[[7]][,2],3),
                              lower.95 = round(age_sig_cox[[8]][,3],3),
                              upper.95 = round(age_sig_cox[[8]][,4],3),
                              p.value=round(age_sig_cox[[7]][,5],3))
age_sig_cox_dat

Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                                data=tcga_cox_dat))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_sig_cox[[8]]),
                                 HR = round(Gender_sig_cox[[7]][,2],3),
                                 lower.95 = round(Gender_sig_cox[[8]][,3],3),
                                 upper.95 = round(Gender_sig_cox[[8]][,4],3),
                                 p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat

A3_T_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~A3_T,
                              data=tcga_cox_dat))
A3_T_sig_cox_dat <- data.frame(Names=rownames(A3_T_sig_cox[[8]]),
                               HR = round(A3_T_sig_cox[[7]][,2],3),
                               lower.95 = round(A3_T_sig_cox[[8]][,3],3),
                               upper.95 = round(A3_T_sig_cox[[8]][,4],3),
                               p.value=round(A3_T_sig_cox[[7]][,5],3))
A3_T_sig_cox_dat

A4_N_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~A4_N,
                              data=tcga_cox_dat))
A4_N_sig_cox_dat <- data.frame(Names=rownames(A4_N_sig_cox[[8]]),
                               HR = round(A4_N_sig_cox[[7]][,2],3),
                               lower.95 = round(A4_N_sig_cox[[8]][,3],3),
                               upper.95 = round(A4_N_sig_cox[[8]][,4],3),
                               p.value=round(A4_N_sig_cox[[7]][,5],3))
A4_N_sig_cox_dat

A5_M_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~A5_M,
                              data=tcga_cox_dat))
A5_M_sig_cox_dat <- data.frame(Names=rownames(A5_M_sig_cox[[8]]),
                               HR = round(A5_M_sig_cox[[7]][,2],3),
                               lower.95 = round(A5_M_sig_cox[[8]][,3],3),
                               upper.95 = round(A5_M_sig_cox[[8]][,4],3),
                               p.value=round(A5_M_sig_cox[[7]][,5],3))
A5_M_sig_cox_dat

Smoking_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Smoking,
                                 data=tcga_cox_dat))
Smoking_sig_cox_dat <- data.frame(Names=rownames(Smoking_sig_cox[[8]]),
                                  HR = round(Smoking_sig_cox[[7]][,2],3),
                                  lower.95 = round(Smoking_sig_cox[[8]][,3],3),
                                  upper.95 = round(Smoking_sig_cox[[8]][,4],3),
                                  p.value=round(Smoking_sig_cox[[7]][,5],3))
Smoking_sig_cox_dat

A6_Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~A6_Stage,
                                  data=tcga_cox_dat))
A6_Stage_sig_cox_dat <- data.frame(Names=rownames(A6_Stage_sig_cox[[8]]),
                                   HR = round(A6_Stage_sig_cox[[7]][,2],3),
                                   lower.95 = round(A6_Stage_sig_cox[[8]][,3],3),
                                   upper.95 = round(A6_Stage_sig_cox[[8]][,4],3),
                                   p.value=round(A6_Stage_sig_cox[[7]][,5],3))
A6_Stage_sig_cox_dat


RiskType_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~RiskType,
                                data=tcga_cox_dat))
RiskType_sig_cox_dat <- data.frame(Names=rownames(RiskType_sig_cox[[8]]),
                                 HR = round(RiskType_sig_cox[[7]][,2],3),
                                 lower.95 = round(RiskType_sig_cox[[8]][,3],3),
                                 upper.95 = round(RiskType_sig_cox[[8]][,4],3),
                                 p.value=round(RiskType_sig_cox[[7]][,5],3))
RiskType_sig_cox_dat

sig_cox_dat <- rbind(age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     A3_T_sig_cox_dat,
                     A4_N_sig_cox_dat,
                     A5_M_sig_cox_dat,
                     Smoking_sig_cox_dat,
                     A6_Stage_sig_cox_dat,
                     RiskType_sig_cox_dat)
class(sig_cox_dat)
rownames(sig_cox_dat) <- c('Age (>65 vs ≤65)',
                           'Gender (Male vs Female)',
                           'T.Stage (T3/T4 vs T1/T2',
                           'N.Stage (N1/N2/N3 vs N0',
                           'M.Stage (M1 vs M0)',
                           'Smoking (2/3/4 vs 1)',
                           'Stage (III/IV vs I/II)',
                           'RiskType (High vs Low)')
sig_cox_dat$Names <- rownames(sig_cox_dat)
sig_cox_dat <- sig_cox_dat[, -1]
colnames(sig_cox_dat)
colnames(sig_cox_dat) <- c('mean', 'lower', 'upper', 'P.Value')
plot_forest <- function(dat, type = 1) {
  # dat 为数据框：行为基因，列名称为mean、upper、lower、P.Value
  # type 标题类型：1为单因素绘图，2为多因素绘图
  library(forestplot)
  library(dplyr)
  dat <- as.data.frame(dat)
  # dat <- arrange(dat, desc(mean))
  dat$Feature <- rownames(dat)
  dat$HR <- paste0(dat$mean, ' (', dat$lower, '~', dat$upper, ')')
  # dat$P.Value[dat$P.Value < 0.01] <- '<0.01'
  lab_text <- dat[, c("Feature", "HR", "P.Value")]
  print(lab_text)
  
  hrzl <- list('1' = gpar(lwd=2, lty = 1, col = "#000044"),
               '2' = gpar(lwd=1, lty = 1, col = "#000044"),
               '3' = gpar(lwd=2, col = "#000044"))
  names(hrzl)=c(1, 2, nrow(lab_text)+2)
  
  if (type) {
    ti <- 'Univariate Analysis'
  } else {
    ti <- 'Multivariate Analysis'
  }
  
  forestplot(rbind(colnames(lab_text),
                   lab_text),  #显示的文本
             mean = c(NA, dat$mean), #误差条的均值(此处为差值的中值)
             lower = c(NA, dat$lower), #误差条的下界(此处为差值的25%分位数)
             upper = c(NA, dat$upper), #误差条的上界(此处为差值的75%分位数)
             hrzl_lines = hrzl,
             title= paste0("Hazard Ratio(", ti, ")"),
             col=fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray"),
             graph.pos=2,
             line.margin = c(0.2, 0.3, 0.2, 0.3),
             # clip=c(0.1,2.5), 
             zero = 1, #显示y=0的垂直线
             xlog=F, #x轴的坐标不取对数
             fn.ci_norm = fpDrawCircleCI, #误差条显示方式
             boxsize = 0.15, ##误差条中的圆心点大小
             lty.ci = 6,   # 误差条的线的线型
             lwd.ci = 1.5,   # 误差条的线的宽度
             ci.vertices.height = 0.1, # # 误差条末端的长度
             txt_gp = fpTxtGp(ticks = gpar(cex = 0.75), xlab = gpar(cex = 1), cex = 1), #文本大小设置
             lineheight = "auto", #线的高度
  )
}

pdf('PDFs/sig_cox_dat.pdf', width = 8, height = 5, onefile = FALSE)
plot_forest(sig_cox_dat)
dev.off()


muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age1+Gender+A3_T+A4_N+A5_M+Smoking+A6_Stage+RiskType, data=tcga_cox_dat))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
rownames(muti_cox_dat) <- c('Age (>65 vs ≤65)',
                           'Gender (Male vs Female)',
                           'T.Stage (T3/T4 vs T1/T2',
                           'N.Stage (N1/N2/N3 vs N0',
                           'M.Stage (M1 vs M0)',
                           'Smoking (2/3/4 vs 1)',
                           'Stage (III/IV vs I/II)',
                           'RiskType (High vs Low)')
muti_cox_dat$Names <- rownames(muti_cox_dat)
muti_cox_dat <- muti_cox_dat[, -1]
colnames(muti_cox_dat)
colnames(muti_cox_dat) <- c('mean', 'lower', 'upper', 'P.Value')

pdf('PDFs/muti_cox_dat.pdf', width = 8, height = 5, onefile = FALSE)
plot_forest(muti_cox_dat)
dev.off()

# 列线图与列线图校正图 ##########
library(rms)
ddist <- datadist(tcga_cox_dat)
options(datadist='ddist')
dat.cph <- cph(formula=Surv(OS.time, OS)~A3_T+A4_N+RiskScore,
               data=tcga_cox_dat,surv=T)
surv=Survival(dat.cph)
surv1 <- function(x)surv(365*1,lp=x) 
surv2 <- function(x)surv(365*3,lp=x) 
surv3 <- function(x)surv(365*5,lp=x)
dat.nomogram <- nomogram(dat.cph,lp= F,
                         fun=list(surv1,surv2,surv3),
                         funlabel=c('1 Year Surv',
                                    '3 Year Surv',
                                    '5 Year Surv'),
                         maxscale=100,
                         fun.at=c('0.9','0.8','0.7','0.5','0.3','0.1'))
pdf('nomogram.pdf', width = 10, height = 6)
plot(dat.nomogram,
     xfrac = .25, varname.label = T, 
     varname.label.sep = '=', ia.space = .2,
     tck=NA, tcl=-0.20, lmgp = 0.3,
     points.label = 'Points', total.points.label = 'Total Points',
     total.sep.page = F,
     cap.labels = F, cex.var = 1.1, cex.axis = 1.05, lwd = 1,
     label.every = 1, col.grid = gray(c(0.8, 0.95)))
dev.off()



f1<-cph(formula = as.formula(paste0("Surv(OS.time, OS) ~",
                                    paste0(c('A3_T','A4_N','RiskScore'),
                                           collapse = '+')))
        ,data=tcga_cox_dat,x=T,y=T,surv = T,
        na.action=na.delete,time.inc = 365*1) 
cal1<-calibrate(f1, cmethod="KM", method="boot",u=365*1,m=120,B=500)
f2<-cph(formula = as.formula(paste0("Surv(OS.time, OS) ~",
                                    paste0(c('A3_T','A4_N','RiskScore'),
                                           collapse = '+')))
        ,data=tcga_cox_dat,x=T,y=T,surv = T,
        na.action=na.delete,time.inc = 365*3) 
cal2<-calibrate(f2, cmethod="KM", method="boot",u=365*3,m=120,B=500)
f3<-cph(formula = as.formula(paste0("Surv(OS.time, OS) ~",
                                    paste0(c('A3_T','A4_N','RiskScore'),
                                           collapse = '+')))
        ,data=tcga_cox_dat,x=T,y=T,surv = T,
        na.action=na.delete,time.inc = 365*5) 
cal3<-calibrate(f3, cmethod="KM", method="boot",u=365*5,m=120,B=500)

pdf('nomogram_cal.pdf', width = 5, height = 5)
plot(cal1,
     lwd = 1,
     # lty = 0,
     errbar.col = color8[1],
     xlim = c(0,1),
     ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",
     ylab = "Observed OS (%)",
     col = color8[1],
     cex.lab=1.2,
     cex.axis=1, 
     cex.main=1.2, 
     cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = color8[1], pch = 16)
mtext("")

plot(cal2,lwd = 1,#lty = 0,
     errbar.col = color8[2],
     xlim = c(0,1),ylim= c(0,1),col = color8[2],add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = color8[2], pch = 16)

plot(cal3,lwd = 1,#lty = 0,
     errbar.col = color8[3],
     xlim = c(0,1),ylim= c(0,1),col = color8[3],add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = color8[3], pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", 
       legend = c("1 Year","3 Year", "5 Year"),
       col =color8[1:3], 
       lwd = 2,
       cex = 1.2,
       bty = "n")
dev.off()


# DCA 图 ################
tcga_roc_dca_dat <- tcga_cli_RS
tcga_roc_dca_dat$A3_T <- as.numeric(factor(tcga_roc_dca_dat$A3_T))
tcga_roc_dca_dat$A4_N <- as.numeric(factor(tcga_roc_dca_dat$A4_N))
tcga_roc_dca_dat$A5_M <- as.numeric(factor(tcga_roc_dca_dat$A5_M))
tcga_roc_dca_dat$A6_Stage <- as.numeric(factor(tcga_roc_dca_dat$A6_Stage))

library(rmda)
set.seed(1984)
T.Stage.model <- decision_curve(OS~A3_T,
                                data = tcga_roc_dca_dat, 
                                bootstraps = 500)
set.seed(1984)
N.Stage.model <- decision_curve(OS~A4_N,
                                data = tcga_roc_dca_dat,
                                bootstraps = 500)
                              bootstraps = 500)
set.seed(1984)
Riskscore.model <- decision_curve(OS~RiskScore,
                                  data = tcga_roc_dca_dat, 
                                  bootstraps = 500)
set.seed(1984)
nomo.model <- decision_curve(OS~A3_T+A4_N+RiskScore,
                             data = tcga_roc_dca_dat,
                             bootstraps = 500)
pdf('PDFs/tcga_cli_DCA.pdf', width = 5, height = 5)
plot_decision_curve( list(T.Stage.model, 
                          N.Stage.model,
                          Riskscore.model, 
                          nomo.model), 
                     curve.names = c("T.Stage model", 
                                     "N.Stage model", 
                                     "RiskScore model", 
                                     'Nomogram model'), 
                     xlim = c(0, 1), 
                     legend.position = "topright",
                     confidence.intervals = FALSE,
                     col = color9)
dev.off() 

# 8.2、高低风险分组的免疫评分与免疫治疗 ########
dim(tcga_exp_log2)
dim(tcga_cluster_cli)
dim(tcga_model_dat)

tcga_estimate_bar <- ggplotViolin_group(tcga_estimate[tcga_cli_RS$A0_Samples, ], 
                                        tcga_cli_RS$RiskType, 
                                        cols = c('#00468BFF', '#ED0000FF'), 
                                        xlab = '', 
                                        ylab = 'Score', 
                                        title = 'RiskType',
                                        angle = 60)
tcga_estimate_bar
ggsave(plot = tcga_estimate_bar,
       filename = 'PDFs/tcga_estimate_bar.pdf',
       width = 4, height = 5)


tcga_MCPcounter_score <- MCPcounter_score(log2(luad_tpm[, tcga_cli_RS$A0_Samples] + 1))
tcga_MCPcounter_bar <- ggplotViolin_group(tcga_MCPcounter_score[tcga_cli_RS$A0_Samples, ], 
                                          tcga_cli_RS$RiskType, 
                                          cols = c('#00468BFF', '#ED0000FF'), 
                                          xlab = '', 
                                          ylab = 'Score', 
                                          title = 'RiskType',
                                          angle = 60)
tcga_MCPcounter_bar
ggsave(plot = tcga_MCPcounter_bar,
       filename = 'PDFs/tcga_MCPcounter_bar.pdf',
       width = 8, height = 5)


tcga_ssGSEA_score <- ssGSEA_score(log2(luad_tpm[, tcga_cli_RS$A0_Samples] + 1))
tcga_ssGSEA_bar <- ggplotViolin_group(tcga_ssGSEA_score[tcga_cli_RS$A0_Samples, ], 
                                      tcga_cli_RS$RiskType, 
                                      cols = c('#00468BFF', '#ED0000FF'), 
                                      xlab = '', 
                                      ylab = 'Score', 
                                      title = 'RiskType',
                                      angle = 60)
tcga_ssGSEA_bar
ggsave(plot = tcga_ssGSEA_bar,
       filename = 'PDFs/tcga_ssGSEA_bar.pdf',
       width = 12, height = 5)


tcga_Timer_score <- Timer_score(log2(luad_tpm[, tcga_cli_RS$A0_Samples] + 1), code = 'LUAD')
tcga_Timer_bar <- ggplotViolin_group(tcga_Timer_score[tcga_cli_RS$A0_Samples, ], 
                                     tcga_cli_RS$RiskType, 
                                     cols = c('#00468BFF', '#ED0000FF'), 
                                     xlab = '', 
                                     ylab = 'Score', 
                                     title = 'RiskType',
                                     angle = 60)
tcga_Timer_bar
ggsave(plot = tcga_Timer_bar,
       filename = 'PDFs/tcga_Timer_bar.pdf',
       width = 5, height = 5)

tcga_cibersort_score <- cibersort_score(log2(luad_tpm[, tcga_cli_RS$A0_Samples] + 1))
tcga_cibersort_score1 <- tcga_cibersort_score[, 1:22]
tcga_cibersort_bar <- ggplotViolin_group(tcga_cibersort_score1[tcga_cli_RS$A0_Samples, ], 
                                         tcga_cli_RS$RiskType, 
                                         cols = c('#00468BFF', '#ED0000FF'), 
                                         xlab = '', 
                                         ylab = 'Score', 
                                         title = 'RiskType',
                                         angle = 60)
tcga_cibersort_bar
ggsave(plot = tcga_cibersort_bar,
       filename = 'PDFs/tcga_cibersort_bar.pdf',
       width = 10, height = 5)



intersect(icg_genes, rownames(tcga_exp_log2))
icg_genes <- read.delim('e:/public/icg_genes_pmid31043417.txt', header = F)[, 1]

tcga_icg <- as.data.frame(t(log2(luad_tpm[icg_genes, tcga_cli_RS$A0_Samples] + 1)))

tcga_icg_bar <- ggplotViolin_group(tcga_icg[tcga_cli_RS$A0_Samples, c("PDCD1", "PDCD1LG2", "CTLA4", 
                                                                   "CD27", "CD274", "TIGIT")], 
                                   tcga_cli_RS$RiskType, 
                                   cols = c('#00468BFF', '#ED0000FF'), 
                                   xlab = '', 
                                   ylab = 'Score', 
                                   title = 'RiskType',
                                   angle = 60)
tcga_icg_bar
ggsave(plot = tcga_icg_bar,
       filename = 'PDFs/tcga_icg_bar.pdf',
       width = 6, height = 5)


# 免疫评分热图 ####################
library(pheatmap)
library(dplyr)

intersect(colnames(tcga_Timer_score), colnames(tcga_cibersort_score1))

colnames(tcga_MCPcounter_score) <- gsub("Neutrophils", "Neutrophils ", colnames(tcga_MCPcounter_score))
colnames(tcga_Timer_score) <- gsub("Neutrophil", "Neutrophil ", colnames(tcga_Timer_score))
colnames(tcga_Timer_score) <- gsub("Macrophage", "Macrophage ", colnames(tcga_Timer_score))

Immune_score_dat <- cbind(tcga_estimate[tcga_cli_RS$A0_Samples, ],
                          tcga_MCPcounter_score[tcga_cli_RS$A0_Samples, ],
                          tcga_Timer_score[tcga_cli_RS$A0_Samples, ],
                          tcga_cibersort_score1[tcga_cli_RS$A0_Samples, ],
                          tcga_ssGSEA_score[tcga_cli_RS$A0_Samples, ])
Immune_score_group <- data.frame(Cells = c(colnames(tcga_estimate),
                                           colnames(tcga_MCPcounter_score),
                                           colnames(tcga_Timer_score),
                                           colnames(tcga_cibersort_score1),
                                           colnames(tcga_ssGSEA_score)),
                                 Methods = c(rep('ESTIMATE', 3),
                                             rep('MCPcounter', 10),
                                             rep('TIMER', 6),
                                             rep('CIBERSORT', 22),
                                             rep('ssGSEA', 28)))
table(Immune_score_group$Methods)



tcga_cli_RS <- arrange(tcga_cli_RS, RiskScore)
table(tcga_cli_RS$RiskType)
annotation_col  <- data.frame(RiskType = factor(tcga_cli_RS$RiskType),
                              RiskScore = as.numeric(tcga_cli_RS$RiskScore))
rownames(annotation_col) <- tcga_cli_RS$A0_Samples


annotation_row  <- data.frame(Methods = factor(Immune_score_group$Methods))
rownames(annotation_row) <- Immune_score_group$Cells


tcga_Immunedat_heatmap <- Immune_score_dat[tcga_cli_RS$A0_Samples,
                                           Immune_score_group$Cells]
bk=unique(c(seq(-1.2, 1.2, length=100)))
ann_colors = list(
  ImmuneGroups = c(high = '#00468BFF', low = '#ED0000FF')
)


# 免疫治疗数据集的验证 ##########
IMvigor210_cli <- read.delim('e:/datas/IMvigor210/IMvigor210_cli.txt', header = T, row.names = 1)

IMvigor210_cli1 <- IMvigor210_cli[, c("os", "censOS", "Best.Confirmed.Overall.Response", 
                                      "binaryResponse", "IC.Level", "TC.Level", "Immune.phenotype",
                                      "FMOne.mutation.burden.per.MB", "Neoantigen.burden.per.MB")]
colnames(IMvigor210_cli1) <-  c("OS.time", "OS", "Response1", "Response", "IC.Level", "TC.Level", "Immune.phenotype",
                                "TMB", "NEO")
IMvigor210_cli1 <- IMvigor210_cli1[!is.na(IMvigor210_cli1$Response), ]


IMvigor210_exp <- read.csv('e:/datas/IMvigor210/IMvigor210_Counts2TPM_symbol.csv', header = T, row.names = 1)
boxplot(IMvigor210_exp[, 1:5])



intersect(rownames(IMvigor210_exp), genes)

IMvigor210_model <- cbind(IMvigor210_cli1[, c("OS.time", "OS")],
                          t(IMvigor210_exp[genes, rownames(IMvigor210_cli1)]))

library(survcomp)
cox.IMvigor210 <- coxph(as.formula(paste0("Surv(OS.time, OS) ~" ,paste0(genes,collapse = '+'))),
                        data =as.data.frame(IMvigor210_model))
lan.IMvigor210 <- coef(cox.IMvigor210)
round(lan.IMvigor210, 3)

library(ggthemes)
risk.IMvigor210 <- as.numeric(lan.IMvigor210%*%as.matrix(t(IMvigor210_model[,genes])))
IMvigor210.dat.plot <- IMvigor210_model[, c("OS.time", "OS", genes)]
IMvigor210.dat.plot$RiskScore <- risk.IMvigor210


library(survcomp)
library(cutoff)
library(survminer)
IMvigor210.dat.plot$RiskType <- ifelse(IMvigor210.dat.plot$RiskScore > median(IMvigor210.dat.plot$RiskScore), 'High', 'Low')


IMvigor210.km <- ggsurvplotKM(IMvigor210.dat.plot[, c("OS.time", "OS", 'RiskType')],
                              title = 'RiskType',
                              lables = c('High', 'Low'),
                              risk.table = T)
IMvigor210.km
ggsave(plot = IMvigor210.km,
       filename = 'PDFs/IMvigor210.km.pdf',
       width = 5, height = 5)

IMvigor210_cli1 <- cbind(IMvigor210_cli1,
                         IMvigor210.dat.plot[rownames(IMvigor210_cli1), 3:9])

dim(IMvigor210_cli1)

Response_Violin <- ggplotViolin_muti_group(IMvigor210_cli1[, c("Response1", "RiskScore")],
                                           cols = color8,
                                           title = 'Response',
                                           ylab = 'RiskScore')
Response_Violin


IC.Level_Violin_dat <- IMvigor210_cli1[, c("IC.Level", "RiskScore")]
IC.Level_Violin_dat <- na.omit(IC.Level_Violin_dat)
IC.Level_Violin <- ggplotViolin_muti_group(IC.Level_Violin_dat,
                                           cols = color8,
                                           title = 'IC.Level',
                                           ylab = 'RiskScore')
IC.Level_Violin


TC.Level_Violin_dat <- IMvigor210_cli1[, c("TC.Level", "RiskScore")]
TC.Level_Violin_dat <- na.omit(TC.Level_Violin_dat)
TC.Level_Violin <- ggplotViolin_muti_group(TC.Level_Violin_dat,
                                           cols = color8,
                                           title = 'TC.Level',
                                           ylab = 'RiskScore')
TC.Level_Violin


Immune.phenotype_Violin_dat <- IMvigor210_cli1[, c("Immune.phenotype", "RiskScore")]
Immune.phenotype_Violin_dat <- na.omit(Immune.phenotype_Violin_dat)
Immune.phenotype_Violin <- ggplotViolin_muti_group(Immune.phenotype_Violin_dat,
                                                   cols = color8,
                                                   title = 'Immune.phenotype',
                                                   ylab = 'RiskScore')
Immune.phenotype_Violin


RS_Violin <- cowplot::plot_grid(Response_Violin,
                                IC.Level_Violin,
                                TC.Level_Violin,
                                Immune.phenotype_Violin,
                                align = 'hv',
                                ncol = 2)
RS_Violin
ggsave(plot = RS_Violin,
       filename = 'PDFs/IMvigor210_RS_Violin.pdf',
       width = 10, height = 10)


ggplot_bar <- function(data, width = 0.6, cols = color8,
                       xlab = '', ylab = '', name = 'Groups') {
  # data：数据框，行和列为不同的分类特征，数据框的数据是特征的统计量
  for (i in 1:ncol(data)) {
    data[, i] = data[, i] / sum(data[, i])
  }
  data=reshape2::melt(data)
  colnames(data)=c('Type','variable','Perc')
  data[, 1] = as.character(data[, 1])
  data[, 2] = as.character(data[, 2])
  bar = ggplot(data, aes(x=variable, y=Perc, fill=Type))+geom_bar(stat = "identity", width=width) + theme_get()
  bar = bar+geom_text(data=data,aes(label=sprintf("%0.2f", round(Perc, digits = 2))),position=position_stack(vjust=0.5)) 
  bar = bar+labs(x=xlab, y=ylab)+scale_fill_manual(name = name, values = cols)+theme(legend.position = "top")# + scale_fill_discrete(name = name)
  return(bar)
}

bar_plot <- function(data, width = 0.6, cols = color9,
                     xlab = '', ylab = '', name = 'Groups',
                     pvalue = T) {
  oridat <- data
  for (i in 1:ncol(data)) {
    data[, i] = data[, i] / sum(data[, i])
  }
  data=reshape2::melt(data)
  colnames(data)=c('Type','variable','Perc')
  data[, 1] = as.character(data[, 1])
  data[, 2] = as.character(data[, 2])
  bar = ggplot(data, aes(x=variable, y=Perc, fill=Type))+geom_bar(stat = "identity", width=width) + theme_get()
  bar = bar+geom_text(data=data,aes(label=sprintf("%0.2f", round(Perc, digits = 2))),position=position_stack(vjust=0.5))
  if (pvalue) {
    pvalue <- round(chisq.p.value(oridat), 4)
    bar = bar+labs(x=paste0('Chi-squared test P=', pvalue), y=ylab)+scale_fill_manual(name = name, values = cols)+theme(legend.position = "top")# + scale_fill_discrete(name = name)
  } else {
    bar = bar+labs(x=xlab, y=ylab)+scale_fill_manual(name = name, values = cols)+theme(legend.position = "top")# + scale_fill_discrete(name = name)
  }
  
  return(bar)
}

bar_dat <- table(IMvigor210_cli1$Response1, IMvigor210_cli1$RiskType)
pdf('PDFs/IMvigor210_bar_plot.pdf', width = 5, height = 5)
bar_plot(bar_dat, cols = color8, name = 'Response')
dev.off()


IMvigor210.mcpcounter <- MCPcounter_score(IMvigor210_exp)


library(Hmisc)
library(corrplot)
tcga_mcpcounter_Cor_dat <- IMvigor210.mcpcounter[rownames(IMvigor210_cli1), ]
tcga_mcpcounter_Cor_dat$RiskScore <- IMvigor210_cli1$RiskScore

mcpcounter.Cor <- rcorr(as.matrix(tcga_mcpcounter_Cor_dat))
col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
# ?corrplot
pdf('IMvigor210/mcpcounter.Cor.pdf', width = 10, height = 10)
corrplot(mcpcounter.Cor$r, method="color", col=col(200),   
         addCoef.col = "black",
         tl.col="black", tl.srt=45,
         diag=T,
         pch = 1,
         pch.cex = 1)
dev.off()




IMvigor210_mcpcounter_Cor_dat <- cbind(IMvigor210.mcpcounter[rownames(IMvigor210_cli1),],
                                      IMvigor210_cli1[, c("NEO","TMB","RiskScore")])
IMvigor210_mcpcounter_cor <- rcorr(as.matrix(IMvigor210_mcpcounter_Cor_dat))
IMvigor210_mcpcounter_cor <- IMvigor210_mcpcounter_cor$r


library(circlize)
library(ComplexHeatmap)
grid_col <-  c(S1 = "red", S2 = "green", S3 = "blue",
               E1 = "grey", E2 = "grey", E3 = "grey", 
               E4 = "grey", E5 = "grey", E6 = "grey")
IMvigor210_mcpcounter_cor[IMvigor210_mcpcounter_cor == 1] <- 0

IMvigor210_mcpcounter_cor1 <- IMvigor210_mcpcounter_cor[c("RiskScore", "NEO", "TMB", "CD8 T cells", "T cells", "NK cells"), 
                                                        c("RiskScore", "NEO", "TMB", "CD8 T cells", "T cells", "NK cells")]
max(IMvigor210_mcpcounter_cor1)
min(IMvigor210_mcpcounter_cor1)
pdf('PDFs/IMvigor210_mcpcounter.Cor.pdf', width = 5, height = 5)
col_fun = colorRamp2(c(0.85, 0, -0.85), c("red", 'white', "blue"))

lgd_links = Legend(at = c(-0.85, 0, 0.85), col_fun = col_fun, 
                   title_position = "topleft", title = "Cor")

lgd_list_vertical = packLegend(lgd_links)
lgd_list_vertical

chordDiagram(IMvigor210_mcpcounter_cor1, grid.col = grid_col, 
             annotationTrack = c("name", "grid", "axis"),
             col = col_fun,
             annotationTrackHeight = c(0.05, 0.1)) 
draw(lgd_list_vertical, x = unit(20, "mm"), y = unit(20, "mm"))
circos.clear()
dev.off()


# ROC 曲线
IMvigor210_roc_dat <- IMvigor210_cli1[, c("OS.time", "OS", "TMB", "NEO", "RiskScore")]
IMvigor210_roc_dat <- na.omit(IMvigor210_roc_dat)
fmla1 <- as.formula(paste0("Surv(OS.time, OS) ~" , 'TMB+NEO+RiskScore'))
coxsig1 <- coxph(fmla1, data =IMvigor210_roc_dat)
IMvigor210_roc_dat$Combine <- predict(coxsig1)

library(timeROC)
library(survival)
library(RColorBrewer)
display.brewer.all(type = "seq")
display.brewer.all(type = "qual")
display.brewer.all(type = "div")
mypal <- brewer.pal(8, 'Dark2')
dat <- IMvigor210_roc_dat

auc <- c()
for (i in seq(0.1, 1.9, 0.01)) {
  print(i)
  fit.rs <- timeROC(T=dat$OS.time,
                    delta=dat$OS,
                    marker=as.numeric(dat$RiskScore),
                    cause=1,weighting="marginal",
                    times=c(365*i),
                    iid=TRUE)
  print(fit.rs$AUC[2])
  auc <- c(auc, fit.rs$AUC[2])
}
max(auc)

fit.rs <- timeROC(T=dat$OS.time,
                  delta=dat$OS,
                  marker=as.numeric(dat$RiskScore),
                  cause=1,weighting="marginal",
                  times=c(365*0.3),
                  iid=TRUE)
fit.rs$AUC[2]


datNEO=data.frame(dat$OS.time,dat$OS,dat$NEO)
colnames(datNEO)=c('OS.time','OS','NEO')
fmlaNEO <- as.formula(paste0("Surv(OS.time, OS) ~NEO"))
coxNEO <- coxph(fmlaNEO, data =datNEO)
NEO.score=predict(coxNEO,datNEO)

fit.NEO <- timeROC(T=dat$OS.time,
                   delta=dat$OS,
                   marker=as.numeric(NEO.score),
                   cause=1,weighting="marginal",
                   times=c(365*1),
                   iid=TRUE)
fit.NEO$AUC[2]

datTMB=data.frame(dat$OS.time,dat$OS,dat$TMB)
colnames(datTMB)=c('OS.time','OS','TMB')
fmlaTMB <- as.formula(paste0("Surv(OS.time, OS) ~TMB"))
coxTMB <- coxph(fmlaTMB, data =datTMB)
TMB.score=predict(coxTMB,datTMB)
fit.TMB <- timeROC(T=dat$OS.time,
                   delta=dat$OS,
                   marker=as.numeric(TMB.score),
                   cause=1,weighting="marginal",
                   times=c(365*1),
                   iid=TRUE)
fit.TMB$AUC[2]


library(timeROC)
auc <- c()
for (i in seq(0.1, 1.9, 0.01)) {
  print(i)
  fit.Complex <- timeROC(T=dat$OS.time,
                         delta=dat$OS,
                         marker=as.numeric(dat$Combine),
                         cause=1,weighting="marginal",
                         times=c(365*i),
                         iid=TRUE)
  print(fit.Complex$AUC[2])
  auc <- c(auc, fit.Complex$AUC[2])
}
max(auc)

fit.Complex <- timeROC(T=dat$OS.time,
                       delta=dat$OS,
                       marker=as.numeric(dat$Combine),
                       cause=1,weighting="marginal",
                       times=c(365*1.26),
                       iid=TRUE)
fit.Complex$AUC[2]

pdf('PDFs/IMvigor210_ROC.pdf',width = 6,height = 6)
plot(c(0,1),c(0,1),type='l',xlab='False positive rate',ylab='Ture positive rate',xaxs = "i",yaxs = "i")
los=lowess(fit.rs$FP[,2], y=fit.rs$TP[,2], f = 1/3, iter = 100)
los$x=c(0,los$x,1)
los$y=c(0,los$y,1)
lines(los,col=mypal[1],lwd = 2)
los=lowess(fit.NEO$FP[,2], y=fit.NEO$TP[,2], f = 1/3, iter = 100)
los$x=c(0,los$x,1)
los$y=c(0,los$y,1)
lines(los,col=mypal[2],lwd = 2)
los=lowess(fit.TMB$FP[,2], y=fit.TMB$TP[,2], f = 2/3, iter = 5)
los$x=c(0,los$x,1)
los$y=c(0,los$y,1)
lines(los,col=mypal[3],lwd = 2)
los=lowess(fit.Complex$FP[,2], y=fit.Complex$TP[,2], f = 2/3, iter = 5)
los$x=c(0,los$x,1)
los$y=c(0,los$y,1)
lines(los,col=mypal[4],lwd = 2)
legend('bottomright',
       c(paste0('Riskscore, AUC=',round(fit.rs$AUC[2],2),' 95%CI(',round(confint(fit.rs,level = 0.9)$CI_AUC[,1]/100,2),'-',round(confint(fit.rs,level = 0.9)$CI_AUC[,2]/100,2),')')
         ,paste0('NEO, AUC=',round(fit.NEO$AUC[2],2),' 95%CI(',round(confint(fit.NEO,level = 0.9)$CI_AUC[,1]/100,2),'-',round(confint(fit.NEO,level = 0.9)$CI_AUC[,2]/100,2),')')
         ,paste0('TMB, AUC=',round(fit.TMB$AUC[2],2),' 95%CI(',round(confint(fit.TMB,level = 0.9)$CI_AUC[,1]/100,2),'-',round(confint(fit.TMB,level = 0.9)$CI_AUC[,2]/100,2),')')
         ,paste0('Combine, AUC=',round(fit.Complex$AUC[2],2),' 95%CI(',round(confint(fit.Complex,level = 0.9)$CI_AUC[,1]/100,2),'-',round(confint(fit.Complex,level = 0.9)$CI_AUC[,2]/100,2),')'))
       ,col = mypal
       ,lty = rep(1,length(fit.rs$times)),lwd=rep(1,length(fit.rs$times))
       ,merge = TRUE,cex = 0.8
       ,title = paste0(paste0(rep(' ',ceiling(strwidth(paste0(fit.rs$times[1]))/strwidth(' '))),collapse = ''),'AUC 95%CI'))

dev.off()

# 保存数据 ########
save.image('20210328_Hypoxia_LUAD_immune_model.Rdata')


# ROC 曲线比较 ########################
library(survcomp)
library(timeROC)
library(survival)

luad_tpm_temp <- t(log2(luad_tpm[, rownames(tcga_model_dat)] + 1))

luad_other_model_genes <- read.delim('raw_dat/luad_models_genes.txt')


# Sun model #################
Sun <- luad_other_model_genes$Sun
Sun
Sun <- intersect(Sun, colnames(luad_tpm_temp))
Sun

Sun.dat <- data.frame(OS.time=tcga_model_dat$OS.time,
                      OS=tcga_model_dat$OS,
                      luad_tpm_temp[, Sun])
fmla1 <- as.formula(paste0("Surv(OS.time, OS) ~"
                           ,paste0(Sun, collapse = '+')))
Sun.cox <- coxph(fmla1, data =as.data.frame(Sun.dat))
Sun_cindex <- concordance.index(predict(Sun.cox),
                                surv.time = Sun.dat$OS.time / 365, 
                                surv.event = Sun.dat$OS,
                                method = "noether")
Sun_cindex$c.index
Sun_cindex$lower
Sun_cindex$upper
Sun_risk <- as.numeric(coef(Sun.cox)%*%as.matrix(t(Sun.dat[,Sun])))
Sun.dat$RiskScore <- Sun_risk
Sun_KM <- ggsurvplotKM(data.frame(Sun.dat$OS.time,
                                  Sun.dat$OS,
                                  ifelse(Sun_risk> median(Sun_risk),'High','Low')),
                       title = 'Sun',
                       lables = c('High', 'Low'),
                       risk.table = T,
                       col = color9)
Sun_KM
Sun_ROC <- ggplotROC(Sun.dat[, c("OS.time", "OS", "RiskScore")])
Sun_ROC

# Xue model #################
Xue <- luad_other_model_genes$Xue
Xue
Xue <- intersect(Xue, colnames(luad_tpm_temp))
Xue

Xue.dat <- data.frame(OS.time=tcga_model_dat$OS.time,
                      OS=tcga_model_dat$OS,
                      luad_tpm_temp[, Xue])
fmla1 <- as.formula(paste0("Surv(OS.time, OS) ~"
                           ,paste0(Xue, collapse = '+')))
Xue.cox <- coxph(fmla1, data =as.data.frame(Xue.dat))
Xue_cindex <- concordance.index(predict(Xue.cox),
                                surv.time = Xue.dat$OS.time / 365, 
                                surv.event = Xue.dat$OS,
                                method = "noether")
Xue_cindex$c.index
Xue_cindex$lower
Xue_cindex$upper
Xue_risk <- as.numeric(coef(Xue.cox)%*%as.matrix(t(Xue.dat[,Xue])))
Xue.dat$RiskScore <- Xue_risk
Xue_KM <- ggsurvplotKM(data.frame(Xue.dat$OS.time,
                                  Xue.dat$OS,
                                  ifelse(Xue_risk> median(Xue_risk),'High','Low')),
                       title = 'Xue',
                       lables = c('High', 'Low'),
                       risk.table = T,
                       col = color9)
Xue_KM
Xue_ROC <- ggplotROC(Xue.dat[, c("OS.time", "OS", "RiskScore")])
Xue_ROC

# Yu model #################
Yu <- luad_other_model_genes$Yu
Yu
Yu <- intersect(Yu, colnames(luad_tpm_temp))
Yu

Yu.dat <- data.frame(OS.time=tcga_model_dat$OS.time,
                      OS=tcga_model_dat$OS,
                      luad_tpm_temp[, Yu])
fmla1 <- as.formula(paste0("Surv(OS.time, OS) ~"
                           ,paste0(Yu, collapse = '+')))
Yu.cox <- coxph(fmla1, data =as.data.frame(Yu.dat))
Yu_cindex <- concordance.index(predict(Yu.cox),
                                surv.time = Yu.dat$OS.time / 365, 
                                surv.event = Yu.dat$OS,
                                method = "noether")
Yu_cindex$c.index
Yu_cindex$lower
Yu_cindex$upper
Yu_risk <- as.numeric(coef(Yu.cox)%*%as.matrix(t(Yu.dat[,Yu])))
Yu.dat$RiskScore <- Yu_risk
Yu_KM <- ggsurvplotKM(data.frame(Yu.dat$OS.time,
                                  Yu.dat$OS,
                                  ifelse(Yu_risk> median(Yu_risk),'High','Low')),
                       title = 'Yu',
                       lables = c('High', 'Low'),
                       risk.table = T,
                       col = color9)
Yu_KM
Yu_ROC <- ggplotROC(Yu.dat[, c("OS.time", "OS", "RiskScore")])
Yu_ROC

# Yue model #################
Yue <- luad_other_model_genes$Yue
Yue
Yue <- intersect(Yue, colnames(luad_tpm_temp))
Yue

Yue.dat <- data.frame(OS.time=tcga_model_dat$OS.time,
                      OS=tcga_model_dat$OS,
                      luad_tpm_temp[, Yue])
fmla1 <- as.formula(paste0("Surv(OS.time, OS) ~"
                           ,paste0(Yue, collapse = '+')))
Yue.cox <- coxph(fmla1, data =as.data.frame(Yue.dat))
Yue_cindex <- concordance.index(predict(Yue.cox),
                                surv.time = Yue.dat$OS.time / 365, 
                                surv.event = Yue.dat$OS,
                                method = "noether")
Yue_cindex$c.index
Yue_cindex$lower
Yue_cindex$upper
Yue_risk <- as.numeric(coef(Yue.cox)%*%as.matrix(t(Yue.dat[,Yue])))
Yue.dat$RiskScore <- Yue_risk
Yue_KM <- ggsurvplotKM(data.frame(Yue.dat$OS.time,
                                  Yue.dat$OS,
                                  ifelse(Yue_risk> median(Yue_risk),'High','Low')),
                       title = 'Yue',
                       lables = c('High', 'Low'),
                       risk.table = T,
                       col = color9)
Yue_KM
Yue_ROC <- ggplotROC(Yue.dat[, c("OS.time", "OS", "RiskScore")])
Yue_ROC


# Zou model #################
Zou <- luad_other_model_genes$Zou
Zou
Zou <- intersect(Zou, colnames(luad_tpm_temp))
Zou

Zou.dat <- data.frame(OS.time=tcga_model_dat$OS.time,
                      OS=tcga_model_dat$OS,
                      luad_tpm_temp[, Zou])
fmla1 <- as.formula(paste0("Surv(OS.time, OS) ~"
                           ,paste0(Zou, collapse = '+')))
Zou.cox <- coxph(fmla1, data =as.data.frame(Zou.dat))
Zou_cindex <- concordance.index(predict(Zou.cox),
                                surv.time = Zou.dat$OS.time / 365, 
                                surv.event = Zou.dat$OS,
                                method = "noether")
Zou_cindex$c.index
Zou_cindex$lower
Zou_cindex$upper
Zou_risk <- as.numeric(coef(Zou.cox)%*%as.matrix(t(Zou.dat[,Zou])))
Zou.dat$RiskScore <- Zou_risk
Zou_KM <- ggsurvplotKM(data.frame(Zou.dat$OS.time,
                                  Zou.dat$OS,
                                  ifelse(Zou_risk> median(Zou_risk),'High','Low')),
                       title = 'Zou',
                       lables = c('High', 'Low'),
                       risk.table = T,
                       col = color9)
Zou_KM
Zou_ROC <- ggplotROC(Zou.dat[, c("OS.time", "OS", "RiskScore")])
Zou_ROC



# Li model #################
Li <- luad_other_model_genes$Li
Li
Li <- intersect(Li, colnames(luad_tpm_temp))
Li

Li.dat <- data.frame(OS.time=tcga_model_dat$OS.time,
                      OS=tcga_model_dat$OS,
                      luad_tpm_temp[, Li])
fmla1 <- as.formula(paste0("Surv(OS.time, OS) ~"
                           ,paste0(Li, collapse = '+')))
Li.cox <- coxph(fmla1, data =as.data.frame(Li.dat))
Li_cindex <- concordance.index(predict(Li.cox),
                                surv.time = Li.dat$OS.time / 365, 
                                surv.event = Li.dat$OS,
                                method = "noether")
Li_cindex$c.index
Li_cindex$lower
Li_cindex$upper
Li_risk <- as.numeric(coef(Li.cox)%*%as.matrix(t(Li.dat[,Li])))
Li.dat$RiskScore <- Li_risk
Li_KM <- ggsurvplotKM(data.frame(Li.dat$OS.time,
                                  Li.dat$OS,
                                  ifelse(Li_risk> median(Li_risk),'High','Low')),
                       title = 'Li',
                       lables = c('High', 'Low'),
                       risk.table = T,
                       col = color9)
Li_KM
Li_ROC <- ggplotROC(Li.dat[, c("OS.time", "OS", "RiskScore")])
Li_ROC

# #################################
library(cowplot)
other_model_km_ROC <- cowplot::plot_grid(Sun_KM, Xue_KM,Yu_KM, Yue_KM, 
                                         Sun_ROC,Xue_ROC,Yu_ROC,Yue_ROC,
                                         align = 'hv',
                                         ncol = 4,
                                         labels = LETTERS[1:8])
other_model_km_ROC
ggsave(plot = other_model_km_ROC,
       filename = 'PDFs/other_model_km_ROC.pdf',
       width = 18, height = 9)

other_model_ROC_dat <- all.data.plot
other_model_ROC_dat$Sun <- Sun_risk
other_model_ROC_dat$Xue <- Xue_risk
other_model_ROC_dat$Zou <- Zou_risk
other_model_ROC_dat$Yu <- Yu_risk
other_model_ROC_dat$Yue <- Yue_risk
other_model_ROC_dat$Li <- Li_risk
other_model_ROC_dat$OS.time <- other_model_ROC_dat$OS.time / 365
other_model_ROC_dat <- other_model_ROC_dat[, -c(3:7, 9)]
other_model_ROC_dat <- other_model_ROC_dat[, -c(6,8)]

get_plot_ROC_dat <- function(dat, year) {
  ROC_dat <- data.frame()
  for (model in colnames(dat)[3:ncol(dat)]) {
    tmp <- timeROC(T=as.numeric(dat[, 1]),
                   delta=as.numeric(dat[, 2]),
                   marker=as.numeric(dat[, model]),
                   cause=1,
                   weighting="marginal",
                   times=year,
                   iid=TRUE)
    tmp.dat <- data.frame(TP = tmp$TP[, 2],
                          FP = tmp$FP[, 2],
                          Model =paste0(model, ' AUC: ', round(tmp$AUC[2], 2)))
    ROC_dat <- rbind(ROC_dat, tmp.dat)
  }
  return(ROC_dat)
}

colnames(other_model_ROC_dat)[3:ncol(other_model_ROC_dat)]
other_ROC1 <- get_plot_ROC_dat(other_model_ROC_dat, 1)
ROC1.plot <- ggplot(other_ROC1, aes(x=FP,y=TP, fill=Model))+
  geom_line(aes(colour=Model),lwd=0.75)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  xlab('False positive rate')+ylab('True positive rate')+ 
  scale_colour_manual(values = color9) +
  theme_few() + 
  theme(legend.position=c(1,0),legend.justification=c(1,0),
        legend.background = element_rect(fill = NA, colour = NA))
ROC1.plot


other_ROC3 <- get_plot_ROC_dat(other_model_ROC_dat, 3)
ROC3.plot <- ggplot(other_ROC3, aes(x=FP,y=TP, fill=Model))+
  geom_line(aes(colour=Model),lwd=0.75)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  xlab('False positive rate')+ylab('True positive rate')+ 
  scale_colour_manual(values = color9) +
  theme_few() + 
  theme(legend.position=c(1,0),legend.justification=c(1,0),
        legend.background = element_rect(fill = NA, colour = NA))
ROC3.plot


other_ROC5 <- get_plot_ROC_dat(other_model_ROC_dat, 5)
ROC5.plot <- ggplot(other_ROC5, aes(x=FP,y=TP, fill=Model))+
  geom_line(aes(colour=Model),lwd=0.75)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  xlab('False positive rate')+ylab('True positive rate')+ 
  scale_colour_manual(values = color9) +
  theme_few() + 
  theme(legend.position=c(1,0),legend.justification=c(1,0),
        legend.background = element_rect(fill = NA, colour = NA))
ROC5.plot


other_model_ROC <- cowplot::plot_grid(ROC1.plot,
                                      ROC3.plot,
                                      ROC5.plot,
                                      align = 'hv',
                                      ncol = 3,
                                      labels = LETTERS[1:3])
other_model_ROC
ggsave(plot = other_model_ROC,
       filename = 'PDFs/other_model_ROC.pdf',
       width = 13.5, height = 4.5)
