library(dplyr)
library(Seurat)
library(patchwork)

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)

load("~/rawdata/scRNA_virus/virus_5D/virus_5D_epi.RData")
sce_5 <- sce
load("~/rawdata/scRNA_virus/virus_3D/virus_3D_epi.RData")
sce_3 <- sce
rm(sce)
sce <- merge(sce_3, y = sce_5, add.cell.ids = c("day3", "day5"), project = "compared")
sce@meta.data$time <- ifelse (grepl("3D",sce@meta.data$orig.ident),'Day3','Day5')

pdf(paste0(outdir,"/","09-",gene,"-compared_MDK.pdf"),height=8,width=12)
VlnPlot(sce, features = c("MDK","NCL"), split.by = "infected",group.by = "group",pt.size = 0)
VlnPlot(sce, features = c("MDK","NCL"), group.by = "group",pt.size = 0)
VlnPlot(sce, 
        group.by = "infected", split.by = 'time',
        cols = c("grey","#A83C1E"),
        features = c("MDK","NCL"), pt.size = 0) &
  stat_summary(fun = median, geom='point', size = 1, colour = "black", position = position_dodge(width = 0.9)) &
  ggpubr::geom_pwc(aes(group = split), 
                   method = "wilcox_test", tip.length = 0, 
                   hide.ns = F, label = "p.adj.format", 
                   p.adjust.method = "fdr") & 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) &
  labs(x='')
dev.off()

save(sce,file = "~/rawdata/scRNA_virus/virus_5D/virus_combined.RData")



##### 富集分析 #####
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)


load("~/rawdata/scRNA_virus/virus_5D/virus_combined.RData")
library(monocle)
library(clusterProfiler)
library(org.Hs.eg.db)

sce <- subset(sce, group %in% c("VG161","VG161_L","VG161_R"))
sce@meta.data$time <- ifelse (grepl("3D",sce@meta.data$orig.ident),'Day3','Day5')
sce <- JoinLayers(sce)


#比较cluster0和cluster1的差异表达基因
dge.cluster <- FindMarkers(sce,ident.1 = "Day5",ident.2 = "Day3",group.by = 'time') #ident1是测试组
write.csv(dge.cluster, file = paste0(outdir,"/09-",gene,"-significant_gene.csv"), row.names = T)



###### 单组火山图 ######
library(ggplot2)
library(dplyr) # 用于数据处理
library(gt) # 制作表格
library(ggrepel)
library(openxlsx)

results <- read.csv(file = paste0(outdir,"/09-",gene,"-significant_gene.csv"),row.names = 1)
results <- results[,c(1,2,5)]
colnames(results) <- c("pvalue","logFC","adj.P.Val")

# 计算显著性
cut_off_FDR =0.05 #设置FDR的阈值
cut_off_log2FC =log2(1.5) #设置log2FC的阈值
results$Sig = ifelse(results$adj.P.Val < cut_off_FDR &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                       abs(results$logFC) >= cut_off_log2FC,  #abs绝对值
                     ifelse(results$logFC > cut_off_log2FC ,'Up','Down'),'no')
results = data.frame(results)
results$gene_name <- row.names(results)


#数据预处理——将上下调基因分开绘制各自的标签框类型#
Up_top_10 =(     #筛选差异显著上调的前10个Gene
  results %>%
    filter(Sig == 'Up') %>%
    filter(!grepl('VG161', gene_name)) %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    head(10)
)
Up_top_10 %>% gt() #数据制成表
Down_top_10 = (       #筛选差异显著下调的前10个Gene
  results %>%
    filter(Sig == 'Down') %>%
    filter(!grepl('VG161', gene_name)) %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    head(10)
)

# Down_top_10= (       #筛选差异显著下调的前10个Gene
#   results %>%
#     filter(gene_name %in% c("Cd274","L2hgdh")) %>%
#     arrange(adj.P.Val, desc(abs(logFC)))
# )
Down_top_10 %>% gt() #数据制成表


pdf(paste0(outdir,"/","09-",gene,"-epi_DEG.pdf"),width = 8,height = 10)

ggplot(results, aes(x =logFC, y= -log10(adj.P.Val), colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + xlim(c(-2, 2)) +  #调整点的颜色和x轴的取值范围
  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = -log10(cut_off_FDR), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  ggtitle("Day 5 vs Day 3") + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank()
  )+geom_label_repel(data = Up_top_10,
                     aes(logFC, -log10(adj.P.Val), label = gene_name),
                     size = 3, fill="#CCFFFF",
                     alpha = 0.65, color = "black")+
  geom_label_repel(data = Down_top_10,
                   aes(logFC, -log10(adj.P.Val), label = gene_name),
                   size = 3, fill="#FFCCCC",
                   alpha = 0.65, color = "black")+
  ylim(0, 300)  # 设置y轴的范围



dev.off()


###### KEGG&GO ######
sig_dge.cluster <- subset(dge.cluster, p_val_adj<0.05 & abs(avg_log2FC)>0.15)
#GO分析(注意是human的)，count表示改变的基因数
ego_ALL <- enrichGO(gene          = row.names(sig_dge.cluster),
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Hs.eg.db', #Mm
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
write.csv(ego_ALL,file = paste0(outdir,"/09-",gene,"-GO.csv"), quote=F, row.names = F)
plotc3 <- barplot(ego_ALL, x = "GeneRatio", color = "p.adjust",font.size = 15, #默认参数（x和color可以根据eG里面的内容更改）
                  showCategory =10, #只显示前10
                  split="ONTOLOGY") + #以ONTOLOGY类型分开
  facet_grid(ONTOLOGY~., scale='free') #以ONTOLOGY类型分开绘图
ego_CC <- enrichGO(gene          = row.names(sig_dge.cluster),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_MF <- enrichGO(gene          = row.names(sig_dge.cluster),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_BP <- enrichGO(gene          = row.names(sig_dge.cluster),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) 
#截取每个description的前70个字符，方便后面作图排版
ego_CC@result$Description <- substring(ego_CC@result$Description,1,100)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,100)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,100)
p_BP <- barplot(ego_BP,showCategory = 10,font.size = 15) + ggtitle("barplot for Biological process")
p_CC <- barplot(ego_CC,showCategory = 10,font.size = 15) + ggtitle("barplot for Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10,font.size = 15) + ggtitle("barplot for Molecular function")
plotc1 <- p_BP/p_MF/p_CC

#KEGG GeneRatio表示差异基因所占比例
genelist <- bitr(row.names(sig_dge.cluster), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')#Hs是人类
# kegg分析的基因名必须要是ENTREZID
genelist <- pull(genelist,ENTREZID)               
options(timeout = 99999)
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa') #hsa是人类，mmu是小鼠

ekegg_df <- as.data.frame(ekegg)
ekegg_df$SYMBOL <- lapply(strsplit(ekegg_df$geneID, "/"), function(entrez_ids) {
  mapIds(org.Hs.eg.db, 
         keys = entrez_ids, 
         keytype = "ENTREZID", 
         column = "SYMBOL",
         multiVals = "first")  # 如果有多个匹配，取第一个
})
ekegg_df$SYMBOL <- sapply(ekegg_df$SYMBOL, paste, collapse = "/")
write.csv(ekegg_df,file = paste0(outdir,"/09-",gene,"-KEGG.csv"),  quote=F, row.names = F)

p1 <- barplot(ekegg, showCategory=10,font.size = 15)+ scale_y_discrete(labels = function(x) str_wrap(x, width = 35))+ggtitle("KEGG")
p2 <- dotplot(ekegg, showCategory=10,font.size = 15)+ scale_y_discrete(labels = function(x) str_wrap(x, width = 35))+ggtitle("KEGG")
plotc2 = p1/p2

genelist <- dge.cluster$avg_log2FC
names(genelist) <- rownames(dge.cluster)
ego <- setReadable(ego_BP, 'org.Hs.eg.db', 'SYMBOL')

pdf(paste0(outdir,"/","09-",gene,"-epi_GO&KEGG.pdf"),height=20,width=20)
plotc1
plotc2
plotc3
barplot(ego_BP,showCategory = 15,font.size = 30) + ggtitle("barplot for Biological process")+
  theme(plot.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 20),       # 图例文字大小
        legend.title = element_text(size = 22, face = "bold") # 图例标题大小（可选）
  )
# cnetplot(ego,color_edge = 'category',foldChange=genelist,showCategory = c("leukocyte migration","leukocyte cell-cell adhesion","regulation of T cell activation"))#igraph::layout_in_circle
dev.off()



##### 单基因比较 #####
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)

source("~/code/new_function.R")
load("~/rawdata/scRNA_virus/virus_5D/virus_combined.RData")
sce <- JoinLayers(sce)
A <- singlecell_gene_test(sce, 
                          genes.use = c('MDK','NCL'),
                          group.by = 'time', 
                          comp = c("Day3", "Day5"))

anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig

plots_violins <- VlnPlot(sce, 
                         cols = c("limegreen", "navy"),
                         pt.size = 0,
                         group.by = "orig.ident",
                         features = c('MDK','NCL'), 
                         ncol = 3, 
                         log = FALSE,
                         combine = FALSE)

for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_discrete(labels = c("Female","Male"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}

CombinePlots(plots_violins)