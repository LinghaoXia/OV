##### 亚群分析 #####
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)


load("~/rawdata/scRNA_virus/virus_3D/virus_3D_anno.RData")
sce <- subset(sce,celltype=="Epithelial")
save(sce,file = "~/rawdata/scRNA_virus/virus_3D/virus_3D_epi.RData")



##### 感染细胞比例 #####
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)


load("~/rawdata/scRNA_virus/virus_3D/virus_3D_epi.RData")

table(sce$group)#查看各组细胞数
prop.table(table(sce$infected))
table(sce$infected, sce$group)#各组不同细胞群细胞数

Cellratio <- prop.table(table(sce$infected, sce$group), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))

##柱状图
library(ggplot2)
library(ggbreak)
pdf(paste0(outdir,"/","04-",gene,"-epi_infected.pdf"),height=8,width=12)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  geom_text(aes(x = Var2, y = Freq, label = scales::percent(Freq)), 
            size = 4, colour = "white")+
  scale_fill_manual(values = c("Bystander" = "#4DBBD5B2", "Infected" = "#DC0000B2", "Naive" = "#91D1C2B2"))  # 设置颜色


VlnPlot(sce, features = "VG161",group.by = "group")
DimPlot(sce, group.by = "infected",split.by = "group", label = TRUE,cols = c("#4DBBD5B2", "#DC0000B2","#91D1C2B2"),pt.size=1.2)
DimPlot(sce, group.by = "infected",label = TRUE,cols = c("#4DBBD5B2", "#DC0000B2","#91D1C2B2"),pt.size=1.2)+DimPlot(sce, group.by = "celltype", label = TRUE,pt.size=1.2)
FeaturePlot(sce, features = "VG161_transcript", cols = c("lightgrey", "red"),pt.size = 0.05)

dev.off()


##### 富集分析 #####
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)


load("~/rawdata/scRNA_virus/virus_3D/virus_3D_epi.RData")
library(monocle)
library(clusterProfiler)
library(org.Hs.eg.db)

sce <- subset(sce,group=="VG161")

#比较cluster0和cluster1的差异表达基因
dge.cluster <- FindMarkers(sce,ident.1 = "Bystander",ident.2 = "Infected",group.by = 'infected')
sig_dge.cluster <- subset(dge.cluster, p_val_adj<0.01&abs(avg_log2FC)>1)


#GO分析(注意是human的)，count表示改变的基因数
ego_ALL <- enrichGO(gene          = row.names(sig_dge.cluster),
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Hs.eg.db', #Mm
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)
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
ego_CC@result$Description <- substring(ego_CC@result$Description,1,30)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,30)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,30)
p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("barplot for Biological process")
p_CC <- barplot(ego_CC,showCategory = 10) + ggtitle("barplot for Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("barplot for Molecular function")
plotc1 <- p_BP/p_CC/p_MF

#KEGG GeneRatio表示差异基因所占比例
genelist <- bitr(row.names(sig_dge.cluster), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')#Hs是人类
# kegg分析的基因名必须要是ENTREZID
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa') #hsa是人类，mmu是小鼠
p1 <- barplot(ekegg, showCategory=10)+ scale_y_discrete(labels = function(x) str_wrap(x, width = 35))
p2 <- dotplot(ekegg, showCategory=10)+ scale_y_discrete(labels = function(x) str_wrap(x, width = 35))
plotc2 = p1/p2

pdf(paste0(outdir,"/","06-",gene,"-epi_GO&KEGG.pdf"),height=8,width=12)
plotc1
plotc2
dev.off()

