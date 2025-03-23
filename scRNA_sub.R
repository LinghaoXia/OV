##### 亚群分析 #####
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)


load("~/rawdata/scRNA_virus/virus_5D/virus_5D_anno.RData")
sce <- subset(sce,celltype=="Epithelial")
save(sce,file = "~/rawdata/scRNA_virus/virus_5D/virus_5D_epi.RData")



##### 感染细胞比例 #####
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)


load("~/rawdata/scRNA_virus/virus_5D/virus_5D_epi.RData")

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

sce <- subset(sce, group %in% c( "VG161"))


#比较cluster0和cluster1的差异表达基因
dge.cluster <- FindMarkers(sce,ident.1 = "Bystander",ident.2 = "Infected",group.by = 'infected')
write.csv(dge.cluster, file = paste0(outdir,"/06-",gene,"-significant_gene.csv"), row.names = T)


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
write.csv(ego_ALL,file = paste0(outdir,"/06-",gene,"-GO.csv"), quote=F, row.names = F)
plotc3 <- barplot(ego_ALL, x = "GeneRatio", color = "p.adjust", #默认参数（x和color可以根据eG里面的内容更改）
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
p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("barplot for Biological process")
p_CC <- barplot(ego_CC,showCategory = 10) + ggtitle("barplot for Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("barplot for Molecular function")
plotc1 <- p_BP/p_MF/p_CC

#KEGG GeneRatio表示差异基因所占比例
genelist <- bitr(row.names(sig_dge.cluster), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')#Hs是人类
# kegg分析的基因名必须要是ENTREZID
genelist <- pull(genelist,ENTREZID)               
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
write.csv(ekegg_df,file = paste0(outdir,"/06-",gene,"-KEGG.csv"),  quote=F, row.names = F)

p1 <- barplot(ekegg, showCategory=10)+ scale_y_discrete(labels = function(x) str_wrap(x, width = 35))+ggtitle("KEGG")
p2 <- dotplot(ekegg, showCategory=10)+ scale_y_discrete(labels = function(x) str_wrap(x, width = 35))+ggtitle("KEGG")
plotc2 = p1/p2

pdf(paste0(outdir,"/","06-",gene,"-epi_GO&KEGG.pdf"),height=20,width=16)
plotc1
plotc2
plotc3
dev.off()

###### GSEA #####
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggstatsplot)
library(GseaVis)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)

sig_dge <- read.csv(file =paste0(outdir,"/06-",gene,"-significant_gene.csv"))
row.names(sig_dge) <- sig_dge[,1]
sig_dge <- sig_dge[,c(2,3)] #选择log2FoldChange和pvalue（凑成数据框）
colnames(sig_dge) <- c('pvalue','log2FoldChange')
sig_dge$SYMBOL <- rownames(sig_dge)

###创建gsea分析的geneList（包含从大到小排列的log2FoldChange和ENTREZID信息）
df <- bitr(rownames(sig_dge), 
           fromType = "SYMBOL",
           toType =  "ENTREZID",
           OrgDb = "org.Hs.eg.db") #人数据库org.Hs.eg.db 小鼠org.Mm.eg.db
sig_dge <- merge(sig_dge, df, by='SYMBOL')  #按照SYMBOL合并注释信息
geneList <- sig_dge$log2FoldChange
names(geneList) <- sig_dge$ENTREZID
geneList <- sort(geneList, decreasing = T)   #从大到小排序

###gsea富集
KEGG_kk_entrez <- gseKEGG(geneList     = geneList,
                          organism     = "hsa", #人hsa 鼠mmu
                          pvalueCutoff = 0.25)  #实际为padj阈值,可调整 
KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez, 
                             OrgDb="org.Hs.eg.db",
                             keyType='ENTREZID')#转化id             

GO_kk_entrez <- gseGO(geneList     = geneList,
                      ont          = "ALL",  # "BP"、"MF"和"CC"或"ALL"
                      OrgDb        = "org.Hs.eg.db",#人类org.Hs.eg.db 鼠org.Mm.eg.db
                      keyType      = "ENTREZID",
                      pvalueCutoff = 0.25)   #实际为padj阈值可调整
GO_kk <- DOSE::setReadable(GO_kk_entrez, 
                           OrgDb= "org.Hs.eg.db",
                           keyType='ENTREZID')#转化id 

###选取富集结果
kk_gse <- GO_kk
kk_gse_entrez <- GO_kk_entrez

###单独的gseaplot
terms <- c("GO:1905517","GO:0071674")

gseaplot_list <- lapply(terms, function(x){
  gseaNb(object = kk_gse_entrez,
         geneSetID = x,
         termWidth = 30,
         addPval = T,
         pvalX = 0.75,
         pvalY = 0.6
  )
})

pdf(paste0(outdir,"/","06-",gene,"-epi_GSEA.pdf"),height=10,width=16)
cowplot::plot_grid(plotlist=gseaplot_list, ncol = 2)
dev.off()

###合并的gseaplot（未改）
#一般认为|NES|>1，NOM pvalue<0.05，FDR（padj）<0.25的通路是显著富集的
kk_gse_cut <- kk_gse[kk_gse$pvalue<0.05 & kk_gse$p.adjust<0.25 & abs(kk_gse$NES)>1]
kk_gse_cut_down <- kk_gse_cut[kk_gse_cut$NES < 0,]
kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0,]
#选择展现NES前几个通路 
down_gsea <- kk_gse_cut_down[tail(order(kk_gse_cut_down$NES,decreasing = T),10),]
up_gsea <- kk_gse_cut_up[head(order(kk_gse_cut_up$NES,decreasing = T),10),]
diff_gsea <- kk_gse_cut[head(order(abs(kk_gse_cut$NES),decreasing = T),10),]
# 合并 GSEA通路 
gseap2 <- gseaplot2(kk_gse,
                    up_gsea$ID,#富集的ID编号
                    title = "UP_GSEA_all",#标题
                    color = "red",#GSEA线条颜色
                    base_size = 20,#基础字体大小
                    rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                    subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                    ES_geom = "line",#enrichment score用线还是用点"dot"
                    pvalue_table = T) #显示pvalue等信息
ggsave(gseap2, filename = "GSEA_up_all.pdf",width =12,height =12)

