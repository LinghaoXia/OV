
##### 代谢组 #####
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggplot2)
library(MNet)
library(data.table)
library(do)

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

metab <- fread("~/rawdata/L-2h_rnaseq/L2HGDH-metabolism.csv")
##处理数据
metab <- as.data.frame(metab)
rownames(metab) <- metab$`Sample Name`
metab <- metab[,-1]
metab <- as.data.frame(t(metab))

metab <- metab %>% arrange(desc(Group))
group <- metab$Group
pair_id <- metab$Pair
metab <- metab[,-c(1:2)]
metab %>% mutate(across(where(is.character), as.numeric))  -> metab
metab <- as.data.frame(t(metab))

group <- ifelse(group == "T", "tumor", "normal")

# ##计算差异——不配对
# result_mlimma_all <- mlimma(metab,group)

##计算差异——配对
library(limma)
design <- model.matrix(~ 0+group + pair_id)
fit <- lmFit(metab,design)
fit <- eBayes(fit)
result_mlimma_all <- topTable(fit,
                           adjust = 'fdr',
                           coef = "groupnormal",
                           n = Inf)
result_mlimma_all$name <- rownames(result_mlimma_all)

write.table(result_mlimma_all,paste0(outdir,"/result_mlimma_all.txt"),quote=F,row.names=F,sep="\t")   

##绘制热图
result_mlimma_filter <- result_mlimma_all %>%
  dplyr::filter(abs(logFC) > 1) %>%
  dplyr::filter(`adj.P.Val` < 0.05)%>% 
  arrange(desc(logFC))%>%
  head(n=10)


if ("L-2-Hydroxyglutaric acid" %in% result_mlimma_filter$name) {
  dat_filter <- metab %>%
    tibble::rownames_to_column(var="label") %>%
    dplyr::filter(label %in% result_mlimma_filter$name) %>%
    tibble::column_to_rownames("label") 
}else{
  dat_filter <- metab %>%
    tibble::rownames_to_column(var="label") %>%
    dplyr::filter(label %in% c("L-2-Hydroxyglutaric acid",result_mlimma_filter$name)) %>%
    tibble::column_to_rownames("label") 
}

##不能出现0
dat_filter[apply(dat_filter, MARGIN = c(1, 2), FUN = function(x) x == 0)] <- 1
pdf(file = paste0(outdir,"/Metabolism-",gene,"-heatmap.pdf"),width = 20,height = 10)
p_heatmap <- MNet::pHeatmap(dat_filter,group,clustering_distance_cols ="manhattan",
                            clustering_method="ward.D",fontsize_row=12)
dev.off()



##### 转录组ifn-y #####
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggplot2)
library(limma)
library(pheatmap)

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)


data_tpm <- as.matrix(read.csv("~/rawdata/L2HGDH_rnaseq/RNA-seq-ifny.csv", row.names = 1))
# data_tpm <- as.matrix(read.csv("~/rawdata/L2HGDH_rnaseq/RNA-seq-ifny.csv", row.names = 1))
sample_info <- data.frame(
  group = factor(c(rep("Control", 3), rep("Kncokdown", 6)))
)

# 创建设计矩阵
design <- model.matrix(~ group, data = sample_info)
# 线性模型拟合
fit <- lmFit(data_tpm, design)
# 进行差异分析
fit <- eBayes(fit)
# 提取差异表达结果
results <- topTable(fit, coef = 2, adjust = "BH", sort.by = "P", number = Inf)
# 查看结果
head(results)

# 提取显著基因并保存
significant_genes <- results[abs(results$logFC)>1.5 & results$adj.P.Val < 0.2, ]
write.csv(significant_genes, file = paste0(outdir,"/IFN-y-",gene,"-significant_genes.csv"), row.names = T)
# 提取 logFC 绝对值前 100 的基因
top_genes <- significant_genes[order(abs(significant_genes$logFC), decreasing = TRUE), ]
# 选择前 100 个基因
top_100_genes <- head(top_genes, 150)


###### 热图 ######
# 提取显著基因的表达数据
heatmap_data <- data_tpm[rownames(top_100_genes), ]
# 标准化数据（可选）
heatmap_data <- t(scale(t(heatmap_data)))
# 定义颜色
my_palette <-  colorRampPalette(colors = c('#11427C','white','#C31E1F'))(100)
# 定义组别
sample_meta <- data.frame(
  row.names = colnames(data_tpm),group = factor(c(rep("ctrl", 3), rep("sh1", 3),rep("sh2", 3)))
)
# 绘制热图
pdf(file = paste0(outdir,"/IFN-y-",gene,"-heatmap_2.pdf"),width = 10,height = 20)
pheatmap(heatmap_data, 
         cluster_rows = T, 
         cluster_cols = F, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         annotation_col = sample_meta,
         angle_col = 45)
         # = list(group = c(case = "#1B9E77",control = "#D95F02")),
         #color = my_palette

dev.off()


###### 单组火山图 ######
library(ggplot2)
library(dplyr) # 用于数据处理
library(gt) # 制作表格
library(ggrepel)
library(openxlsx)

# 计算显著性
cut_off_FDR =0.2 #设置FDR的阈值
cut_off_log2FC =1 #设置log2FC的阈值
results$Sig = ifelse(results$adj.P.Val < cut_off_FDR &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                    abs(results$logFC) >= cut_off_log2FC,  #abs绝对值
                  ifelse(results$logFC > cut_off_log2FC ,'Up','Down'),'no')
results = data.frame(results)
results$gene_name <- row.names(results)


#数据预处理——将上下调基因分开绘制各自的标签框类型#
Up_top_10 =(     #筛选差异显著上调的前10个Gene
  results %>%
    filter(Sig == 'Up') %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    head(10)
)
Up_top_10 %>% gt() #数据制成表
Down_top_10 = (       #筛选差异显著下调的前10个Gene
  results %>%
    filter(Sig == 'Down') %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    head(10)
)

# Down_top_10= (       #筛选差异显著下调的前10个Gene
#   results %>%
#     filter(gene_name %in% c("Cd274","L2hgdh")) %>%
#     arrange(adj.P.Val, desc(abs(logFC)))
# )
Down_top_10 %>% gt() #数据制成表


pdf(file = paste0(outdir,"/IFN-y-",gene,"-volcano_1.pdf"),width = 8,height = 10)

ggplot(results, aes(x =logFC, y= -log10(adj.P.Val), colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + xlim(c(-30, 30)) +  #调整点的颜色和x轴的取值范围
  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = -log10(cut_off_FDR), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  ggtitle("IFN-y") + #标题
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
  ylim(0, 4)  # 设置y轴的范围



dev.off()

###### 多组火山图 ######
library(ggplot2)
library(tidyverse)
library(ggrepel)


# https://zhuanlan.zhihu.com/p/516955474


###### KEGG&GO ######
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
rm(list=ls())
gc()

getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

sig_dge <- read.csv(file = paste0(outdir,"/IFN-y-",gene,"-significant_genes.csv"))
rownames(sig_dge) <- sig_dge[,1]
sig_dge <- sig_dge[,-1]

#GO分析(注意是小鼠的)，count表示改变的基因数
ego_ALL <- enrichGO(gene          = row.names(sig_dge),
                    OrgDb         = 'org.Mm.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
write.csv(ego_ALL,file = paste0(outdir,"/IFN-y-",gene,"-GO.csv"), quote=F, row.names = F)
ego_all <- data.frame(ego_ALL)
plotc3 <- barplot(ego_ALL, x = "GeneRatio", color = "p.adjust", #默认参数（x和color可以根据eG里面的内容更改）
        showCategory =10, #只显示前10
        split="ONTOLOGY") + #以ONTOLOGY类型分开
  facet_grid(ONTOLOGY~., scale='free') #以ONTOLOGY类型分开绘图

ego_CC <- enrichGO(gene          = row.names(sig_dge),
                   OrgDb         = 'org.Mm.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_MF <- enrichGO(gene          = row.names(sig_dge),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Mm.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_BP <- enrichGO(gene          = row.names(sig_dge),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Mm.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) 
#截取每个description的前70个字符，方便后面作图排版
ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)
p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("barplot for Biological process")
p_CC <- barplot(ego_CC,showCategory = 10) + ggtitle("barplot for Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("barplot for Molecular function")
plotc1 <- p_BP/p_CC/p_MF
#KEGG GeneRatio表示差异基因所占比例
genelist <- bitr(row.names(sig_dge), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Mm.eg.db')#Hs是人类
# kegg分析的基因名必须要是ENTREZID
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'mmu',qvalueCutoff = 0.2,pvalueCutoff = 0.2) #hsa是人类
ekegg_df <- as.data.frame(ekegg)
ekegg_df$SYMBOL <- lapply(strsplit(ekegg_df$geneID, "/"), function(entrez_ids) {
  mapIds(org.Mm.eg.db, 
         keys = entrez_ids, 
         keytype = "ENTREZID", 
         column = "SYMBOL",
         multiVals = "first")  # 如果有多个匹配，取第一个
})
ekegg_df$SYMBOL <- sapply(ekegg_df$SYMBOL, paste, collapse = "/")
write.csv(ekegg_df,file = paste0(outdir,"/IFN-y-",gene,"-KEGG.csv"),  quote=F, row.names = F)

p1 <- barplot(ekegg_df, showCategory=10)+ scale_y_discrete(labels = function(x) str_sub(x,1,nchar(x)-28))
p2 <- dotplot(ekegg, showCategory=10)+ scale_y_discrete(labels = function(x) str_sub(x,1,nchar(x)-28))
plotc2 = p1/p2

pdf(paste0(outdir,"/IFN-y-",gene,"-go&kegg_2.pdf"),height=16,width=12)
plotc1
plotc2
plotc3
dev.off()


##### 蛋白组 #####
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggplot2)
library(data.table)
library(do)
library(limma)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gt) # 制作表格

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

pro <- fread("~/rawdata/L-2h_rnaseq/LLC-proteome-ubi.csv")
group <- factor(c(rep("NC", 3), rep("SH", 6)))


##顺转/稳转
expr_data <- as.data.frame(pro[,c(4,5:13)])
rownames(expr_data) <- expr_data[,1]
expr_data <- expr_data[,-1]
# 对数转换
expr_log2 <- log2(expr_data + 1)
# z-score标准化（按蛋白）
expr_scaled <- t(scale(t(expr_log2)))  # 行标准化（每个蛋白）
# 构建设计矩阵
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
# 构建对比矩阵
contrast.matrix <- makeContrasts(SH_vs_NC = SH - NC, levels = design)
# 使用 limma 分析
fit <- lmFit(expr_scaled, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
# 获取结果
results <- topTable(fit2, coef = "SH_vs_NC", number = Inf, adjust.method = "BH")

##绘制热图
# 筛选显著差异蛋白（可自定义阈值）
sig_prots <- results[which(results$P.Value < 0.05 & abs(results$logFC) > 1), ]
write.csv(sig_prots,file=paste0(outdir,"/Proteome-",gene,"-pheatmaap-virus.csv"),row.names = TRUE)
# 选取前100的蛋白
sig_prots_sorted <- sig_prots[order(abs(sig_prots$logFC), decreasing = TRUE), ]
top100_sig_prots <- sig_prots_sorted[1:50, ]
sig_ids <- rownames(top100_sig_prots)
#加上要看的
sig_ids <- c(sig_ids,
"Ubiquitin-like protein ISG15",
"Ubiquitin-like modifier-activating enzyme 7",
"E3 ubiquitin-protein ligase TM129",
"E3 ubiquitin-protein ligase MARCHF6",
"E3 ubiquitin-protein ligase RNF149",
"Ubiquitin domain-containing protein 1",
"Ubiquitin-like protein 3")

# 从表达矩阵中提取这些蛋白的表达值
heat_data <- expr_scaled[sig_ids, ]
# 构建分组信息
group <- c(rep("NC", 3), rep("SH", 6))
subgroup <- c(rep("NC", 3), rep("SH1", 3),rep("SH2", 3))
annotation_col <- data.frame(
  Group = factor(group),
  Subgroup = factor(subgroup)
)
rownames(annotation_col) <- colnames(heat_data)

# 绘制热图
pdf(paste0(outdir,"/Proteome-",gene,"-pheatmaap-virus.pdf"),height=16,width=12)
pheatmap(
  heat_data,
  scale = "row",  # 每行（蛋白）标准化
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  fontsize_row = 6,
  fontsize_col = 10,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Significant Differential Proteins_Stable"
)
dev.off()

# 绘制火山图
special_genes <- read.csv("~/rawdata/L-2h_rnaseq/ubi.csv", stringsAsFactors = FALSE)[[1]]
# 添加注释列
results$gene_name <- rownames(results)
# 添加标识字段
# 计算显著性
cut_off_FDR =0.05 #设置FDR的阈值
cut_off_log2FC =1 #设置log2FC的阈值
results$Significance <- "Not Sig"
results$Significance[results$P.Value < 0.05 & results$logFC > 1]  <- "Up"
results$Significance[results$P.Value  < 0.05 & results$logFC < -1] <- "Down"
# 添加是否是 special gene 的列
results$Special <- ifelse(results$gene_name %in% special_genes, "Special", "Other")
Up_top_10 =(     #筛选差异显著上调的前10个Gene
  results %>%
    filter(Significance == 'Up') %>%
    arrange(P.Value, desc(abs(logFC))) %>%
    filter(Special == "Special")
)
Up_top_10 %>% gt() #数据制成表
Down_top_10 = (       #筛选差异显著下调的前10个Gene
  results %>%
    filter(Significance == 'Down') %>%
    arrange(P.Value, desc(abs(logFC))) %>%
    filter(Special == "Special")
)
Down_top_10 %>% gt() #数据制成表



pdf(paste0(outdir,"/Proteome-",gene,"-volcano-sz.pdf"),height=14,width=12)
ggplot(results, aes(x =logFC, y= -log10(P.Value), colour=Significance)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + xlim(c(-2, 2)) +  #调整点的颜色和x轴的取值范围
  geom_hline(yintercept=-log10(cut_off_FDR),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_vline(xintercept = c(-cut_off_log2FC,cut_off_log2FC), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  ggtitle("Significant Differential Proteins_Transient") + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank()
  )+geom_label_repel(data = Up_top_10,
                     aes(logFC, -log10(P.Value), label = gene_name),
                     size = 5, fill="#CCFFFF",
                     alpha = 0.65, color = "black")+
  geom_label_repel(data = Down_top_10,
                   aes(logFC, -log10(P.Value), label = gene_name),
                   size = 5, fill="#FFCCCC",
                   alpha = 0.65, color = "black")+
  coord_cartesian(ylim = c(1, max(-log10(results$P.Value), na.rm = TRUE))) # 设置y轴的范围

dev.off()


##### 转录组 #####
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggplot2)
library(limma)
library(pheatmap)
library(DESeq2)
library("BiocParallel") #启用多核计算

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

##设定 实验组exp / 对照组ctr
exp="SH1"
ctr="CTRL"

data_tpm <- as.matrix(read.csv("~/rawdata/L2HGDH_rnaseq/counts.csv", row.names = 1))

sample_info <- data.frame(
  group = factor(c(rep("CTRL", 3), rep("SH1", 3)))
)

##构建dds DESeqDataSet
if(T){
  dds <- DESeqDataSetFromMatrix(countData = data_tpm,
                                colData = sample_info,
                                design = ~ group)
}
if(F){  #若上游为salmon
  dds <- DESeqDataSetFromTximport(txi, 
                                  colData = sample_info,
                                  design = ~ group)
}


dds$group <- relevel(dds$group, ref = ctr)   #指定 control group

keep <- rowSums(counts(dds) > 0) >= 2  #Pre-filtering ，过滤低表达基因
dds <- dds[keep,] 
dds <- DESeq(dds,quiet = F) 
res <- results(dds,contrast=c("group", exp, ctr))  #指定提取为exp/ctr结果
resOrdered <- res[order(res$padj),]  #order根据padj从小到大排序结果
tempDEG <- as.data.frame(resOrdered)
DEG_DEseq2 <- na.omit(tempDEG)

write.csv(DEG_DEseq2, file = paste0(outdir,"/",gene,"-significant_gene.csv"), row.names = T)

###### GSEA #####
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggstatsplot)
library(GseaVis)
rm(list=ls())
gc()

getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

sig_dge <- read.csv(file =paste0(outdir,"/",gene,"-significant_gene.csv"))
row.names(sig_dge) <- sig_dge[,1]
sig_dge <- sig_dge[,c(3,6)] #选择log2FoldChange和pvalue（凑成数据框）
colnames(sig_dge) <- c('log2FoldChange','pvalue')
sig_dge$SYMBOL <- rownames(sig_dge)

###创建gsea分析的geneList（包含从大到小排列的log2FoldChange和ENTREZID信息）
df <- bitr(rownames(sig_dge), 
           fromType = "SYMBOL",
           toType =  "ENTREZID",
           OrgDb = "org.Mm.eg.db") #人数据库org.Hs.eg.db 小鼠org.Mm.eg.db
sig_dge <- merge(sig_dge, df, by='SYMBOL')  #按照SYMBOL合并注释信息

geneList <- sig_dge$log2FoldChange
names(geneList) <- sig_dge$ENTREZID
geneList <- sort(geneList, decreasing = T)   #从大到小排序

###gsea富集
options(timeout = 120)  # 设置超时时间为 120 秒
KEGG_kk_entrez <- gseKEGG(geneList     = geneList,
                          organism     = "mmu", #人hsa 鼠mmu
                          pvalueCutoff = 0.05)  #实际为padj阈值,可调整 
KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez, 
                             OrgDb="org.Mm.eg.db",
                             keyType='ENTREZID')#转化id             

GO_kk_entrez <- gseGO(geneList     = geneList,
                      ont          = "ALL",  # "BP"、"MF"和"CC"或"ALL"
                      OrgDb        = "org.Mm.eg.db",#人类org.Hs.eg.db 鼠org.Mm.eg.db
                      keyType      = "ENTREZID",
                      pvalueCutoff = 0.25)   #实际为padj阈值可调整
GO_kk <- DOSE::setReadable(GO_kk_entrez, 
                           OrgDb= "org.Mm.eg.db",
                           keyType='ENTREZID')#转化id 

###选取富集结果
kk_gse <- GO_kk
kk_gse_entrez <- GO_kk_entrez

###单独的gseaplot
terms <- c("GO:0002703","GO:0002443","GO:0002252")  

gseaplot_list <- lapply(terms, function(x){
  gseaNb(object = kk_gse,
         geneSetID = x,
         termWidth = 30,
         addPval = T,
         pvalX = 0.75,
         pvalY = 0.6,
         addGene = "Cd274",
         geneCol= '#4d4d4d',
         kegg = T
  )
})
# addGene= T, #是否添加基因
# markTopgene= T, #是否标注Top基因
# topGeneN= 25, #标注前多少个gene

pdf(paste0(outdir,"/",gene,"-epi_GSEA.pdf"),height=10,width=16)
cowplot::plot_grid(plotlist=gseaplot_list, ncol = 3)
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