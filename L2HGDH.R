
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
write.csv(ego_ALL,file = paste0(outdir,"/IFN-y-",gene,"-GO.csv"), sep="\t", quote=F, row.names = F)
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
write.csv(ekegg,file = paste0(outdir,"/IFN-y-",gene,"-KEGG.csv"), sep="\t", quote=F, row.names = F)
p1 <- barplot(ekegg, showCategory=10)+ scale_y_discrete(labels = function(x) str_sub(x,1,nchar(x)-28))
p2 <- dotplot(ekegg, showCategory=10)+ scale_y_discrete(labels = function(x) str_sub(x,1,nchar(x)-28))
plotc2 = p1/p2

pdf(paste0(outdir,"/IFN-y-",gene,"-go&kegg_2.pdf"),height=16,width=12)
plotc1
plotc2
plotc3
dev.off()
