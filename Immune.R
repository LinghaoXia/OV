library(tidyverse)
library(patchwork)

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)


##### 数据处理——easier输入格式 #####
library(data.table)
library("easier")

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

load("~/rawdata/LUAD/LUAD_tpm.RData")
# 考虑到部分基因有多个对应关系，需要进一步处理（保留作者给定的gene symbol）
genes_info <- easier:::reannotate_genes(cur_genes = colnames(data_tpm))

## 去除不支持的基因symbol
non_na <- !is.na(genes_info$new_names)
data_tpm <- data_tpm[,non_na]
genes_info <- genes_info[non_na, ]

## 去除 entries that are withdrawn
if(length(which(genes_info$new_names == "entry withdrawn"))!=0){
  data_tpm <- data_tpm[,-which(genes_info$new_names == "entry withdrawn")]
  genes_info <- genes_info[-which(genes_info$new_names == "entry withdrawn"),]
}
## 找出重复基因
newnames_dup <- unique(genes_info$new_names[duplicated(genes_info$new_names)])
newnames_dup_ind <- do.call(c, lapply(newnames_dup, function(X) which(genes_info$new_names == X)))
newnames_dup <- genes_info$new_names[newnames_dup_ind]

if(length(newnames_dup)!=0){
  ## 检索重复基因的数据
  tmp <- data_tpm[,genes_info$old_names[genes_info$new_names %in% newnames_dup]]
  ## 删除重复基因的数据
  data_tpm <- data_tpm[,-which(colnames(data_tpm) %in% colnames(tmp))]
  ## 整合重复基因的数据
  dup_genes <- genes_info$new_names[which(genes_info$new_names %in% newnames_dup)]
  names(dup_genes) <- rownames(tmp)
  if (anyDuplicated(newnames_dup)){
    tmp2 <- stats::aggregate(tmp, by = list(dup_genes), FUN = "mean")
    rownames(tmp2) <- tmp2$Group.1
    tmp2$Group.1 <- NULL
    # 整理归纳到一个表达矩阵
    data_tpm <- rbind(data_tpm, tmp2)
  }
}

data_tpm <- as.data.frame(t(data_tpm))
save(data_tpm,file="~/rawdata/LUAD/LUAD_tpm_immune.RData")



##### 计算免疫治疗反应的Hallmarks #####
library(data.table)
library("easier")

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

load("~/rawdata/LUAD/LUAD_tpm_immune.RData")
##hallmarkers
hallmarks_of_immune_response <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
immune_response_scores <- compute_scores_immune_response(RNA_tpm = data_tpm, 
                                                         selected_scores = hallmarks_of_immune_response)

head(immune_response_scores) 



##### 免疫细胞浸润  #####
library(data.table)
library("easier")
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(corrplot) 
library(cowplot)

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

load("~/rawdata/LUAD/LUAD_tpm_immune.RData")
##quanTIseq
cell_fractions <- compute_cell_fractions(RNA_tpm = data_tpm)
colnames(cell_fractions)[8] <- "CD8 T"

cell <- colnames(cell_fractions)
data_tpm <- as.data.frame(t(data_tpm))

##对免疫检查点基因进行循环，找出与目标基因具有相关的免疫检查点
x=as.numeric(data_tpm[,gene])
outTab=data.frame()
pFilter=0.001             #相关性检验pvalue的过滤条件
#循环遍历每个与目标基因相关的基因（sameGene 中的每个基因）
for(i in cell){
  if(i==gene){next}
  # 提取当前相关基因的表达量数据，并转换为数字向量
  y=as.numeric(cell_fractions[,i])
  # 使用cor.test函数计算目标基因和当前相关基因的皮尔逊相关系数和p-value
  corT=cor.test(x, y, method = 'pearson')
  cor=corT$estimate
  # 如果p-value小于设定的过滤条件（pFilter），则将结果添加到结果表格（outTab）
  pvalue=corT$p.value
  if(pvalue < pFilter) {
    outTab=rbind(outTab, cbind(Query=gene, Cell=i, cor, pvalue))
  }
}
#输出相关性结果文件
write.table(file=paste0(outdir,"/LUAD-","CELL-corResult.txt"), outTab, sep="\t", quote=F, row.names=F) 


##相关性矩阵
# 重新组织数据矩阵，提取目标基因和之前筛选出来的相关基因的表达数据
data <- data_frame(data_tpm[,gene],cell_fractions)
colnames(data)[1] <- gene
# 计算新的数据矩阵的相关性矩阵
data <- na.omit(data)
M=cor(data)

#绘制相关性图形
pdf(file=paste0(outdir,"/LUAD-","CELL-corpot.pdf"),width=7,height=7) 
# 绘制相关性矩阵的热图
corrplot(M,           # 相关性矩阵
         method = "circle",          # 使用圆圈的方式来绘制热图
         order = "original",
         type = "upper",
         col=colorRampPalette(c("#104E8B", "white", "red"))(50)
)
dev.off()

##相关性散点图
Cox_Plot <- function(cell){  
  ggscatter(data, x = gene, y = cell, 
            color = "red3",fill = "lightgray",
            add = "reg.line", conf.int = TRUE, 
            add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
            cor.coef = T,
            cor.coeff.args = list(),
            cor.method = "pearson",
            cor.coef.size = 4)
}

if("CD8 T" %in% as.vector(outTab[,2])){
  cox_list <- lapply(as.vector(outTab[,2]), Cox_Plot)
}else {
  
  cox_list <- lapply(c("CD8 T",as.vector(outTab[,2])), Cox_Plot)
}
pdf(file = paste0(outdir,"/LUAD-","CELL-corpot2.pdf"),width = 17,height = 10)
plot_grid(plotlist = cox_list, align = "h",ncol = 3)
dev.off()



##### 通路评分 #####
library(data.table)
library("easier")

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

load("~/rawdata/LUAD/LUAD_tpm_immune.RData")
# 考虑到部分通路相关基因可能被纳入计算免疫反应评分，
# 因此参数remove_sig_genes_immune_response  设置为True， 去除这部分重复基因进行计算，
# 初次使用可以尝试使用效果。
pathway_activities <- compute_pathway_activity(RNA_counts = data_tpm,
                                               remove_sig_genes_immune_response = TRUE)

##### 转录因子分析 #####
library(data.table)
library("easier")

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

load("~/rawdata/LUAD/LUAD_tpm_immune.RData")
tf_activities <- compute_TF_activity(RNA_tpm = RNA_tpm)

##### 配体-受体分析与细胞互作 #####
library(data.table)
library("easier")

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

load("~/rawdata/LUAD/LUAD_tpm_immune.RData")
#> lrpair_weights 是个data.frame
lrpair_weights <- compute_LR_pairs(RNA_tpm = RNA_tpm,
                                   cancer_type = "LUAD")
#> LR signature genes found in data set: 629/644 (97.7%)
#> Ligand-Receptor pair weights computed 
head(lrpair_weights)[,1:5]
#根据上一步量化细胞间的关联强度
ccpair_scores <- compute_CC_pairs(lrpairs = lrpair_weights, 
                                  cancer_type = "LUAD")
head(ccpair_scores) 



##### 免疫检查点 #####
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(corrplot) 
library(cowplot)

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

load("~/rawdata/LUAD/LUAD_tpm_immune.RData")
checkpoint <- fread("~/rawdata/Immune/checkpoint.csv")

data_tpm <- as.data.frame(t(data_tpm))

##对免疫检查点基因进行循环，找出与目标基因具有相关的免疫检查点
x=as.numeric(data_tpm[,gene])
outTab=data.frame()
pFilter=0.001             #相关性检验pvalue的过滤条件
#循环遍历每个与目标基因相关的基因（sameGene 中的每个基因）
for(i in checkpoint$`Check-point`){
  if(i==gene){next}
  # 提取当前相关基因的表达量数据，并转换为数字向量
  y=as.numeric(data_tpm[,i])
  # 使用cor.test函数计算目标基因和当前相关基因的皮尔逊相关系数和p-value
  corT=cor.test(x, y, method = 'pearson')
  cor=corT$estimate
  # 如果p-value小于设定的过滤条件（pFilter），则将结果添加到结果表格（outTab）
  pvalue=corT$p.value
  if(pvalue < pFilter) {
    outTab=rbind(outTab, cbind(Query=gene, Gene=i, cor, pvalue))
  }
}
#输出相关性结果文件
write.table(file=paste0(outdir,"/LUAD-","Checkpoint-corResult.txt"), outTab, sep="\t", quote=F, row.names=F) 


##相关性矩阵
# 重新组织数据矩阵，提取目标基因和之前筛选出来的相关基因的表达数据
if("CD274" %in% outTab$Gene){
  data=data_tpm[,c(gene, outTab$Gene)]
}else {
  data=data_tpm[,c(gene, "CD274",outTab$Gene)]
}
# 计算新的数据矩阵的相关性矩阵
data <- na.omit(data)
M=cor(data)

#绘制相关性图形
pdf(file=paste0(outdir,"/LUAD-","Checkpoint-corpot.pdf"),width=7,height=7) 
# 绘制相关性矩阵的热图
corrplot(M,           # 相关性矩阵
         method = "circle",          # 使用圆圈的方式来绘制热图
         order = "original",
         type = "upper",
         col=colorRampPalette(c("#104E8B", "white", "red"))(50)
)
dev.off()

##相关性散点图
Cox_Plot <- function(checkpoint){  
  ggscatter(data, x = gene, y = checkpoint, 
            color = "red3",fill = "lightgray",
            add = "reg.line", conf.int = TRUE, 
            add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
            cor.coef = T,
            cor.coeff.args = list(),
            cor.method = "pearson",
            cor.coef.size = 4)
}

if("CD274" %in% outTab$Gene){
  cox_list <- lapply(outTab$Gene, Cox_Plot)
}else {
  cox_list <- lapply(c("CD274",outTab$Gene), Cox_Plot)
}
pdf(file = paste0(outdir,"/LUAD-","Checkpoint-corpot2.pdf"),width = 17,height = 10)

plots_per_page <- 6
total_plots <- length(cox_list)
pages <- ceiling(total_plots / plots_per_page)

for (i in seq_len(pages)) {
  start_index <- (i - 1) * plots_per_page + 1
  end_index <- min(i * plots_per_page, total_plots)
  # 从图列表中提取对应的图
  current_plots <- cox_list[start_index:end_index]
  # 填充空图以保证每页完整
  while (length(current_plots) < plots_per_page) {
    current_plots <- c(current_plots, ggplot() + theme_void())
  }
  # 绘制当前页面
  print(plot_grid(plotlist = current_plots, align = "h", nrow = 2, ncol = 3))
}

dev.off()


