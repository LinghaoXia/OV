##### T细胞亚群聚类 #####
library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)
library(harmony)
library(sctransform)
library(future)
library(glmGamPoi)
plan("multisession", workers = 16)
options(future.globals.maxSize= 1024^4)
plan()
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)


load("~/rawdata/scRNA_virus/virus_3D/virus_3D_anno.RData")
###提取亚群
sce <- subset(sce, celltype=="Tcell")

###### NormalizeData ######
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4)  #对数据进行标准化
sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000) #寻找高变基因
# 找出前10高可变基因用于后续可视化
top10 <- head(VariableFeatures(sce), 10)
# 高变基因可视化
plot1 <- VariableFeaturePlot(sce)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
# PCA前的scale和PCA
all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) #默认最大PC数为50，可查阅函数help自行修改参数
# 线性降维（PCA），默认用高变基因集，但也可通过features参数自己指定；
sce <- RunPCA(sce,features=VariableFeatures(object=sce))
# 检查PCA分群结果，这里只展示前12个PC，每个PC只显示3个基因；
print(sce[["pca"]],dims=1:12,nfeatures = 3)
# # 方法1：Jackstraw置换检验算法：重复取样（原数据的1%），重跑PCA，鉴定p-value较小的PC；计算’null distribution‘（即零假设成立时）时的基因scores；
# sce <- JackStraw(sce,num.replicate = 100)
# sce <- ScoreJackStraw(sce,dims=1:20)
# JackStrawPlot(sce,dims=1:30)
# # 方法2：肘部图（碎石图），基于每个主成分对方差解释率的排名；
# ElbowPlot(sce)
# 方法3：生信技能树
pct <- sce [["pca"]]@stdev / sum( sce [["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
# 基于PCA空间中的欧式距离计算nearest neighbor graph，优化任意两个细胞间的距离权重（输入上一步得到的PC维数）；
sce <- FindNeighbors(sce, dims = 1:16) # 前10个PC
# 用umap的方法，并可视化
sce <- RunUMAP(sce, dims = 1:16)
# 用tsne的方法，并可视化
sce <- RunTSNE(sce,dims=1:16)

###### 合适的分辨率 ######
#接着优化模型，resolution参数决定下游聚类分析得到的分群数，对于3k左右的细胞，设为0.4-1.2能得到较好的结果（官方说明）；如果数据量增大，该参数也应该适当增大；
library(clustree)
library(patchwork)
library(cluster)
sce <- FindClusters(sce, resolution = c(seq(.1,1.5,.1))) # 多个分辨率
clustree(sce, prefix = 'RNA_snn_res.') + coord_flip()
clustree_plt <- clustree(sce, prefix = 'RNA_snn_res.')
cell_dists <- dist(sce@reductions$pca@cell.embeddings,method = "euclidean")
head(cell_dists)
cluster_info <- sce@meta.data[,grepl(paste0(DefaultAssay(sce),"_snn_res"),
                                     colnames(sce@meta.data))] %>%
  dplyr::mutate_all(as.character) %>%
  dplyr::mutate_all(as.numeric)
head(cluster_info)[,1:8]
si= silhouette(cluster_info[,1], cell_dists) %>%head()
si
silhouette_res <- apply(cluster_info, 2, function(x){
  si <- silhouette(x, cell_dists)
  if(!any(is.na(si))) {
    mean(si[, 'sil_width'])
  } else {
    NA
  }
})
silhouette_res#峰顶最优
sce[["opt_clust_integrated"]] <- sce[["RNA_snn_res.0.8"]] #"RNA_snn_res.0.6"#names(which.max(silhouette_res))
Idents(sce) = "opt_clust_integrated"
# 去除多余分辨率
spam_cols <- grepl(paste0(DefaultAssay(sce), "_snn_res"),
                   colnames(sce@meta.data)) |
  grepl("seurat_clusters",colnames(sce@meta.data))
sce@meta.data <- sce@meta.data[,!spam_cols]
save(sce,file = "~/rawdata/scRNA_virus/virus_3D/virus_3D_tcell.RData")

pdf(paste0(outdir,"/","08-",gene,"-tcell.pdf"),height=10,width=6)
DimPlot(sce, reduction = 'umap', group.by = 'opt_clust_integrated',label = TRUE, pt.size = 0.5)
dev.off()


##### T细胞亚群注释 #####
trace(scRNAtoolVis:::jjDotPlot, edit = T)
library(Seurat) 
library(ggplot2)
library(dplyr)
library(scRNAtoolVis)
library(ggdendro)
library(legendry)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)

load("~/rawdata/scRNA_virus/virus_3D/virus_3D_tcell.RData")

###### genes to check ######
markers <- c("CD4", #cd4+ Tcells            
             "CD8A","CD8B",  #cd8+ Tcells 
             "CCR7", "SELL", #T naive CD44-CCR7+SELL(CD62L)+ 
             "FASLG","FAS",#T effecctor             
             "IL7R")  #T memory   CCR7和CD62在Tcm和Tem上不同  
                

markers_plot <- data.frame(cluster = c(rep("CD4 Tcell",1),                                  
                                       rep("CD8 Tcell",2),                                  
                                       rep("naive Tcell",2),                                  
                                       rep("effector Tcell",2),                                  
                                       rep("memory Tcell",1)),                                
                           gene = markers)

pdf(paste0(outdir,"/","08-",gene,"-Tcell_markers.pdf"),height=10,width=6)
CellDimPlot(
  srt = sce, group.by = c("opt_clust_integrated"),
  reduction = "UMAP", theme_use = "theme_blank"
)
scRNAtoolVis::jjDotPlot(object = sce,          
                        markerGene = markers_plot,   
                        anno = TRUE,
                        id = 'opt_clust_integrated', 
                        ytree=T,
                        tree.pos = 'right',
                        textSize = 15,          
                        base_size= 15,          
                        plot.margin = c(8,1,3,1))
FeaturePlot(sce,features = markers,cols = c("lightgrey" ,"red"),
            combine = TRUE,raster=FALSE)
dev.off()


###### annotation ######
celltype <- c(    "0"="CD4+ Tcell",
                  "1"="other Tcell",
                  "2"="other Tcell",
                  "3"="CD8+ Tcell",
                  "4"="CD8+ Tcell",
                  "5"="other Tcell",
                  "6"="CD8+ Tcell",
                  "7"="CD8+ Tcell",
                  "8"="CD8+ Tcell",
                  "9"="other Tcell",
                  "10" = "CD4+ Tcell")
sce@meta.data$t_type <- celltype[sce@meta.data$seurat_clusters]

tcelltype <- c(   "0"="naive Tcell",
                  "1"="naive Tcell",
                  "2"="naive Tcell",
                  "3"="effector Tcell",
                  "4"="effector Tcell",
                  "5"="memory Tcell",
                  "6"="effector Tcell",
                  "7"="memory Tcell",
                  "8"="naive Tcell",
                  "9"="memory Tcell",
                  "10" = "memory Tcell")
sce@meta.data$t_function <- tcelltype[sce@meta.data$seurat_clusters]
save(sce,file = "~/rawdata/scRNA_virus/virus_3D/virus_3D_tcell.RData")

###### plotting ######
pdf(paste0(outdir,"/","08-",gene,"-tcell_ratio.pdf"),height=8,width=12)

####### *CD4 CD8 比例 #######
Cellratio <- prop.table(table(celltype=sce$t_type, group=sce$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- data.frame(Cellratio)

library(reshape2)
cellper <- dcast(Cellratio,group~celltype, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]
#添加分组信息
sample <- unique(sce$orig.ident)
group <- sub("(_1|_2|_3)$", "", sample)
samples <- data.frame(sample, group)#创建数据框

rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$group <- samples[rownames(cellper),'group']#R添加列
pplist = list()
sce_groups = unique(sce$t_type)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
for(group_ in sce_groups){
  cellper_  = cellper[,c('sample','group',group_)]
  colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
  cellper_$group <- factor(cellper_$group , levels =c("Vehicle_3D","VG161_3D"))
  cellper_$percent = as.numeric(cellper_$percent)#数值型数据
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#上下分位数
  print(group_)
  print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
    scale_y_continuous(expand =  expansion(mult = c(0, 0.05)),limits = c(0, NA))+
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10,angle = 45, hjust = 1),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain')) + 
    labs(title = group_,y=paste0(group_,"/T cells")) +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  pplist[[group_]] = pp1
}

library(cowplot)
plot_grid(pplist[[1]],
          pplist[[2]],
          pplist[[3]],nrow = 1)


###### *Te Tm Tn比例 ######

Cellratio <- prop.table(table(celltype=sce$t_function, group=sce$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- data.frame(Cellratio)

library(reshape2)
cellper <- dcast(Cellratio,group~celltype, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]
#添加分组信息
sample <- unique(sce$orig.ident)
group <- sub("(_1|_2|_3)$", "", sample)
samples <- data.frame(sample, group)#创建数据框

rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$group <- samples[rownames(cellper),'group']#R添加列
pplist = list()
sce_groups = unique(sce$t_function)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
for(group_ in sce_groups){
  cellper_  = cellper[,c('sample','group',group_)]
  colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
  cellper_$group <- factor(cellper_$group , levels =c("Vehicle_3D","VG161_3D"))
  cellper_$percent = as.numeric(cellper_$percent)#数值型数据
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#上下分位数
  print(group_)
  print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
    scale_y_continuous(expand =  expansion(mult = c(0, 0.05)),limits = c(0, NA))+
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10,angle = 45, hjust = 1),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain')) + 
    labs(title = group_,y=paste0(group_,"/T cells")) +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  pplist[[group_]] = pp1
}

library(cowplot)
plot_grid(pplist[[1]],
          pplist[[2]],
          pplist[[3]],nrow = 1)

###### *Te Tm Tn占 CD4/CD8比例 ######

Cellratio <- prop.table(table(celltype=sce$t_function, group=sce$orig.ident,tcell=sce$t_type), margin = c(2, 3))#计算各组样本不同细胞群比例
Cellratio <- data.frame(Cellratio)

library(reshape2)
cellper <- dcast(Cellratio,group+tcell~celltype, value.var = "Freq")#长数据转为宽数据
sce_tcell <- unique(sce$t_type) 
pplist = list()

for(tcell_ in sce_tcell){
  cellper_1 <- cellper[cellper$tcell == tcell_,]
  rownames(cellper_1) <- cellper_1[,1]
  cellper_1 <- cellper_1[,-1]
  #添加分组信息
  sample <- unique(sce$orig.ident)
  group <- sub("(_1|_2|_3)$", "", sample)
  samples <- data.frame(sample, group)#创建数据框
  
  rownames(samples)=samples$sample
  cellper_1$sample <- samples[rownames(cellper_1),'sample']#R添加列
  cellper_1$group <- samples[rownames(cellper_1),'group']#R添加列
  
  sce_groups = unique(sce$t_function)
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  library(cowplot)
  library(grid)
  for(group_ in sce_groups){
    cellper_  = cellper_1[,c('sample','group',group_)]
    colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
    cellper_$group <- factor(cellper_$group , levels =c("Vehicle_3D","VG161_3D"))
    cellper_$percent = as.numeric(cellper_$percent)#数值型数据
    cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                        lower = quantile(percent, 0.25),
                                                        mean = mean(percent),
                                                        median = median(percent))#上下分位数
    print(group_)
    print(cellper_$median)
    
    
    pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
      scale_y_continuous(expand =  expansion(mult = c(0, 0.05)),limits = c(0, NA))+
      geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
      stat_summary(fun=mean, geom="point", color="grey60") +
      theme_cowplot() +
      theme(axis.text = element_text(size = 10,angle = 45, hjust = 1),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain')) + 
      labs(title = group_,y=paste0(group_,"/",tcell_)) +
      geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
    
    pplist[[paste0(tcell_,group_)]] = pp1
  }
  
}

library(cowplot)
plot_grid(pplist[[1]],
          pplist[[2]],
          pplist[[3]],
          pplist[[4]],
          pplist[[5]],
          pplist[[6]])

dev.off()
