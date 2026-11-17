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
sce <- FindClusters(sce, resolution = c(seq(.1,1,.1))) # 多个分辨率
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
sce[["opt_clust_integrated"]] <- sce[["RNA_snn_res.0.3"]] #"RNA_snn_res.0.6"#names(which.max(silhouette_res))
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
jjDotPlot(object = sce,          
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


##### 运行前提 #####
# 续接现有T细胞分析：当前会话须有sce、outdir及已计算的UMAP。
# 不清理环境、不重新聚类；只读取已有RNA log-normalized data，不改变sce。

##### 加载依赖 #####
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
# scRNAtoolVis存在且兼容时优先沿用jjDotPlot；否则用同口径ggplot2绘制。

##### 路径与绘图参数 #####
cluster_col <- "opt_clust_integrated"
sample_col <- "orig.ident"
marker_assay <- "RNA"
marker_reduction <- "umap"
marker_dot_engine <- "jjDotPlot" # 也可显式设为"ggplot2"，不需安装scRNAtoolVis。
marker_show_tree <- FALSE       # 默认所有模块固定cluster顺序；TRUE时jj按各模块重排树序。
marker_feature_per_page <- 9L
marker_feature_ncol <- 3L
marker_raster <- TRUE           # 只栅格化UMAP点层，文字和图例仍为矢量。
marker_point_size <- 0.3
marker_z_limits <- c(-2.5, 2.5)
marker_feature_colors <- c("#E6E9EC", "#2D5375")
marker_z_colors <- c("#426B91", "#F7F7F5", "#C35B3F")
marker_na_color <- "#B8BDC3"
marker_include_sample_summary <- TRUE
# 同一现有outdir下新增子目录，不覆盖原08-marker PDF。
if (!exists("outdir") || !is.character(outdir) || length(outdir) != 1L || is.na(outdir) || !nzchar(outdir)) stop("请先定义当前分析的outdir。")
marker_outdir <- file.path(path.expand(outdir), "Tcell_markers")

##### Marker清单：分组是核对用途，不是自动注释结果 #####
marker_sets <- list(
  "T lineage" = c("TRAC", "CD3D", "CD3E"),
  "CD4 axis" = c("CD4", "ZBTB7B"),
  "CD8 axis" = c("CD8A", "CD8B", "RUNX3"),
  "Naive/Tcm" = c("CCR7", "SELL", "TCF7", "LEF1", "MAL"),
  "Memory" = c("IL7R", "LTB", "CD27", "CD28", "GZMK", "KLRB1", "S100A4", "GPR183"),
  "Cytotoxic" = c("NKG7", "PRF1", "GZMB", "GZMH", "GNLY", "CTSW", "CCL5", "FGFBP2", "KLRG1"),
  "Activation" = c("FASLG", "FAS", "CD69", "ICOS"),
  "Treg" = c("FOXP3", "IL2RA", "CTLA4", "IKZF2"),
  "NK check" = c("KLRD1", "NCR1", "FCGR3A"),
  "gdT check" = c("TRDC", "TRGC1", "TRGC2"),
  "Cycling" = c("MKI67", "TOP2A"),
  "Checkpoint" = c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "TOX", "ENTPD1"))
marker_panels <- list(
  "T cell identity" = c("T lineage", "CD4 axis", "CD8 axis", "Treg", "NK check", "gdT check"),
  "T cell differentiation" = c("Naive/Tcm", "Memory", "Cytotoxic"),
  "T cell state" = c("Activation", "Checkpoint", "Cycling"))
markers <- unname(unlist(marker_sets, use.names = FALSE))
markers_plot <- data.frame(cluster = rep(names(marker_sets), lengths(marker_sets)), gene = markers,
                           stringsAsFactors = FALSE)
if (anyDuplicated(markers) || anyDuplicated(unlist(marker_panels)) ||
    !setequal(unlist(marker_panels), names(marker_sets))) stop("marker须唯一，所有组须恰好分配一个总览页。")

##### 校验当前对象、聚类及表达层 #####
if (!exists("sce") || !inherits(sce, "Seurat")) stop("请先载入当前T细胞Seurat对象sce。")
if (!marker_assay %in% names(sce@assays) || !marker_reduction %in% Reductions(sce)) stop("缺少指定RNA assay或已有UMAP；本脚本不会重新计算嵌入。")
marker_meta <- sce[[]]
if (!all(c(cluster_col, sample_col) %in% names(marker_meta))) stop("缺少真实聚类/样本列；不要改用Idents猜分组。")
marker_cells <- colnames(sce)
marker_cluster_id <- as.character(marker_meta[marker_cells, cluster_col])
marker_samples <- as.character(marker_meta[marker_cells, sample_col])
if (anyNA(marker_cluster_id) || any(!nzchar(trimws(marker_cluster_id))) ||
    anyNA(marker_samples) || any(!nzchar(trimws(marker_samples)))) stop("cluster或样本身份有缺失。")
marker_cluster_levels <- unique(marker_cluster_id)
if (all(grepl("^[0-9]+$", marker_cluster_levels))) {
  marker_cluster_levels <- marker_cluster_levels[order(as.integer(marker_cluster_levels))]
} else marker_cluster_levels <- sort(marker_cluster_levels)
if (!marker_dot_engine %in% c("jjDotPlot", "ggplot2") || marker_feature_per_page < 1L || marker_feature_ncol < 1L) stop("绘图参数有误。")
marker_coordinates <- Embeddings(sce, reduction = marker_reduction)
if (!all(marker_cells %in% rownames(marker_coordinates)) || ncol(marker_coordinates) < 2L ||
    any(!is.finite(marker_coordinates[marker_cells, 1:2]))) stop("已有UMAP坐标缺失或不完整。")

# 只读取现有RNA data。不能在仅52个marker上重新NormalizeData，也不以SCT/integrated代替RNA。
marker_rna <- sce[[marker_assay]]
if (inherits(marker_rna, "Assay5")) {
  marker_data_layers <- Layers(marker_rna, search = "^data($|\\.)")
  if (!length(marker_data_layers)) stop("RNA没有data层；请先在完整RNA assay上完成LogNormalize。")
  marker_layer_matrices <- lapply(marker_data_layers, function(k) LayerData(marker_rna, layer = k, fast = FALSE))
} else {
  marker_data_layers <- "data"
  marker_layer_matrices <- list(if (packageVersion("SeuratObject") >= "5.0.0")
    GetAssayData(sce, assay = marker_assay, layer = "data") else GetAssayData(sce, assay = marker_assay, slot = "data"))
}
names(marker_layer_matrices) <- marker_data_layers
marker_layer_matrices <- lapply(marker_layer_matrices, function(z) z[, intersect(colnames(z), marker_cells), drop = FALSE])
marker_layer_matrices <- Filter(function(z) ncol(z) > 0L, marker_layer_matrices)
if (!length(marker_layer_matrices)) stop("RNA data层为空。")
marker_layer_cells <- unlist(lapply(marker_layer_matrices, colnames), use.names = FALSE)
if (anyDuplicated(marker_layer_cells) || !setequal(marker_layer_cells, marker_cells)) stop("RNA data层细胞重叠或没有覆盖当前sce；先核对表达层。")
marker_common_features <- Reduce(intersect, lapply(marker_layer_matrices, rownames))
marker_any_features <- Reduce(union, lapply(marker_layer_matrices, rownames))
marker_measured <- intersect(markers, marker_common_features)
if (!length(marker_measured)) stop("没有可共同评估的marker，请核对gene symbol、物种前缀和RNA表达层。")
marker_expression <- do.call(cbind, lapply(marker_layer_matrices, function(z) z[marker_measured, , drop = FALSE]))
marker_expression <- as(Matrix(marker_expression[, marker_cells, drop = FALSE], sparse = TRUE), "dgCMatrix")
if (any(!is.finite(marker_expression@x)) || any(marker_expression@x < 0) || any(!is.finite(expm1(marker_expression@x)))) stop("RNA data不是有效的非负log-normalized表达；不能使用scale.data或integrated值。")
marker_layer_audit <- data.frame(layer = names(marker_layer_matrices),
                                 n_cells = vapply(marker_layer_matrices, ncol, integer(1)), n_features = vapply(marker_layer_matrices, nrow, integer(1)))
rm(marker_layer_matrices, marker_rna)

##### 建立只用于绘图的副本，保留原sce与原嵌入 #####
marker_sce <- sce
DefaultAssay(marker_sce) <- marker_assay
marker_sce <- subset(marker_sce, features = marker_measured)
if (inherits(marker_sce[[marker_assay]], "Assay5")) {
  marker_sce <- JoinLayers(marker_sce, assay = marker_assay, layers = "data", new = "data")
}
marker_sce[[cluster_col]] <- factor(marker_cluster_id, levels = marker_cluster_levels)
# 不调用NormalizeData、ScaleData、RunUMAP、FindClusters，不写回sce或原对象文件。

##### 按cluster及样本记录marker表达与缺失状态 #####
marker_status <- markers_plot
marker_status$feature_status <- ifelse(markers %in% marker_common_features, "measured",
                                       ifelse(markers %in% marker_any_features, "incomplete_layer_coverage", "not_measured"))
marker_status$feature_status[marker_status$gene %in% marker_measured[
  Matrix::rowSums(marker_expression > 0) == 0]] <- "measured_zero"
marker_summary_rows <- list()
marker_summarize <- function(cells, cluster, scope, sample_id = NA_character_) {
  value <- pct <- mean_log <- rep(NA_real_, length(markers))
  if (length(cells)) for (j in which(markers %in% marker_measured)) {
    x <- as.numeric(marker_expression[markers[j], cells])
    value[j] <- mean(expm1(x))
    mean_log[j] <- mean(x)
    pct[j] <- 100 * mean(x > 0)
  }
  data.frame(scope = scope, sample_id = sample_id, cluster_id = cluster, marker_group = markers_plot$cluster,
             gene = markers, n_cells = length(cells), feature_status = marker_status$feature_status,
             mean_normalized = value, mean_log_normalized = mean_log, detected_percent = pct,
             stringsAsFactors = FALSE)
}
for (id in marker_cluster_levels) {
  cc <- marker_cells[marker_cluster_id == id]
  marker_summary_rows[[length(marker_summary_rows) + 1L]] <- marker_summarize(cc, id, "cluster")
  if (marker_include_sample_summary) for (sm in unique(marker_samples)) {
    marker_summary_rows[[length(marker_summary_rows) + 1L]] <- marker_summarize(
      cc[marker_samples[match(cc, marker_cells)] == sm], id, "sample", sm)
  }
}
marker_summary <- do.call(rbind, marker_summary_rows)
marker_summary$dot_z <- NA_real_
for (g in markers) {
  ii <- which(marker_summary$scope == "cluster" & marker_summary$gene == g)
  x <- marker_summary$mean_normalized[ii]
  if (sum(is.finite(x)) >= 2L && sd(x, na.rm = TRUE) > 1e-10) {
    marker_summary$dot_z[ii] <- pmax(marker_z_limits[1], pmin(marker_z_limits[2], (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
  }
}
marker_dot_data <- marker_summary[marker_summary$scope == "cluster", , drop = FALSE]

##### 绘图函数：每图一次打印，常数/缺失marker保留状态 #####
marker_theme <- theme_bw(base_size = 10, base_family = "sans") + theme(
  panel.grid = element_blank(), panel.border = element_rect(colour = "#454545", linewidth = 0.4),
  strip.background = element_rect(fill = "#E7EEF4", colour = "#454545", linewidth = 0.35),
  strip.text = element_text(face = "bold", size = 9), plot.title = element_text(face = "bold", size = 12))
marker_dot_scale <- function() scale_fill_gradient2(low = marker_z_colors[1], mid = marker_z_colors[2],
                                                    high = marker_z_colors[3], midpoint = 0, limits = marker_z_limits, na.value = marker_na_color,
                                                    name = "Mean expression\nz score")
marker_dot_ggplot <- function(groups, title) {
  z <- marker_dot_data[marker_dot_data$marker_group %in% groups, ]
  z$marker_group <- factor(z$marker_group, levels = groups)
  z$gene <- factor(z$gene, levels = markers_plot$gene[markers_plot$cluster %in% groups])
  z$cluster_id <- factor(z$cluster_id, levels = rev(marker_cluster_levels))
  p <- ggplot(z, aes(gene, cluster_id)) +
    geom_blank() + # 先训练完整marker次序，缺失基因不会被挪到分组末尾。
    geom_point(data = z[is.finite(z$detected_percent) & is.finite(z$dot_z), ], aes(size = detected_percent, fill = dot_z),
               shape = 21, colour = "#454545", stroke = 0.25) +
    facet_grid(~ marker_group, scales = "free_x", space = "free_x") +
    marker_dot_scale() + scale_size_continuous(range = c(0, 5), limits = c(0, 100), name = "Detected (%)") +
    scale_x_discrete(drop = TRUE) + scale_y_discrete(drop = FALSE) +
    labs(title = NULL, x = NULL, y = "Cluster") + marker_theme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "right")
  missing <- z[!is.finite(z$detected_percent), ]
  constant <- z[is.finite(z$detected_percent) & !is.finite(z$dot_z), ]
  if (nrow(missing)) p <- p + geom_point(data = missing, aes(shape = "Not measured"), size = 1.7, colour = marker_na_color, stroke = 0.6)
  if (nrow(constant)) p <- p + geom_point(data = constant, aes(size = detected_percent, shape = "z score NA"),
                                          fill = marker_na_color, colour = "#454545", stroke = 0.25)
  if (nrow(missing) || nrow(constant)) p <- p + scale_shape_manual(values = c("Not measured" = 4, "z score NA" = 21), name = "Status")
  p
}
marker_build_jj <- function(groups, title) {
  # jjAnno内部可能print；在空设备中构造，避免同一张图在目标PDF中重复。
  grDevices::pdf(NULL, width = 12, height = 9)
  marker_build_device <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(which = marker_build_device), add = TRUE)
  z <- markers_plot[markers_plot$cluster %in% groups, , drop = FALSE]
  z <- z[order(match(z$cluster, groups)), , drop = FALSE]
  z$cluster <- factor(z$cluster, levels = groups)
  p <- scRNAtoolVis::jjDotPlot(object = marker_sce, assay = marker_assay, slot = "data",
                               markerGene = z, anno = TRUE, id = cluster_col, gene.order = z$gene,
                               cluster.order = rev(marker_cluster_levels), ytree = marker_show_tree, xtree = FALSE,
                               tree.pos = "right", textSize = 11, base_size = 11, dot.min = 0, dot.max = 5,
                               dot.col = marker_z_colors, midpoint = 0, scale = TRUE,
                               bar.legendTitle = "Mean expression\nz score", point.legendTitle = "Detected (%)",
                               col.min = marker_z_limits[1], col.max = marker_z_limits[2], plot.margin = c(2.8, 0.4, 0.4, 0.4))
  suppressMessages(p + marker_dot_scale() +
                     scale_size_continuous(range = c(0, 5), limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), name = "Detected (%)") +
                     guides(fill = guide_colorbar(direction = "vertical", barheight = grid::unit(3, "cm"), barwidth = grid::unit(0.35, "cm"))) +
                     labs(title = NULL) + theme(panel.grid = element_blank(), legend.direction = "vertical"))
}
marker_strip <- function(p, gene) {
  p$data$marker_title <- gene
  p + facet_wrap(~ marker_title) + labs(title = NULL) +
    theme(strip.background = element_rect(fill = "#E7EEF4", colour = "#454545", linewidth = 0.35),
          strip.text = element_text(face = "bold", size = 10), legend.text = element_text(size = 8),
          legend.title = element_text(size = 8), plot.margin = margin(3, 3, 3, 3))
}
marker_placeholder <- function(gene) {
  dd <- data.frame(x = marker_coordinates[marker_cells, 1], y = marker_coordinates[marker_cells, 2])
  p <- ggplot(dd, aes(x, y)) + geom_point(colour = "#E6E9EC", size = 0.2) +
    annotate("label", x = mean(range(dd$x)), y = mean(range(dd$y)), label = "NA\nnot measured", size = 3) +
    coord_fixed(xlim = range(dd$x), ylim = range(dd$y)) + theme_void(base_family = "sans")
  marker_strip(p, gene)
}

##### 保存PDF：cluster参考、三张DotPlot及全部marker分页FeaturePlot #####
dir.create(marker_outdir, recursive = TRUE, showWarnings = FALSE)
marker_panel_audit <- list()
marker_write_pdf <- function() {
  grDevices::pdf(file.path(marker_outdir, "plots.pdf"), width = 12, height = 9, onefile = TRUE, useDingbats = FALSE)
  marker_output_device <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(which = marker_output_device), add = TRUE)
  print(DimPlot(marker_sce, reduction = marker_reduction, group.by = cluster_col,
                label = TRUE, repel = TRUE, raster = marker_raster, pt.size = marker_point_size) +
          labs(title = "T cell clusters") + theme_classic(base_size = 12) + NoLegend())
  for (title in names(marker_panels)) {
    groups <- marker_panels[[title]]
    genes <- markers_plot$gene[markers_plot$cluster %in% groups]
    z <- marker_dot_data[marker_dot_data$gene %in% genes, ]
    suitable <- all(is.finite(z$dot_z)) && length(marker_cluster_levels) >= 2L
    reason <- "requested ggplot2"
    p <- NULL
    if (marker_dot_engine == "jjDotPlot" && suitable && requireNamespace("scRNAtoolVis", quietly = TRUE)) {
      p <- tryCatch(marker_build_jj(groups, title), error = function(e) {
        reason <<- paste("jjDotPlot failed:", conditionMessage(e)); NULL
      })
    } else if (marker_dot_engine == "jjDotPlot") reason <- if (!suitable) "missing or constant marker; explicit NA display" else "scRNAtoolVis unavailable"
    used <- if (is.null(p)) "ggplot2" else "jjDotPlot"
    if (is.null(p)) p <- marker_dot_ggplot(groups, title)
    marker_panel_audit[[title]] <<- data.frame(panel = title, engine = used, status = if (used == "jjDotPlot") "ok" else reason)
    if (used != marker_dot_engine) message(title, "：", reason, "；使用同口径ggplot2。")
    print(wrap_plots(list(p), ncol = 1) + plot_annotation(title = title))
  }
  pages <- split(markers, ceiling(seq_along(markers) / marker_feature_per_page))
  for (page in seq_along(pages)) {
    gg <- pages[[page]]
    available <- intersect(gg, marker_measured)
    plots <- if (length(available)) FeaturePlot(marker_sce, features = available, reduction = marker_reduction,
                                                cols = marker_feature_colors, order = TRUE, min.cutoff = 0, keep.scale = "feature",
                                                pt.size = marker_point_size, raster = marker_raster, combine = FALSE, slot = "data") else list()
    names(plots) <- available
    shown <- lapply(gg, function(g) {
      if (!g %in% available) return(marker_placeholder(g))
      p <- plots[[g]] + NoAxes() + coord_fixed(xlim = range(marker_coordinates[, 1]), ylim = range(marker_coordinates[, 2]))
      # Seurat对全零基因可能使用色标中点；显式保持浅灰及真实0刻度。
      if (marker_status$feature_status[match(g, marker_status$gene)] == "measured_zero") {
        p <- suppressMessages(p + scale_colour_gradient(low = marker_feature_colors[1],
                                                        high = marker_feature_colors[1], breaks = 0, name = NULL))
      }
      marker_strip(p, g)
    })
    print(wrap_plots(shown, ncol = marker_feature_ncol) + plot_annotation(title = "T cell markers"))
  }
}
marker_write_pdf()

##### 保存表达表和运行记录 #####
write.csv(marker_summary, file.path(marker_outdir, "summary.csv"), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
marker_panel_audit <- do.call(rbind, marker_panel_audit)
marker_info <- capture.output({
  cat("Recorded:", as.character(Sys.time()), "\n")
  cat("Object: existing sce; original object and annotations are unchanged.\n")
  cat("Assay:", marker_assay, "data; inherited log-normalization, no marker-only normalization.\n")
  cat("Cluster column:", cluster_col, "Reduction:", marker_reduction, "\n")
  cat("Clusters:", paste(marker_cluster_levels, collapse = ", "), "Cells:", length(marker_cells), "Markers:", length(markers), "\n")
  cat("Dot color: per-gene z score of mean(expm1(RNA data)), clipped to", paste(marker_z_limits, collapse = ", "), "\n")
  cat("Dot size: percentage RNA data > 0; constant gene z is NA, observed zero is retained.\n")
  cat("FeaturePlot: existing RNA data, each gene has its own displayed color scale.\n")
  cat("Requested engine:", marker_dot_engine, "Requested tree:", marker_show_tree, "Raster points:", marker_raster, "\n")
  print(marker_layer_audit)
  print(marker_panel_audit, row.names = FALSE)
  print(marker_status[marker_status$feature_status != "measured", ], row.names = FALSE)
  cat("Cells by sample and cluster:\n"); print(table(marker_samples, marker_cluster_id))
  print(vapply(c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork"), function(p) as.character(packageVersion(p)), character(1)))
  if (requireNamespace("scRNAtoolVis", quietly = TRUE)) cat("scRNAtoolVis:", as.character(packageVersion("scRNAtoolVis")), "\n")
  print(sessionInfo())
})
writeLines(marker_info, file.path(marker_outdir, "info.txt"), useBytes = TRUE)
message("完成：", marker_outdir, "；请将plots.pdf发回，并保留summary.csv/info.txt供核对。")

