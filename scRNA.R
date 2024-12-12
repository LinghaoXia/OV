##直接对存储文件夹读取seurat包会自动识别三个文件并进行构建
library(dplyr)
library(Seurat)
library(patchwork)

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)

##### 整理数据 #####
##virus_3D
setwd('~/rawdata/scRNA_virus/virus_3D')
folders=list.files('./')
folders
library(Seurat)
scList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder,
                     min.cells = 3, min.features = 200)
})

virus_3D <- merge(scList[[1]], 
                  y = c(scList[[2]],scList[[3]],scList[[4]]),
                  add.cell.ids = c("Vehicle_3D_1","Vehicle_3D_2","VG161_3D_1","VG161_3D_2"), 
                  project = "virus_3D")

save(virus_3D,file="virus_3D.RDS")

##virus_5D
setwd('~/rawdata/scRNA_virus/virus_5D')
folders=list.files('./')
folders
library(Seurat)
scList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder,
                     min.cells = 3, min.features = 200)
})

virus_5D <- merge(scList[[1]], 
                  y = c(scList[[2]],scList[[3]],scList[[4]],scList[[5]],scList[[6]]),
                  add.cell.ids = c("Vehicle_5D_1","Vehicle_5D_2","VG161_5D_L1","VG161_5D_L2","VG161_5D_R1","VG161_5D_R2"), 
                  project = "virus_5D")

save(virus_5D,file="virus_5D.RDS")

##virus_PD1
setwd('~/rawdata/scRNA_virus/virus_PD1')
folders=list.files('./')
folders
library(Seurat)
scList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder,
                     min.cells = 3, min.features = 200)
})

virus_PD1 <- merge(scList[[1]], 
                   y = c(scList[[2]],scList[[3]],scList[[4]],scList[[5]],scList[[6]]),
                   add.cell.ids = c("PD1_5D_1","PD1_5D_2","PD1_VG161_5D_L1","PD1_VG161_5D_L2","PD1_VG161_5D_R1","PD1_VG161_5D_R2"), 
                   project = "virus_PD1")

save(virus_PD1,file="virus_PD1.RDS")


##### 质控降维 #####
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
##virus_3D
load("~/rawdata/scRNA_virus/virus_PD1/virus_PD1.RDS")
sce <- virus_PD1
###########SCTransform V2标准化质控降维
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
preQC <- VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,
                 group.by = "orig.ident", 
                 pt.size = 0)
#细胞基因数量与mRNA、核糖体基因数量的相关性
plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

sce <- subset(sce, subset = nFeature_RNA > 200  & nFeature_RNA < 5000 & percent.mt < 15)
postQC <- VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,
                  group.by = "orig.ident", 
                  pt.size = 0)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
#top10 <- head(VariableFeatures(scRNA), 10)
#plot1 <- VariableFeaturePlot(scRNA) 
#LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) 
#如下图，横坐标是某基因在所有细胞中的平均表达值，纵坐标是此基因的方差;红点即为高变基因（2000个）
###########Harmony去批次后降维
sce <- SCTransform(sce,vst.flavor = "v2", verbose = FALSE, method = "glmGamPoi",vars.to.regress = "percent.mt")
sce=RunPCA(sce,assay="SCT",verbose = FALSE)
##去批次
sce=RunHarmony(sce,group.by.vars="orig.ident",assay.use="SCT", plot_convergence = TRUE,max.iter.harmony =50 )
##最佳PC数量
pct <- sce [["harmony"]]@stdev / sum( sce [["harmony"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2
pcs <- min(co1, co2)
pcs
#获取了最佳PC用于UMAP和FindNeighbors
bestpc=1:pcs
sce<- sce %>% RunUMAP(reduction = "harmony", dims = bestpc) %>% 
  FindNeighbors(reduction = "harmony", dims = bestpc)
sce=FindClusters(sce,resolution = 0.5)#需要对粒度进行调整
#保存结果
#sce@meta.data$group <- ifelse (grepl("VG161",sce@meta.data$orig.ident),'VG161','Vehicle')
#sce@meta.data$group <- ifelse (grepl("VG161",sce@meta.data$orig.ident),ifelse(grepl("_L",sce@meta.data$orig.ident),"VG161_L","VG161_R"),'Vehicle')
sce@meta.data$group <- ifelse (grepl("VG161",sce@meta.data$orig.ident),ifelse(grepl("_L",sce@meta.data$orig.ident),"PD1_VG161_L","PD1_VG161_R"),'PD1')

pdf(paste0(outdir,"/","01-",gene,"-cluster_sample.pdf"),height=10,width=6)
DimPlot(sce, reduction = "umap", group.by = "group",label = TRUE,repel = T, pt.size = .1)
dev.off()

pdf(paste0(outdir,"/","01-",gene,"-cluster.pdf"),height=10,width=6)
DimPlot(sce, reduction = "umap",label = TRUE,repel = T, pt.size = .1)
dev.off()

save(sce,file = "~/rawdata/scRNA_virus/virus_PD1/virus_PD1_cluster.RData")

##### 自动注释-singleR #####
library(SingleR)
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

load("~/rawdata/scRNA_virus/virus_3D/virus_3D_cluster.RData")
load("~/rawdata/scRNA_virus/refdata.RData")
# refdata <- HumanPrimaryCellAtlasData()
# save(refdata,file = "~/rawdata/scRNA_virus/refdata.RData")
#参考数据库，等待时间较长。建议下载成功后，储存为Rdata，以后方便使用。
testdata <- GetAssayData(sce, slot="data")
clusters <- sce@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    # 注意此时labels参数为 refdata$label.main，与下一节亚类再注释时的设置不同
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
#如下为singleR的细胞cluster鉴定结果。
celltype
#结合上述结果，给scRNA增添celltype注释信息
sce@meta.data$celltype = "NA"
#先新增列celltype，值均为NA，然后利用下一行代码循环填充
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
DimPlot(sce, group.by="celltype", label=T, label.size=5, reduction='umap')


##### 自动注释-cellassign #####
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleCellExperiment)
library(cellassign)
library(scran)
library(dplyr)
rm(list=ls())
gc()

load("~/rawdata/scRNA_virus/virus_3D/virus_3D_cluster.RData")
data(example_TME_markers)
#自己定义marker
# marker_gene_list <- list(
#   "Naive CD4+ T" = c("IL7R", "CCR7"),
#   "CD14+ Mono" = c("CD14", "LYZ"),
#   "Memory CD4+" = c("IL7R", "S100A4"),
#   "B" = c("MS4A1"),
#   "CD8+ T" = c("CD8A"),
#   "FCGR3A+ Mono" = c("FCGR3A","MS4A7"),
#   "NK" = c("GNLY","NKG7"),
#   "DC" = c("FCER1A","CST3"),
#   "Platelet" = c("PPBP")
# )
# marker_gene_mat <- marker_list_to_mat(marker_gene_list)
# head(marker_gene_mat)

marker_gene_mat <- marker_list_to_mat(example_TME_markers$symbol)#include_other=FALSE
#marker_gene_mat <- marker_gene_mat[rownames(marker_gene_mat) !="CLEC14A",]

count<-as.matrix(sce@assays$SCT@counts)
# sizeFactors的准备
cell_anns <- data.frame(Barcode =  names(Idents(sce)),celltype=Idents(sce),samples=sce@meta.data$orig.ident)
rownames(cell_anns) <- names(Idents(sce))
sceset <- SingleCellExperiment(assays = list(counts = as.matrix(count)),colData=cell_anns)
str(sceset)
rm(count)
gc()

qclust <- quickCluster(sceset, min.size = 30)
sceset <- computeSumFactors(sceset, sizes = 15, clusters = qclust)
s <- sizeFactors(sceset)


fit  <- cellassign(exprs_obj = sceset[as.vector(rownames(marker_gene_mat)),], 
                   marker_gene_info = marker_gene_mat,
                   s = s,
                   learning_rate = 1e-2,
                   shrinkage = TRUE,
                   verbose = FALSE)
print(head(celltypes(fit)))

library(ggplot2)
sce.filter <- sce[,colnames(mat)]
sce.filter@meta.data$cellassign <-celltypes(fit)
p1 <- DimPlot(sce.filter, label = T,group.by = "cellassign")
##点图，看细胞的基因特征
feature <- c("MS4A1","FCGR3A","MS4A7","FCER1A","CST3","PPBP","CD8A","CCR7","IL7R","S100A4","CD14", "LYZ","GNLY","NKG7")
p2 <- DotPlot(obj.filter,features = feature,group.by = "cellassign")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1+p2

##### 半自动注释-SCINA #####
library(SCINA)
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

load("~/rawdata/scRNA_virus/virus_PD1/virus_PD1_cluster.RData")
##详细
# read.gmt=function(filename){
#   if(! file.exists(filename)) stop('File ',filename,' not available\n')
#   dat=readLines(filename)
#   n=length(dat)
#   res=list(genesets=vector(mode = "list", length = n),geneset.names=vector(mode = "character", length = n),geneset.descriptions=vector(mode = "character", length = n))
#   for(i in 1:n){
#     s=strsplit(dat[i],'\t')[[1]]
#     res$genesets[[i]]=s[-c(1:2)]
#     res$geneset.names[i]=s[1]
#     res$geneset.descriptions[i]=s[2]
#   }
#   names(res$genesets)=res$geneset.names
#   res
# }
# markers <- read.gmt('~/rawdata/scRNA_virus/30_ImmuneCells_StromalCells.gmt')

##模糊
#细胞类型marker列表
Epithelial <- c("EPCAM","CDH1","KRT19")
Endothelial <- c("CDH5","VWF")
Fibroblast <- c("COL1A1","DCN","LUM")
Mast <- c("TPSAB1","MS4A2")
Bcell <- c("CD79A","CD19","MS4A1")
Tcell <- c("CD3D","CD3E")
Myeloid <- c("FCER1G","CD68")

markers <- list(Tcell= Tcell,
                Myeloid= Myeloid,
                Mast=Mast,
                Bcell=Bcell,
                Endothelial=Endothelial,
                Epithelial=Epithelial,
                Fibroblast=Fibroblast)

expr.data <- as.matrix(GetAssayData(sce))

predictions.scina = SCINA(exp = expr.data, signatures = markers, #详细markers$genesets
                          allow_unknown = FALSE,
                          rm_overlap = FALSE)

rm(expr.data)
gc()
sce$celltype <- predictions.scina$cell_labels

pdf(paste0(outdir,"/","02-",gene,"-anno.pdf"),height=20,width=12)
DimPlot(sce, reduction = "umap",
                   group.by = "celltype",
                   label = TRUE, label.size = 3, repel = TRUE)
dev.off()
save(sce,file = "~/rawdata/scRNA_virus/virus_PD1/virus_PD1_anno.RData")



##### 病毒感染 #####
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)

load("~/rawdata/scRNA_virus/virus_PD1/virus_PD1_anno.RData")
##病毒早期基因比例
early_viral_genes <- c("VG161-UL4","VG161-UL5","VG161-UL8","VG161-UL9","VG161-UL11","VG161-UL12","VG161-UL13","VG161-UL14",
                       "VG161-UL22","VG161-UL23","VG161-UL29",
                       "VG161-UL30","VG161-UL39","VG161-UL40","VG161-UL42","VG161-UL50","VG161-UL52","VG161-UL53")
viral_expr <- FetchData(sce, vars = early_viral_genes)
viral_expr_total <- rowSums(viral_expr)
# 计算病毒基因总表达量占该细胞所有基因表达量的比例
total_expr <- rowSums(FetchData(sce,vars = rownames(sce)))  # 获取每个细胞所有基因的表达量
viral_expression_ratio <- viral_expr_total / total_expr * 100  # 计算百分比
# 如果病毒基因表达超过 0.01%，则标记为感染细胞
threshold <- 0.01  # 定义阈值为 0.01%
infected_cells <- viral_expression_ratio > threshold
sce@meta.data$infected <- ifelse(grepl("VG161",sce@meta.data$orig.ident),ifelse(infected_cells, "Infected", "Bystander"),"Naive")

##病毒总转录本
VirTranscript <- function(obj,viral_genes){
  viral_counts <- FetchData(obj, vars = viral_genes, slot = 'counts')
  viral_counts_total <- rowSums(viral_counts)
  total_counts <- rowSums(FetchData(obj,slot = 'counts',vars = rownames(obj))) 
  scaled_viral_counts <- log2((viral_counts_total / total_counts) * 10000+1)
  obj$VG161_transcript <- scaled_viral_counts
  return(obj)
}

viral_genes <-  rownames(sce@assays$SCT)[grep("VG161", rownames(sce@assays$SCT))]
sce <- VirTranscript(sce,viral_genes)

pdf(paste0(outdir,"/","03-",gene,"-infected.pdf"),height=20,width=12)
DimPlot(sce, group.by = "infected",split.by = "group", label = TRUE)
DimPlot(sce, group.by = "infected", label = TRUE)+DimPlot(sce, group.by = "celltype", label = TRUE)
FeaturePlot(sce, features = "VG161_transcript", cols = c("lightgrey", "red"))
dev.off()

save(sce,file = "~/rawdata/scRNA_virus/virus_PD1/virus_PD1_anno.RData")


##### 细胞比例 #####
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)


load("~/rawdata/scRNA_virus/virus_3D/virus_3D_anno.RData")
table(sce$group)#查看各组细胞数
prop.table(table(sce$infected))
table(sce$infected, sce$group)#各组不同细胞群细胞数

Cellratio <- prop.table(table(sce$infected, sce$group), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))

##柱状图
library(ggplot2)
library(ggbreak)
pdf(paste0(outdir,"/","03-",gene,"-infected_ratio.pdf"),height=8,width=12)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))#+
  # geom_text(aes(x = Var2, y = Freq, label = scales::percent(Freq)), 
  #           position = position_stack(vjust = 0.5), size = 4, colour = "white")

group <- unique(Cellratio$Var2)
Cellratio_1 <- Cellratio[which(Cellratio$Var2==group[1]),]
Cellratio_2 <- Cellratio[which(Cellratio$Var2==group[3]),]
top<-function(x){
  return(mean(x)+sd(x)/sqrt(length(x)))
}
bottom<-function(x){
  return(mean(x)-sd(x)/sqrt(length(x)))
}
p1 <- ggplot(data=Cellratio_1,aes(x=Var1,y=Freq,fill=Var1))+
  stat_summary(geom = "bar",fun = "mean",
               position = position_dodge(0.9))+
  stat_summary(geom = "errorbar",
               fun.min = bottom,
               fun.max = top,
               position = position_dodge(0.9),
               width=0.2)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Celltype",y="Proportion")+
  geom_point(data=Cellratio_1,aes(Var1,Freq),size=3,pch=19)+
  theme(
    axis.text.x.bottom = element_text(size = 14,hjust = 1, vjust = 1, angle = 45)
  ) +ggtitle(group[1])

p2 <- ggplot(data=Cellratio_2,aes(x=Var1,y=Freq,fill=Var1))+
  stat_summary(geom = "bar",fun = "mean",
               position = position_dodge(0.9))+
  stat_summary(geom = "errorbar",
               fun.min = bottom,
               fun.max = top,
               position = position_dodge(0.9),
               width=0.2)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  scale_y_break(c(0.005,0.1),#截断位置及范围
                space = 0.3,#间距大小
                scales = 1.5)+#上下显示比例，大于1上面比例大，小于1下面比例大
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Celltype",y="Proportion")+
  geom_point(data=Cellratio_2,aes(Var1,Freq),size=3,pch=19)+
  ggtitle(group[3])+
  theme(axis.text.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.title.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())
p1+p2
dev.off()


##比较
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]

#添加分组信息
samples <- data.frame(sce$orig.ident, sce$group)#创建数据框

rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$group <- samples[rownames(cellper),'group']#R添加列
pplist = list()
sce_groups = unique(sce$infected)
library(ggplot2)
library(dplyr)
library(ggpubr)
for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group',group_)))#选择一组数据
  colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
  cellper_$percent = as.numeric(cellper_$percent)#数值型数据
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#上下分位数
  print(group_)
  print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  ###组间t检验分析
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list( c("GM", "BM") )
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}


