##直接对存储文件夹读取seurat包会自动识别三个文件并进行构建
library(dplyr)
library(Seurat)
library(patchwork)

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

##### 整理数据 #####
##injected
setwd('~/rawdata/L2HGDH_scRNA/injected')
folders=list.files('./')
folders
library(Seurat)
scList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder,
                     min.cells = 3, min.features = 200)
})

sce <- merge(scList[[1]], 
                  y = c(scList[[2]],scList[[3]],scList[[4]],scList[[5]],scList[[6]],scList[[7]],scList[[8]],scList[[9]]),
                  add.cell.ids = c("PBS_1","PBS_2","PBS_3","S_1","S_2","S_3","VG21_1","VG21_2","VG21_3"), 
                  project = "injected")

save(sce,file="injected.RDS")

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
load("~/rawdata/L2HGDH_scRNA/injected/injected.RDS")
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
sce@meta.data$group <- ifelse (grepl("VG161",sce@meta.data$orig.ident),"VG21",ifelse(grepl("PBS",sce@meta.data$orig.ident),"PBS","S"))

pdf(paste0(outdir,"/","01-",gene,"-cluster_sample.pdf"),height=10,width=6)
DimPlot(sce, reduction = "umap", group.by = "group",label = TRUE,repel = T, pt.size = .1)
dev.off()

pdf(paste0(outdir,"/","01-",gene,"-cluster.pdf"),height=10,width=6)
DimPlot(sce, reduction = "umap",label = TRUE,repel = T, pt.size = .1)
dev.off()

save(sce,file = "~/rawdata/L2HGDH_scRNA/injected/injected_cluster.RData")

##### 自动注释-singleR #####
library(SingleR)
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

load("~/rawdata/L2HGDH_scRNA/injected/injected_cluster.RData")
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
