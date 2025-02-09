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

##uninjected
setwd('~/rawdata/L2HGDH_scRNA/uninjected')
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

save(sce,file="uninjected.RDS")

##blood
setwd('~/rawdata/L2HGDH_scRNA/blood')
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

save(sce,file="blood.RDS")


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

setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)
##virus_3D
load("~/rawdata/L2HGDH_scRNA/uninjected/uninjected.RDS")
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
sce@meta.data$group <- ifelse (grepl("VG21",sce@meta.data$orig.ident),'VG21',ifelse(grepl("PBS",sce@meta.data$orig.ident),"PBS","S"))

pdf(paste0(outdir,"/","01-",gene,"-cluster_sample.pdf"),height=10,width=6)
DimPlot(sce, reduction = "umap", group.by = "group",label = TRUE,repel = T, pt.size = .1)
dev.off()

pdf(paste0(outdir,"/","01-",gene,"-cluster.pdf"),height=10,width=6)
DimPlot(sce, reduction = "umap",label = TRUE,repel = T, pt.size = .1)
dev.off()

save(sce,file = "~/rawdata/L2HGDH_scRNA/uninjected/uninjected_cluster.RData")


##### 半自动注释-SCINA #####
library(SCINA)
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

load("~/rawdata/L2HGDH_scRNA/injected/injected_cluster.RData")
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
DC <- c("Batf3","Irf8","Itgax","H2-DMb1","H2-Ab1")
Macrophage <- c("Itgam","Adgre1","Cd68")
Bcell <- c("Ms4a1","Cd19","Igalpha")
Tcell <- c("Cd3d","Cd3e","Cd81","Rora","Stat3","Tgfb1")
Monocyte <- c("Itgam","Cd14","Cd16")


markers <- list(Tcell= Tcell,
                DC= DC,
                Macrophage=Macrophage,
                Bcell=Bcell,
                Monocyte=Monocyte)

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
save(sce,file = "~/rawdata/L2HGDH_scRNA/injected/injected_anno.RData")

##### 细胞比例 #####
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)


load("~/rawdata/L2HGDH_scRNA/injected/injected_anno.RData")
table(sce$group)#查看各组细胞数
prop.table(table(sce$celltype))
table(sce$celltype, sce$group)#各组不同细胞群细胞数


pdf(paste0(outdir,"/","02-",gene,"-anno_ratio.pdf"),height=8,width=12)
##柱状图
Cellratio <- prop.table(table(celltype=sce$celltype, group=sce$group), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$celltype))
library(ggplot2)
library(ggbreak)
ggplot(Cellratio) + 
  geom_bar(aes(x =group, y= Freq, fill = celltype),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))


##比较散点图
Cellratio <- prop.table(table(celltype=sce$celltype, group=sce$orig.ident), margin = 2)#计算各组样本不同细胞群比例
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
sce_groups = unique(sce$celltype)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
for(group_ in sce_groups){
  cellper_  = cellper[,c('sample','group',group_)]
  colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
  cellper_$group <- factor(cellper_$group , levels =c("PBS","VG21","S"))
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
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  pplist[[group_]] = pp1
}

library(cowplot)
plot_grid(pplist[[1]],
          pplist[[2]],
          pplist[[3]],
          pplist[[4]],
          pplist[[5]])

dev.off()



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
##virus_3D
load("~/rawdata/L2HGDH_scRNA/injected/injected_anno.RData")

###提取亚群
sce <- subset(sce, celltype=="Tcell")
###降维
DefaultAssay(sce) <- "RNA"
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 

sce <- FindNeighbors(sce, dims = 1:10)
sce <- FindClusters(sce, resolution = 0.5)
sce <- RunUMAP(sce, dims = 1:10)
pdf(paste0(outdir,"/","03-",gene,"-Tcell_cluster.pdf"),height=10,width=6)
DimPlot(sce, reduction = 'umap', group.by = 'seurat_clusters',label = TRUE, pt.size = 0.5)
dev.off()

save(sce,file = "~/rawdata/L2HGDH_scRNA/injected/injected_tcell.RData")


##### T细胞亚群注释 #####
library(Seurat) 
library(ggplot2)
library(dplyr)
library(scRNAtoolVis)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

load("~/rawdata/L2HGDH_scRNA/uninjected/uninjected_tcell.RData")

###### genes to check ######
# DefaultAssay(sce) <- "RNA"
markers <- c("Cd4",  #cd4+ Tcells            
             "Cd8a","Cd8b1",  #cd8+ Tcells 
             "Ccr7", "Sell","Cd44", #T naive CD44-CCR7+SELL(CD62L)+ 
             "Fasl","Fas",#T effecctor             
             "Il7r")#T memory   CCR7和CD62在Tcm和Tem上不同        
             
markers_plot <- data.frame(cluster = c(rep("CD4 Tcell",1),                                  
                                       rep("CD8 Tcell",2),                                  
                                       rep("naive Tcell",3),                                  
                                       rep("effector Tcell",2),                                  
                                       rep("memory Tcell",1)),                                
                           gene = markers)

pdf(paste0(outdir,"/","03-",gene,"-Tcell_markers.pdf"),height=10,width=6)
jjDotPlot(object = sce,          
          markerGene = markers_plot,          
          anno = T,          
          id = 'seurat_clusters',          
          textSize = 10,          
          base_size= 10,          
          plot.margin = c(4,1.5,1.5,1.5))
dev.off()


###### annotation ######
# celltype <- c(    "0"="CD8+ Tcell",
#                   "1"="CD8+ Tcell",
#                   "2"="CD4+ Tcell",
#                   "3"="other Tcell",
#                   "4"="CD8+ Tcell",
#                   "5"="CD4+ Tcell",
#                   "6"="CD8+ Tcell",
#                   "7"="CD4+ Tcell",
#                   "8"="CD4+ Tcell",
#                   "9"="CD8+ Tcell",
#                   "10" = "CD8+ Tcell",
#                   "11" = "CD8+ Tcell",
#                   "12" = "CD4+ Tcell",
#                   "13" = "CD8+ Tcell",
#                   "14" = "other Tcell",
#                   "15" = "CD8+ Tcell",
#                   "16" = "CD4+ Tcell",
#                   "17" = "CD8+ Tcell") ##injected
celltype <- c(    "0"="CD8+ Tcell",
                  "1"="CD4+ Tcell",
                  "2"="CD4+ Tcell",
                  "3"="CD8+ Tcell",
                  "4"="other Tcell",
                  "5"="CD4+ Tcell",
                  "6"="CD8+ Tcell",
                  "7"="CD8+ Tcell",
                  "8"="CD4+ Tcell",
                  "9"="CD8+ Tcell",
                  "10" = "CD4+ Tcell",
                  "11" = "CD8+ Tcell",
                  "12" = "CD8+ Tcell",
                  "13" = "other Tcell",
                  "14" = "other Tcell",
                  "15" = "CD8+ Tcell",
                  "16" = "CD4+ Tcell")##uninjected
sce@meta.data$t_type <- celltype[sce@meta.data$seurat_clusters]

# tcelltype <- c(   "0"="effector Tcell",
#                   "1"="naive Tcell",
#                   "2"="memory Tcell",
#                   "3"="naive Tcell",
#                   "4"="effector Tcell",
#                   "5"="effector Tcell",
#                   "6"="effector Tcell",
#                   "7"="naive Tcell",
#                   "8"="memory Tcell",
#                   "9"="memory Tcell",
#                   "10" = "naive Tcell",
#                   "11" = "memory Tcell",
#                   "12" = "effector Tcell",
#                   "13" = "effector Tcell",
#                   "14" = "memory Tcell",
#                   "15" = "effector Tcell",
#                   "16" = "naive Tcell",
#                   "17" = "memory Tcell")##injected
tcelltype <- c(   "0"="naive Tcell",
                  "1"="memory Tcell",
                  "2"="effector Tcell",
                  "3"="memory Tcell",
                  "4"="naive Tcell",
                  "5"="memory Tcell",
                  "6"="effector Tcell",
                  "7"="effector Tcell",
                  "8"="memory Tcell",
                  "9"="naive Tcell",
                  "10" = "effector Tcell",
                  "11" = "effector Tcell",
                  "12" = "effector Tcell",
                  "13" = "memory Tcell",
                  "14" = "memory Tcell",
                  "15" = "effector Tcell",
                  "16" = "memory Tcell")##uninjjected
sce@meta.data$t_function <- tcelltype[sce@meta.data$seurat_clusters]
save(sce,file = "~/rawdata/L2HGDH_scRNA/uninjected/uninjected_tcell.RData")

###### plotting ######
pdf(paste0(outdir,"/","03-",gene,"-Tcell_ratio.pdf"),height=8,width=12)

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
  cellper_$group <- factor(cellper_$group , levels =c("PBS","VG21","S"))
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
  cellper_$group <- factor(cellper_$group , levels =c("PBS","VG21","S"))
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
  cellper_$group <- factor(cellper_$group , levels =c("PBS","VG21","S"))
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


##### 巨噬细胞细胞亚群聚类 #####
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
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)
##virus_3D
load("~/rawdata/L2HGDH_scRNA/uninjected/uninjected_anno.RData")

###提取亚群
sce <- subset(sce, celltype=="Macrophage")
###降维
DefaultAssay(sce) <- "RNA"
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 

sce <- FindNeighbors(sce, dims = 1:10)
sce <- FindClusters(sce, resolution = 0.5)
sce <- RunUMAP(sce, dims = 1:10)
pdf(paste0(outdir,"/","04-",gene,"-M_cluster.pdf"),height=10,width=6)
DimPlot(sce, reduction = 'umap', group.by = 'seurat_clusters',label = TRUE, pt.size = 0.5)
dev.off()

save(sce,file = "~/rawdata/L2HGDH_scRNA/uninjected/uninjected_macrophage.RData")


##### 巨噬细胞亚群注释 #####
library(Seurat) 
library(ggplot2)
library(dplyr)
library(scRNAtoolVis)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

load("~/rawdata/L2HGDH_scRNA/injected/injected_macrophage.RData")
###### genes to check ######
# DefaultAssay(sce) <- "RNA"
markers <- c("Nos2","Cd86","Cd80","Il12b",  #M1            
             "Arg1","Mrc1","Cd163","Il10")#M2      

markers_plot <- data.frame(cluster = c(rep("M1 macrophage",4),                                  
                                       rep("M2 macrophage",4)),                                
                           gene = markers)

pdf(paste0(outdir,"/","04-",gene,"-M_markers.pdf"),height=10,width=6)
jjDotPlot(object = sce,          
          markerGene = markers_plot,          
          anno = T,          
          id = 'seurat_clusters',          
          textSize = 10,          
          base_size= 10,          
          plot.margin = c(5,1.5,1.5,1.5))
dev.off()


###### annotation ######
celltype <- c(    "0"="CD8+ Tcell",
                  "1"="CD8+ Tcell",
                  "2"="CD4+ Tcell",
                  "3"="other Tcell",
                  "4"="CD8+ Tcell",
                  "5"="CD4+ Tcell",
                  "6"="CD8+ Tcell",
                  "7"="CD4+ Tcell",
                  "8"="CD4+ Tcell",
                  "9"="CD8+ Tcell",
                  "10" = "CD8+ Tcell",
                  "11" = "CD8+ Tcell",
                  "12" = "CD4+ Tcell",
                  "13" = "CD8+ Tcell",
                  "14" = "other Tcell",
                  "15" = "CD8+ Tcell",
                  "16" = "CD4+ Tcell",
                  "17" = "CD8+ Tcell") ##injected
# celltype <- c(    "0"="CD8+ Tcell", 
#                   "1"="CD8+ Tcell", 
#                   "2"="CD4+ Tcell", 
#                   "3"="other Tcell",
#                   "4"="CD8+ Tcell",
#                   "5"="CD4+ Tcell",
#                   "6"="CD8+ Tcell",
#                   "7"="CD4+ Tcell",
#                   "8"="CD4+ Tcell",
#                   "9"="CD8+ Tcell",
#                   "10" = "CD8+ Tcell",
#                   "11" = "CD8+ Tcell",
#                   "12" = "CD4+ Tcell",
#                   "13" = "CD8+ Tcell",
#                   "14" = "other Tcell",
#                   "15" = "CD8+ Tcell",
#                   "16" = "CD4+ Tcell")##uninjected
sce@meta.data$m_type <- celltype[sce@meta.data$seurat_clusters]

save(sce,file = "~/rawdata/L2HGDH_scRNA/injected/injected_macrophage.RData")

###### plotting--M1&M2 ######
pdf(paste0(outdir,"/","03-",gene,"-Tcell_ratio.pdf"),height=8,width=12)


Cellratio <- prop.table(table(celltype=sce$m_type, group=sce$orig.ident), margin = 2)#计算各组样本不同细胞群比例
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
  cellper_$group <- factor(cellper_$group , levels =c("PBS","VG21","S"))
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
          pplist[[2]])


dev.off()


##### T细胞富集分析 #####
library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle)
library(clusterProfiler)
library(org.Hs.eg.db)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

load("~/rawdata/L2HGDH_scRNA/injected/injected_tcell.RData")

sce <- JoinLayers(sce)
dge <- FindMarkers(sce,ident.1 = "VG21",ident.2 = "S",group.by = "group")##1是ctrl，2是treat
sig_dge <- subset(dge, p_val_adj<0.01&abs(avg_log2FC)>1)

#GO分析(注意是小鼠的)，count表示改变的基因数
ego_ALL <- enrichGO(gene          = row.names(sig_dge),
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Mm.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)
ego_CC <- enrichGO(gene          = row.names(sig_dge),
                   #universe     = row.names(dge.celltype),
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
p1 <- barplot(ekegg, showCategory=10)+ scale_y_discrete(labels = function(x) str_sub(x,1,nchar(x)-28))
p2 <- dotplot(ekegg, showCategory=10)+ scale_y_discrete(labels = function(x) str_sub(x,1,nchar(x)-28))
plotc2 = p1/p2

pdf(paste0(outdir,"/","05-",gene,"-Tcell_GO&KEGG.pdf"),height=16,width=12)
plotc1
plotc2
dev.off()


##### 免疫检查点表达量 #####
library(Seurat)
library(ggpubr) # 继承ggplot语法
library(patchwork) # 拼图包
library(ggsci) #配色包
library(wesanderson)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(data.table)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

load("~/rawdata/L2HGDH_scRNA/injected/injected_tcell.RData")
checkpoint <- fread("~/rawdata/Immune/checkpoint.csv")


##鼠源-人源转换
#使用homologene函进行转换
#@genelist是要转换的基因列表
#@inTax是输入的基因列表所属的物种号，10090是小鼠
#@outTax是要转换成的物种号，9606是人
#install.packages('homologene')
library(homologene)
#homologene(genelist, inTax = 10090, outTax = 9606)
checkpoint <- human2mouse(checkpoint$`Check-point`)
checkpoint$`Check-point` <- checkpoint[,2]


##挑出checkpoint
sce <- JoinLayers(sce)
data_count <- as.data.frame(GetAssayData(object = sce, slot = "data"))
data_count <- as.data.frame(t(data_count))
data_count <- data_count[,colnames(data_count) %in% checkpoint$`Check-point`]

gene_meta <- as.data.frame(sce$group)
colnames(gene_meta) <- "gene_meta"



##批量箱线图
col <- wes_palette("Darjeeling1" , n = 3)

Exp_plot <- merge(gene_meta,data_count,by="row.names")
rownames(Exp_plot) <- Exp_plot$Row.names
Exp_plot <- Exp_plot[,-1]

library(dplyr)
library(rstatix)
library(reshape2)
result=melt(Exp_plot) %>%
  group_by(variable) %>%
  t_test(value ~ gene_meta) %>%
  adjust_pvalue(method = "fdr")#计算p值
result <- result[order(result$p),]
Exp_plot <- Exp_plot[, c(1,na.omit(match(unique(result$variable),colnames(Exp_plot))))]

gene_list <- colnames(Exp_plot[,-1])
sample <- "gene_meta"


plist2<-list()
for (i in 1:length(gene_list)){
  bar_tmp<-Exp_plot[,c(gene_list[i],sample)]
  bar_tmp <- na.omit(bar_tmp)
  colnames(bar_tmp)<-c("Expression","sam")
  my_comparisons1 <- list(c("PBS", "VG21"),c("PBS","S"),c("VG21","S")) 
  bar_tmp$sam <- factor(bar_tmp$sam, levels = c("PBS", "VG21", "S"))
  pb1<-ggviolin(bar_tmp,
                 x="sam",
                 y="Expression",
                 color="sam",
                 fill=NULL,
                 add = "jitter",
                 bxp.errorbar.width = 0.6,
                 width = 0.6,
                 size=0.01,
                 font.label = list(size=30), 
                 palette = col)+theme(panel.background =element_blank())+ylab("Expression") 
  pb1<-pb1+theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())
  pb1<-pb1+theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))
  pb1<-pb1+theme(axis.text.y = element_text(size = 15))+ggtitle(gene_list[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
  pb1<-pb1+theme(legend.position = "NA")#
  pb1<-pb1+stat_compare_means(method="t.test",
                              comparisons =c(my_comparisons1),
                              label="p.signif")
  pb1 <- pb1 + geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA)
  plist2[[i]]<-pb1 
}


n_pdf=ceiling(length(plist2)/8)   ##判断一页4张图，需要多少页pdf
pdf(paste0(outdir,"/","05-",gene,"-Tcell_checkpoint_boxplot.pdf"),height=12,width=16)
plist1 <- plist2
for (i in 1:n_pdf){ 
  if (length(plist1) > 8){
    new_list=plist1[1:8]
    plist1=plist1[-c(1:8)]
    print(plot_grid(plotlist=new_list,ncol=4,nrow=2,hjust=1.5))
  }   ##当list中存在的图数量>4时，把前4个给到画图的newlist中并在一页pdf画出来；最后在原始list删除这4个图
  else {
    new_list=plist1
    print(plot_grid(plotlist=new_list,ncol=4,nrow=2,hjust=1.5))
  }##当list中存在的图数量≤4时，整个list给到现在的newlist中并在一页pdf画出来
}
dev.off()


##热图
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(pheatmap)
library(gtable)
library(grid)

feature <- intersect(rownames(sce[["RNA"]]$scale.data),checkpoint$`Check-point`)
Idents(sce) <- sce$group


pdf(paste0(outdir,"/","05-",gene,"-Tcell_checkpoint_heatmap.pdf"),height=30,width=20)
DoHeatmap(sce,features = feature,label = F,slot = "scale.data")
dev.off()
