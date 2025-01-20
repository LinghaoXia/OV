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
sce@meta.data$group <- ifelse (grepl("VG161",sce@meta.data$orig.ident),ifelse(grepl("_L",sce@meta.data$orig.ident),"PD1_VG161_L","PD1_VG161_R"),'PD1')

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

load("~/rawdata/L2HGDH_scRNA/uninjected/uninjected_cluster.RData")
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
DC <- c("Batf3","Irf8")
Macrophage <- c("Itgam","Adgre1")
Bcell <- c("Blk","Cd19","Casp3","Cd27","Cd5","Cd69","Cd83")
Tcell <- c("Cd3d","Cd3e","Cd81","Rora","Stat3","Tgfb1")
Monocyte <- c("Cd14","Cd16")


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
save(sce,file = "~/rawdata/L2HGDH_scRNA/uninjected/uninjected_anno.RData")

##### 细胞比例 #####
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)


load("~/rawdata/L2HGDH_scRNA/uninjected/uninjected_anno.RData")
table(sce$group)#查看各组细胞数
prop.table(table(sce$celltype))
table(sce$infected, sce$group)#各组不同细胞群细胞数


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