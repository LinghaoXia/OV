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


###降维
sce <- SCTransform(sce,vst.flavor = "v2", verbose = FALSE, method = "glmGamPoi",vars.to.regress = "percent.mt")
sce=RunPCA(sce,assay="SCT",verbose = FALSE)
##去批次
sce=RunHarmony(sce,group.by.vars="orig.ident",assay.use="SCT", plot_convergence = TRUE,max.iter.harmony =50 )
##最佳PC数量
pct <- sce [["harmony"]]@stdev / sum( sce [["harmony"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
#获取了最佳PC用于UMAP和FindNeighbors
bestpc=1:pcs
sce<- sce %>% RunUMAP(reduction = "harmony", dims = bestpc) %>% 
  FindNeighbors(reduction = "harmony", dims = bestpc)
sce=FindClusters(sce,resolution = 0.2)#需要对粒度进行调整
save(sce,file = "~/rawdata/scRNA_virus/virus_3D/virus_3D_tcell.RData")

pdf(paste0(outdir,"/","08-",gene,"-tcell.pdf"),height=10,width=6)
DimPlot(sce, reduction = 'umap', group.by = 'seurat_clusters',label = TRUE, pt.size = 0.5)
dev.off()


##### T细胞亚群注释 #####
library(Seurat) 
library(ggplot2)
library(dplyr)
library(scRNAtoolVis)
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
          anno = T,          
          id = 'seurat_clusters',          
          textSize = 10,          
          base_size= 10,          
          plot.margin = c(4,1.5,1.5,1.5))
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
