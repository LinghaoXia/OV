##### 髓系细胞细胞亚群聚类 #####
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


load("~/rawdata/scRNA_virus/virus_5D/virus_5D_anno.RData")
sce <- subset(sce,celltype=="Myeloid")


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
sce=FindClusters(sce,resolution = 0.3)#需要对粒度进行调整
save(sce,file = "~/rawdata/scRNA_virus/virus_5D/virus_5D_mye.RData")

pdf(paste0(outdir,"/","07-",gene,"-mye.pdf"),height=10,width=6)
DimPlot(sce, reduction = 'umap', group.by = 'seurat_clusters',label = TRUE, pt.size = 0.5)
dev.off()

##### 髓系细胞细胞亚群注释 #####
library(Seurat) 
library(ggplot2)
library(dplyr)
library(scRNAtoolVis)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)


load("~/rawdata/scRNA_virus/virus_5D/virus_5D_mye.RData")

###### genes to check ######
markers <- c("ADGRE1","CD68","CD86","CD163",  #macrophage           
             "FCN1","CD14","FCGR3A",  #monocyte 
             "BATF3","IRF8","ITGAX","HLA-DMB","HLA-DRB1" #DC 
             )        

markers_plot <- data.frame(cluster = c(rep("Macrophage",4),                                  
                                       rep("monocyte",3),                                  
                                       rep("DC cell",5)),                                  
                           gene = markers)

# sce_X <- subset(sce,idents=c(0,1,2,4))
pdf(paste0(outdir,"/","07-",gene,"-mye_markers.pdf"),height=10,width=6)
jjDotPlot(object = sce,          
          markerGene = markers_plot,          
          anno = T,          
          id = 'seurat_clusters',          
          textSize = 10,          
          base_size= 10,          
          plot.margin = c(4,1.5,1.5,1.5))
dev.off()


###### annotation ######
celltype <- c(    "0"="DC cell",
                  "1"="DC cell",
                  "2"="Macrophage",
                  "3"="Macrophage",
                  "4"="Macrophage",
                  "5"="Macrophage",
                  "6"="DC cell",
                  "7"="Macrophage",
                  "8"="Monocyte",
                  "9"="DC cell",
                  "10" = "DC cell",
                  "11" = "DC cell")
sce@meta.data$mye_type <- celltype[sce@meta.data$seurat_clusters]
save(sce,file = "~/rawdata/scRNA_virus/virus_5D/virus_5D_mye.RData")

###### plotting ######
pdf(paste0(outdir,"/","07-",gene,"-mye_ratio.pdf"),height=8,width=12)

####### Mye 比例 #######
Cellratio <- prop.table(table(celltype=sce$mye_type, group=sce$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- data.frame(Cellratio)

library(reshape2)
cellper <- dcast(Cellratio,group~celltype, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]
#添加分组信息
sample <- unique(sce$orig.ident)
group <- sub("(_1|_2|_3|1|2)$", "", sample)
samples <- data.frame(sample, group)#创建数据框

rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$group <- samples[rownames(cellper),'group']#R添加列
pplist = list()
sce_groups = unique(sce$mye_type)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
for(group_ in sce_groups){
  cellper_  = cellper[,c('sample','group',group_)]
  colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
  cellper_$group <- factor(cellper_$group , levels =c("Vehicle_5D","VG161_5D_R","VG161_5D_L"))
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
    labs(title = group_,y=paste0(group_)) +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  pplist[[group_]] = pp1
}

library(cowplot)
plot_grid(pplist[[1]],
          pplist[[2]],
          pplist[[3]],
          pplist[[4]],
          nrow = 1)

dev.off()



##### 细胞互作 #####
###### 建立实验对象 ######
library(Seurat)
library(tidyverse)
library(patchwork)
library(CellChat)
library(ggalluvial)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)

load("~/rawdata/scRNA_virus/virus_5D/virus_5D_mye.RData") 

##合并
sce.mye <- sce   # 细分后的亚群
load("~/rawdata/scRNA_virus/virus_5D/virus_5D_anno.RData")  # 原数据集
sce <- subset(sce,celltype %in% c("Epithelial","Myeloid"))
Idents(sce.mye) <- "mye_type" # 设置亚群标识
Idents(sce) <- "celltype"
Idents(sce, cells = colnames(sce.mye)) <- Idents(sce.mye)
sce$celltype_new <- Idents(sce)
rm(sce.mye)


##提取表达矩阵和细胞类别创建cellchat对象
group <- unique(sce$group)
sce <- subset(sce,group==group[length(group)])#VG161
sce$celltype_new <- ifelse(sce$celltype == "Epithelial", sce$infected ,as.character(sce$celltype_new))
data.input <- GetAssayData(sce, assay = "SCT", slot = "data")
identity <- subset(sce@meta.data, select = "celltype_new")
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "celltype_new")

rm(sce)
rm(data.input)
gc()

#导入配受体数据库（人或者鼠）
CellChatDB <- CellChatDB.human
#CellChatDB <- CellChatDB.mouse

#查看可以选择的侧面（选择特定的信息描述细胞间的作用）
unique(CellChatDB$interaction$annotation)
# "Secreted Signaling" ，"ECM-Receptor"， "Cell-Cell Contact" 三种可选
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
#对数据进行子集化，节省计算成本
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)

#计算通信概率推断细胞互作的通信网络
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
#过滤掉低质量的细胞间的通信
cellchat <- filterCommunication(cellchat, min.cells = 3)

#在信号通路水平上推断细胞间的通讯
cellchat <- computeCommunProbPathway(cellchat)
##汇总细胞间的通讯
cellchat <- aggregateNet(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")

save(cellchat,file = "~/rawdata/scRNA_virus/virus_5D/virus_5D_chat_mye.rds")


###### 细胞互作--实验组分析 ######
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)


load("~/rawdata/scRNA_virus/virus_5D/virus_5D_chat_mye.rds")
library(CellChat)
library(ggalluvial)

###总的信号通路分析
pdf(paste0(outdir,"/","07-",gene,"-mye_chat.pdf"),height=8,width=12)
#提取所有推断的配体/受体级别的细胞-细胞通信
df.net <- subsetCommunication(cellchat)
#计算聚合细胞互作通信网络
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
#互作的数量
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
#互作的权重
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#箭头发出的细胞为表达配体，箭头指向的细胞表达受体。细胞圆圈的大小代表细胞数量。
##每个细胞如何跟别的细胞互作（互作的强度或概率图）
mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {
  mat1 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat1[i, ] <- mat[i, ]
  netVisual_circle(mat1, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
##每个细胞如何跟别的细胞互作（number+of+interaction图）
mat <- cellchat@net$count
#par(mfrow = c(2,1), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
##气泡图
levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:6), remove.isolate = FALSE, thresh = 0.01)
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:6), remove.isolate = FALSE, thresh = 0.01)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:6), remove.isolate = FALSE, thresh = 0.01)
dev.off()


###单个信号通路分析
pdf(paste0(outdir,"/","07-",gene,"-mye_chat_.pdf"),height=8,width=12)
cellchat@netP$pathways  #查看都有哪些信号通路
##层次图,要找到一个明显的通路先，下面同理
vertex.receiver = c(2,3,8)
netVisual_aggregate(cellchat, signaling = "PROS",  vertex.receiver = vertex.receiver,layout="hierarchy")
##圈图
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="MK", layout = "circle")
##热图
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "MK", color.heatmap = "Reds")
#配体-受体层级的可视化
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()


###### 细胞互作--建立对照对象&实验组对象 ######
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)


load("~/rawdata/scRNA_virus/virus_5D/virus_5D_mye.RData")
library(CellChat)
library(ggalluvial)

##合并
sce.mye <- sce   # 细分后的亚群
load("~/rawdata/scRNA_virus/virus_5D/virus_5D_anno.RData")  # 原数据集
sce <- subset(sce,celltype %in% c("Epithelial","Myeloid"))
Idents(sce.mye) <- "mye_type" # 设置亚群标识
Idents(sce) <- "celltype"
Idents(sce, cells = colnames(sce.mye)) <- Idents(sce.mye)
sce$celltype_new <- Idents(sce)
rm(sce.mye)


##提取表达矩阵和细胞类别创建cellchat对象
group <- unique(sce$group)
sce <- subset(sce,group==group[1])#Vehicle
# sce <- subset(sce,group==group[length(group)])#VG161
data.input <- GetAssayData(sce, assay = "SCT", slot = "data")
identity <- subset(sce@meta.data, select = "celltype_new")
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "celltype_new")

rm(sce)
rm(data.input)
gc()

#导入配受体数据库（人或者鼠）
CellChatDB <- CellChatDB.human
#CellChatDB <- CellChatDB.mouse

#查看可以选择的侧面（选择特定的信息描述细胞间的作用）
unique(CellChatDB$interaction$annotation)
# "Secreted Signaling" ，"ECM-Receptor"， "Cell-Cell Contact" 三种可选
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
#对数据进行子集化，节省计算成本
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)

#计算通信概率推断细胞互作的通信网络
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
#过滤掉低质量的细胞间的通信
cellchat <- filterCommunication(cellchat, min.cells = 3)

#在信号通路水平上推断细胞间的通讯
cellchat <- computeCommunProbPathway(cellchat)
##汇总细胞间的通讯
cellchat <- aggregateNet(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")

cellchat_0 <- cellchat
save(cellchat_0,file = "~/rawdata/scRNA_virus/virus_5D/virus_5D_chat0_mye.rds")
# cellchat_1 <- cellchat
# save(cellchat_1,file = "~/rawdata/scRNA_virus/virus_5D/virus_5D_chat1_mye.rds")


###### 细胞互作--差异分析 ######
library(Seurat)
library(tidyverse)
library(patchwork)
library(CellChat)
library(ggalluvial)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)


load("~/rawdata/scRNA_virus/virus_5D/virus_5D_chat1_mye.rds")
load("~/rawdata/scRNA_virus/virus_5D/virus_5D_chat0_mye.rds")


##合并对象
cco.list <- list(Vehicle=cellchat_0,VG161=cellchat_1)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)
rm(cellchat_0)
rm(cellchat_1)
gc()

pdf(paste0(outdir,"/","07-",gene,"-mye_chat_differential.pdf"),height=12,width=16)
##柱状图-所有细胞群总体观：通讯数量与强度对比
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#数量与强度差异网络图
par(mfrow = c(1,2),xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


#数量与强度差异热图
par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2
dev.off()


###### Cellcall ######
library(devtools)
library(Seurat)
library(tidyverse)
library(patchwork)
library(cellcall)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)

load("~/rawdata/scRNA_virus/virus_5D/virus_5D_mye.RData") 

##合并
sce.mye <- sce   # 细分后的亚群
load("~/rawdata/scRNA_virus/virus_5D/virus_5D_anno.RData")  # 原数据集
sce <- subset(sce,celltype %in% c("Epithelial","Myeloid"))
Idents(sce.mye) <- "mye_type" # 设置亚群标识
Idents(sce) <- "celltype"
Idents(sce, cells = colnames(sce.mye)) <- Idents(sce.mye)
sce$celltype_new <- Idents(sce)
rm(sce.mye)


##提取VG161组，比较Bystander和Infected
group <- unique(sce$group)
sce <- subset(sce,group==group[length(group)])#VG161
sce$celltype_new <- ifelse(sce$celltype == "Epithelial", sce$infected ,as.character(sce$celltype_new))

##构建对象
sce_obj <- CreateObject_fromSeurat(Seurat.object=sce, #seurat对象
                                slot="counts", 
                                cell_type="celltype_new", #细胞类型
                                data_source="UMI",
                                scale.factor = 10^6, 
                                Org = "Homo sapiens") #物种信息

##计算得分
mt <- TransCommuProfile(object = sce_obj,
                        pValueCor = 0.05,
                        CorValue = 0.1,
                        topTargetCor=1,
                        p.adjust = 0.2,
                        use.type="median",
                        probs = 0.9,
                        method="max",
                        IS_core = TRUE,
                        Org = 'Homo sapiens')

pdf(paste0(outdir,"/","07-",gene,"-mye_cellcall.pdf"),height=12,width=16)

##涉及的主要通路
n <- mt@data$expr_l_r_log2_scale
pathway.hyper.list <- lapply(colnames(n), function(i){
  print(i)
  tmp <- getHyperPathway(data = n, object = mt, cella_cellb = i, Org="Homo sapiens")
  return(tmp)
})
myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=colnames(n))
plotBubble(myPub.df)

##细胞互作可视化(圈图)
#有多少细胞类型就设置多少个颜色
cell_color <- data.frame(color=c("#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F"), stringsAsFactors = FALSE) #,"#FF6A6A","#7CFC00"
rownames(cell_color) <- c("Bystander","DC cell","Infected","Macrophage","Monocyte","other Myeloid")
ViewInterCircos(object = mt, font = 2, cellColor = cell_color, 
                lrColor = c("#F16B6F", "#84B1ED"),
                arr.type = "big.arrow",arr.length = 0.04,
                trackhight1 = 0.05, slot="expr_l_r_log2_scale",
                linkcolor.from.sender = TRUE,
                linkcolor = NULL, gap.degree = 0.16, #细胞类型多的话设置小点，不然图太大画不出来
                order.vector=c("Bystander","DC cell","Infected","Macrophage","Monocyte","Naive","other Myeloid"),
                trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = FALSE)

##细胞互作可视化(热图)
viewPheatmap(object = mt, slot="expr_l_r_log2_scale", show_rownames = T,
             show_colnames = T,treeheight_row=0, treeheight_col=10,
             cluster_rows = T,cluster_cols = F,fontsize = 12,angle_col = "45",  
             main="score")

##转录因子富集图（某一细胞）
tf <- names(mt@data$gsea.list$Macrophage@geneSets)
tf ##根据tf结果选择对应的转录因子
getGSEAplot(gsea.list=mt@data$gsea.list, geneSetID = tf, 
            myCelltype="Macrophage", fc.list=mt@data$fc.list,  
            selectedGeneID = mt@data$gsea.list$Macrophage@geneSets$NR1H2[1:10],##选择对应的转录因子
            mycol = NULL)


dev.off()

#### 基因表达量 ####
library(tidyverse)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(gghalves)
library(ggpubr)
library(paletteer)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "scRNA_virus"
outdir = paste0("~/OV/",gene)

load("~/rawdata/scRNA_virus/virus_5D/virus_5D_mye.RData") 

##合并
sce.mye <- sce   # 细分后的亚群
load("~/rawdata/scRNA_virus/virus_5D/virus_5D_anno.RData")  # 原数据集
Idents(sce.mye) <- "mye_type" # 设置亚群标识
Idents(sce) <- "celltype"
Idents(sce, cells = colnames(sce.mye)) <- Idents(sce.mye)
sce$celltype_new <- Idents(sce)
rm(sce.mye)


##小提琴图--感染前后
gene_list <- c("MDK","SDC1","SDC2","SDC4","LRP1","NCL","ITGA6","ITGB1")
exprs <- data.frame(FetchData(object = sce, vars = c("celltype_new",gene_list,"group")))
exprs$Proj <- "Seurat"
exprs$Cell <- rownames(exprs)
exprs.melt <- reshape2::melt(exprs,                              
                             id.vars = c("Cell","celltype_new","group"),                              
                             measure.vars = gene_list,                             
                             variable.name = "gene",                              
                             value.name = "Expr")

p1 <- ggplot()+  
  geom_half_violin(data = exprs.melt[exprs.melt$group == 'Vehicle',],                   
                   aes(x = celltype_new, y = Expr, fill = group),
                   color = 'black',                   
                   scale = 'width') +   
  facet_grid(rows = vars(gene), scales = 'free_y') +   
  geom_half_violin(data = exprs.melt[exprs.melt$group == 'VG161',],                   
                   aes(x = celltype_new, y = Expr, fill = group),                   
                   color = 'black',                   
                   scale = 'width',                   
                   side = 'r') +   
  facet_grid(rows = vars(gene), scales = 'free_y')+
  theme_bw() +  
  theme(panel.grid = element_blank()) +  
  scale_fill_manual(values = c("#E39A35","#68A180")) +  
  labs(x = "", y = 'Expression Level') #y轴标题本文内容修改

##小提琴图--旁观者和感染
group <- unique(sce$group)
sce2 <- subset(sce,celltype=="Epithelial")
sce2 <- subset(sce2,group==group[length(group)])
sce2$celltype_new <- ifelse(sce2$celltype == "Epithelial", sce2$infected ,as.character(sce2$celltype_new))
gene_list <- c("MDK","SDC1","SDC2","SDC4","LRP1","NCL","ITGA6","ITGB1")
exprs <- data.frame(FetchData(object = sce2, vars = c("celltype_new",gene_list)))
exprs$Proj <- "Epithelial"
exprs$Cell <- rownames(exprs)
exprs.melt <- reshape2::melt(exprs,                              
                             id.vars = c("Cell","celltype_new","Proj"),                              
                             measure.vars = gene_list,                             
                             variable.name = "gene",                              
                             value.name = "Expr")
  
p2 <-  ggplot()+  
  geom_half_violin(data = exprs.melt[exprs.melt$celltype_new == 'Bystander',],                   
                   aes(x = Proj, y = Expr, fill = celltype_new),
                   color = 'black',                   
                   scale = 'width') +   
  facet_grid(cols = vars(gene), scales = 'free') +   
  geom_half_violin(data = exprs.melt[exprs.melt$celltype_new == 'Infected',],                   
                   aes(x = Proj, y = Expr, fill = celltype_new),                   
                   color = 'black',                   
                   scale = 'width',                   
                   side = 'r') +   
  facet_grid(cols = vars(gene), scales = 'free')+
  theme_bw() +  
  theme(panel.grid = element_blank()) +  
  scale_fill_manual(values = c("#E39A35","#68A180")) +  
  labs(x = "", y = 'Expression Level') #y轴标题本文内容修改



pdf(paste0(outdir,"/","09-",gene,"-MDK-NCL.pdf"),height=12,width=12)
pal <- viridis(n = 10, option = "D")
p1
p2
##MDK
DimPlot(sce, reduction = 'umap', group.by = 'celltype_new',label = TRUE, pt.size = 0.5)+FeaturePlot_scCustom(seurat_object = sce, features = "MDK", colors_use = pal)
FeaturePlot_scCustom(seurat_object = sce, features = "MDK", split.by = "group",num_columns = 2,colors_use = pal)
##NCL
DimPlot(sce, reduction = 'umap', group.by = 'celltype_new',label = TRUE, pt.size = 0.5)+FeaturePlot_scCustom(seurat_object = sce, features = "NCL", colors_use = pal)
FeaturePlot_scCustom(seurat_object = sce, features = "NCL", split.by = "group",num_columns = 2,colors_use = pal)
dev.off()




