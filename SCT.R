##### 数据加载 #####
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "SCT"
outdir = paste0("~/OV/",gene)

sct <- Load10X_Spatial(data.dir ="./SCT/result/3D_Vehicle_2/outs", 
                               filename = "filtered_feature_bc_matrix.h5")

pdf(paste0(outdir,"/","01-",gene,"-pre.pdf"))
plot1 <- VlnPlot(sct, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(sct, features = "nCount_Spatial",pt.size.factor = 3) + theme(legend.position = "right")
plot_grid(plot1, plot2)
dev.off()

##标准化基因可视化
sct <- SCTransform(sct, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

pdf(paste0(outdir,"/","01-",gene,"-MDK_NCL.pdf"))
DefaultAssay(sct) <- "SCT"
SpatialFeaturePlot(sct, features = c("MDK", "NCL"), pt.size.factor = 3)
dev.off()


##降维、聚类和可视化
sct <- RunPCA(sct, assay = "SCT", verbose = FALSE)
sct <- FindNeighbors(sct, reduction = "pca", dims = 1:30)
sct <- FindClusters(sct, verbose = FALSE)
sct <- RunUMAP(sct, reduction = "pca", dims = 1:30)


pdf(paste0(outdir,"/","01-",gene,"-cluster.pdf"))
DimPlot(sct, reduction = "umap", label = TRUE)
SpatialDimPlot(sct, label = TRUE, label.size = 3,pt.size.factor = 3)
dev.off()

save(sct,file="~/rawdata/SCT/analysis/vehicle_3D_1/vehicle_3D_cluster.RData")

##### 结合单细胞数据注释细胞类型 #####
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "SCT"
outdir = paste0("~/OV/",gene)

load("~/rawdata/SCT/analysis/virus_3D_1/virus_3D_cluster.RData")
load("~/rawdata/scRNA_virus/virus_3D/virus_3D_mye.RData") 


###### 合并 ######
sce.mye <- sce   # 细分后的亚群
load("~/rawdata/scRNA_virus/virus_3D/virus_3D_anno.RData")  # 原数据集
Idents(sce.mye) <- "mye_type" # 设置亚群标识
Idents(sce) <- "celltype"
Idents(sce, cells = colnames(sce.mye)) <- Idents(sce.mye)
sce$celltype_new <- Idents(sce)
rm(sce.mye)
sce <- subset(sce,orig.ident=="VG161_3D_2")
sce$celltype_new <- ifelse(sce$celltype == "Epithelial", sce$infected ,as.character(sce$celltype_new))

###### Seurat-Mapping ######
##平衡采样，每种取300个
sce_sub <- subset(sce, cells = unlist(lapply(split(Cells(sce), sce$celltype_new), function(x) head(x, 300))))
DefaultAssay(sct) <- "SCT"
anchors <- FindTransferAnchors(reference = sce_sub, 
                               query = sct, 
                               normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = sce_sub$celltype_new, 
                                  prediction.assay = TRUE,
                                  weight.reduction = sct[["pca"]], 
                                  dims = 1:30)
predictions.res <- predictions.assay@data

sct[["predictions"]] <- predictions.assay
save(sct,file="~/rawdata/SCT/analysis/virus_5D_1/virus_5D_anno.RData")

pdf(paste0(outdir,"/","02-",gene,"-anno_infected.pdf"))
DefaultAssay(sct) <- "predictions"
SpatialFeaturePlot(sct, features = c("Macrophage","Tcell","Epithelial","Infected"),
                   pt.size.factor = 3, ncol = 3, crop = TRUE)
dev.off()

##细胞类型占比矩阵
library(SPOTlight)

anno <- sct@assays[["predictions"]]@data
mat <- t(anno[rownames(anno)!="max", ])
ct <- colnames(mat)
mat[mat < 0.1] <- 0

# 颜色定义
paletteMartin <- c(
 "#004949","#009292","#ff6db6","#ffb6db",
  "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
  "#920000","#924900","#db6d00","#24ff24","#ffff6d")

pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct

# plot
spe <- sct@images$slice1$centroids@coords
rownames(spe) <- rownames(mat)
head(spe)
pdf(paste0(outdir,"/","02-",gene,"-anno_spotlight.pdf"))
plotSpatialScatterpie(x = spe, y = mat, cell_types = colnames(mat), img = FALSE, scatterpie_alpha = 1, pie_scale = 0.4) +
  scale_fill_manual(values = pal, breaks = names(pal))
dev.off()


###### SPOTlight ######
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)


### singlecellexperiment
scm <- as.SingleCellExperiment(sce)
### spatialexperiment
counts <- as.matrix(GetAssayData(sct, assay = "Spatial", slot = "counts"))
coords <- GetTissueCoordinates(sct) 
spe <- SpatialExperiment(
  assays = list(counts = counts),
  spatialCoords = as.matrix(coords[, c("x", "y")])  
)
spe <- addImg(spe, sample_id="sample01",#和coldata sample_id里的内容一致              
              image_id = "xx",
              imageSource="~/rawdata/SCT/3D_VG161_2/image.png",#课上文件夹里有但其实是我随便截个图              
              scaleFactor = 1,               
              load = TRUE)


### Feature selection
scm <- logNormCounts(scm)

### Variance modelling
# 去掉核糖体和线粒体基因
genes <- !grepl(pattern = "^RP[L|S]|MT", x = rownames(scm))
dec <- modelGeneVar(scm , subset.row = genes)
# 计算高变基因
hvg <- getTopHVGs(dec, n = 3000)
# 加上细胞注释信息
colLabels(scm) <- colData(scm)$celltype_new
# Compute marker genes
mgs <- scoreMarkers(scm, subset.row = genes)
# 保留最相关的marker基因
mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.6, ]
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)

res <- SPOTlight(
  x = scm, 
  y = spe,
  groups = as.character(scm$celltype_new), # 也可以是cluster，
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")

sct[["SPOTlight"]] <- CreateAssayObject(t(res$mat))
pdf(paste0(outdir,"/","02-",gene,"-spotlight_infected.pdf"))
DefaultAssay(sct) <- "SPOTlight"
SpatialFeaturePlot(sct, features = c("Macrophage","Tcell","Epithelial","Infected"),
                   pt.size.factor = 3, ncol = 3, crop = TRUE)
dev.off()
### 结果可视化
head(mat <- res$mat)[, seq_len(length(unique(scm$celltype_new)))]
mod <- res$NMF
res.data <- (mat <- res$mat)[, seq_len(length(unique(scm$celltype_new)))]

ct <- colnames(mat)
# 占比小于0.1的不展示
mat[mat < 0.1] <- 0
# 颜色设置
paletteMartin <- c(
  "#000000","#004949","#009292","#ff6db6","#ffb6db",
  "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
  "#920000","#924900","#db6d00","#24ff24","#ffff6d")
pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct
pal

plotSpatialScatterpie(
  x = spe,
  y = mat,
  cell_types = colnames(mat),
  img = F,
  scatterpie_alpha = 1,
  pie_scale = 0.4) +
  scale_fill_manual(
    values = pal,
    breaks = names(pal))
pdf(paste0(outdir,"/","02-",gene,"-anno_spotlight.pdf"))
plotSpatialScatterpie(x = spe, y = mat, cell_types = colnames(mat), img = FALSE, scatterpie_alpha = 1, pie_scale = 0.4) +
  scale_fill_manual(values = pal, breaks = names(pal))
dev.off()

####### CRAD ######
library(CARD)
library(MuSiC)
library(Seurat)
library(patchwork)
library(tidyverse)

###获取空转的counts表达矩阵
spatial_count <-  GetAssayData(sct, assay = "Spatial", slot = "counts")
spatial_count[1:4,1:4]
###获取空转的空间位置矩阵
spatial_loca <- GetTissueCoordinates(sct) 
spatial_location <- spatial_loca[,1:2]
#名字必须是x y ，否则后面CARD_deconvolution会报错
colnames(spatial_location) <- c("x","y")
spatial_location[1:3,]

###获取单细胞counts矩阵
sc_count <- sce@assays$RNA@counts
###获取单细胞细胞注释矩阵
sc_meta <- sce@meta.data %>% 
  rownames_to_column("cellID") %>%
  dplyr::select(cellID,orig.ident,celltype_new) %>% 
  mutate(CB = cellID) %>% 
  column_to_rownames("CB")
head(sc_meta)

###构建CARD对象，并进行空间细胞成分反卷积
CARD_obj = createCARDObject( 
  sc_count = sc_count, 
  sc_meta = sc_meta, 
  spatial_count = spatial_count, 
  spatial_location = spatial_location, 
  ct.varname = "celltype_new", 
  ct.select = unique(sc_meta$celltype_new), #细胞类型列名
  sample.varname = "orig.ident")
#CARD 解卷积
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
#CARD-spot 可视化spot的细胞类型分布饼图
colors = c("#4DAF4A","#F0027F","#377EB8","#FDC086","#A6761D","#FFFF00","#BEAED4",
           "#BF5B17","#666666",'#AA0000')
p1<-CARD.visualize.pie(proportion = CARD_obj@Proportion_CARD,
                       spatial_location = CARD_obj@spatial_location, 
                       colors = colors)
ct.visualize = c("Macrophage","Tcell","Infected")

p2 <- CARD.visualize.prop(proportion = CARD_obj@Proportion_CARD,        
                          spatial_location = CARD_obj@spatial_location, 
                          ct.visualize = ct.visualize,                
                          colors = c("lightblue","lightyellow","red"), 
                          NumCols = 3,pointSize = 1)#图中spot大小


##### 病毒感染 #####

DefaultAssay(sct) <- "SCT"
VirTranscript <- function(obj,viral_genes){
  viral_counts <- FetchData(obj, vars = viral_genes, slot = 'counts')
  viral_counts_total <- rowSums(viral_counts)
  total_counts <- rowSums(FetchData(obj,slot = 'counts',vars = rownames(obj))) 
  scaled_viral_counts <- log2((viral_counts_total / total_counts) * 10000+1)
  obj$VG161_transcript <- scaled_viral_counts
  return(obj)
}

pdf(paste0(outdir,"/","03-",gene,"-infected.pdf"))
viral_genes <-  rownames(sct@assays$SCT)[grep("VG161", rownames(sct@assays$SCT))]
sct <- VirTranscript(sct,viral_genes)
SpatialFeaturePlot(sct, features = "VG161_transcript",pt.size.factor = 3)


viral_genes <-  rownames(sce@assays$SCT)[grep("VG161", rownames(sce@assays$SCT))]
sct[["VG161"]] <- PercentageFeatureSet(sct,pattern = c("^VG161-UL","^VG161-ICP"))
SpatialFeaturePlot(sct, features = "VG161",pt.size.factor = 3)

dev.off()



##### 细胞互作 ####
library(CellChat)
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "SCT"
outdir = paste0("~/OV/",gene)
load("~/rawdata/SCT/analysis/virus_3D_1/virus_3D_anno.RData")

SpatialFeaturePlot(sct, features = c("Macrophage","Tcell","Epithelial","Infected","Bystander"),
                   pt.size.factor = 3, ncol = 2, crop = TRUE)

###区域定义
sct@meta.data$Region<-NA
sct@meta.data$Region[sct@meta.data$seurat_clusters %in% c('2','6')] <- "Macrophage"
sct@meta.data$Region[sct@meta.data$seurat_clusters %in% c('0','4')] <- "Tcell"
sct@meta.data$Region[sct@meta.data$seurat_clusters %in% c('1','5')] <- "Infected"
sct@meta.data$Region[sct@meta.data$seurat_clusters %in% c('3')] <- "Bystander"
SpatialPlot(sct, label = TRUE, label.size = 3,pt.size.factor = 3,group.by = 'Region',cols = c('Bystander'='#4b5cc4','Macrophage'='#FE8D3C','Infected'='#AA0000','Tcell'='#4DAF4A'))
Idents(sct) <- "Region" 

###获取空转矩阵信息
data.input = Seurat::GetAssayData(sct, slot = "data", assay = "SCT") 
###获取meta信息
meta = data.frame(labels = Idents(sct),
                  row.names = names(Idents(sct)))
###获取空间位置信息
spatial.locs = Seurat::GetTissueCoordinates(sct, scale = NULL,cols = c("imagerow", "imagecol")) 
spatial.locs <- spatial.locs[,1:2]
#名字必须是x y ，否则后面CARD_deconvolution会报错
colnames(spatial.locs) <- c("x","y")
#scalefactors_json存于文件夹下
scalefactors = jsonlite::fromJSON(txt = file.path("~/rawdata/SCT/result/3D_VG161_2/outs/spatial", 'scalefactors_json.json')) 
scalefactors = list(spot.diameter = 65, spot = scalefactors$spot_diameter_fullres, # these two information are required
                     fiducial = scalefactors$fiducial_diameter_fullres, hires = scalefactors$tissue_hires_scalef, lowres = scalefactors$tissue_lowres_scalef # these three information are not required
)


###创建CellChat对象
cellchat <- createCellChat(object = data.input, 
                           meta = meta, 
                           group.by = "labels", #定义的名字是labels
                           datatype = "spatial", #数据类型：空转
                           coordinates = spatial.locs,
                           scale.factors = scalefactors)


###CellChat分析
#设置参考数据库
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
#使用CellChatDB的子集进行细胞间通信分析
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") #选择Secreted Signaling
cellchat@DB <- CellChatDB.use
#CellChat预处理
cellchat <- subsetData(cellchat) #即使使用整个数据库，此步骤也是必要的
future::plan("multisession", workers = 4) #多线程
#识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
#识别过表达配体受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
#细胞间通信网络的推断
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE,scale.distance = 0.01)
#默认情况下，每个细胞组中用于细胞间通信所需的最小细胞数为10
cellchat <- filterCommunication(cellchat, min.cells = 10)
#在信号通路水平上推断细胞间通讯
cellchat <- computeCommunProbPathway(cellchat)
#计算聚合的 cell-cell 通信网络
cellchat <- aggregateNet(cellchat)

###可视化
pdf(paste0(outdir,"/","03-",gene,"-cellchat.pdf"))
#可视化交互次数或总交互次数 
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#热图显示celltype间的通讯次数（左）或总通讯强度(右)
p1 <- netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
p2 <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")
p1 + p2

#展示显著通路结果
cellchat@netP$pathways
par(mfrow=c(1,1), xpd = TRUE)# xpd = TRUE以显示标题
pathways.show <- c("MK")
#可视化 'PTN' 信号网络
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
#在空间转录组上显示'MDK'信号网络
netVisual_aggregate(cellchat, 
                    signaling = pathways.show, 
                    layout = "spatial", 
                    edge.width.max = 3, 
                    alpha.image = 0.2, 
                    vertex.weight = "outgoing", #以更大的圆圈表示更大的传出信号
                    vertex.size.max = 6, 
                    vertex.label.cex = 4.5)
#在热图上显示'MDK'信号网络
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                  width = 8, height = 2.5, font.size = 10)

#取配体-受体对的输入，并以气泡图显示表达
netVisual_bubble(cellchat, sources.use = c(3,4), 
                 targets.use = c(1,2,3,4), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(1,2,3,4),                  
                 signaling = c("MK"), remove.isolate = FALSE)##指定通路
#取配体-受体对的输入，并以二进制形式显示表达
spatialFeaturePlot(cellchat, pairLR.use = "MDK_NCL", point.size = 1, do.binary = TRUE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)
dev.off()