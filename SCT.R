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

sct <- Load10X_Spatial(data.dir ="./SCT/result/3D_VG161_2/outs", 
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

save(sct,file="~/rawdata/SCT/analysis/virus_3D_1/virus_3D_cluster.RData")

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


##合并
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
anchors <- FindTransferAnchors(reference = sce, 
                               query = sct, 
                               normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = sce$celltype_new, 
                                  prediction.assay = TRUE,
                                  weight.reduction = sct[["pca"]], 
                                  dims = 1:30)
predictions.res <- predictions.assay@data

sct[["predictions"]] <- predictions.assay
save(sct,file="~/rawdata/SCT/analysis/virus_3D_1/virus_3D_anno.RData")

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
DefaultAssay(sct) <- "SPOTlight"

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

