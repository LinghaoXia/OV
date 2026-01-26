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

sct <- Load10X_Spatial(data.dir ="./SCT/virus_3D/VG161_3D_1", 
                               filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial",
                       slice = "slice1",
                       filter.matrix = TRUE)

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
#load("~/rawdata/scRNA_virus/virus_3D/virus_3D_mye.RData") 


###### 合并 ######
#sce.mye <- sce   # 细分后的亚群
load("~/rawdata/scRNA_virus/virus_3D/virus_3D_anno.RData")  # 原数据集
#Idents(sce.mye) <- "mye_type" # 设置亚群标识
Idents(sce) <- "celltype"
#Idents(sce, cells = colnames(sce.mye)) <- Idents(sce.mye)
sce$celltype_new <- Idents(sce)
#rm(sce.mye)
sce <- subset(sce,orig.ident=="VG161_3D_1")


###### Seurat-Mapping ######
##平衡采样，每种取300个
sce <- SCTransform(sce, ncells = 3000, verbose = FALSE) %>% 
  RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30)
sct <- SCTransform(sct, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
anchors <- FindTransferAnchors(reference = sce, 
                               query = sct, 
                               dims = 1:50, 
                               normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = sce$celltype_new, 
                                  prediction.assay = TRUE,
                                  weight.reduction = sct[["pca"]], 
                                  dims = 1:30)
predictions <- TransferData(anchorset = anchors, 
                                  refdata = sce$celltype_new,
                                  weight.reduction = sct[["pca"]], 
                                  dims = 1:30)
predictions.res <- predictions.assay@data
predictions.id <- predictions$predicted.id

sct[["predictions"]] <- predictions.assay
sct[["predictions_type"]] <- predictions.id
save(sct,file="~/rawdata/SCT/analysis/virus_3D_1/virus_3D_anno.RData")

pdf(paste0(outdir,"/","02-",gene,"-anno_infected.pdf"))
DefaultAssay(sct) <- "predictions"
SpatialFeaturePlot(sct, features = c("Tcell","Epithelial","Fibroblast"),
                   pt.size.factor = 3, ncol = 3, crop = TRUE)
SpatialDimPlot(sct, group.by = "predictions_type", label = TRUE)
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
              imageSource="~/rawdata/SCT/virus_3D/VG161_3D_1/spatial/tissue_hires_image.png",#课上文件夹里有但其实是我随便截个图              
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
  x <- x[x$mean.AUC > 0.8, ]
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
SpatialFeaturePlot(sct, features = c("Fibroblast","Tcell","Epithelial","Myeloid"),
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
plotCorrelationMatrix(mat)
# 共定位（Co-localization）
plotInteractions(mat, which = "heatmap", metric = "prop")
plotInteractions(mat, which = "heatmap", metric = "jaccard")
plotInteractions(mat, which = "network")
main_type <- apply(sct$SPOTlight@data, 2, function(x) names(x)[which.max(x)])
sct$predicted.id <- main_type
SpatialDimPlot(sct, group.by = "predicted.id", label = TRUE,pt.size.factor = 3)
dev.off()
save(sct,file="~/rawdata/SCT/analysis/virus_3D_1/virus_3D_anno.RData")


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
sct[["VG161"]] <- PercentageFeatureSet(sct,pattern = "^VG161-UL|^VG161-ICP|^VG161-US|^VG161-IRL")
SpatialFeaturePlot(sct, features = "VG161",pt.size.factor = 3)

dev.off()

##### CellTrek #####
library(reshape2)
library(CMAP) 
library(Seurat) 
library(e1071)
library(purrr)  
library(dplyr)
library(preprocessCore)
library(reticulate)
library(smfishHmrf)
library(Giotto)
library(ggplot2)
rm(list=ls())
gc()

setwd('~/rawdata')
gene = "SCT"
outdir = paste0("~/OV/",gene,"/CMAP")
python_path <- '/home/xialinghao/anaconda3/envs/CMAP_Env/bin/python'
use_condaenv(python_path)
save_directory <- outdir
if(!file.exists(save_directory)) dir.create(save_directory, recursive = T)


load("~/rawdata/SCT/analysis/virus_3D_1/virus_3D_anno.RData")
load("~/rawdata/scRNA_virus/virus_3D/virus_3D_anno.RData")
source("~/code/new_function.R")
Idents(sce) <- "celltype"
sce$celltype_new <- Idents(sce)
sce <- subset(sce,orig.ident=="VG161_3D_1")
sce$celltype_new <- ifelse(sce$celltype_new=="Epithelial",sce$infected,sce$celltype)
unique(sce$celltype_new)

sc_counts <- GetAssayData(sce, assay = "RNA", layer = "counts")
sc_meta <- data.frame(sce@meta.data,row.names=rownames(sce@meta.data))

spatial_count <- as.matrix(GetAssayData(sct, assay = "Spatial", layer = "counts"))
spatial_location <- GetTissueCoordinates(sct) [,1:2]

sc_counts <- sc_counts[rowSums(sc_counts)>0,]
sc_norm = as.matrix(log1p(sweep(sc_counts,2,Matrix::colSums(sc_counts),FUN = '/') * 1e4))

spatial_count <- spatial_count[rowSums(spatial_count)>0,]
st_norm = log1p(sweep(spatial_count,2,Matrix::colSums(spatial_count),FUN = '/') * 1e4)

cluster_k <- 3
# Create specific instructions for Giotto analysis workflow
instrs <- createGiottoInstructions(save_plot = TRUE,
                                   show_plot = TRUE,
                                   return_plot = TRUE,
                                   python_path = python_path,
                                   save_dir = save_directory)
spatial_obj <- createGiottoObject(raw_exprs = spatial_count,
                                  spatial_locs = spatial_location[,c('x','y')],
                                  instructions = instrs,
                                  cell_metadata = spatial_location)
# Filter genes and cells. If you have filtered some low quality spots before, you can skip this step
spatial_obj <- filterGiotto(gobject = spatial_obj,
                            expression_threshold = 1,
                            gene_det_in_min_cells = 50,
                            min_det_genes_per_cell = 250,
                            expression_values = c('raw'),
                            verbose = T)
spatial_obj <- normalizeGiotto(gobject = spatial_obj, scalefactor = 6000, verbose = T)

# Create spatial network
#@ maximum_distance_knn: Visium data, tissue_hires_scalef set ceiling(24.8/tissue_hires_scalef) or as maximum_distance_knn,  tissue_hires_scalef is saved in scalefactors_json.json; slide-seq/ST: set 1.5
spatial_obj <- createSpatialNetwork(gobject = spatial_obj,
                                    method = 'kNN',
                                    k = 6, # this k represents the number of neighbors
                                    maximum_distance_knn = 370, 
                                    minimum_k = 1,
                                    name = 'KNN_network')
kmtest  <- binSpect(spatial_obj, calc_hub = T, hub_min_int = 5,spatial_network_name = 'KNN_network')

hmrf_folder = paste0(save_directory,'/11_HMRF')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)
spatial_genes_selected <- hmrf_spatial_gene(spatial_obj,
                                            kmtest,
                                            k = cluster_k) # k: Number of spatial domains; set according to your data.

#@ betas: For detailed settings, see https://search.r-project.org/CRAN/refmans/smfishHmrf/html/smfishHmrf.hmrfem.multi.it.min.html
# For quick results, we recomoned setting betas to 45(non-tumor) or 0(tumor sample).
# If you don't mind taking more time and want the best results, you can iteratively test values between 0 and 100 and select the best one.
HMRF_spatial_genes = doHMRF(gobject = spatial_obj,
                            expression_values = 'scaled',
                            spatial_genes = spatial_genes_selected,
                            k = cluster_k, # This value should match the number of spatial domains (k).
                            spatial_network_name="KNN_network",
                            betas = c(0, 45, 2), 
                            python_path = python_path,
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_topgenes_elbow_k_scaled'))
#@betas_to_add: Results from different betas that you want to add
# Recommendations: Tumor sample: beta=0; Non-tumor: beta=45.
beta = 0
spatial_obj = addHMRF(gobject = spatial_obj,
                      HMRFoutput = HMRF_spatial_genes,
                      k = cluster_k,
                      betas_to_add = beta,  # according to the above beta settings
                      hmrf_name = 'HMRF')
# Add spatial domain to spatial metadata. You can also save the spatial_location as an intermediate file, which must include spatial genes and spatial cluster labels.
spatial_location = spatial_location[as.data.frame(spatial_obj@cell_metadata)[,'cell_ID'],]
column <- paste0('HMRF_k',cluster_k,'_b.',beta)
spatial_location$HMRF_cluster <- spatial_obj@cell_metadata[[column]]# this coloumn needs to be set as described above (the number of domains and beta)
st_norm = st_norm[,rownames(spatial_location)]

matrix <- data_to_transform(sc_norm,st_norm,spatial_genes_selected,batch=TRUE,pca_method='prcomp_irlba')
train_set <- cbind(as.data.frame(t(matrix[,colnames(st_norm)])),label=spatial_location$HMRF_cluster)
test_set <- as.data.frame(t(matrix[,colnames(sc_norm)]))
train_set$label = as.factor(train_set$label)
# Predict spatial domain of individual cells
# This tuning step requires some time. You can adjust the cross-validation proportion using `cross_para` parameter in the `tune_parameter()` function.
parameters <- tune_parameter(train_set, test_set, kernel = "radial", scale = TRUE, class.weight = TRUE, verbose = TRUE, cross_para=4)
pred_st_svm <- PredictDomain(train_set, test_set, cost=parameters[['cross_4']][['cost']],
                             gamma=parameters[['cross_4']][['gamma']], st_svm=TRUE,verbose = FALSE)
pred_sc_svm <- PredictDomain(train_set, test_set, cost=parameters[['cross_4']][['cost']],
                             gamma=parameters[['cross_4']][['gamma']], scale = TRUE, verbose = TRUE)

sc_meta <- sc_meta[apply(attr(pred_sc_svm, "probabilities"),1,max)>0.8,] 
pred_sc_svm <- pred_sc_svm[apply(attr(pred_sc_svm, "probabilities"),1,max)>0.8]
sc_norm <- sc_norm[,rownames(sc_meta)]

# 确保基因名唯一、非空、顺序一致
rownames(sc_norm) <- make.unique(rownames(sc_norm))
rownames(st_norm) <- make.unique(rownames(st_norm))

common_genes <- intersect(rownames(sc_norm), rownames(st_norm))
sc_norm <- sc_norm[common_genes, ]
st_norm <- st_norm[common_genes, ]


cell_spot_map <- map_cell_to_spot(sc_norm=sc_norm,sc_meta=sc_meta,
                                  st_norm=st_norm,spatial_location=spatial_location,
                                  pred_sc_svm=pred_sc_svm, pred_st_svm=pred_st_svm,
                                  python_path=python_path,
                                  batch=TRUE,
                                  num_epochs=2000L,
                                  para_distance=2.0,
                                  para_density=0.5)
matched_spots <- unique(cell_spot_map$Spot)
coverage <- length(intersect(matched_spots, rownames(sct@meta.data))) / nrow(sct@meta.data)
cat("新映射覆盖率:", round(coverage * 100, 2), "%\n")

spot_neigh_list <- spatial_relation_all(spatial_location,
                                        spatial_data_type=c('honeycomb'))

sc_meta_coord <- calculate_cell_location(cell_spot_map=cell_spot_map,
                                         st_meta =spatial_location,
                                         sc_meta=sc_meta,
                                         sc_norm=sc_norm,
                                         st_norm=st_norm,
                                         parallel = TRUE,                   
                                         batch = TRUE,
                                         spot_neigh_list=spot_neigh_list,
                                         radius = 1.5)
pdf(paste0(outdir,"/","02-",gene,"-CMAP.pdf"))
color_use <- c("Infected" = "#CE4D4C",
               "Fibroblast" = "#EBC948", 
               "Mast" = "#DDBEAD",
               "Myeloid" = "#8C564B", 
               "Tcell" = "#5954A4",
               "Bystander" = "#5279BB")    
ggplot(sc_meta_coord,aes(pred_loc_x,pred_loc_y,color=celltype_new))+
  geom_point(size=2)+ 
  theme_bw()+
  scale_color_manual(values = color_use)+
  theme(panel.grid.major=element_line(colour=NA), 
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_fixed()+
  labs(x=NULL,y=NULL,color= 'Cell type')+
  theme(legend.position='right',
        legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  guides(color = guide_legend(override.aes = list(size = 4)))

cell_spot_map_annot <- cell_spot_map %>%
  left_join(
    sc_meta_coord %>%
      dplyr::select(celltype_new) %>%
      tibble::rownames_to_column("Single_cell"),
    by = "Single_cell"
  )

spot_best <- cell_spot_map_annot %>%
  group_by(Spot) %>%
  slice_max(Probability, n = 1, with_ties = FALSE) %>%
  ungroup()

spot_celltype <- spot_best %>%
  dplyr::select(Spot, celltype_new) %>%
  tibble::column_to_rownames("Spot")

sct <- AddMetaData(sct, metadata = spot_celltype)
sct@meta.data$celltype_new[is.na(sct@meta.data$celltype_new)] <- "Unassigned"
SpatialDimPlot(sct, group.by = "celltype_new", label = TRUE,pt.size.factor = 3,
               cols =c("Infected" = "#CE4D4C","Fibroblast" = "#EBC948", "Mast" = "#DDBEAD","Myeloid" = "#8C564B", "Tcell" = "#5954A4","Bystander" = "#5279BB"))+ 
  theme(
    legend.text = element_text(size = 14),         # 图例文字大小
    legend.title = element_text(size = 16, face = "bold"),  # 图例标题大小
    legend.key.size = unit(1.2, "cm"),             # 图例小圆点大小
    legend.box.spacing = unit(0.5, "cm")           # 图例项间距
  )
dev.off()

save(sct,file="~/rawdata/SCT/analysis/virus_3D_1/virus_3D_anno.RData")


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


###获取空转矩阵信息
sct <- subset(sct, subset = celltype_new != "Unassigned")
data.input = Seurat::GetAssayData(sct, slot = "data", assay = "SCT") 
meta        <- sct@meta.data
meta$labels <- meta$celltype_new
###获取meta信息
Idents(sct) <- "predicted.id" 
meta = data.frame(labels = Idents(sct),
                  row.names = names(Idents(sct)))
###获取空间位置信息
spatial.locs = Seurat::GetTissueCoordinates(sct, scale = NULL,cols = c("imagerow", "imagecol")) 
spatial.locs <- spatial.locs[,1:2]
#名字必须是x y ，否则后面CARD_deconvolution会报错
colnames(spatial.locs) <- c("x","y")
#scalefactors_json存于文件夹下
scalefactors = jsonlite::fromJSON(txt = file.path("~/rawdata/SCT/virus_3D/VG161_3D_1/spatial", 'scalefactors_json.json')) 
scalefactors = list(spot.diameter = 65, spot = scalefactors$spot_diameter_fullres, # these two information are required
                     fiducial = scalefactors$fiducial_diameter_fullres, hires = scalefactors$tissue_hires_scalef, lowres = scalefactors$tissue_lowres_scalef # these three information are not required
)


###创建CellChat对象
cellchat <- createCellChat(object = data.input, 
                           meta = meta, 
                           group.by = "labels", #定义的名字是labels
                           datatype = "spatial", #数据类型：空转
                           coordinates = data.matrix(spatial.locs),
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
save(cellchat,file = "~/rawdata/SCT/analysis/virus_3D_1/virus_3D_chat.rds")

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
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")
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
levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use = c(1,3,6), 
                 targets.use = c(1,3,6), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1,3,6), targets.use = c(1,3,6),                  
                 signaling = c("MK"), remove.isolate = FALSE)##指定通路
#取配体-受体对的输入，并以二进制形式显示表达
spatialFeaturePlot(cellchat, pairLR.use = "MDK_NCL", point.size = 1, do.binary = TRUE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)
dev.off()
