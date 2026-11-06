##### 清理运行环境 #####
rm(list = ls(all.names = TRUE))
graphics.off()
invisible(gc())

##### D3 空间 TLS 候选信号筛查 #####
# 第一轮：原始 Spatial counts、B/T 支持、组织化/GC信号及作者公开TLS基因集。
# 不自动诊断 TLS 或成熟度；不重新判定病毒感染；不保存大型 Seurat 对象。
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)

##### 只需优先修改这里 #####
input_type <- "object"  # object：读RData/rds；10x：读原始Visium目录。清理后重新加载。
input_file <- "~/rawdata/SCT/analysis/virus_3D_1/virus_3D_cluster.RData"
input_object <- "sct"
data_dir <- "~/rawdata/SCT/virus_3D/VG161_3D_1" # 原SCT.R已使用的D3目录。
# D5路径模板（请按服务器实际目录修改后使用）：
# input_file <- "~/rawdata/SCT/analysis/virus_5D_1/virus_5D_cluster.RData"
# data_dir <- "~/rawdata/SCT/virus_5D/VG161_5D_1"

sample_id <- "VG161_3D_1"
day <- "D3"
treatment <- "VG161"
side <- "injected"       # 本次按用户指定为注射侧，不根据L/R自动推断。
mouse_id <- NA_character_ # 尚未核实，不从sample_id猜小鼠身份。
humanization <- "possibly PBMC; unconfirmed"
outdir <- file.path(path.expand("~/OV/SCT"), paste0("TLS_", day), sample_id)
save_details <- FALSE     # TRUE时仅额外保存一个压缩details.rds。
image_name <- NULL        # 只有一张图时自动选择；多图对象请填写Images(sct)中的名称。
gene_prefix <- ""         # 人源symbol无前缀时留空；有前缀时例如填写GRCh38-。
mt_pattern <- "^MT-"
viral_pattern <- "^VG161[-_]"

# 第一轮只排除没有宿主UMI的spot，其余阈值先不凭经验一刀切。
min_host_UMI <- 0
min_features <- 0
max_percent_mt <- Inf
min_B_genes <- 2          # 仅为“标记支持”的探索规则，不是TLS诊断阈值。
min_T_genes <- 2

##### 加载包含原始counts、图像及坐标的空间对象 #####
if (input_type == "object") {
  if (!file.exists(path.expand(input_file))) stop("输入文件不存在，请修改input_file。")
  input_env <- new.env()
  loaded <- tryCatch(suppressWarnings(load(input_file, envir = input_env)),
                     error = function(e) NULL)
  if (is.null(loaded)) sct <- readRDS(input_file) else {
    if (!input_object %in% loaded) stop("RData没有input_object指定的对象。")
    sct <- get(input_object, envir = input_env)
  }
  rm(input_env, loaded)
} else if (input_type == "10x") {
  sct <- Load10X_Spatial(data.dir = path.expand(data_dir), assay = "Spatial",
                         filename = "filtered_feature_bc_matrix.h5",
                         slice = "slice1", filter.matrix = TRUE)
} else stop("input_type必须为object或10x。")
if (!exists("sct", inherits = FALSE) || !inherits(sct, "Seurat")) stop("没有Seurat对象sct。")
# 更新内存中的对象后再取坐标，避免旧图像对象与当前Seurat不兼容。
object_version_before <- as.character(sct@version)
sct <- SeuratObject::UpdateSeuratObject(sct)
if (!"Spatial" %in% names(sct@assays)) stop("没有Spatial assay，请勿用SCT或predictions替代。")
if (!length(Images(sct))) stop("对象没有空间图像/坐标，请改用完整空间对象或input_type='10x'。")
if (is.null(image_name)) {
  if (length(Images(sct)) != 1) stop("有多张切片，请明确填写image_name后逐片分析。")
  image_name <- Images(sct)[1]
}
if (!image_name %in% Images(sct)) stop("image_name不在对象中。")

sp <- sct[["Spatial"]]
if (inherits(sp, "Assay5")) {
  layers_use <- Layers(sp, search = "^counts($|\\.)")
  if (!length(layers_use)) stop("Spatial中没有counts层。")
  mats <- lapply(layers_use, function(x) LayerData(sp, layer = x, fast = FALSE))
  if (anyDuplicated(unlist(lapply(mats, colnames), use.names = FALSE))) {
    stop("counts层之间条码重复，不能自动合并。")
  }
  feature_by_layer <- setNames(lapply(mats, rownames), layers_use)
  layer_by_spot <- setNames(rep(layers_use, vapply(mats, ncol, integer(1))),
                            unlist(lapply(mats, colnames), use.names = FALSE))
  if (length(layers_use) > 1 || layers_use != "counts") {
    sp <- JoinLayers(sp, layers = "^counts($|\\.)", new = "counts")
  }
  ct <- LayerData(sp, layer = "counts", fast = FALSE)
  rm(mats)
} else {
  ct <- if (packageVersion("SeuratObject") >= "5.0.0") {
    GetAssayData(sct, assay = "Spatial", layer = "counts")
  } else GetAssayData(sct, assay = "Spatial", slot = "counts")
  feature_by_layer <- list(counts = rownames(ct))
  layer_by_spot <- setNames(rep("counts", ncol(ct)), colnames(ct))
}
if (is.null(rownames(ct)) || is.null(colnames(ct)) ||
    anyDuplicated(rownames(ct)) || anyDuplicated(colnames(ct))) stop("counts命名缺失或重复。")
ct <- as(Matrix(ct, sparse = TRUE), "dgCMatrix")
if (any(!is.finite(ct@x)) || any(ct@x < 0) || any(abs(ct@x-round(ct@x)) > 1e-6)) {
  stop("Spatial counts不是非负近整数，请核对原始数据来源。")
}

coords <- GetTissueCoordinates(sct, image = image_name, scale = "lowres")
if (all(c("imagerow", "imagecol") %in% names(coords))) {
  coords <- data.frame(x = coords$imagecol, y = coords$imagerow, row.names = rownames(coords))
} else if (all(c("x", "y") %in% names(coords))) {
  if ("cell" %in% names(coords)) rownames(coords) <- coords$cell
  coords <- coords[, c("x", "y"), drop = FALSE]
} else stop("未识别坐标列，请检查GetTissueCoordinates(sct)结果。")
if (anyDuplicated(rownames(coords)) || any(!is.finite(as.matrix(coords)))) stop("坐标重复或无效。")
spots_use <- intersect(colnames(ct), rownames(coords))
if (!length(spots_use)) stop("counts与切片坐标没有共同条码。")
if (!all(rownames(coords) %in% colnames(ct))) warning("部分图像spot缺少Spatial counts，未当成零。")
ct <- ct[, spots_use, drop = FALSE]
coords <- coords[spots_use, , drop = FALSE]
stopifnot(identical(colnames(ct), rownames(coords)))
selected_layers <- unique(unname(layer_by_spot[spots_use]))
selected_features <- unique(unlist(feature_by_layer[selected_layers], use.names = FALSE))

##### 固定基因集，不按每张切片的高分位数强制寻找TLS #####
# Atlas22是Science 2026作者公开markers.tsv中tls=1的22个symbol。
# 此处只复用作者基因清单，均值评分不等于其iSTAR/分割/成熟分类模型。
atlas_source <- "https://raw.githubusercontent.com/Coolgenome/Pan-Cancer_Spatial_Atlas_TLS/main/data/TLS_signature/markers.tsv"
marker_sets <- list(
  Bcell = c("MS4A1", "CD79A", "CD79B", "CD19", "CD22"),
  Tcell = c("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"),
  Organization = c("CXCL13", "CCL19", "CCL21"),
  GC_related = c("AICDA", "BCL6", "RGS13", "MEF2B", "SERPINA9"),
  Atlas22 = c("CD19", "CD1D", "CD22", "CD3D", "CD3E", "CD52", "CD74", "CD79A",
              "CD79B", "CD8A", "CD8B", "CR2", "CXCL13", "CXCR5", "FCER2", "IL7R",
              "ITGAE", "MS4A1", "PDCD1", "PTGDS", "CD4", "TRBC2"))
marker_sets$Lymphoid17 <- setdiff(marker_sets$Atlas22,
                                  c("CR2", "CXCL13", "CXCR5", "FCER2", "PDCD1"))
# MKI67仅作定位参考，不加入TLS/GC评分，避免把任意肿瘤增殖灶视为GC。
all_symbols <- unique(c(unlist(marker_sets, use.names = FALSE), "MKI67"))
features_use <- paste0(gene_prefix, all_symbols)
marker_audit <- data.frame(symbol = all_symbols, feature = features_use,
                           in_counts = features_use %in% selected_features)
marker_audit$n_spots_assessable <- vapply(features_use, function(g) {
  available <- vapply(feature_by_layer, function(x) g %in% x, logical(1))
  sum(available[layer_by_spot[spots_use]])
}, integer(1))
marker_audit$total_UMI <- NA_real_
marker_audit$n_spots_detected <- NA_real_
idx <- which(marker_audit$in_counts)
marker_audit$total_UMI[idx] <- Matrix::rowSums(ct[features_use[idx], , drop = FALSE])
marker_audit$n_spots_detected[idx] <- Matrix::rowSums(ct[features_use[idx], , drop = FALSE] > 0)

##### 从Spatial原始counts计算QC和宿主归一化表达 #####
host_features <- rownames(ct)[!grepl(viral_pattern, rownames(ct))]
mt_features <- grep(mt_pattern, rownames(ct), value = TRUE)
spotper <- data.frame(spot = colnames(ct), x = coords$x, y = coords$y,
                      counts_layer = unname(layer_by_spot[colnames(ct)]),
                      total_UMI = Matrix::colSums(ct),
                      host_UMI = Matrix::colSums(ct[host_features, , drop = FALSE]),
                      nFeature = Matrix::colSums(ct > 0))
spotper$percent_mt_host <- if (length(mt_features)) {
  100 * Matrix::colSums(ct[mt_features, , drop = FALSE]) / pmax(spotper$host_UMI, 1)
} else NA_real_
spotper$QC_pass <- spotper$host_UMI > min_host_UMI & spotper$nFeature > min_features
if (is.finite(max_percent_mt)) {
  if (anyNA(spotper$percent_mt_host)) stop("线粒体指标缺失，不能使用max_percent_mt过滤。")
  spotper$QC_pass <- spotper$QC_pass & spotper$percent_mt_host < max_percent_mt
}
if (!any(spotper$QC_pass)) stop("没有通过QC的spot，请核对counts和阈值。")

present_markers <- intersect(features_use, selected_features)
norm_marker <- ct[present_markers, , drop = FALSE] %*%
  Diagonal(x = 1e4 / pmax(spotper$host_UMI, 1))
norm_marker <- log1p(norm_marker)
dimnames(norm_marker) <- list(present_markers, colnames(ct))
coverage <- coverage_by_layer <- list()
for (module in names(marker_sets)) {
  symbols <- marker_sets[[module]]
  features <- paste0(gene_prefix, symbols)
  present <- intersect(features, selected_features)
  complete_layer <- vapply(feature_by_layer, function(x) all(features %in% x), logical(1))
  assessable <- spotper$QC_pass & complete_layer[spotper$counts_layer]
  score <- detected <- marker_UMI <- rep(NA_real_, ncol(ct))
  if (length(present)) {
    score[assessable] <- Matrix::colSums(norm_marker[present, , drop = FALSE])[assessable] / length(features)
    detected[assessable] <- Matrix::colSums(ct[present, , drop = FALSE] > 0)[assessable]
    marker_UMI[assessable] <- Matrix::colSums(ct[present, , drop = FALSE])[assessable]
  }
  spotper[[paste0(module, "_meanlog")]] <- score
  spotper[[paste0(module, "_genes")]] <- detected
  spotper[[paste0(module, "_UMI")]] <- marker_UMI
  coverage[[module]] <- data.frame(module = module, expected_n = length(features),
                                   present_n = length(present),
                                   n_complete_layers = sum(complete_layer[selected_layers]),
                                   n_selected_layers = length(selected_layers),
                                   n_assessable_spots = sum(assessable),
                                   missing = paste(setdiff(features, present), collapse = ";"))
  coverage_by_layer[[module]] <- bind_rows(lapply(selected_layers, function(l) {
    data.frame(module = module, layer = l, expected_n = length(features),
               present_n = sum(features %in% feature_by_layer[[l]]),
               missing = paste(setdiff(features, feature_by_layer[[l]]), collapse = ";"))
  }))
}
coverage <- bind_rows(coverage)
coverage_by_layer <- bind_rows(coverage_by_layer)
if (any(coverage$n_complete_layers < coverage$n_selected_layers)) {
  warning("有模块缺少基因，完整固定面板不可评估的spot记NA；请查看info.txt，不要把NA解释为无TLS。")
}

##### B/T标记支持只用于筛查，不自动圈定TLS区域 #####
spotper$BT_support <- "Not_assessable"
ok <- spotper$QC_pass & !is.na(spotper$Bcell_genes) & !is.na(spotper$Tcell_genes)
spotper$BT_support[ok] <- "No_BT_support"
spotper$BT_support[ok & spotper$Bcell_genes >= min_B_genes] <- "B_support"
spotper$BT_support[ok & spotper$Tcell_genes >= min_T_genes] <- "T_support"
spotper$BT_support[ok & spotper$Bcell_genes >= min_B_genes &
                     spotper$Tcell_genes >= min_T_genes] <- "BT_co_detection"
spotper$BT_support[!spotper$QC_pass] <- "QC_excluded"
# 同spot共检出不能证明B/T分区；B/T分处相邻spot时也可能没有共检出。

##### 精简汇总：一张小CSV #####
summary_table <- data.frame(sample_id = sample_id, day = day, treatment = treatment,
                            side = side, mouse_id = mouse_id,
                            n_spots_image = nrow(coords), n_spots_QC = sum(spotper$QC_pass),
                            host_UMI_median = median(spotper$host_UMI[spotper$QC_pass]),
                            n_BT_assessable = sum(ok),
                            n_B_support = if (all(is.na(spotper$Bcell_genes))) NA_integer_ else
                              sum(spotper$Bcell_genes >= min_B_genes, na.rm = TRUE),
                            n_T_support = if (all(is.na(spotper$Tcell_genes))) NA_integer_ else
                              sum(spotper$Tcell_genes >= min_T_genes, na.rm = TRUE),
                            n_BT_co_detection = if (!any(ok)) NA_integer_ else
                              sum(spotper$BT_support == "BT_co_detection"),
                            treatment_comparison_performed = FALSE)
for (module in names(marker_sets)) {
  v <- spotper[[paste0(module, "_meanlog")]]
  summary_table[[paste0(module, "_genes_present")]] <- coverage$present_n[coverage$module == module]
  summary_table[[paste0(module, "_genes_expected")]] <- length(marker_sets[[module]])
  summary_table[[paste0(module, "_spots_assessable")]] <- sum(!is.na(v))
  summary_table[[paste0(module, "_score_median")]] <- if (all(is.na(v))) NA_real_ else median(v, na.rm = TRUE)
  summary_table[[paste0(module, "_spots_any_marker")]] <- if (all(is.na(v))) NA_integer_ else
    sum(spotper[[paste0(module, "_genes")]] > 0, na.rm = TRUE)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE)
print(summary_table)

##### 空间图：仅在首图嵌入H&E，控制PDF体积 #####
rownames(spotper) <- spotper$spot
support_cols <- c(No_BT_support = "#D9D9D9", B_support = "#4575B4", T_support = "#1B9E77",
                  BT_co_detection = "#D73027", Not_assessable = "#756BB1", QC_excluded = "#333333")
# 直接使用公开接口返回的原始图像和lowres坐标叠图。
# 兼容旧Seurat对VisiumV1误报UpdateSeuratObject的情况，不手改图像slot。
he_image <- GetImage(sct[[image_name]], mode = "raster")
if (length(dim(he_image)) != 2) stop("H&E没有返回二维raster，请核对图像对象。")
he_height <- nrow(he_image)
he_width <- ncol(he_image)
support_data <- data.frame(x = spotper$x, y = -spotper$y,
                           support = factor(spotper$BT_support, levels = names(support_cols)))
p_support <- ggplot(support_data, aes(x, y)) +
  annotation_raster(he_image, xmin = 0, xmax = he_width, ymin = -he_height, ymax = 0) +
  geom_point(aes(color = support), size = 1.2, alpha = 0.8) +
  scale_color_manual(values = support_cols, name = "B/T support") +
  coord_fixed(xlim = c(0, he_width), ylim = c(-he_height, 0), expand = FALSE) +
  theme_void() + labs(title = paste(sample_id, "B/T marker support: inspect with H&E"))
map_value <- function(values, title) {
  dat <- data.frame(x = spotper$x, y = spotper$y, value = values)
  upper <- if (all(is.na(values))) 1 else max(values, na.rm = TRUE)
  if (!is.finite(upper) || upper <= 0) upper <- 1
  ggplot(dat, aes(x, y, color = value)) + geom_point(size = 1.1) +
    scale_color_gradient(low = "#EEEEEE", high = "#B2182B", na.value = "#808080",
                         limits = c(0, upper)) +
    scale_y_reverse() + coord_fixed() + theme_void() +
    labs(title = title, color = NULL) + theme(plot.title = element_text(size = 11))
}
pdf(file.path(outdir, "plots.pdf"), width = 11, height = 8, compress = TRUE)
tryCatch({
  print(p_support)
  p1 <- map_value(log1p(spotper$host_UMI), "Host UMI: log1p")
  p2 <- map_value(spotper$nFeature, "Detected features")
  print(p1 + p2)
  module_plots <- lapply(names(marker_sets), function(m) {
    map_value(spotper[[paste0(m, "_meanlog")]], paste(m, "mean logCP10k (host)"))
  })
  print(wrap_plots(module_plots, ncol = 2))
  key_markers <- c("CD3D", "MS4A1", "CR2", "FCER2", "BCL6", "MKI67",
                   "CXCL13", "CCL19", "CCL21", "AICDA")
  marker_plots <- lapply(key_markers, function(g) {
    feature <- paste0(gene_prefix, g)
    v <- rep(NA_real_, ncol(ct))
    if (feature %in% rownames(norm_marker)) {
      v <- as.numeric(norm_marker[feature, ])
      available <- vapply(feature_by_layer, function(x) feature %in% x, logical(1))
      v[!available[spotper$counts_layer] | !spotper$QC_pass] <- NA_real_
    }
    map_value(v, paste(g, "logCP10k (host)"))
  })
  print(wrap_plots(marker_plots[1:6], ncol = 2))
  print(wrap_plots(marker_plots[7:10], ncol = 2))
}, finally = dev.off()) # 即使某页绘图失败也关闭PDF，便于修正后重跑。

##### 参数、时间、覆盖及参考集中放进一个info #####
capture.output({
  cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
  print(list(input_type = input_type, input_file = input_file, data_dir = data_dir,
             sample_id = sample_id, day = day, treatment = treatment, side = side,
             mouse_id = mouse_id, humanization = humanization,
             image_name = image_name, gene_prefix = gene_prefix,
             object_version_before = object_version_before, object_version_after = as.character(sct@version),
             HE_plot = "GetImage raster + lowres coordinates; x=imagecol, y=imagerow",
             min_host_UMI = min_host_UMI, min_features = min_features,
             max_percent_mt = max_percent_mt, min_B_genes = min_B_genes, min_T_genes = min_T_genes))
  cat("\nPanel coverage: a score requires every predefined gene in that counts layer\n")
  print(coverage, row.names = FALSE)
  cat("\nPanel coverage by original selected counts layer\n")
  print(coverage_by_layer, row.names = FALSE)
  cat("\nMarker audit: selected image, before optional spot QC\n")
  print(marker_audit, row.names = FALSE)
  cat("\nScore = mean(log1p(raw marker UMI / host UMI * 10000)).\n")
  cat("This is a descriptive marker score, not calibrated TLS probability or cell proportion.\n")
  cat("Atlas22 gene list source:", atlas_source, "\n")
  cat("Atlas paper: https://doi.org/10.1126/science.adz2742\n")
  cat("Atlas uses inferred high-resolution expression plus spatial segmentation; that model is not reproduced here.\n")
  cat("No automatic TLS count, maturity label, tissue-area estimate, or treatment P value was generated.\n")
  cat("Single-sample screening; add matched Vehicle and independent animals before evaluating treatment change.\n")
  cat("Same-spot B/T support does not establish structure; lack of co-detection does not exclude adjacent B/T zones.\n")
  cat("Species reference is unconfirmed. Missing human stromal/GC markers do not prove biological absence.\n")
  cat("Panel score color scales are feature-specific. Do not compare colors across different modules as cell fractions.\n\n")
  sessionInfo()
}, file = file.path(outdir, "info.txt"))
if (save_details) saveRDS(list(spots = spotper, markers = marker_audit, coverage = coverage,
                               coverage_by_layer = coverage_by_layer,
                               marker_sets = marker_sets, summary = summary_table),
                          file.path(outdir, "details.rds"), compress = "gzip")
message("输出目录：", normalizePath(outdir, winslash = "/", mustWork = FALSE))
message("先检查B/T是否有支持，再结合H&E查看区域；不能把红色spot直接称作TLS。")
