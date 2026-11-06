##### 清理运行环境 #####
rm(list = ls(all.names = TRUE))
graphics.off()
invisible(gc())

##### D3注射侧：样本整理、统一TLS评分、待复核ROI #####
# 本版只做步骤1—3；候选ROI不等于确诊TLS，不自动判断成熟度。
# 不保存Seurat对象或全基因表达矩阵；默认输出summary.csv、plots.pdf、results.rds。
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)
if (!requireNamespace("AUCell", quietly = TRUE)) {
  stop("请先安装AUCell：BiocManager::install('AUCell', ask = FALSE, update = FALSE)")
}

##### 输入集中修改：目前Vehicle、VG161各一张D3注射侧切片 #####
data_root <- "~/rawdata/SCT/virus_3D"
outdir <- "~/OV/SCT/TLS_D3"
samples <- data.frame(
  sample_id = c("Vehicle_3D_1", "VG161_3D_1"),
  treatment = c("Vehicle", "VG161"),
  day = "D3", side = "injected",
  mouse_id = NA_character_, donor_id = NA_character_, batch = NA_character_,
  input_type = "10x", # 推荐直接读取各自原始Visium目录；也支持object。
  data_dir = file.path(data_root, c("Vehicle_3D_1", "VG161_3D_1")),
  object_file = c(NA_character_, "~/rawdata/SCT/analysis/virus_3D_1/virus_3D_cluster.RData"),
  object_name = "sct", image_name = NA_character_, barcode_prefix = "",
  stringsAsFactors = FALSE)
# 如改用object，填写每行自己的object_file；data_dir仍用于原始坐标、比例尺和H&E。
# D5/其他侧别：以后修改data_root、samples和outdir；本版不混入D5或对侧。

gene_prefix <- ""                    # 人源symbol无前缀时留空；不自动改变大小写或合并物种。
viral_pattern <- "^VG161[-_]"
mt_pattern <- "^MT-"
min_host_UMI <- 0                    # 与前一轮一致，先排除空spot。
min_features <- 0
max_percent_mt <- Inf
auc_fraction <- 0.05                 # 所有样本同一个共有宿主基因全集、同一个aucMaxRank。
auc_chunk <- 128L                    # 分块计算排名，降低峰值内存。
seed <- 2026L

# 统一的探索性筛查参数；不是论文验证的TLS诊断阈值。
# 背景为两组共用、按表达均值/检出率匹配的随机基因集，不强制每片产生候选。
n_null <- 199L
null_match_pool <- 200L
null_tail_cutoff <- 0.01
min_TLS_genes <- 3L
min_ROI_spots <- 3L                  # 直接六邻接；小于此值的离散点仍保存在spot表。
imprint_min_coverage <- 0.80          # 辅助面板允许部分覆盖，但固定使用各样本共有的同一批基因。

# 现有H&E人工ROI可选：CSV列为sample_id,roi_id,x,y，每个多边形至少3个顶点。
# 也可填上次results.rds路径，复用其中HE_polygons及已完成复核的记录。
# x/y必须为原始fullres像素坐标，原点左上；不要填写PDF页坐标。
HE_roi_file <- ""
draw_HE_ROI <- FALSE                 # TRUE：先在交互式R/RStudio的匿名H&E画面圈图，再看表达图。
image_max_px <- 1200L                # 只压缩PDF预览；不改动原始H&E。

##### 固定基因集及来源 #####
marker_sets <- list(
  TLS22 = c("CD19", "CD1D", "CD22", "CD3D", "CD3E", "CD52", "CD74", "CD79A",
            "CD79B", "CD8A", "CD8B", "CR2", "CXCL13", "CXCR5", "FCER2", "IL7R",
            "ITGAE", "MS4A1", "PDCD1", "PTGDS", "CD4", "TRBC2"),
  Imprint29 = c("IGHA1", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHGP", "IGHM",
                "IGKC", "IGLC1", "IGLC2", "IGLC3", "JCHAIN", "CD79A", "FCRL5",
                "MZB1", "SSR4", "XBP1", "TRBC2", "IL7R", "CXCL12", "LUM",
                "C1QA", "C7", "CD52", "APOE", "PTLP", "PTGDS", "PIM2", "DERL3"),
  Bcell = c("MS4A1", "CD79A", "CD79B", "CD19", "CD22"),
  Tcell = c("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"),
  Organization = c("CXCL13", "CCL19", "CCL21"),
  GC_related = c("AICDA", "BCL6", "RGS13", "MEF2B", "SERPINA9"))
# HNSCC原文写PTLP；作者代码未核实更正，不猜测替换成PLTP。
# 如PTLP不在矩阵，Imprint29明确按实际覆盖率报告为部分面板，不称完整29基因复现。
sources <- c(
  TLS22 = "https://raw.githubusercontent.com/Coolgenome/Pan-Cancer_Spatial_Atlas_TLS/main/data/TLS_signature/markers.tsv",
  Science = "https://doi.org/10.1126/science.adz2742",
  Imprint29 = "https://doi.org/10.1038/s41467-025-59341-w",
  AUCell = "https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html")

##### 小函数：加载一个样本，只取它自己的原始counts和配套空间信息 #####
read_sample <- function(i) {
  s <- samples[i, ]
  pos_dir <- file.path(path.expand(s$data_dir), "spatial")
  if (s$input_type == "10x") {
    sct <- Load10X_Spatial(data.dir = path.expand(s$data_dir), assay = "Spatial",
                           filename = "filtered_feature_bc_matrix.h5", filter.matrix = TRUE)
  } else {
    input_env <- new.env()
    loaded <- tryCatch(suppressWarnings(load(path.expand(s$object_file), envir = input_env)),
                       error = function(e) NULL)
    if (is.null(loaded)) sct <- readRDS(path.expand(s$object_file)) else {
      if (!s$object_name %in% loaded) stop(s$sample_id, "：RData中没有", s$object_name)
      sct <- get(s$object_name, envir = input_env)
    }
    rm(input_env, loaded)
  }
  if (!inherits(sct, "Seurat") || !"Spatial" %in% names(sct@assays)) {
    stop(s$sample_id, "：需要含原始Spatial counts的Seurat对象。")
  }
  sct <- SeuratObject::UpdateSeuratObject(sct)
  image_names <- Images(sct)
  image_name <- s$image_name
  if (is.na(image_name) || !nzchar(image_name)) {
    if (length(image_names) != 1) stop(s$sample_id, "：请明确指定image_name。")
    image_name <- image_names[1]
  }
  if (!image_name %in% image_names) stop(s$sample_id, "：image_name不存在。")
  image_cells <- Cells(sct[[image_name]])
  sp <- sct[["Spatial"]]
  if (inherits(sp, "Assay5")) {
    ln <- Layers(sp, search = "^counts($|\\.)")
    if (!length(ln)) stop(s$sample_id, "：没有counts层。")
    mats <- lapply(ln, function(l) LayerData(sp, layer = l, fast = FALSE))
    if (anyDuplicated(unlist(lapply(mats, colnames), use.names = FALSE))) {
      stop(s$sample_id, "：counts层条码重复，不能自动合并。")
    }
    use <- vapply(mats, function(x) any(colnames(x) %in% image_cells), logical(1))
    mats <- mats[use]
    if (!length(mats)) stop(s$sample_id, "：切片没有对应counts。")
    # 层间缺失基因不能补零：仅使用所选切片各层都保留的基因。
    feature_common <- Reduce(intersect, lapply(mats, rownames))
    mats <- lapply(mats, function(x) x[feature_common, intersect(colnames(x), image_cells), drop = FALSE])
    ct <- do.call(cbind, mats)
    layers_used <- ln[use]
  } else {
    ct <- if (packageVersion("SeuratObject") >= "5.0.0") {
      GetAssayData(sct, assay = "Spatial", layer = "counts")
    } else GetAssayData(sct, assay = "Spatial", slot = "counts")
    ct <- ct[, intersect(colnames(ct), image_cells), drop = FALSE]
    layers_used <- "counts"
  }
  ct <- as(Matrix(ct, sparse = TRUE), "dgCMatrix")
  if (is.null(rownames(ct)) || is.null(colnames(ct)) || anyDuplicated(rownames(ct)) ||
      anyDuplicated(colnames(ct)) || any(!is.finite(ct@x)) || any(ct@x < 0) ||
      any(abs(ct@x-round(ct@x)) > 1e-6)) stop(s$sample_id, "：counts命名/数值无效。")
  
  # 使用原始Visium位置表的阵列坐标建邻接，不用k近邻跨接孔洞。
  pos_file <- file.path(pos_dir, "tissue_positions.csv")
  if (file.exists(pos_file)) pos <- read.csv(pos_file, stringsAsFactors = FALSE) else {
    pos_file <- file.path(pos_dir, "tissue_positions_list.csv")
    pos <- read.csv(pos_file, header = FALSE, stringsAsFactors = FALSE)
  }
  if (ncol(pos) != 6) stop(s$sample_id, "：需要标准Visium六列位置表；HD请单独适配。")
  names(pos) <- c("spot", "in_tissue", "array_row", "array_col", "y", "x")
  pos$spot <- paste0(s$barcode_prefix, pos$spot)
  if (anyDuplicated(pos$spot) || any(!is.finite(as.matrix(pos[, 2:6])))) stop("位置表无效。")
  pos <- pos[pos$in_tissue == 1 & pos$spot %in% colnames(ct), ]
  if (!nrow(pos)) stop(s$sample_id, "：位置表与counts没有共同的组织内条码。")
  pos <- pos[order(pos$array_row, pos$array_col), ]
  ct <- ct[, pos$spot, drop = FALSE]
  # 额外核对更新后的Seurat坐标，防止对象与data_dir指向不同切片。
  obj_xy <- GetTissueCoordinates(sct, image = image_name, scale = NULL)
  if (all(c("imagecol", "imagerow") %in% names(obj_xy))) {
    obj_xy <- data.frame(x = obj_xy$imagecol, y = obj_xy$imagerow, row.names = rownames(obj_xy))
  } else {
    if ("cell" %in% names(obj_xy)) rownames(obj_xy) <- obj_xy$cell
    obj_xy <- obj_xy[, c("x", "y"), drop = FALSE]
  }
  if (!all(pos$spot %in% rownames(obj_xy)) ||
      any(abs(as.matrix(obj_xy[pos$spot, c("x", "y")])-as.matrix(pos[, c("x", "y")])) > 1)) {
    stop(s$sample_id, "：对象与原始位置表坐标不一致，请检查data_dir/image_name。")
  }
  sf <- jsonlite::fromJSON(file.path(pos_dir, "scalefactors_json.json"))
  he_file <- file.path(pos_dir, "tissue_hires_image.png")
  he_scale <- sf$tissue_hires_scalef
  if (!file.exists(he_file)) {
    he_file <- file.path(pos_dir, "tissue_lowres_image.png")
    he_scale <- sf$tissue_lowres_scalef
  }
  if (length(he_scale) != 1 || !is.finite(he_scale) || he_scale <= 0) stop("H&E比例尺无效。")
  he <- png::readPNG(he_file)
  he_dim <- dim(he)[1:2]
  if (any(pos$x * he_scale > he_dim[2] + 2) || any(pos$y * he_scale > he_dim[1] + 2) ||
      any(pos$x < 0 | pos$y < 0)) stop(s$sample_id, "：坐标超出H&E范围。")
  rm(he, sct, sp)
  
  host <- !grepl(viral_pattern, rownames(ct))
  mt <- grepl(mt_pattern, rownames(ct)) & host
  met <- pos
  met$host_UMI <- Matrix::colSums(ct[host, , drop = FALSE])
  met$nFeature_host <- Matrix::colSums(ct[host, , drop = FALSE] > 0)
  met$percent_mt <- if (any(mt)) 100 * Matrix::colSums(ct[mt, , drop = FALSE]) / pmax(met$host_UMI, 1) else NA_real_
  met$QC_pass <- met$host_UMI > min_host_UMI & met$nFeature_host > min_features
  if (is.finite(max_percent_mt)) {
    if (anyNA(met$percent_mt)) stop(s$sample_id, "：线粒体基因缺失，不能使用该过滤阈值。")
    met$QC_pass <- met$QC_pass & met$percent_mt < max_percent_mt
  }
  list(counts = ct, spots = met, layers = layers_used,
       image = list(file = he_file, scale = he_scale, height = he_dim[1], width = he_dim[2],
                    position_md5 = unname(tools::md5sum(pos_file)), image_md5 = unname(tools::md5sum(he_file))))
}

##### 小函数：真实六邻域、连通候选和保留孔洞的网格轮廓 #####
make_graph <- function(pos) {
  if (anyDuplicated(paste(pos$array_row, pos$array_col))) stop("阵列位置重复。")
  keys <- paste(pos$array_row, pos$array_col, sep = ":")
  offsets <- rbind(c(0, 2), c(1, 1), c(1, -1))
  edges <- bind_rows(lapply(seq_len(nrow(offsets)), function(k) {
    to <- match(paste(pos$array_row + offsets[k, 1], pos$array_col + offsets[k, 2], sep = ":"), keys)
    data.frame(from = pos$spot[!is.na(to)], to = pos$spot[to[!is.na(to)]])
  }))
  igraph::graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = pos$spot))
}
find_regions <- function(g, flag, min_size = min_ROI_spots) {
  out <- rep(NA_integer_, igraph::vcount(g))
  if (!any(flag)) return(out)
  subg <- igraph::induced_subgraph(g, vids = which(flag))
  comp <- igraph::components(subg)
  keep <- which(comp$csize >= min_size)
  out[match(names(comp$membership), igraph::V(g)$name)] <- match(comp$membership, keep)
  out
}
grid_edges <- function(pos) {
  if (!nrow(pos)) return(data.frame())
  # 阵列坐标中的整数顶点；共享边删除后保留凹陷和内部孔洞，不画凸包。
  off <- rbind(c(3, 1), c(0, 2), c(-3, 1), c(-3, -1), c(0, -2), c(3, -1))
  verts <- lapply(seq_len(nrow(pos)), function(i) {
    z <- sweep(off, 2, c(3*pos$array_col[i], 3*pos$array_row[i]), "+")
    a <- paste(z[, 1], z[, 2], sep = ":")
    b <- c(a[-1], a[1])
    data.frame(c1 = z[, 1]/3, r1 = z[, 2]/3,
               c2 = c(z[-1, 1], z[1, 1])/3, r2 = c(z[-1, 2], z[1, 2])/3,
               key = ifelse(a < b, paste(a, b), paste(b, a)))
  })
  e <- bind_rows(verts)
  e[!duplicated(e$key) & !duplicated(e$key, fromLast = TRUE), ]
}

##### 第1步：输入预检查及两组共同的基因全集 #####
stopifnot(!anyDuplicated(samples$sample_id), all(samples$day == "D3"),
          all(samples$side == "injected"), setequal(samples$treatment, c("Vehicle", "VG161")))
if (any(!samples$input_type %in% c("10x", "object"))) stop("input_type须为10x或object。")
if (any(!dir.exists(path.expand(samples$data_dir)))) stop("以下目录不存在：", paste(samples$data_dir[!dir.exists(path.expand(samples$data_dir))], collapse = "; "))
dir_keys <- vapply(path.expand(samples$data_dir), normalizePath, "", winslash = "/", mustWork = TRUE)
if (.Platform$OS.type == "windows") dir_keys <- tolower(dir_keys)
if (anyDuplicated(dir_keys)) stop("各样本data_dir重复，请勿将同一片改名当成另一组。")
obj_rows <- samples$input_type == "object"
if (any(obj_rows & (is.na(samples$object_file) | !file.exists(path.expand(samples$object_file))))) stop("请填写每个object样本自己的object_file。")
if (any(obj_rows)) {
  obj_keys <- vapply(path.expand(samples$object_file[obj_rows]), normalizePath, "", winslash = "/", mustWork = TRUE)
  if (.Platform$OS.type == "windows") obj_keys <- tolower(obj_keys)
  obj_images <- samples$image_name[obj_rows]
  obj_images[is.na(obj_images) | !nzchar(obj_images)] <- "AUTO"
  # 同一多图对象仅在每行明确指定不同切片时才允许复用。
  for (key in unique(obj_keys)) {
    ii <- which(obj_keys == key)
    if (length(ii) > 1 && ("AUTO" %in% obj_images[ii] || anyDuplicated(obj_images[ii]))) {
      stop("重复使用同一个object文件/图像，不能将同一片标成两组。")
    }
  }
}
if (n_null < 99 || null_tail_cutoff < 1/(n_null+1)) stop("背景集数量不足以支持当前尾部阈值。")
if (auc_fraction <= 0 || auc_fraction >= 1 || auc_chunk < 1 || min_ROI_spots < 1) stop("评分/ROI参数无效。")
set.seed(seed)
samples$alias <- sample(sprintf("S%02d", seq_len(nrow(samples))))
print(samples[, c("sample_id", "treatment", "input_type", "data_dir")])
feature_lists <- stat_lists <- precheck <- images <- vector("list", nrow(samples))
for (i in seq_len(nrow(samples))) {
  message("预检查：", samples$sample_id[i])
  x <- read_sample(i)
  host <- rownames(x$counts)[!grepl(viral_pattern, rownames(x$counts)) & startsWith(rownames(x$counts), gene_prefix)]
  feature_lists[[i]] <- host
  good <- which(x$spots$QC_pass)
  means <- detected <- setNames(rep(0, length(host)), host)
  if (length(good)) {
    h <- x$counts[host, good, drop = FALSE]
    means <- Matrix::rowMeans(h %*% Diagonal(x = 1e4 / x$spots$host_UMI[good]))
    detected <- Matrix::rowMeans(h > 0)
    names(means) <- names(detected) <- host
  }
  stat_lists[[i]] <- list(mean = means, detected = detected, eligible = length(good) > 0)
  images[[i]] <- x$image
  precheck[[i]] <- data.frame(sample_id = samples$sample_id[i], n_spots = ncol(x$counts),
                              n_QC = length(good), n_host_features = length(host),
                              counts_layers = paste(x$layers, collapse = ";"))
  rm(x)
  invisible(gc())
}
precheck <- bind_rows(precheck)
common_genes <- sort(Reduce(intersect, feature_lists))
aucMaxRank <- ceiling(length(common_genes) * auc_fraction)
if (aucMaxRank < 2 || length(common_genes) < 100) stop("共同基因全集过小，请检查基因命名/物种前缀。")
full_sets <- lapply(marker_sets, function(g) paste0(gene_prefix, g))
gene_sets <- lapply(full_sets, intersect, y = common_genes)
coverage <- data.frame(module = names(full_sets), expected = lengths(full_sets), used = lengths(gene_sets),
                       missing = vapply(full_sets, function(g) paste(setdiff(g, common_genes), collapse = ";"), ""))
if (length(gene_sets$TLS22) != 22) stop("TLS22未完整保留：", coverage$missing[coverage$module == "TLS22"], "。请使用未筛基因的原始counts。")
score_ok <- coverage$used == coverage$expected
score_ok[coverage$module == "Imprint29"] <- length(gene_sets$Imprint29)/29 >= imprint_min_coverage
names(score_ok) <- coverage$module
if (any(coverage$used < coverage$expected)) warning("部分辅助面板覆盖不足；实际使用基因和不可评价模块记录在results.rds$coverage。")
print(coverage)
message("共同宿主基因：", length(common_genes), "；统一aucMaxRank：", aucMaxRank)

##### 与分组无关的共同背景：匹配基因表达均值和检出率 #####
usable <- which(vapply(stat_lists, function(x) x$eligible, logical(1)))
if (!length(usable)) stop("两组均无通过基础QC的spot。")
mean_cp <- Reduce(`+`, lapply(stat_lists[usable], function(x) x$mean[common_genes])) / length(usable)
det_rate <- Reduce(`+`, lapply(stat_lists[usable], function(x) x$detected[common_genes])) / length(usable)
mean_rank <- rank(mean_cp, ties.method = "average")/length(common_genes)
det_rank <- rank(det_rate, ties.method = "average")/length(common_genes)
background <- setdiff(common_genes, unique(unlist(full_sets, use.names = FALSE)))
if (length(background) < max(50L, length(gene_sets$TLS22)*2)) stop("可匹配的背景基因不足。")
match_pools <- lapply(gene_sets$TLS22, function(g) {
  z <- match(g, common_genes)
  b <- match(background, common_genes)
  distance <- (mean_rank[b]-mean_rank[z])^2 + (det_rank[b]-det_rank[z])^2
  background[order(distance, background)][seq_len(min(null_match_pool, length(background))) ]
})
set.seed(seed)
null_sets <- lapply(seq_len(n_null), function(k) {
  chosen <- character()
  for (pool in match_pools) {
    available <- setdiff(pool, chosen)
    if (!length(available)) stop("匹配池太小，背景基因集无法保持不重复。")
    chosen <- c(chosen, available[sample.int(length(available), 1)])
  }
  chosen
})
names(null_sets) <- sprintf("Null%03d", seq_len(n_null))
null_matrix <- do.call(rbind, null_sets)
background_matching <- bind_rows(lapply(seq_along(gene_sets$TLS22), function(j) {
  target <- match(gene_sets$TLS22[j], common_genes)
  chosen <- match(null_matrix[, j], common_genes)
  data.frame(feature = gene_sets$TLS22[j], target_mean_CP10k = mean_cp[target],
             target_detection = det_rate[target], null_mean_CP10k_median = median(mean_cp[chosen]),
             null_detection_median = median(det_rate[chosen]),
             abs_mean_rank_difference_median = median(abs(mean_rank[chosen]-mean_rank[target])),
             abs_detection_rank_difference_median = median(abs(det_rank[chosen]-det_rank[target])))
}))
rownames(background_matching) <- NULL

##### H&E独立候选：可读入多边形，也可先在匿名H&E画面圈图 #####
HE_polygons <- data.frame(sample_id = character(), roi_id = character(), x = numeric(), y = numeric())
HE_reviewed <- character()
if (nzchar(HE_roi_file)) {
  if (grepl("\\.rds$", HE_roi_file, ignore.case = TRUE)) {
    previous <- readRDS(path.expand(HE_roi_file))
    if (is.null(previous$HE_polygons) || is.null(previous$samples)) stop("该RDS不是本流程的results.rds。")
    old_idx <- match(samples$sample_id, previous$samples$sample_id)
    if (anyNA(old_idx) || any(path.expand(samples$data_dir) != path.expand(previous$samples$data_dir[old_idx]))) {
      stop("旧ROI的样本/目录与当前输入不一致，不能自动复用。")
    }
    for (i in seq_len(nrow(samples))) {
      old_image <- previous$images[[old_idx[i]]]
      if (is.null(old_image$position_md5) || is.null(old_image$image_md5) ||
          !identical(old_image$position_md5, images[[i]]$position_md5) ||
          !identical(old_image$image_md5, images[[i]]$image_md5)) {
        stop("旧ROI的坐标表/H&E指纹与当前输入不一致，不能自动复用。")
      }
    }
    HE_polygons <- previous$HE_polygons
    HE_reviewed <- previous$summary$sample_id[which(previous$summary$HE_review_done)]
    rm(previous)
  } else {
    HE_polygons <- read.csv(path.expand(HE_roi_file), stringsAsFactors = FALSE)
    HE_reviewed <- unique(HE_polygons$sample_id)
  }
  if (!all(c("sample_id", "roi_id", "x", "y") %in% names(HE_polygons))) stop("H&E ROI文件列名不完整。")
  if (any(!HE_polygons$sample_id %in% samples$sample_id) || anyNA(HE_polygons[, c("sample_id", "roi_id", "x", "y")]) ||
      any(!is.finite(as.matrix(HE_polygons[, c("x", "y")])))) stop("H&E ROI样本名/坐标无效。")
}
if (draw_HE_ROI) {
  if (!interactive()) stop("draw_HE_ROI=TRUE需要交互式R/RStudio；批处理请使用HE_roi_file。")
  for (i in order(samples$alias)) {
    im <- images[[i]]
    he <- png::readPNG(im$file)
    w <- im$width/im$scale; h <- im$height/im$scale
    plot(NA, xlim = c(0, w), ylim = c(h, 0), asp = 1, xaxs = "i", yaxs = "i",
         xlab = "fullres x", ylab = "fullres y", main = paste(samples$alias[i], "H&E only"))
    rasterImage(he, 0, h, w, 0)
    k <- 1L
    while (paste0("draw", k) %in% HE_polygons$roi_id[HE_polygons$sample_id == samples$sample_id[i]]) k <- k+1L
    while (tolower(readline("是否圈定一个H&E候选区域？y/n：")) == "y") {
      while (paste0("draw", k) %in% HE_polygons$roi_id[HE_polygons$sample_id == samples$sample_id[i]]) k <- k+1L
      message("依次左键圈顶点；右键或Esc结束当前多边形。")
      z <- locator(type = "l")
      if (!is.null(z) && length(z$x) >= 3) {
        HE_polygons <- bind_rows(HE_polygons, data.frame(sample_id = samples$sample_id[i],
                                                         roi_id = paste0("draw", k), x = z$x, y = z$y))
        polygon(z$x, z$y, border = "#0072B2", lwd = 2)
        k <- k+1L
      }
    }
    HE_reviewed <- union(HE_reviewed, samples$sample_id[i])
    rm(he)
  }
}
if (nrow(HE_polygons)) {
  polygons <- split(HE_polygons, interaction(HE_polygons$sample_id, HE_polygons$roi_id, drop = TRUE))
  if (any(vapply(polygons, nrow, integer(1)) < 3)) stop("每个H&E多边形至少需要3个顶点。")
  for (p in polygons) {
    q <- unique(p[, c("x", "y")])
    next_row <- c(seq.int(2L, nrow(p)), 1L)
    area2 <- abs(sum(p$x * p$y[next_row] - p$y * p$x[next_row]))
    if (nrow(q) < 3 || area2 < 1e-6) stop("H&E多边形顶点重复或共线，不能形成有效区域。")
  }
  for (i in seq_len(nrow(samples))) {
    p <- HE_polygons[HE_polygons$sample_id == samples$sample_id[i], ]
    im <- images[[i]]
    if (any(p$x < 0 | p$y < 0 | p$x > im$width/im$scale | p$y > im$height/im$scale)) {
      stop(samples$sample_id[i], "：人工ROI顶点超出H&E，检查是否误用了hires/lowres或PDF坐标。")
    }
  }
}

##### 第2步：各样本用同一基因全集、同一截断值计算AUCell #####
spots_all <- markers_all <- expr_edges <- graphs <- marker_audit <- vector("list", nrow(samples))
score_sets <- c(gene_sets[score_ok], null_sets)
for (i in seq_len(nrow(samples))) {
  message("统一评分：", samples$sample_id[i])
  x <- read_sample(i)                 # 分两遍读取，避免同时把两组大对象留在内存。
  ct <- x$counts[common_genes, , drop = FALSE]
  met <- x$spots
  met$nFeature_common <- Matrix::colSums(ct > 0)
  met$AUC_assessable <- met$QC_pass & met$nFeature_common >= aucMaxRank
  met$TLS22_null_tail <- NA_real_
  for (m in names(marker_sets)) {
    g <- gene_sets[[m]]
    met[[paste0(m, "_AUC")]] <- NA_real_
    met[[paste0(m, "_genes")]] <- if (score_ok[m]) Matrix::colSums(ct[g, , drop = FALSE] > 0) else NA_real_
    met[[paste0(m, "_UMI")]] <- if (score_ok[m]) Matrix::colSums(ct[g, , drop = FALSE]) else NA_real_
  }
  eligible <- which(met$AUC_assessable)
  chunks <- split(eligible, ceiling(seq_along(eligible)/auc_chunk))
  set.seed(seed)
  for (ids in chunks) {
    rankings <- AUCell::AUCell_buildRankings(ct[, ids, drop = FALSE], plotStats = FALSE,
                                             splitByBlocks = TRUE, keepZeroesAsNA = FALSE, BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
    # 要求检出基因数>=K，确保计分范围内没有零表达随机排名；另设全零面板断言。
    scores <- AUCell::getAUC(AUCell::AUCell_calcAUC(score_sets, rankings,
                                                    aucMaxRank = aucMaxRank, normAUC = TRUE, nCores = 1, verbose = FALSE))
    scores <- scores[, met$spot[ids], drop = FALSE]
    for (m in names(gene_sets)[score_ok]) {
      a <- as.numeric(scores[m, ])
      a[met[[paste0(m, "_genes")]][ids] == 0] <- 0
      met[[paste0(m, "_AUC")]][ids] <- a
    }
    actual <- met$TLS22_AUC[ids]
    met$TLS22_null_tail[ids] <- (1 + colSums(sweep(scores[names(null_sets), , drop = FALSE], 2, actual, ">="))) / (n_null+1)
    rm(rankings, scores)
  }
  stopifnot(all(met$TLS22_AUC[met$AUC_assessable & met$TLS22_genes == 0] == 0))
  met$expression_seed <- met$AUC_assessable & met$TLS22_genes >= min_TLS_genes &
    !is.na(met$TLS22_null_tail) & met$TLS22_null_tail <= null_tail_cutoff & met$TLS22_AUC > 0
  graphs[[i]] <- make_graph(met)
  met$expression_region <- find_regions(graphs[[i]], met$expression_seed)
  met$sample_id <- samples$sample_id[i]
  met$alias <- samples$alias[i]
  spots_all[[i]] <- met
  markers <- unique(paste0(gene_prefix, c(unlist(marker_sets, use.names = FALSE), "MKI67")))
  present <- intersect(markers, rownames(x$counts))
  markers_all[[i]] <- x$counts[present, , drop = FALSE]
  marker_audit[[i]] <- data.frame(sample_id = samples$sample_id[i], feature = markers,
                                  in_counts = markers %in% rownames(x$counts), total_UMI = NA_real_, n_spots_detected = NA_real_)
  j <- match(present, markers)
  marker_audit[[i]]$total_UMI[j] <- Matrix::rowSums(markers_all[[i]])
  marker_audit[[i]]$n_spots_detected[j] <- Matrix::rowSums(markers_all[[i]] > 0)
  rm(x, ct)
  invisible(gc())
}

##### 第3步：表达候选与H&E候选取并集，保留来源及待复核状态 #####
roi_list <- members_list <- polygons_list <- boundaries_list <- sensitivity <- list()
for (i in seq_len(nrow(samples))) {
  met <- spots_all[[i]]
  regions <- list()
  er <- sort(unique(na.omit(met$expression_region)))
  for (k in er) regions[[paste0("E", k)]] <- met$spot[which(met$expression_region == k)]
  hp <- HE_polygons[HE_polygons$sample_id == samples$sample_id[i], ]
  he_ids <- unique(hp$roi_id)
  for (k in he_ids) {
    p <- hp[hp$roi_id == k, ]
    inside <- sp::point.in.polygon(met$x, met$y, p$x, p$y) > 0
    regions[[paste0("H_", k)]] <- met$spot[inside]
  }
  # 只有共享捕获spot的候选才自动归并；不以凸包覆盖或距离补洞合并。
  node_names <- names(regions)
  if (length(regions)) {
    pairs <- if (length(regions) > 1) t(combn(seq_along(regions), 2)) else matrix(integer(), ncol = 2)
    linked <- if (nrow(pairs)) vapply(seq_len(nrow(pairs)), function(k) {
      length(intersect(regions[[pairs[k, 1]]], regions[[pairs[k, 2]]])) > 0
    }, logical(1)) else logical()
    edges <- data.frame(from = node_names[pairs[linked, 1]], to = node_names[pairs[linked, 2]])
    g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = node_names))
    membership <- igraph::components(g)$membership
    # 用全片阵列位置拟合原始像素坐标，供六边网格轮廓显示。
    if (length(er)) {
      xy_fit <- lm(cbind(x, y) ~ array_col + array_row, data = met)
      if (xy_fit$rank != 3) stop("阵列位置不足以定位表达ROI轮廓。")
    }
    for (k in sort(unique(membership))) {
      nodes <- names(membership)[membership == k]
      ss <- unique(unlist(regions[nodes], use.names = FALSE))
      rid <- paste0(samples$alias[i], "_R", sprintf("%02d", k))
      is_e <- any(startsWith(nodes, "E")); is_h <- any(startsWith(nodes, "H_"))
      source <- if (is_e && is_h) "both" else if (is_e) "expression" else "HE_only"
      ix <- match(ss, met$spot)
      row <- data.frame(sample_id = samples$sample_id[i], alias = samples$alias[i], roi_id = rid,
                        source = source, status = "pending_review", n_spots = length(ss), n_AUC_assessable = sum(met$AUC_assessable[ix]),
                        n_B_marker_spots = if (score_ok["Bcell"]) sum(met$Bcell_genes[ix] >= 2, na.rm = TRUE) else NA_integer_,
                        n_T_marker_spots = if (score_ok["Tcell"]) sum(met$Tcell_genes[ix] >= 2, na.rm = TRUE) else NA_integer_)
      for (m in names(marker_sets)) {
        v <- met[[paste0(m, "_AUC")]][ix]
        row[[paste0(m, "_AUC_median")]] <- if (!length(v) || all(is.na(v))) NA_real_ else median(v, na.rm = TRUE)
      }
      roi_list[[length(roi_list)+1L]] <- row
      if (length(ss)) members_list[[length(members_list)+1L]] <- data.frame(sample_id = samples$sample_id[i], roi_id = rid, spot = ss)
      ep <- met[which(met$expression_region %in% as.integer(sub("E", "", nodes[startsWith(nodes, "E")]))), ]
      if (nrow(ep)) {
        e <- grid_edges(ep)
        p1 <- cbind(1, e$c1, e$r1) %*% coef(xy_fit)
        p2 <- cbind(1, e$c2, e$r2) %*% coef(xy_fit)
        boundaries_list[[length(boundaries_list)+1L]] <- data.frame(sample_id = samples$sample_id[i], roi_id = rid,
                                                                    x = p1[, 1], y = p1[, 2], xend = p2[, 1], yend = p2[, 2])
      }
      for (h in sub("H_", "", nodes[startsWith(nodes, "H_")])) {
        p <- hp[hp$roi_id == h, ]; p$roi_id <- rid
        p$polygon_id <- paste(rid, h, sep = "_")
        polygons_list[[length(polygons_list)+1L]] <- p
      }
    }
  }
  # 固定范围敏感性检查，不据此挑选最漂亮的一组参数。
  for (tail in unique(c(0.005, null_tail_cutoff, 0.025))) for (size in unique(c(2L, min_ROI_spots, 4L))) {
    flag <- met$AUC_assessable & met$TLS22_genes >= min_TLS_genes & !is.na(met$TLS22_null_tail) &
      met$TLS22_null_tail <= tail & met$TLS22_AUC > 0
    ids <- find_regions(graphs[[i]], flag, size)
    sensitivity[[length(sensitivity)+1L]] <- data.frame(sample_id = samples$sample_id[i],
                                                        null_tail = tail, min_spots = size, n_expression_ROI = length(unique(na.omit(ids))))
  }
}
rois <- if (length(roi_list)) bind_rows(roi_list) else data.frame(sample_id = character(), alias = character(), roi_id = character(), source = character(), status = character(), n_spots = integer(), n_AUC_assessable = integer())
members <- if (length(members_list)) bind_rows(members_list) else data.frame(sample_id = character(), roi_id = character(), spot = character())
boundaries <- if (length(boundaries_list)) bind_rows(boundaries_list) else data.frame(sample_id = character(), roi_id = character(), x = numeric(), y = numeric(), xend = numeric(), yend = numeric())
polygons_plot <- if (length(polygons_list)) bind_rows(polygons_list) else data.frame(sample_id = character(), roi_id = character(), x = numeric(), y = numeric(), polygon_id = character())
spots <- bind_rows(spots_all)
sensitivity <- bind_rows(sensitivity)

##### 精简汇总：每组目前一个样本，只做描述性展示 #####
summary_table <- samples[, c("sample_id", "treatment", "day", "side", "mouse_id", "donor_id", "batch", "alias")]
for (i in seq_len(nrow(samples))) {
  m <- spots_all[[i]]
  summary_table$n_spots[i] <- nrow(m)
  summary_table$n_QC[i] <- sum(m$QC_pass)
  summary_table$n_AUC_assessable[i] <- sum(m$AUC_assessable)
  summary_table$AUC_assessable_fraction_QC[i] <- if (any(m$QC_pass)) sum(m$AUC_assessable)/sum(m$QC_pass) else NA_real_
  summary_table$host_UMI_median[i] <- if (any(m$QC_pass)) median(m$host_UMI[m$QC_pass]) else NA_real_
  summary_table$n_expression_seed[i] <- sum(m$expression_seed)
  summary_table$n_expression_ROI[i] <- if (any(m$AUC_assessable)) length(unique(na.omit(m$expression_region))) else NA_integer_
  summary_table$n_HE_ROI[i] <- length(unique(HE_polygons$roi_id[HE_polygons$sample_id == samples$sample_id[i]]))
  summary_table$n_union_candidate_ROI[i] <- if (any(m$AUC_assessable) || summary_table$n_HE_ROI[i] > 0) sum(rois$sample_id == samples$sample_id[i]) else NA_integer_
  summary_table$HE_review_done[i] <- samples$sample_id[i] %in% HE_reviewed
  for (module in names(marker_sets)) {
    v <- m[[paste0(module, "_AUC")]]
    summary_table[[paste0(module, "_AUC_median")]][i] <- if (all(is.na(v))) NA_real_ else median(v, na.rm = TRUE)
    summary_table[[paste0(module, "_AUC_mean")]][i] <- if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
    summary_table[[paste0(module, "_AUC_nonzero_fraction")]][i] <- if (all(is.na(v))) NA_real_ else mean(v > 0, na.rm = TRUE)
    summary_table[[paste0(module, "_genes_used")]][i] <- length(gene_sets[[module]])
  }
}
print(summary_table)

##### 图像函数：匿名H&E预览、统一色阶评分图、候选轮廓 #####
he_plot <- function(i, candidate = FALSE) {
  im <- images[[i]]
  he <- png::readPNG(im$file)
  step <- max(1L, ceiling(max(dim(he)[1:2])/image_max_px))
  rr <- seq(1, dim(he)[1], by = step); cc <- seq(1, dim(he)[2], by = step)
  he <- if (length(dim(he)) == 2) he[rr, cc, drop = FALSE] else he[rr, cc, , drop = FALSE]
  p <- ggplot() + annotation_raster(as.raster(he), xmin = 0, xmax = im$width/im$scale,
                                    ymin = -im$height/im$scale, ymax = 0) + coord_fixed(expand = FALSE,
                                                                                        xlim = c(0, im$width/im$scale), ylim = c(-im$height/im$scale, 0)) +
    theme_void() + labs(title = paste(samples$alias[i], if (candidate) "candidate ROI: pending review" else "H&E only: review before expression"))
  if (candidate) {
    b <- boundaries[boundaries$sample_id == samples$sample_id[i], ]
    h <- polygons_plot[polygons_plot$sample_id == samples$sample_id[i], ]
    p <- p + geom_segment(data = b, aes(x = x, y = -y, xend = xend, yend = -yend), color = "#D55E00", linewidth = 0.5) +
      geom_polygon(data = h, aes(x, -y, group = polygon_id), fill = NA, color = "#0072B2", linewidth = 0.6)
    rr <- rois[rois$sample_id == samples$sample_id[i], ]
    labs <- bind_rows(lapply(rr$roi_id, function(r) {
      ss <- members$spot[members$sample_id == samples$sample_id[i] & members$roi_id == r]
      mm <- spots_all[[i]][match(ss, spots_all[[i]]$spot), ]
      if (nrow(mm)) data.frame(roi_id = r, x = mean(mm$x), y = mean(mm$y)) else {
        hp <- h[h$roi_id == r, ]; data.frame(roi_id = r, x = mean(hp$x), y = mean(hp$y))
      }
    }))
    if (nrow(labs)) p <- p + geom_label(data = labs, aes(x, -y, label = roi_id),
                                        size = 3, fill = "white", color = "#111111", linewidth = 0.2, alpha = 0.9)
    p <- p + labs(caption = "Orange: expression grid outline; blue: H&E polygon. Neither is a confirmed TLS boundary.")
  }
  p
}
score_map <- function(i, column, title, upper) {
  m <- spots_all[[i]]
  p <- ggplot(m, aes(x, -y, color = .data[[column]])) + geom_point(size = 1) +
    scale_color_gradient(low = "#EEEEEE", high = "#B2182B", na.value = "#777777", limits = c(0, upper)) +
    coord_fixed() + theme_void() + labs(title = paste(samples$alias[i], title), color = NULL)
  p
}

##### 输出：三份文件，运行信息及小型中间结果集中放入RDS #####
dir.create(path.expand(outdir), recursive = TRUE, showWarnings = FALSE)
write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE)
pdf(file.path(outdir, "plots.pdf"), width = 11, height = 8, compress = TRUE)
tryCatch({
  for (i in order(samples$alias)) print(he_plot(i))
  print(ggplot(spots, aes(nFeature_common)) + geom_histogram(bins = 40, fill = "#999999", color = "white") +
          geom_vline(xintercept = aucMaxRank, color = "#D55E00", linetype = 2) + facet_wrap(~alias, ncol = 2) +
          theme_classic() + labs(title = "Common-gene detection and the same AUCell rank cutoff", x = "Detected common genes", y = "Spots"))
  columns <- c("host_UMI", "nFeature_common", paste0(names(marker_sets), "_AUC"))
  for (column in columns) {
    upper <- if (all(is.na(spots[[column]]))) 1 else max(spots[[column]], na.rm = TRUE)
    if (!is.finite(upper) || upper <= 0) upper <- 1
    print(wrap_plots(lapply(order(samples$alias), function(i) score_map(i, column, column, upper)), ncol = 2))
  }
  for (i in order(samples$alias)) print(he_plot(i, candidate = TRUE))
  # TLS可能仅占小区域，全片中位数常为0；同时展示均值和候选数，不只比较中位数。
  plot_data <- summary_table
  plot_data$treatment <- factor(plot_data$treatment, levels = c("Vehicle", "VG161"))
  p1 <- ggplot(plot_data, aes(treatment, TLS22_AUC_mean, color = treatment)) + geom_point(size = 3) +
    theme_classic() + labs(title = "TLS22 mean AUC", subtitle = "Evaluable spots only", y = "Mean AUCell AUC", x = NULL)
  p2 <- ggplot(plot_data, aes(treatment, Imprint29_AUC_mean, color = treatment)) + geom_point(size = 3) +
    theme_classic() + labs(title = "Imprint29 mean AUC", subtitle = "See measured gene coverage", y = "Mean AUCell AUC", x = NULL)
  p3 <- ggplot(plot_data, aes(treatment, n_expression_ROI, color = treatment)) + geom_point(size = 3) +
    theme_classic() + labs(title = "Expression candidate ROI", subtitle = "Pending review; not confirmed TLS", y = "Candidates", x = NULL)
  p4 <- ggplot(plot_data, aes(treatment, AUC_assessable_fraction_QC, color = treatment)) + geom_point(size = 3) +
    scale_y_continuous(limits = c(0, 1)) + theme_classic() +
    labs(title = "AUCell assessable fraction", subtitle = "Among spots passing basic QC", y = "Fraction", x = NULL)
  print(wrap_plots(p1, p2, p3, p4, ncol = 2) + plot_annotation(title = "D3 injected: one slice per group, descriptive only"))
}, finally = dev.off())
results <- list(samples = samples, summary = summary_table, spots = spots, rois = rois,
                members = members, expression_boundaries = boundaries, HE_polygons = HE_polygons,
                HE_polygons_merged = polygons_plot, marker_counts = setNames(markers_all, samples$sample_id),
                marker_audit = bind_rows(marker_audit), coverage = coverage, gene_sets_used = gene_sets,
                common_genes = common_genes, null_gene_sets = null_sets, background_matching = background_matching,
                precheck = precheck, images = images,
                sensitivity = sensitivity,
                info = list(run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), sources = sources,
                            seed = seed, aucMaxRank = aucMaxRank, auc_fraction = auc_fraction, auc_chunk = auc_chunk,
                            min_host_UMI = min_host_UMI, min_features = min_features, max_percent_mt = max_percent_mt,
                            n_null = n_null, null_match_pool = null_match_pool, null_tail_cutoff = null_tail_cutoff,
                            min_TLS_genes = min_TLS_genes, min_ROI_spots = min_ROI_spots, gene_prefix = gene_prefix,
                            viral_pattern = viral_pattern, mt_pattern = mt_pattern, imprint_min_coverage = imprint_min_coverage,
                            HE_roi_file = HE_roi_file, draw_HE_ROI = draw_HE_ROI, image_max_px = image_max_px,
                            interpretation = c("TLS22 AUCell is an adapted screening score, not the Science iSTAR model or TLS probability.",
                                               "Null tails compare matched random gene sets, not treatment P values or validated TLS diagnostic thresholds.",
                                               "Null tails are not multiplicity-adjusted and do not establish ROI-level FDR.",
                                               "Screening can include single-lineage lymphoid aggregates; B/T structure is assessed during ROI review.",
                                               "H&E panels are anonymous and precede expression display, but the script is not a formally blinded assessment.",
                                               "All ROIs are pending review; H&E absence of annotation does not mean absence of TLS.",
                                               "Expression grid outlines preserve holes but are not histological TLS area; no physical area is inferred.",
                                               "With one sample per group, comparison is descriptive; no biological replicate P value is calculated.",
                                               "Imprint29 PTLP is kept literally; missing genes are reported, not silently replaced by PLTP."),
                            sessionInfo = capture.output(sessionInfo())))
saveRDS(results, file.path(outdir, "results.rds"), compress = "gzip")
message("完成。输出目录：", normalizePath(outdir, winslash = "/", mustWork = FALSE))
message("候选ROI均待复核；results.rds中保存坐标、基因覆盖、阈值敏感性和运行信息。")
