##### 清理运行环境 #####
rm(list = ls(all.names = TRUE))
graphics.off()
invisible(gc())

##### D3注射侧：样本整理、统一TLS评分、待复核ROI #####
# 前半做步骤1—3，后半续接步骤4—6；候选ROI不等于确诊TLS，不自动判断成熟度。
# 前半输出summary.csv、plots.pdf、results.rds；后半只新增review.csv和review.pdf。
# 不保存Seurat对象或全基因表达矩阵。
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
##### 图版样式 #####
tls_group_colors <- c(Vehicle = "#426B91", VG161 = "#C35B3F")
tls_expression_colors <- c("#EDF2F6", "#91AEC7", "#24537D")
tls_unscored <- "#B8BDC3"
tls_strip_theme <- theme(strip.background = element_rect(fill = "#EDF2F6", color = NA),
                         strip.text = element_text(face = "bold", size = 9, color = "#22313F", hjust = 0.5,
                                                   margin = margin(3, 4, 3, 4)))
tls_map_theme <- theme_void(base_size = 10, base_family = "Helvetica") + tls_strip_theme +
  theme(plot.title = element_text(face = "bold", size = 10, color = "#22313F"),
        legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        plot.margin = margin(6, 8, 6, 8))
tls_chart_theme <- theme_classic(base_size = 10, base_family = "Helvetica") + tls_strip_theme +
  theme(plot.title = element_text(face = "bold", size = 10, color = "#22313F"),
        panel.border = element_rect(fill = NA, color = "#526170", linewidth = 0.35),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.35, color = "#526170"),
        strip.background = element_rect(fill = "#EDF2F6", color = "#526170", linewidth = 0.35),
        plot.margin = margin(8, 12, 8, 12))

he_plot <- function(i, candidate = FALSE) {
  im <- images[[i]]
  he <- png::readPNG(im$file)
  step <- max(1L, ceiling(max(dim(he)[1:2])/image_max_px))
  rr <- seq(1, dim(he)[1], by = step); cc <- seq(1, dim(he)[2], by = step)
  he <- if (length(dim(he)) == 2) he[rr, cc, drop = FALSE] else he[rr, cc, , drop = FALSE]
  p <- ggplot() + annotation_raster(as.raster(he), xmin = 0, xmax = im$width/im$scale,
                                    ymin = -im$height/im$scale, ymax = 0) + coord_fixed(expand = FALSE,
                                                                                        xlim = c(0, im$width/im$scale), ylim = c(-im$height/im$scale, 0)) +
    tls_map_theme + labs(title = paste(samples$alias[i], if (candidate) "ROI candidates" else "H&E"))
  if (candidate) {
    b <- boundaries[boundaries$sample_id == samples$sample_id[i], ]
    h <- polygons_plot[polygons_plot$sample_id == samples$sample_id[i], ]
    p <- p + geom_segment(data = b, aes(x = x, y = -y, xend = xend, yend = -yend,
                                        color = "Expression candidate"), linewidth = 0.45) +
      geom_polygon(data = h, aes(x, -y, group = polygon_id, color = "H&E candidate"),
                   fill = NA, linewidth = 0.5, key_glyph = "path") +
      scale_color_manual(name = NULL, values = c(`Expression candidate` = "#BB521E", `H&E candidate` = "#24537D"))
    rr <- rois[rois$sample_id == samples$sample_id[i], ]
    labs <- bind_rows(lapply(rr$roi_id, function(r) {
      ss <- members$spot[members$sample_id == samples$sample_id[i] & members$roi_id == r]
      mm <- spots_all[[i]][match(ss, spots_all[[i]]$spot), ]
      if (nrow(mm)) data.frame(roi_id = r, x = mean(mm$x), y = mean(mm$y)) else {
        hp <- h[h$roi_id == r, ]; data.frame(roi_id = r, x = mean(hp$x), y = mean(hp$y))
      }
    }))
    if (nrow(labs)) p <- p + geom_label(data = labs, aes(x, -y, label = roi_id),
                                        size = 3, family = "Helvetica", fill = "white", color = "#22313F", linewidth = 0.15)
  }
  p
}
score_map <- function(i, column, title, upper, unit, show_unscored) {
  m <- spots_all[[i]]
  label <- samples$treatment[i]
  if (sum(samples$treatment == label) > 1) label <- paste(label, samples$alias[i])
  unscored_guide <- if (show_unscored) {
    guide_legend(order = 2, override.aes = list(colour = tls_unscored, size = 2))
  } else "none"
  p <- ggplot(m, aes(x, -y, color = .data[[column]])) +
    geom_point(aes(shape = is.na(.data[[column]])), size = 1.2, show.legend = c(shape = TRUE)) +
    scale_color_gradientn(colours = tls_expression_colors, na.value = tls_unscored,
                          limits = c(0, upper), name = unit, guide = guide_colorbar(order = 1)) +
    scale_shape_manual(name = NULL, values = c(`FALSE` = 16, `TRUE` = 16),
                       limits = c(FALSE, TRUE), breaks = TRUE, labels = "Unscored", guide = unscored_guide) +
    coord_fixed() + facet_wrap(vars(!!paste(label, title))) + tls_map_theme
  p
}

##### 输出：三份文件，运行信息及小型中间结果集中放入RDS #####
dir.create(path.expand(outdir), recursive = TRUE, showWarnings = FALSE)
write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE)
pdf(file.path(outdir, "plots.pdf"), width = 9, height = 6, family = "Helvetica",
    useDingbats = FALSE, compress = TRUE)
tryCatch({
  for (i in order(samples$alias)) print(he_plot(i))
  print(ggplot(spots, aes(nFeature_common)) + geom_histogram(bins = 40, fill = "#91AEC7", color = "white", linewidth = 0.2) +
          geom_vline(xintercept = aucMaxRank, color = "#BB521E", linetype = 2, linewidth = 0.5) +
          facet_wrap(~alias, ncol = 2) + tls_chart_theme +
          labs(title = "AUCell rank cutoff", x = "Detected common genes", y = "Spots"))
  map_titles <- c(host_UMI = "Host UMI", nFeature_common = "Detected genes", TLS22_AUC = "TLS22",
                  Imprint29_AUC = "Imprint29", Bcell_AUC = "B-cell markers", Tcell_AUC = "T-cell markers",
                  Organization_AUC = "Organization", GC_related_AUC = "GC markers")
  spatial_panels <- list(c("host_UMI", "nFeature_common"), c("TLS22_AUC", "Imprint29_AUC"),
                         c("Bcell_AUC", "Tcell_AUC"), c("Organization_AUC", "GC_related_AUC"))
  map_order <- order(match(samples$treatment, names(tls_group_colors)), samples$alias)
  for (pair in spatial_panels) {
    rows <- lapply(pair, function(column) {
      upper <- if (all(is.na(spots[[column]]))) 1 else max(spots[[column]], na.rm = TRUE)
      if (!is.finite(upper) || upper <= 0) upper <- 1
      unit <- if (endsWith(column, "_AUC")) "AUC" else if (column == "host_UMI") "UMI" else "Genes"
      show_unscored <- anyNA(spots[[column]])
      maps <- lapply(map_order, function(i) score_map(i, column, map_titles[[column]], upper, unit, show_unscored))
      for (j in seq_len(length(maps) - 1L)) maps[[j]] <- maps[[j]] + theme(legend.position = "none")
      wrap_plots(maps, ncol = 2)
    })
    print(wrap_plots(rows, ncol = 1))
  }
  for (i in order(samples$alias)) print(he_plot(i, candidate = TRUE))
  # TLS可能仅占小区域，全片中位数常为0；同时展示均值和候选数，不只比较中位数。
  plot_data <- summary_table
  plot_data$treatment <- factor(plot_data$treatment, levels = c("Vehicle", "VG161"))
  p1 <- ggplot(plot_data, aes(treatment, TLS22_AUC_mean, color = treatment)) + geom_point(size = 3) +
    tls_chart_theme + labs(y = "Mean AUC", x = NULL)
  p2 <- ggplot(plot_data, aes(treatment, Imprint29_AUC_mean, color = treatment)) + geom_point(size = 3) +
    tls_chart_theme + labs(y = "Mean AUC", x = NULL)
  p3 <- ggplot(plot_data, aes(treatment, n_expression_ROI, color = treatment)) + geom_point(size = 3) +
    tls_chart_theme + labs(y = "ROIs", x = NULL)
  p4 <- ggplot(plot_data, aes(treatment, AUC_assessable_fraction_QC, color = treatment)) + geom_point(size = 3) +
    scale_y_continuous(limits = c(0, 1)) + tls_chart_theme +
    labs(y = "Fraction of QC spots", x = NULL)
  summary_strips <- c(TLS22 = "#DCE9F3", Imprint29 = "#DDEEEA",
                      `Expression candidates` = "#F5E9DC", `AUCell coverage` = "#EDF2F6")
  summary_plots <- Map(function(p, label, fill) {
    p + facet_wrap(vars(!!label)) + scale_color_manual(values = tls_group_colors) +
      theme(legend.position = "none",
            strip.background = element_rect(fill = fill, color = "#526170", linewidth = 0.35))
  }, list(p1, p2, p3, p4), names(summary_strips), unname(summary_strips))
  print(wrap_plots(summary_plots, ncol = 2))
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

##### 第4—6步：续接现有结果，不再次清理运行环境 #####
# 前面的清理只在整个脚本最开始执行。后续可以从这里单独运行，不必重算AUCell。
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)

##### 后续输入与参数：默认沿用上一步输出目录 #####
review_dir <- if (exists("outdir", inherits = FALSE)) outdir else "~/OV/SCT/TLS_D3"
review_input <- file.path(review_dir, "results.rds")
review_file <- file.path(review_dir, "review.csv")
context_rings <- 2L                  # 原始六邻域外扩两圈，不跨组织孔洞。
depth_bins <- 3L                     # UMI与检出基因数各分三档，作技术敏感性检查。
min_stratum_spots <- 20L             # 两片均达到此数量的分层才用于共同权重描述。
review_image_max_px <- 1600L
review_image_files <- character()    # 图像搬家时填写命名向量，如c(VG161_3D_1="/path/image.png")。

tls <- readRDS(path.expand(review_input))
needed <- c("samples", "summary", "spots", "rois", "members", "marker_counts", "info")
if (!all(needed %in% names(tls))) stop("review_input不是本流程的results.rds。")
if (context_rings < 0 || depth_bins < 1 || min_stratum_spots < 1) stop("后续参数无效。")
tls_spots <- tls$spots
tls_rois <- tls$rois
prefix <- tls$info$gene_prefix
if (is.null(prefix)) prefix <- ""
review_panels <- list(
  Bcell = c("MS4A1", "CD79A", "CD79B", "CD19", "CD22"),
  Tcell = c("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"),
  Organization = c("CXCL13", "CCL19", "CCL21"),
  FDC_related = c("CR2", "FCER2"),
  GC_B_related = c("AICDA", "RGS13", "SERPINA9"),
  GC_reference = c("BCL6", "MEF2B", "MKI67"))
review_genes <- unique(c(unlist(review_panels, use.names = FALSE), substring(tls$gene_sets_used$TLS22, nchar(prefix)+1)))
for (s in tls$samples$sample_id) {
  m <- tls$marker_counts[[s]]
  ss <- tls_spots$spot[tls_spots$sample_id == s]
  if (is.null(m) || anyDuplicated(colnames(m)) || !setequal(colnames(m), ss)) stop(s, "：marker counts与spot表不一致。")
}
if (anyDuplicated(paste(tls_rois$sample_id, tls_rois$roi_id))) stop("ROI标识重复。")
if (any(!paste(tls$members$sample_id,tls$members$roi_id) %in% paste(tls_rois$sample_id,tls_rois$roi_id))) stop("ROI成员表含未知ROI。")

##### 原始六邻域和ROI邻域：这里独立定义，便于只运行后半段 #####
review_graph <- function(m) {
  key <- paste(m$array_row, m$array_col, sep = ":")
  if (anyDuplicated(key)) stop("阵列坐标重复。")
  off <- rbind(c(0, 2), c(1, 1), c(1, -1))
  e <- bind_rows(lapply(seq_len(nrow(off)), function(k) {
    j <- match(paste(m$array_row+off[k, 1], m$array_col+off[k, 2], sep = ":"), key)
    data.frame(from = m$spot[!is.na(j)], to = m$spot[j[!is.na(j)]])
  }))
  igraph::graph_from_data_frame(e, directed = FALSE, vertices = data.frame(name = m$spot))
}
roi_context <- graphs_review <- vector("list", nrow(tls_rois))
gene_evidence <- roi_evidence <- composition <- list()
all_candidate_keys <- paste(tls$members$sample_id, tls$members$spot, sep = "|")
for (i in seq_len(nrow(tls_rois))) {
  s <- tls_rois$sample_id[i]; rid <- tls_rois$roi_id[i]
  m <- tls_spots[tls_spots$sample_id == s, ]
  core <- unique(tls$members$spot[tls$members$sample_id == s & tls$members$roi_id == rid])
  if (!all(core %in% m$spot)) stop(rid, "：成员条码不存在。")
  g <- review_graph(m[m$QC_pass | m$spot %in% core, ])
  near <- if (length(core)) unique(unlist(igraph::ego(g, order = context_rings, nodes = core), use.names = FALSE)) else integer()
  near_spots <- igraph::V(g)$name[near]
  ring <- setdiff(near_spots, core)
  # 背景只用QC合格、且不属于任何候选的邻域spot。
  ring <- ring[m$QC_pass[match(ring, m$spot)] & !paste(s, ring, sep = "|") %in% all_candidate_keys]
  roi_context[[i]] <- list(sample_id = s, roi_id = rid, core = core, context = union(core, near_spots), ring = ring)
  ct <- tls$marker_counts[[s]]
  for (gene in review_genes) {
    feature <- paste0(prefix, gene)
    for (zone in c("core", "ring")) {
      ss <- if (zone == "core") core else ring
      assessable <- feature %in% rownames(ct) && length(ss) > 0
      raw <- if (assessable) as.numeric(ct[feature, ss]) else numeric()
      qc <- if (length(ss)) m$QC_pass[match(ss, m$spot)] else logical()
      norm <- if (assessable && any(qc)) log1p(raw[qc]/m$host_UMI[match(ss[qc], m$spot)]*1e4) else numeric()
      gene_evidence[[length(gene_evidence)+1L]] <- data.frame(sample_id = s, roi_id = rid,
                                                              zone = zone, gene = gene, in_matrix = feature %in% rownames(ct), n_spots = length(ss),
                                                              UMI = if (assessable) sum(raw) else NA_real_,
                                                              n_detected_spots = if (assessable) sum(raw > 0) else NA_integer_,
                                                              mean_logCP10k_QC = if (length(norm)) mean(norm) else NA_real_)
    }
  }
  row <- data.frame(sample_id = s, roi_id = rid, n_core = length(core), n_ring = length(ring),
                    n_core_QC = sum(m$QC_pass[match(core,m$spot)]), n_core_AUC = sum(m$AUC_assessable[match(core,m$spot)]),
                    core_host_UMI_median = if (length(core)) median(m$host_UMI[match(core, m$spot)]) else NA_real_,
                    ring_host_UMI_median = if (length(ring)) median(m$host_UMI[match(ring, m$spot)]) else NA_real_)
  for (panel in names(review_panels)) {
    features <- paste0(prefix, review_panels[[panel]])
    available <- all(features %in% rownames(ct)) && length(core) > 0
    row[[paste0(panel, "_UMI")]] <- if (available) sum(ct[features, core, drop = FALSE]) else NA_real_
    row[[paste0(panel, "_genes_detected")]] <- if (available) sum(Matrix::rowSums(ct[features, core, drop = FALSE]) > 0) else NA_integer_
  }
  row$RNA_B_support <- if (is.na(row$Bcell_UMI) || row$n_core_QC == 0) "not_assessable" else if (row$Bcell_UMI == 0) "not_detected" else "detected_requires_review"
  row$RNA_GC_B_support <- if (is.na(row$GC_B_related_UMI) || row$n_core_QC == 0) "not_assessable" else if (row$GC_B_related_UMI == 0) "not_detected" else "detected_requires_review"
  roi_evidence[[i]] <- row
  for (panel in intersect(c("TLS22", "Imprint29"), names(tls$gene_sets_used))) {
    gg <- intersect(tls$gene_sets_used[[panel]], rownames(ct))
    if (length(core) && length(gg)) {
      uu <- Matrix::rowSums(ct[gg, core, drop = FALSE])
      composition[[length(composition)+1L]] <- data.frame(sample_id = s, roi_id = rid, panel = panel,
                                                          gene = substring(names(uu), nchar(prefix)+1), raw_UMI = as.numeric(uu),
                                                          raw_UMI_fraction = if (sum(uu) > 0) as.numeric(uu)/sum(uu) else NA_real_)
    }
  }
}
gene_evidence <- if (length(gene_evidence)) bind_rows(gene_evidence) else data.frame(sample_id=character(), roi_id=character(), zone=character(), gene=character(), UMI=numeric())
roi_evidence <- if (length(roi_evidence)) bind_rows(roi_evidence) else data.frame(sample_id=character(), roi_id=character(), n_core=integer(), n_ring=integer(), RNA_B_support=character(), RNA_GC_B_support=character())
composition <- if (length(composition)) bind_rows(composition) else data.frame(sample_id=character(), roi_id=character(), panel=character(), gene=character(), raw_UMI=numeric(), raw_UMI_fraction=numeric())

##### 第4步：生成一张审核表，保留样本整片审核行和逐ROI审核行 #####
# review.csv中只填写人工列；自动证据列不要当作已经确认的TLS结论。
# decision：pending / uncertain / rejected / TLS_supported。
# maturity：not_assessed / undetermined / early_like / follicle_like / GC_like。
sample_keys <- setNames(vapply(tls$samples$sample_id, function(s) {
  j <- match(s, tls$samples$sample_id)
  digest::digest(list(s, tls_spots[tls_spots$sample_id == s, c("spot", "x", "y", "host_UMI")],
                      tls$images[[j]][c("position_md5", "image_md5")], tls$rois[tls$rois$sample_id == s, ], tls$marker_counts[[s]],
                      tls$info[c("aucMaxRank", "null_tail_cutoff", "min_TLS_genes", "min_ROI_spots")] ), algo = "xxhash64")
}, ""), tls$samples$sample_id)
review <- data.frame(record_type = "sample", sample_id = tls$samples$sample_id, roi_id = "",
                     review_key = unname(sample_keys), whole_section_HE_reviewed = FALSE,
                     HE_reviewed = NA, decision = "not_applicable", morphology = "not_applicable", BT_support = NA,
                     maturity = "not_applicable", GC_structure_support = NA, tissue_area_mm2 = NA_real_,
                     n_TLS_structures = NA_real_, TLS_area_mm2 = NA_real_, structures_nonoverlapping = FALSE,
                     areas_nonoverlapping = FALSE, notes = "", stringsAsFactors = FALSE)
if (nrow(tls_rois)) {
  rr <- data.frame(record_type = "roi", sample_id = tls_rois$sample_id, roi_id = tls_rois$roi_id,
                   review_key = vapply(seq_len(nrow(tls_rois)), function(i) {
                     ctx <- roi_context[[i]]
                     digest::digest(list(sample_keys[ctx$sample_id], sort(ctx$core),
                                         tls$HE_polygons_merged[tls$HE_polygons_merged$roi_id == ctx$roi_id, ]), algo = "xxhash64")
                   }, ""), whole_section_HE_reviewed = NA, HE_reviewed = FALSE, decision = "pending",
                   morphology = "not_assessed", BT_support = NA, maturity = "not_assessed", GC_structure_support = NA,
                   tissue_area_mm2 = NA_real_, n_TLS_structures = NA_real_, TLS_area_mm2 = NA_real_,
                   structures_nonoverlapping = NA, areas_nonoverlapping = NA, notes = "", stringsAsFactors = FALSE)
  review <- bind_rows(review, rr)
}
review$n_core_spots <- NA_integer_; review$RNA_B_support <- NA_character_; review$RNA_GC_B_support <- NA_character_
if (nrow(tls_rois)) {
  ii <- which(review$record_type == "roi")
  review$n_core_spots[ii] <- roi_evidence$n_core
  review$RNA_B_support[ii] <- roi_evidence$RNA_B_support
  review$RNA_GC_B_support[ii] <- roi_evidence$RNA_GC_B_support
}
editable <- c("whole_section_HE_reviewed", "HE_reviewed", "decision", "morphology", "BT_support",
              "maturity", "GC_structure_support", "tissue_area_mm2", "n_TLS_structures", "TLS_area_mm2",
              "structures_nonoverlapping", "areas_nonoverlapping", "notes")
if (file.exists(path.expand(review_file))) {
  previous_review <- read.csv(path.expand(review_file), stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("record_type", "sample_id", "roi_id", "review_key", editable) %in% names(previous_review))) stop("review.csv列名不完整。")
  previous_review$roi_id[is.na(previous_review$roi_id) & previous_review$record_type == "sample"] <- ""
  old_key <- paste(previous_review$record_type, previous_review$sample_id, previous_review$roi_id, sep = "|")
  new_key <- paste(review$record_type, review$sample_id, review$roi_id, sep = "|")
  if (anyDuplicated(old_key) || !setequal(old_key, new_key)) stop("候选列表已变化或审核行重复；请保留旧表，并把review_file改为新的文件名后再运行。")
  mm <- match(new_key, old_key)
  if (any(review$review_key != previous_review$review_key[mm])) stop("旧审核表与当前ROI/坐标/筛查版本不一致，未套用旧结论。请保留旧表并更换review_file。")
  review[, editable] <- previous_review[mm, editable]
}
review$notes[is.na(review$notes)] <- ""
as_review_bool <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z[z %in% c("", "NA")] <- NA_character_
  if (any(!is.na(z) & !z %in% c("TRUE", "FALSE", "1", "0"))) stop("审核布尔列请填写TRUE/FALSE；未知留NA。")
  ifelse(is.na(z), NA, z %in% c("TRUE", "1"))
}
for (cc in c("whole_section_HE_reviewed", "HE_reviewed", "BT_support", "GC_structure_support", "structures_nonoverlapping", "areas_nonoverlapping")) review[[cc]] <- as_review_bool(review[[cc]])
is_roi <- review$record_type == "roi"
if (anyNA(review$decision[is_roi]) || any(!review$decision[is_roi] %in% c("pending", "uncertain", "rejected", "TLS_supported"))) stop("ROI decision值无效。")
if (anyNA(review$maturity[is_roi]) || any(!review$maturity[is_roi] %in% c("not_assessed", "undetermined", "early_like", "follicle_like", "GC_like"))) stop("maturity值无效。")
if (anyNA(review$morphology[is_roi]) || any(!review$morphology[is_roi] %in% c("not_assessed","uncertain","no_aggregate","aggregate","follicle_like","GC_like"))) stop("morphology值无效。")
accepted <- is_roi & review$decision == "TLS_supported"
if (any(accepted & (is.na(review$HE_reviewed) | !review$HE_reviewed | is.na(review$BT_support) | !review$BT_support))) stop("TLS_supported需要HE_reviewed和BT_support均为TRUE；未知证据不能自动补齐。")
if (any(accepted & !review$morphology %in% c("aggregate","follicle_like","GC_like"))) stop("TLS_supported还需要明确的聚集/滤泡/GC样形态支持。")
if (any(accepted & review$RNA_B_support %in% c("not_detected","not_assessable") & !nzchar(trimws(review$notes)), na.rm = TRUE)) stop("该ROI缺少可用B标记支持；若人工支持TLS，请在notes说明额外的B/T与形态依据。")
if (any(is_roi & !accepted & review$maturity %in% c("early_like", "follicle_like", "GC_like"))) stop("未支持TLS的ROI不能先赋成熟相关类别。")
if (any(accepted & review$maturity == "GC_like" & (is.na(review$GC_structure_support) | !review$GC_structure_support))) stop("GC_like需要明确的GC结构支持，不能仅凭BCL6/MKI67评分。")
if (any(accepted & review$maturity == "GC_like" & review$morphology != "GC_like")) stop("GC_like成熟相关类别需要一致的GC样形态记录。")
for (cc in c("tissue_area_mm2", "TLS_area_mm2")) {
  raw <- review[[cc]]
  num <- suppressWarnings(as.numeric(raw))
  if (any(!is.na(raw) & nzchar(as.character(raw)) & is.na(num)) || any(!is.na(num) & (!is.finite(num) | num <= 0))) stop(cc, "须为正数或NA，不自动换算网格面积。")
  review[[cc]] <- num
}
nn <- suppressWarnings(as.numeric(review$n_TLS_structures))
if (any(!is.na(review$n_TLS_structures) & is.na(nn)) || any(!is.na(nn) & (!is.finite(nn) | nn < 0 | nn != round(nn)))) stop("n_TLS_structures须为非负整数或NA。")
if (any(accepted & !is.na(nn) & nn < 1)) stop("支持TLS的ROI至少需包含1个独立结构；未数清时保留NA。")
if (any(is_roi & !accepted & !is.na(nn) & nn > 0)) stop("未支持TLS的ROI不能预先填写阳性的TLS结构数量。")
review$n_TLS_structures <- nn
if (any(accepted & review$maturity %in% c("early_like","follicle_like","GC_like") & (is.na(nn) | nn != 1))) stop("成熟相关类别需要逐结构ROI；多结构ROI请先拆分，未明确数量时保留undetermined。")

##### 第5步：样本级TLS结构指标，未审完保持NA，不能把待审当成零 #####
structure_summary <- list()
for (s in tls$samples$sample_id) {
  sr <- review[review$record_type == "sample" & review$sample_id == s, ]
  rr <- review[is_roi & review$sample_id == s, ]
  ss <- tls_spots[tls_spots$sample_id == s, ]
  ids <- rr$roi_id[rr$decision == "TLS_supported"]
  whole <- isTRUE(sr$whole_section_HE_reviewed)
  ROI_complete <- whole && all(rr$decision %in% c("rejected", "TLS_supported")) && all(!is.na(rr$HE_reviewed) & rr$HE_reviewed)
  nn <- rr$n_TLS_structures[rr$decision == "TLS_supported"]
  # 一个表达/形态候选ROI可以包含多个TLS，不能默认一ROI=一结构。
  distinct <- length(ids) <= 1 || isTRUE(sr$structures_nonoverlapping)
  complete <- ROI_complete && all(!is.na(nn)) && distinct
  count <- if (complete) sum(nn) else NA_integer_
  captured <- unique(tls$members$spot[tls$members$sample_id == s & tls$members$roi_id %in% ids])
  RNA_regions <- if (length(ids)) vapply(ids, function(id) any(tls$members$sample_id == s & tls$members$roi_id == id & tls$members$spot %in% ss$spot[ss$QC_pass]), logical(1)) else logical()
  coverage_QC <- if (ROI_complete && any(ss$QC_pass) && (!length(ids) || all(RNA_regions))) {
    sum(ss$spot %in% captured & ss$QC_pass)/sum(ss$QC_pass)
  } else NA_real_
  density <- if (complete && !is.na(sr$tissue_area_mm2)) count/sr$tissue_area_mm2 else NA_real_
  area_fraction <- NA_real_
  if (complete && !is.na(sr$tissue_area_mm2)) {
    aa <- rr$TLS_area_mm2[rr$decision == "TLS_supported"]
    if (!length(ids)) area_fraction <- 0 else if ((length(ids) <= 1 || isTRUE(sr$areas_nonoverlapping)) && all(!is.na(aa))) {
      if (sum(aa) > sr$tissue_area_mm2*(1+1e-6)) stop(s, "：支持ROI面积总和超过可评估组织面积。")
      area_fraction <- sum(aa)/sr$tissue_area_mm2
    }
  }
  stage <- rr$maturity[rr$decision == "TLS_supported"]
  stages_complete <- complete && length(ids) > 0 && all(nn == 1) && all(stage %in% c("early_like", "follicle_like", "GC_like"))
  n_reviewed <- sum(rr$decision != "pending" & !is.na(rr$HE_reviewed) & rr$HE_reviewed)
  n_stage_reviewed <- sum(stage %in% c("early_like","follicle_like","GC_like"))
  structure_summary[[length(structure_summary)+1L]] <- data.frame(sample_id = s,
                                                                  whole_section_HE_reviewed = whole, ROI_review_complete = ROI_complete, TLS_review_complete = complete,
                                                                  n_ROI_reviewed = n_reviewed,
                                                                  n_TLS_supported_ROI_reviewed_subset = if (n_reviewed) length(ids) else NA_integer_,
                                                                  n_TLS_structures_reviewed_subset = if (n_reviewed && all(!is.na(nn)) && distinct) sum(nn) else NA_integer_,
                                                                  n_TLS_supported_ROI_whole_section = if (ROI_complete) length(ids) else NA_integer_,
                                                                  n_TLS_supported_whole_section = count,
                                                                  supported_ROI_capture_fraction_QC = coverage_QC, TLS_density_per_mm2 = density,
                                                                  TLS_area_fraction = area_fraction,
                                                                  n_GC_like_supported_ROI_reviewed = if (n_stage_reviewed) sum(stage == "GC_like") else NA_integer_,
                                                                  GC_like_fraction_supported_TLS = if (stages_complete) mean(stage == "GC_like") else NA_real_)
}
structure_summary <- bind_rows(structure_summary)

##### 第6步：描述性比较，并检查可评估spot的测序复杂度分布影响 #####
# 共同分层加权只是技术敏感性检查，不是因果校正，也不能补回未评分spot。
depth_data <- tls_spots[tls_spots$AUC_assessable, ]
bin_variable <- function(x) {
  br <- unique(as.numeric(quantile(x, probs = seq(0, 1, length.out = depth_bins+1), names = FALSE)))
  if (length(br) < 2) rep(1L, length(x)) else as.integer(cut(x, br, include.lowest = TRUE))
}
depth_data$stratum <- if (nrow(depth_data)) paste(bin_variable(log1p(depth_data$host_UMI)), bin_variable(depth_data$nFeature_common), sep = "_") else character()
tab <- table(factor(depth_data$sample_id, levels = tls$samples$sample_id), depth_data$stratum)
keep <- if (ncol(tab)) colnames(tab)[apply(tab, 2, min) >= min_stratum_spots] else character()
weights <- if (length(keep)) apply(tab[, keep, drop = FALSE], 2, min) else numeric()
if (length(weights)) weights <- weights/sum(weights)
depth_summary <- depth_details <- list()
for (s in tls$samples$sample_id) {
  m <- depth_data[depth_data$sample_id == s, ]
  row <- data.frame(sample_id = s, n_common_depth_strata = length(keep),
                    n_spots_common_depth_strata = sum(m$stratum %in% keep),
                    fraction_AUC_spots_in_common_depth = if (nrow(m)) mean(m$stratum %in% keep) else NA_real_)
  for (module in intersect(c("TLS22", "Imprint29"), names(tls$gene_sets_used))) {
    values <- vapply(keep, function(k) {
      v <- m[[paste0(module, "_AUC")]][m$stratum == k]
      if (!length(v) || all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
    }, numeric(1))
    row[[paste0(module, "_AUC_depth_standardized")]] <- if (length(keep) && all(is.finite(values))) sum(values*weights[keep]) else NA_real_
    if (length(keep)) depth_details[[length(depth_details)+1L]] <- data.frame(sample_id=s, module=module, stratum=keep, weight=as.numeric(weights[keep]), mean_AUC=values)
  }
  depth_summary[[length(depth_summary)+1L]] <- row
}
depth_summary <- bind_rows(depth_summary)
depth_details <- bind_rows(depth_details)
new_cols <- union(setdiff(names(structure_summary), "sample_id"), setdiff(names(depth_summary), "sample_id"))
final_summary <- tls$summary[, setdiff(names(tls$summary), new_cols), drop = FALSE]
final_summary <- left_join(left_join(final_summary, structure_summary, by = "sample_id"), depth_summary, by = "sample_id")
for (s in final_summary$sample_id) for (module in c("TLS22","Imprint29")) {
  v <- tls_spots[[paste0(module,"_AUC")]][tls_spots$sample_id == s & tls_spots$AUC_assessable]
  final_summary[[paste0(module,"_AUC_mean")]][final_summary$sample_id == s] <- if (!length(v) || all(is.na(v))) NA_real_ else mean(v,na.rm=TRUE)
}
descriptive <- list()
for (metric in c("TLS22_AUC_mean", "Imprint29_AUC_mean", "TLS22_AUC_depth_standardized", "Imprint29_AUC_depth_standardized",
                 "n_union_candidate_ROI", "n_TLS_supported_whole_section", "supported_ROI_capture_fraction_QC", "TLS_density_per_mm2", "GC_like_fraction_supported_TLS")) {
  a <- final_summary[[metric]][final_summary$treatment == "Vehicle"]
  b <- final_summary[[metric]][final_summary$treatment == "VG161"]
  if (length(a) == 1 && length(b) == 1) descriptive[[length(descriptive)+1L]] <- data.frame(metric = metric,
                                                                                            Vehicle = a, VG161 = b, difference_VG161_minus_Vehicle = b-a,
                                                                                            ratio_VG161_over_Vehicle = if (is.finite(a) && is.finite(b) && a != 0) b/a else NA_real_)
}
descriptive <- bind_rows(descriptive)
print(final_summary[, c("sample_id", "n_union_candidate_ROI", "TLS_review_complete", "n_TLS_supported_whole_section", "TLS22_AUC_mean", "TLS22_AUC_depth_standardized")])
print(roi_evidence)

##### 复核图：ROI局部H&E、原始marker、组成及描述性比较 #####
# 本段可独立运行，沿用前半段的字体和配色。
review_group_colors <- c(Vehicle = "#426B91", VG161 = "#C35B3F")
review_expression_colors <- c("#EDF2F6", "#91AEC7", "#24537D")
review_unavailable <- "#B8BDC3"
review_strip_theme <- theme(strip.background = element_rect(fill = "#EDF2F6", color = NA),
                            strip.text = element_text(face = "bold", size = 9, color = "#22313F", hjust = 0.5,
                                                      margin = margin(3, 4, 3, 4)))
review_map_theme <- theme_void(base_size = 10, base_family = "Helvetica") + review_strip_theme +
  theme(plot.title = element_text(face = "bold", size = 10, color = "#22313F"),
        legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        plot.margin = margin(6, 8, 6, 8))
review_chart_theme <- theme_classic(base_size = 10, base_family = "Helvetica") + review_strip_theme +
  theme(plot.title = element_text(face = "bold", size = 10, color = "#22313F"),
        panel.border = element_rect(fill = NA, color = "#526170", linewidth = 0.35),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.35, color = "#526170"),
        strip.background = element_rect(fill = "#EDF2F6", color = "#526170", linewidth = 0.35),
        plot.margin = margin(8, 12, 8, 12))
review_score_strips <- c(TLS22 = "#DCE9F3", Imprint29 = "#DDEEEA")
review_marker_panels <- list(c("MS4A1", "CD79A", "CD3D", "TRBC2", "CCL19", "CXCL13"),
                             c("CR2", "FCER2", "AICDA", "BCL6", "MKI67", "SERPINA9"))
marker_limits <- vapply(roi_context, function(ctx) {
  counts <- tls$marker_counts[[ctx$sample_id]]
  genes <- intersect(paste0(prefix, unlist(review_marker_panels)), rownames(counts))
  if (!length(genes) || !length(ctx$context)) return(0)
  m <- tls_spots[tls_spots$sample_id == ctx$sample_id, ]
  host <- pmax(m$host_UMI[match(ctx$context, m$spot)], 1)
  scaled <- sweep(as.matrix(counts[genes, ctx$context, drop = FALSE]), 2, host, "/")
  max(log1p(scaled * 1e4))
}, numeric(1))
review_marker_max <- max(c(0, marker_limits), na.rm = TRUE)
if (review_marker_max == 0) review_marker_max <- 1
review_auc_max <- max(c(0, tls_spots$TLS22_AUC), na.rm = TRUE)
if (review_auc_max == 0) review_auc_max <- 1

review_view <- function(i) {
  ctx <- roi_context[[i]]
  m <- tls_spots[tls_spots$sample_id == ctx$sample_id, ]
  mm <- m[m$spot %in% ctx$context, ]
  j <- match(ctx$sample_id, tls$samples$sample_id)
  im <- tls$images[[j]]
  hp <- tls$HE_polygons_merged[tls$HE_polygons_merged$sample_id == ctx$sample_id & tls$HE_polygons_merged$roi_id == ctx$roi_id, ]
  xx <- c(mm$x, hp$x); yy <- c(mm$y, hp$y)
  if (!length(xx)) { xx <- m$x; yy <- m$y }
  padding <- max(diff(range(xx)), diff(range(yy)), 1)*0.15
  limits <- list(x = range(xx)+c(-padding, padding), y = range(yy)+c(-padding, padding))
  file <- if (ctx$sample_id %in% names(review_image_files)) review_image_files[[ctx$sample_id]] else im$file
  p <- ggplot()
  image_ok <- !is.null(file) && length(file) == 1 && file.exists(path.expand(file))
  if (image_ok) {
    if (!is.null(im$image_md5) && !identical(unname(tools::md5sum(path.expand(file))), im$image_md5)) stop("H&E文件指纹改变；不能沿用旧比例尺/ROI，请检查review_image_files。")
    he <- png::readPNG(path.expand(file))
    if (!identical(as.integer(dim(he)[1:2]),as.integer(c(im$height,im$width)))) stop("H&E尺寸与保存的坐标比例不一致，未自动拉伸。")
    step <- max(1L, ceiling(max(dim(he)[1:2])/review_image_max_px))
    rr <- seq(1, dim(he)[1], by = step); cc <- seq(1, dim(he)[2], by = step)
    he <- if (length(dim(he)) == 2) he[rr, cc, drop=FALSE] else he[rr, cc, , drop=FALSE]
    p <- p + annotation_raster(as.raster(he), xmin=0, xmax=im$width/im$scale, ymin=-im$height/im$scale, ymax=0)
  } else p <- p + annotate("text", x=mean(limits$x), y=-mean(limits$y), label="H&E unavailable", size=3)
  bb <- tls$expression_boundaries[tls$expression_boundaries$sample_id == ctx$sample_id & tls$expression_boundaries$roi_id == ctx$roi_id, ]
  p <- p + geom_segment(data = bb, aes(x = x, y = -y, xend = xend, yend = -yend,
                                       color = "Expression candidate"), linewidth = 0.5) +
    geom_polygon(data = hp, aes(x, -y, group = polygon_id, color = "H&E candidate"),
                 fill = NA, linewidth = 0.5, key_glyph = "path") +
    scale_color_manual(name = NULL, values = c(`Expression candidate` = "#BB521E", `H&E candidate` = "#24537D")) +
    coord_fixed(xlim=limits$x, ylim=-rev(limits$y), expand=FALSE) +
    facet_wrap(vars(!!paste(ctx$roi_id, "H&E"))) + review_map_theme
  list(plot=p, met=mm, limits=limits, core=ctx$core)
}
review_point_map <- function(m, view, title, upper, unit, unavailable = NULL) {
  m$spot_role <- ifelse(m$spot %in% view$core, "ROI spot", "Nearby spot")
  p <- ggplot(m, aes(x, -y)) +
    coord_fixed(xlim = view$limits$x, ylim = -rev(view$limits$y), expand = FALSE) +
    facet_wrap(vars(!!title)) + review_map_theme
  if (!nrow(m)) return(p)
  
  shapes <- c(`ROI spot` = 16, `Nearby spot` = 1)
  key_colors <- rep("#526170", length(shapes))
  # 追加灰点图例，同时保留不可评分点的原有内外圈点形。
  if (!is.null(unavailable)) {
    shapes[unavailable] <- 16
    key_colors <- c(key_colors, review_unavailable)
  }
  p + geom_point(aes(color = value, shape = spot_role), size = 2.6, stroke = 0.5,
                 show.legend = c(shape = TRUE)) +
    scale_color_gradientn(colours = review_expression_colors, na.value = review_unavailable,
                          limits = c(0, upper), name = unit, guide = guide_colorbar(order = 1)) +
    scale_shape_manual(name = NULL, values = shapes, limits = names(shapes),
                       guide = guide_legend(order = 2, override.aes = list(colour = key_colors)))
}
review_marker_map <- function(i, gene, view, unavailable) {
  s <- roi_context[[i]]$sample_id; ct <- tls$marker_counts[[s]]; feature <- paste0(prefix, gene)
  m <- view$met
  m$value <- if (feature %in% rownames(ct) && nrow(m)) log1p(as.numeric(ct[feature,m$spot])/pmax(m$host_UMI,1)*1e4) else rep(NA_real_,nrow(m))
  review_point_map(m, view, gene, review_marker_max, "log1p(CP10K)", unavailable)
}
pdf(file.path(review_dir, "review.pdf"), width = 9, height = 6, family = "Helvetica",
    useDingbats = FALSE, compress = TRUE)
tryCatch({
  for (i in seq_len(nrow(tls_rois))) {
    v <- review_view(i)
    score_spots <- v$met
    score_spots$value <- score_spots$TLS22_AUC
    score_unavailable <- if (anyNA(score_spots$value)) "Unscored" else NULL
    p_score <- review_point_map(score_spots, v, "TLS22", review_auc_max, "AUC", score_unavailable)
    if (!nrow(score_spots)) p_score <- p_score +
      annotate("text", x = mean(v$limits$x), y = -mean(v$limits$y), label = "No captured spots", size = 3)
    print(wrap_plots(v$plot, p_score, widths = c(1, 1)))
    available_features <- rownames(tls$marker_counts[[roi_context[[i]]$sample_id]])
    for (genes in review_marker_panels) {
      marker_unavailable <- if (any(!paste0(prefix, genes) %in% available_features)) "Gene unavailable" else NULL
      maps <- lapply(genes, function(gene) review_marker_map(i, gene, v, marker_unavailable))
      marker_page <- wrap_plots(lapply(maps, function(p) p + theme(legend.position = "none")), ncol = 3)
      if (nrow(v$met)) {
        # 整组marker共用一份图例，放在网格旁而不是某一张小图旁。
        guide_plot <- ggplotGrob(maps[[1L]])
        guide_box <- guide_plot$grobs[[which(guide_plot$layout$name %in% c("guide-box-right", "guide-box"))]]
        marker_page <- wrap_plots(marker_page, wrap_elements(full = guide_box), widths = c(1, 0.18))
      }
      print(marker_page + plot_annotation(title = paste(tls_rois$roi_id[i], "markers")))
    }
    comp <- composition[composition$roi_id == tls_rois$roi_id[i], ]
    pp <- lapply(c("TLS22","Imprint29"),function(panel) {
      d <- comp[comp$panel==panel, ]; d <- head(d[order(d$raw_UMI,decreasing=TRUE), ],10)
      d$gene <- factor(d$gene,levels=rev(d$gene))
      ggplot(d, aes(raw_UMI, gene)) + geom_col(fill = "#426B91", width = 0.7) + review_chart_theme +
        facet_wrap(vars(!!panel)) +
        theme(strip.background = element_rect(fill = review_score_strips[[panel]], color = "#526170", linewidth = 0.35)) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
        labs(x = "Core UMI", y = NULL)
    })
    print(wrap_plots(pp, ncol = 2) + plot_annotation(title = paste(tls_rois$roi_id[i], "marker counts")))
  }
  dd <- bind_rows(data.frame(sample_id=final_summary$sample_id,method="Unadjusted",value=final_summary$TLS22_AUC_mean),
                  data.frame(sample_id=final_summary$sample_id,method="Common depth strata",value=final_summary$TLS22_AUC_depth_standardized))
  dd <- left_join(dd,tls$samples[,c("sample_id","treatment")],by="sample_id")
  dd$treatment <- factor(dd$treatment,levels=c("Vehicle","VG161"))
  # 两组共用的方法图例放在右下空白处。
  score_legend_position <- c(0.98, 0.05)
  score_legend <- if (packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.position = "inside", legend.position.inside = score_legend_position)
  } else {
    theme(legend.position = score_legend_position)
  }
  score_panel <- ggplot(dd, aes(treatment, value, color = treatment, shape = method, group = method)) +
    geom_point(size = 3.1, stroke = 0.85, position = position_dodge(width = 0.24)) +
    scale_color_manual(values = review_group_colors, guide = "none") +
    scale_shape_manual(values = c(Unadjusted = 16, `Common depth strata` = 1),
                       guide = guide_legend(override.aes = list(colour = "#526170"))) +
    facet_wrap(vars("TLS22 scores")) + review_chart_theme + score_legend +
    theme(legend.justification = c(1, 0), legend.direction = "vertical",
          legend.background = element_blank(), legend.text = element_text(size = 8),
          strip.background = element_rect(fill = review_score_strips[["TLS22"]], color = "#526170", linewidth = 0.35)) +
    labs(x = NULL, y = "Mean AUC", shape = NULL)
  status <- final_summary[,c("sample_id","n_union_candidate_ROI","TLS_review_complete","n_TLS_supported_whole_section")]
  status_table <- data.frame(Sample = status$sample_id, Candidates = as.character(status$n_union_candidate_ROI),
                             Review = ifelse(status$TLS_review_complete, "Complete", "Pending"),
                             `Supported TLS` = ifelse(is.na(status$n_TLS_supported_whole_section), "NA", status$n_TLS_supported_whole_section),
                             check.names = FALSE)
  cells <- stack(status_table)
  cells$row <- rep(rev(seq_len(nrow(status_table))), ncol(status_table))
  review_panel <- ggplot(cells, aes(ind, row)) +
    geom_tile(fill = "#F0F4F7", color = "white", linewidth = 0.5, width = 1, height = 0.95) +
    geom_text(aes(label = values), family = "Helvetica", size = 3.2, color = "#22313F") +
    scale_x_discrete(position = "top") + review_map_theme +
    theme(axis.text.x = element_text(face = "bold", size = 9, color = "#22313F", margin = margin(b = 4))) +
    labs(title = "TLS review")
  print(wrap_plots(score_panel, review_panel, ncol = 1, heights = c(2.4, 1)))
},finally=dev.off())

##### 沿用原结果：只新增review.csv和review.pdf #####
write.csv(review,path.expand(review_file),row.names=FALSE,na="NA")
tls$summary <- final_summary
tls$followup <- list(review=review,roi_evidence=roi_evidence,gene_evidence=gene_evidence,
                     raw_marker_composition=composition,ROI_context=roi_context,structure_summary=structure_summary,
                     depth_summary=depth_summary,depth_details=depth_details,descriptive_comparison=descriptive,
                     info=list(run_time=format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z"),context_rings=context_rings,
                               depth_bins=depth_bins,min_stratum_spots=min_stratum_spots,review_file=review_file,
                               interpretation=c("RNA evidence and manual decisions are separate; no ROI is automatically confirmed.",
                                                "Missing or uncaptured marker evidence is NA; measured zeros remain zero.",
                                                "Whole-section supported counts stay NA until full H&E and all candidates are resolved.",
                                                "Grid/capture coverage is not histological area; physical metrics require measured areas.",
                                                "Raw UMI fractions are not exact AUCell contribution fractions.",
                                                "Depth standardization is a descriptive sensitivity analysis, not causal correction or recovery of excluded spots.",
                                                "The D3 comparison uses one slice per group. AUCell summaries use evaluable captured spots."),
                               sessionInfo=capture.output(sessionInfo())))
saveRDS(tls,path.expand(review_input),compress="gzip")
write.csv(final_summary,file.path(review_dir,"summary.csv"),row.names=FALSE)
results <- tls
message("第4—6步完成：查看review.pdf，并在review.csv回填审核列后重跑本段。")
