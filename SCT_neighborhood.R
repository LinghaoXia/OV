##### 清理运行环境 #####
rm(list = ls(all.names = TRUE))
graphics.off()
invisible(gc())

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

##### 输入路径与参数 #####
sc_file <- "~/rawdata/scRNA_virus/virus_3D/virus_3D_anno.RData"
sc_object <- "sce"
plot_sc_violin <- TRUE  # 有单细胞原始counts时在PDF末尾画逐细胞表达。
sample_col <- "orig.ident"
celltype_col <- "celltype"
sc_groups <- c(Vehicle_3D_1 = "Vehicle", Vehicle_3D_2 = "Vehicle",
               VG161_3D_1 = "VG161", VG161_3D_2 = "VG161")
st_root <- "~/rawdata/SCT/virus_3D"
samples <- data.frame(
  sample_id = c("Vehicle_3D_1", "VG161_3D_1"),
  treatment = c("Vehicle", "VG161"),
  data_dir = file.path(st_root, c("Vehicle_3D_1", "VG161_3D_1")),
  composition_csv = ""  # 留空运行RCTD。复用比例时填写CSV，第一列为spot，其余列为各细胞类型。
)
outdir <- "~/OV/SCT/MDK_NCL_D3"
reuse_details_rds <- ""  # 填旧details.rds可复用坐标、计数和RCTD比例，仍需下方的CMAP标签。
mapping_files <- c(VG161_3D_1 = "~/rawdata/SCT/analysis/virus_3D_1/virus_3D_anno.RData")
mapping_object <- "sct"
mapping_label_col <- "celltype_new"  # SCT.R完成CMAP并保存后的metadata列。
matrix_file <- "filtered_feature_bc_matrix.h5"
virus_group <- "VG161"
epithelial_types <- c("Epithelial")  # 完整组成保留这些类别，计算非上皮分母时排除。
sc_host_prefix <- ""  # 若人鼠混合且有物种前缀，填人源前缀。当前D3符号无前缀。
st_host_prefix <- ""
all_viral_pattern <- "^VG161[-_]"
viral_pattern <- "^VG161[-_](UL|ICP|US|IRL)"
viral_genes_use <- NULL  # 病毒UMI只做质控，可用精确病毒基因名替代正则。
min_host_UMI <- 100L
rings <- c(1L, 2L)
min_neighbors <- 2L
rctd_cores <- 1L
reference_min_cells <- 25L
reference_max_cells <- 1500L
seed <- 123L
# RCTD自动识别GitHub版和Bioconductor版spacexr。缺包时在运行脚本的同一个R环境安装：
# install.packages("remotes")  # 已有remotes时跳过。
# options(timeout = max(600, getOption("timeout")))
# remotes::install_github("dmcable/spacexr", build_vignettes = FALSE, upgrade = "never")
# R >= 4.5也可用对应版本的Bioconductor安装：BiocManager::install("spacexr", update = FALSE)
# 组成是参考覆盖的细胞比例估计。表达与共检出用于描述空间邻近证据。

##### 原始计数和配准位置 #####
check_counts <- function(ct, label) {
  if (is.null(rownames(ct)) || is.null(colnames(ct)) || anyNA(rownames(ct)) ||
      anyNA(colnames(ct)) || anyDuplicated(rownames(ct)) || anyDuplicated(colnames(ct))) {
    stop(label, "：counts的基因名或条码缺失、重复。")
  }
  if (!inherits(ct, "sparseMatrix")) ct <- Matrix(ct, sparse = TRUE)
  ct <- as(ct, "dgCMatrix")
  if (any(!is.finite(ct@x) | ct@x < 0 | ct@x != round(ct@x))) {
    stop(label, "：需要非负整数原始counts。")
  }
  ct
}

read_seurat <- function(path, object_name) {
  path <- path.expand(path)
  if (!file.exists(path)) stop("找不到Seurat文件：", path)
  if (grepl("\\.(rdata|rda)$", path, ignore.case = TRUE)) {
    e <- new.env()
    loaded <- load(path, envir = e)
    if (!object_name %in% loaded) stop("RData中没有对象：", object_name)
    object <- e[[object_name]]
  } else {
    object <- readRDS(path)
  }
  if (!inherits(object, "Seurat")) stop("需要Seurat对象：", path)
  object
}

read_sc <- function(path, object_name, assay = "RNA") {
  sce <- read_seurat(path, object_name)
  if (!assay %in% names(sce@assays)) stop("需要带RNA assay的Seurat对象。")
  meta <- sce[[]]
  if (!all(c(sample_col, celltype_col) %in% names(meta))) stop("单细胞metadata缺少样本列或注释列。")
  id <- as.character(meta[[sample_col]])
  if (anyNA(id) || any(!id %in% names(sc_groups))) {
    stop("请用sc_groups映射单细胞样本列的所有取值：", paste(unique(id), collapse = ", "))
  }
  a <- sce[[assay]]
  if (inherits(a, "Assay5")) {
    layers <- Layers(a, search = "^counts($|\\.)")
    if (!length(layers)) stop("RNA assay中没有counts层。")
    cells <- lapply(layers, function(k) colnames(LayerData(a, layer = k, fast = FALSE)))
    features <- lapply(layers, function(k) rownames(LayerData(a, layer = k, fast = FALSE)))
    barcodes <- unlist(cells, use.names = FALSE)
    if (anyDuplicated(barcodes)) stop("RNA counts层含重复细胞，请先确定正确的原始层。")
    cell_layer <- setNames(rep(seq_along(layers), lengths(cells)), barcodes)
    if (length(layers) > 1L) a <- JoinLayers(a, layers = "^counts($|\\.)", new = "counts")
    ct <- LayerData(a, layer = if (length(layers) > 1L) "counts" else layers, fast = FALSE)
  } else {
    ct <- if (packageVersion("SeuratObject") >= "5.0.0") GetAssayData(sce, assay = assay, layer = "counts") else
      GetAssayData(sce, assay = assay, slot = "counts")
    features <- list(rownames(ct))
    cell_layer <- setNames(rep(1L, ncol(ct)), colnames(ct))
  }
  ct <- check_counts(ct, "RNA")
  if (!all(colnames(ct) %in% rownames(meta))) stop("RNA条码无法匹配metadata。")
  list(counts = ct, meta = meta[colnames(ct), , drop = FALSE], features = features,
       layer = unname(cell_layer[colnames(ct)]), missing_cells = nrow(meta) - ncol(ct))
}

read_spatial <- function(s) {
  root <- path.expand(s$data_dir)
  path <- file.path(root, matrix_file)
  if (!file.exists(path)) stop(s$sample_id, "找不到原始矩阵：", path, "。请修改samples$data_dir。")
  ct <- Read10X_h5(path, use.names = TRUE)
  if (is.list(ct)) ct <- ct[["Gene Expression"]]
  ct <- check_counts(ct, s$sample_id)
  path <- file.path(root, "spatial", "tissue_positions.csv")
  pos_cols <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  if (file.exists(path)) {
    pos <- read.csv(path, stringsAsFactors = FALSE)
  } else {
    path <- file.path(root, "spatial", "tissue_positions_list.csv")
    if (!file.exists(path)) stop(s$sample_id, "缺少Space Ranger位置表。")
    pos <- read.csv(path, header = FALSE, col.names = pos_cols, stringsAsFactors = FALSE)
  }
  if (!all(pos_cols %in% names(pos)) || anyNA(pos$barcode) || anyDuplicated(pos$barcode) ||
      !all(colnames(ct) %in% pos$barcode)) stop(s$sample_id, "的条码与位置表不一致。")
  pos <- pos[match(colnames(ct), pos$barcode), pos_cols]
  if (anyNA(pos$in_tissue) || !all(pos$in_tissue %in% c(0, 1))) stop("in_tissue应为0或1。")
  pos <- pos[pos$in_tissue == 1, ]
  if (!nrow(pos) || any(!is.finite(pos$pxl_col_in_fullres)) ||
      any(!is.finite(pos$pxl_row_in_fullres))) stop(s$sample_id, "缺少有效组织坐标。")
  spots <- data.frame(spot = pos$barcode, array_row = pos$array_row, array_col = pos$array_col,
                      x = pos$pxl_col_in_fullres, y = pos$pxl_row_in_fullres)
  list(counts = ct[, spots$spot, drop = FALSE], spots = spots)
}

##### 导入SCT.R的CMAP代表细胞标签 #####
check_mapping_files <- function(ids) {
  if (!length(ids) || is.null(names(mapping_files)) || anyDuplicated(names(mapping_files)) ||
      any(!ids %in% names(mapping_files))) {
    stop("请在mapping_files中按样本名指定每张VG161切片的CMAP结果。")
  }
  paths <- path.expand(unname(mapping_files[ids]))
  if (anyNA(paths) || any(!nzchar(paths)) || any(!file.exists(paths))) {
    stop("找不到CMAP结果，请检查mapping_files：", paste(ids[is.na(paths) | !nzchar(paths) | !file.exists(paths)],
                                              collapse = ", "))
  }
  if (anyDuplicated(normalizePath(paths))) stop("不同切片须使用各自的CMAP结果文件。")
}

read_mapped_labels <- function(file, spots, sample_id) {
  file <- path.expand(file)
  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    z <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
    if ("sample_id" %in% names(z) && !all(!is.na(z$sample_id) & z$sample_id == sample_id)) {
      stop(sample_id, "：映射CSV的sample_id与切片不一致。")
    }
  } else {
    sct <- read_seurat(file, mapping_object)
    z <- sct[[]]
    z$spot <- rownames(z)
  }
  if (!all(c("spot", mapping_label_col) %in% names(z)) || anyDuplicated(names(z)) ||
      anyNA(z$spot) || any(!nzchar(z$spot)) || anyDuplicated(z$spot)) {
    stop(sample_id, "：需要spot及", mapping_label_col, "，且条码须非空、唯一。请核对是否为CMAP保存后的对象。")
  }
  if (!nrow(z) || any(!z$spot %in% spots)) {
    stop(sample_id, "：映射结果有无法匹配的条码，请确认文件来自同一切片。")
  }
  label <- as.character(z[[mapping_label_col]][match(spots, z$spot)])
  label[!is.na(label) & !nzchar(label)] <- NA_character_
  if (!all(c("Infected", "Bystander") %in% label)) {
    stop(sample_id, "：CMAP结果须包含Infected和Bystander。请核对文件版本及标签原始拼写。")
  }
  label
}

host_counts <- function(ct, prefix) {
  keep <- !grepl(all_viral_pattern, rownames(ct)) & startsWith(rownames(ct), prefix)
  m <- ct[keep, , drop = FALSE]
  rownames(m) <- substring(rownames(m), nchar(prefix) + 1L)
  if (!all(c("MDK", "NCL") %in% rownames(m))) stop("缺少人源MDK或NCL，请核对计数矩阵和物种前缀。")
  m
}

viral_features <- function(genes) {
  v <- if (is.null(viral_genes_use)) grep(viral_pattern, genes, value = TRUE) else unique(viral_genes_use)
  if (length(v) && !all(grepl(all_viral_pattern, v))) stop("all_viral_pattern必须覆盖用于判定的病毒基因。")
  v
}

extract_values <- function(ct, host, vg) {
  genes <- intersect(rownames(ct), vg)
  total <- Matrix::colSums(host)
  viral <- if (length(vg) && length(genes) == length(vg)) Matrix::colSums(ct[genes, , drop = FALSE]) else
    rep(NA_real_, ncol(ct))
  d <- data.frame(host_UMI = unname(total), viral_UMI = unname(viral))
  for (gene in c("MDK", "NCL")) {
    raw <- as.numeric(host[gene, ])
    d[[paste0(gene, "_UMI")]] <- raw
    d[[paste0(gene, "_logCP10K")]] <- ifelse(total > 0, log1p(1e4 * raw / pmax(total, 1)), NA_real_)
  }
  d
}

##### 细胞组成与spot邻域 #####
read_composition <- function(file, spots) {
  z <- read.csv(path.expand(file), stringsAsFactors = FALSE, check.names = FALSE)
  if (!"spot" %in% names(z) || anyNA(z$spot) || anyDuplicated(z$spot) || anyDuplicated(names(z))) {
    stop("composition_csv须有唯一spot及各细胞类型比例列。")
  }
  m <- as.matrix(z[, names(z) != "spot", drop = FALSE])
  if (ncol(m) < 2L || !is.numeric(m) || any(!is.finite(m) | m < 0 | m > 1) ||
      any(abs(rowSums(m) - 1) > 0.02)) {
    stop("composition_csv应给出完整组成比例，范围0–1，行和接近1。请删除resid等非细胞类型列。")
  }
  rownames(m) <- z$spot
  m <- m / rowSums(m)
  if (!any(spots %in% rownames(m))) stop("composition_csv没有匹配到该切片的QC条码。")
  m
}

run_rctd <- function(host, spots, reference, genes, backend) {
  ct <- host[genes, spots$spot, drop = FALSE]
  coords <- data.frame(x = spots$x, y = spots$y, row.names = spots$spot)
  if (backend == "GitHub") {
    puck <- spacexr::SpatialRNA(coords, ct, Matrix::colSums(ct))
    fit <- spacexr::create.RCTD(puck, reference, max_cores = rctd_cores, CELL_MIN_INSTANCE = reference_min_cells)
    fit <- spacexr::run.RCTD(fit, doublet_mode = "full")
    w <- as.matrix(fit@results$weights)
  } else {
    puck <- SpatialExperiment::SpatialExperiment(assays = list(counts = ct), spatialCoords = as.matrix(coords))
    fit <- spacexr::createRctd(puck, reference, ref_n_cells_min = reference_min_cells,
                               ref_n_cells_max = reference_max_cells)
    fit <- spacexr::runRctd(fit, rctd_mode = "full", max_cores = rctd_cores)
    # Bioconductor返回细胞类型 × spot，统一成spot × 细胞类型。
    w <- t(as.matrix(SummarizedExperiment::assay(fit, "weights")))
  }
  weight_sum <- rowSums(w)
  valid <- is.finite(weight_sum) & weight_sum > 0
  w[valid, , drop = FALSE] / weight_sum[valid]
}

##### Visium一/二环邻域：用QC后的阵列位置建图 #####
make_rings <- function(spots, rings = c(1L, 2L)) {
  if (!is.numeric(rings) || anyNA(rings) || any(!rings %in% c(1L, 2L))) {
    stop("rings须为1或2。")
  }
  rings <- sort(unique(c(1L, as.integer(rings))))
  if (!all(c("spot", "array_row", "array_col") %in% names(spots))) {
    stop("spots须包含spot、array_row和array_col。")
  }
  row <- spots$array_row
  col <- spots$array_col
  if (!is.numeric(row) || !is.numeric(col) || any(!is.finite(row)) ||
      any(!is.finite(col)) || any(row != round(row)) || any(col != round(col)) ||
      anyNA(spots$spot) || anyDuplicated(spots$spot)) {
    stop("spot和Visium阵列坐标须有效。")
  }
  key <- paste(row, col, sep = ":")
  if (anyDuplicated(key)) stop("单切片内Visium阵列坐标有重复。")
  
  offsets <- rbind(c(0L, -2L), c(0L, 2L), c(-1L, -1L),
                   c(-1L, 1L), c(1L, -1L), c(1L, 1L))
  neighbor_index <- lapply(seq_len(nrow(offsets)), function(k) {
    match(paste(row + offsets[k, 1], col + offsets[k, 2], sep = ":"), key)
  })
  first <- lapply(seq_len(nrow(spots)), function(i) {
    ix <- vapply(neighbor_index, function(z) z[i], integer(1))
    sort(ix[!is.na(ix)])
  })
  names(first) <- spots$spot
  out <- list("1" = first)
  if (2L %in% rings) {
    out[["2"]] <- lapply(seq_len(nrow(spots)), function(i) {
      reached <- as.integer(unlist(first[first[[i]]], use.names = FALSE))
      sort(setdiff(reached, c(i, first[[i]])))
    })
    names(out[["2"]]) <- spots$spot
  }
  out
}

mean_available <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_

ring_metrics <- function(neighbors, spots, fractions) {
  comp <- t(vapply(neighbors, function(ix) {
    ix <- ix[complete.cases(fractions[ix, , drop = FALSE])]
    if (length(ix) < min_neighbors) return(rep(NA_real_, ncol(fractions)))
    colMeans(fractions[ix, , drop = FALSE])
  }, numeric(ncol(fractions))))
  colnames(comp) <- paste0("Fraction.", colnames(fractions))
  mdk <- vapply(neighbors, function(ix) mean_available(spots$MDK_logCP10K[ix]), numeric(1))
  ncl <- vapply(neighbors, function(ix) mean_available(spots$NCL_logCP10K[ix]), numeric(1))
  # 每个中心的共检出是满足中心基因和邻点基因均有UMI的邻点比例，中心基因未检出时为0。
  mdk_ncl <- vapply(seq_along(neighbors), function(i) {
    mean_available(spots$MDK_UMI[i] > 0 & spots$NCL_UMI[neighbors[[i]]] > 0)
  }, numeric(1))
  ncl_mdk <- vapply(seq_along(neighbors), function(i) {
    mean_available(spots$NCL_UMI[i] > 0 & spots$MDK_UMI[neighbors[[i]]] > 0)
  }, numeric(1))
  expr <- cbind(Neighbor_MDK = mdk, Neighbor_NCL = ncl,
                MDK_center_NCL_neighbor = mdk_ncl, NCL_center_MDK_neighbor = ncl_mdk)
  expr[lengths(neighbors) < min_neighbors, ] <- NA_real_
  cbind(comp, expr)
}

center_summary <- function(x, labels, ring, s, analysis = NULL) {
  bind_rows(lapply(colnames(x), function(metric) {
    a <- x[labels %in% "Infected", metric]
    b <- x[labels %in% "Bystander", metric]
    mean_a <- mean_available(a)
    mean_b <- mean_available(b)
    section <- if (startsWith(metric, "Fraction.")) "Composition" else "MDK_NCL"
    data.frame(analysis = if (is.null(analysis)) section else analysis,
               sample_id = s$sample_id, treatment = s$treatment, ring = ring,
               metric = metric, n_Infected = sum(is.finite(a)), n_Bystander = sum(is.finite(b)),
               Infected = mean_a, Bystander = mean_b, difference = mean_a - mean_b)
  }))
}

virus_summary <- function(d, source) {
  d %>% group_by(sample_id, treatment, state) %>%
    summarise(analysis = source, metric = "Viral UMI QC", n = n(),
              n_assessable = sum(is.finite(viral_UMI)), n_with_viral_UMI = sum(viral_UMI > 0, na.rm = TRUE),
              median_viral_UMI = median(viral_UMI, na.rm = TRUE), .groups = "drop")
}

##### 单细胞表达及参考准备 #####
single_cell_values <- function(sc, host, vg) {
  d <- extract_values(sc$counts, host, vg)
  # JoinLayers补齐的结构零仍按原counts层标为不可判定。
  has_virus <- vapply(sc$features, function(x) length(vg) > 0L && all(vg %in% x), logical(1))
  d$viral_UMI[!has_virus[sc$layer]] <- NA_real_
  for (gene in c("MDK", "NCL")) {
    available <- vapply(sc$features, function(x) paste0(sc_host_prefix, gene) %in% x, logical(1))
    d[!available[sc$layer], paste0(gene, c("_UMI", "_logCP10K"))] <- NA_real_
  }
  d$sample_id <- as.character(sc$meta[[sample_col]])
  d$treatment <- unname(sc_groups[d$sample_id])
  d$cell_type <- as.character(sc$meta[[celltype_col]])
  d$state <- if ("infected" %in% names(sc$meta)) as.character(sc$meta$infected) else NA_character_
  d
}

single_cell_plot_values <- function(d) {
  epi <- filter(d, host_UMI > 0, cell_type %in% epithelial_types)
  labeled <- filter(epi, treatment == virus_group, state %in% c("Infected", "Bystander"))
  if (all(c("Infected", "Bystander") %in% labeled$state)) {
    epi <- labeled
  } else {
    epi$state <- epi$treatment
  }
  epi %>% select(sample_id, state, MDK_logCP10K, NCL_logCP10K)
}

read_sc_plot_values <- function() {
  sc <- read_sc(sc_file, sc_object)
  host <- host_counts(sc$counts, sc_host_prefix)
  single_cell_plot_values(single_cell_values(sc, host, viral_features(rownames(sc$counts))))
}

if (!all(rings %in% c(1L, 2L)) || !1L %in% rings || min_neighbors < 1 ||
    min_neighbors > 6 || reference_max_cells < reference_min_cells) stop("请核对邻域及参考参数。")
sc_plot <- NULL
if (!nzchar(reuse_details_rds)) {
  if (anyDuplicated(samples$sample_id) || anyDuplicated(path.expand(samples$data_dir))) stop("切片名及路径须唯一。")
  check_mapping_files(samples$sample_id[samples$treatment == virus_group])
  use_rctd <- !nzchar(samples$composition_csv)
  rctd_backend <- rctd_version <- NA_character_
  if (any(use_rctd)) {
    if (!length(find.package("spacexr", quiet = TRUE))) {
      stop("当前R库中未找到spacexr。请在运行脚本的同一个R环境安装，命令在脚本开头。\n",
           "R库路径：", paste(.libPaths(), collapse = " | "), call. = FALSE)
    }
    ns <- tryCatch(loadNamespace("spacexr"), error = identity)
    if (inherits(ns, "error")) {
      stop("已找到spacexr，但加载失败：\n", conditionMessage(ns), call. = FALSE)
    }
    rctd_version <- as.character(getNamespaceVersion(ns))
    api <- getNamespaceExports(ns)
    if (all(c("Reference", "SpatialRNA", "create.RCTD", "run.RCTD") %in% api)) {
      rctd_backend <- "GitHub"
    } else if (all(c("createRctd", "runRctd") %in% api)) {
      rctd_backend <- "Bioconductor"
    } else {
      stop("当前spacexr ", rctd_version, " 的接口不支持。请按脚本开头的命令安装支持的版本。\n",
           "包路径：", getNamespaceInfo(ns, "path"), call. = FALSE)
    }
    message("RCTD使用 ", rctd_backend, " 版spacexr ", rctd_version)
  }
  set.seed(seed)
  sc <- read_sc(sc_file, sc_object)
  sc_host <- host_counts(sc$counts, sc_host_prefix)
  sc_vg <- viral_features(rownames(sc$counts))
  sc_data <- single_cell_values(sc, sc_host, sc_vg)
  if (plot_sc_violin) sc_plot <- single_cell_plot_values(sc_data)
  sc_ok <- sc_data$host_UMI > 0 & !is.na(sc_data$cell_type) & nzchar(sc_data$cell_type)
  epithelial_sc <- sc_ok & sc_data$cell_type %in% epithelial_types
  sc_all <- sc_data[sc_ok, ]
  sc_all$state <- "All"
  i <- epithelial_sc & sc_data$treatment == virus_group & sc_data$state %in% c("Infected", "Bystander")
  sc_summary <- bind_rows(lapply(c("MDK", "NCL"), function(g) {
    bind_rows(sc_all, sc_data[i, ]) %>% group_by(sample_id, treatment, cell_type, state) %>%
      summarise(analysis = "scRNA expression", metric = g, n = n(),
                n_assessable = sum(!is.na(.data[[paste0(g, "_UMI")]])),
                fraction_detected = mean_available(.data[[paste0(g, "_UMI")]] > 0),
                mean_logCP10K = mean_available(.data[[paste0(g, "_logCP10K")]]), .groups = "drop")
  }))
  out <- list(virus_summary(bind_rows(mutate(sc_data[epithelial_sc, ], state = "All"), sc_data[i, ]), "scRNA virus QC"),
              sc_summary)
  
  ##### 逐切片分析 #####
  sp <- lapply(seq_len(nrow(samples)), function(i) read_spatial(samples[i, ]))
  for (i in which(samples$treatment == virus_group)) {
    id <- samples$sample_id[i]
    sp[[i]]$spots$mapping_label <- read_mapped_labels(mapping_files[[id]], sp[[i]]$spots$spot, id)
  }
  st_vg <- viral_features(unique(unlist(lapply(sp, function(z) rownames(z$counts)))))
  st_hosts <- lapply(sp, function(z) host_counts(z$counts, st_host_prefix))
  reference <- NULL
  genes <- character()
  ref_n <- table(sc_data$cell_type[sc_ok])
  if (any(use_rctd)) {
    # 只用各实际拟合切片及各参考counts层共享的宿主基因。排除MDK/NCL，避免参与分型后再用来验证。
    shared_sc <- Reduce(intersect, sc$features)
    shared_sc <- substring(shared_sc[startsWith(shared_sc, sc_host_prefix)], nchar(sc_host_prefix) + 1L)
    genes <- setdiff(Reduce(intersect, c(list(intersect(rownames(sc_host), shared_sc)),
                                         lapply(st_hosts[use_rctd], rownames))), c("MDK", "NCL"))
    if (length(genes) < 500L) stop("反卷积共有宿主基因过少，请核对基因命名和物种。")
    ct <- sc_host[genes, , drop = FALSE]
    depth <- Matrix::colSums(ct)
    ref_ok <- sc_ok & depth >= 100
    ref_n <- table(sc_data$cell_type[ref_ok])
    types <- names(ref_n)[ref_n >= reference_min_cells]
    if (length(types) < 2L) stop("达到最少参考细胞数的细胞类型不足两个。")
    ix <- which(ref_ok & sc_data$cell_type %in% types)
    labels <- setNames(factor(sc_data$cell_type[ix]), colnames(ct)[ix])
    if (rctd_backend == "GitHub") {
      reference <- spacexr::Reference(ct[, ix, drop = FALSE], labels, depth[ix],
                                      n_max_cells = reference_max_cells)
    } else {
      # 在切片循环前抽样一次，确保各切片复用同一批参考细胞。
      cells <- unlist(lapply(split(names(labels), labels), function(ids) {
        if (length(ids) > reference_max_cells) sample(ids, reference_max_cells) else ids
      }), use.names = FALSE)
      reference <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = ct[, cells, drop = FALSE]),
        colData = S4Vectors::DataFrame(cell_type = unname(labels[cells]), row.names = cells))
    }
  }
  
  details <- setNames(vector("list", nrow(samples)), samples$sample_id)
  for (i in seq_len(nrow(samples))) {
    s <- samples[i, ]
    message("分析切片 ", s$sample_id)
    spots <- cbind(sp[[i]]$spots, extract_values(sp[[i]]$counts, st_hosts[[i]], st_vg))
    spots$sample_id <- s$sample_id
    spots$treatment <- s$treatment
    spots$QC_pass <- spots$host_UMI >= min_host_UMI & spots$host_UMI > 0
    if (!"mapping_label" %in% names(spots)) spots$mapping_label <- NA_character_
    q <- spots[spots$QC_pass, ]
    if (!nrow(q)) stop(s$sample_id, "没有通过宿主UMI质控的spot。")
    m <- if (use_rctd[i]) run_rctd(st_hosts[[i]], q, reference, genes, rctd_backend) else
      read_composition(s$composition_csv, q$spot)
    fractions <- m[match(q$spot, rownames(m)), , drop = FALSE]
    rownames(fractions) <- q$spot
    graph <- make_rings(q, rings)
    details[[i]] <- list(spots = spots, fractions = fractions, rings = lapply(graph, function(z) {
      lapply(z, function(ix) q$spot[ix])
    }), virus_features = intersect(st_vg, rownames(sp[[i]]$counts)))
  }
  info <- list(sample_table = samples, sc_file = sc_file, sc_groups = sc_groups,
               sc_viral_features = sc_vg, sc_feature_layers = lapply(sc$features, intersect, y = sc_vg),
               sc_missing_counts = sc$missing_cells, st_viral_features = st_vg,
               sc_host_prefix = sc_host_prefix, st_host_prefix = st_host_prefix,
               reference_counts = ref_n, n_RCTD_genes = length(genes),
               rctd_backend = rctd_backend, spacexr_version = rctd_version, min_host_UMI = min_host_UMI,
               seed = seed, counts_session = capture.output(sessionInfo()))
} else {
  if (!file.exists(path.expand(reuse_details_rds))) stop("找不到reuse_details_rds。")
  cache <- readRDS(path.expand(reuse_details_rds))
  details <- cache$spatial
  info <- cache$info
  if (!length(details) || is.null(names(details)) || anyDuplicated(names(details)) || !is.list(info)) {
    stop("reuse_details_rds需要本脚本保存的details.rds。")
  }
  treated <- names(details)[vapply(details, function(d) identical(as.character(d$spots$treatment[1L]), virus_group),
                                   logical(1))]
  check_mapping_files(treated)
  for (id in names(details)) {
    spots <- details[[id]]$spots$spot
    details[[id]]$spots$mapping_label <- if (id %in% treated) {
      read_mapped_labels(mapping_files[[id]], spots, id)
    } else rep(NA_character_, length(spots))
  }
  if (plot_sc_violin) {
    if (file.exists(path.expand(sc_file))) {
      sc_plot <- read_sc_plot_values()
    } else {
      message("未找到sc_file，继续输出空间图。")
    }
  }
  out <- list()
  # 如旧汇总在同一目录，保留单细胞表达和病毒QC记录。
  previous_summary <- file.path(dirname(path.expand(reuse_details_rds)), "summary.csv")
  if (file.exists(previous_summary)) {
    old <- read.csv(previous_summary, check.names = FALSE)
    if (all(c("analysis", "state") %in% names(old))) {
      out <- list(filter(old, analysis == "scRNA expression", state %in% c("All", "Infected", "Bystander")) %>%
                    select(any_of(c("sample_id", "treatment", "analysis", "metric", "n", "n_assessable",
                                    "cell_type", "state", "fraction_detected", "mean_logCP10K"))))
    }
  }
}

##### 全部可评估中心的邻域组成和表达 #####
plot_centers <- list()
for (id in names(details)) {
  d <- details[[id]]
  d$pairs <- d$epithelial_pairs <- d$centers <- NULL  # 按本次标签重新生成中心表。
  d$spots$region <- d$spots$tumor_fraction <- d$spots$n_composition_neighbors <- NULL
  details[[id]] <- d
  s <- list(sample_id = id, treatment = d$spots$treatment[1L])
  q <- d$spots[d$spots$QC_pass, , drop = FALSE]
  if (!all(as.character(rings) %in% names(d$rings)) || !all(q$spot %in% names(d$rings[["1"]]))) {
    stop(id, "缺少所需的邻接图，请从原始输入重跑。")
  }
  fractions <- d$fractions[q$spot, , drop = FALSE]
  graph <- lapply(d$rings[as.character(rings)], function(z) {
    lapply(z[q$spot], function(ids) match(ids, q$spot))
  })
  if (anyNA(unlist(graph, use.names = FALSE))) stop(id, "邻接图含无法匹配的条码。")
  is_center <- q$mapping_label %in% c("Infected", "Bystander")
  mapped_q <- if (s$treatment == virus_group) {
    mutate(q[is_center, ], state = mapping_label)
  } else q[0, ]
  out[[length(out) + 1L]] <- virus_summary(bind_rows(mutate(q, state = "All"), mapped_q), "Spatial virus QC")
  out[[length(out) + 1L]] <- data.frame(analysis = "Composition QC", s,
                                        metric = c("QC spots", "Spots with composition"), value = c(nrow(q), sum(complete.cases(fractions))))
  if (s$treatment != virus_group) next
  if (length(epithelial_types) && !any(colnames(fractions) %in% epithelial_types)) {
    warning(id, "的组成里没有匹配到epithelial_types，非上皮分母会等于完整组成。请核对类别名。", call. = FALSE)
  }
  
  out[[length(out) + 1L]] <- d$spots %>%
    mutate(state = coalesce(mapping_label, "Unassigned")) %>% group_by(sample_id, treatment, state) %>%
    summarise(analysis = "Mapping QC", metric = "QC spots", n = n(), value = sum(QC_pass), .groups = "drop")
  centers <- data.frame(spot = q$spot[is_center], state = q$mapping_label[is_center])
  details[[id]]$centers <- centers
  
  # 先计算每个中心自己的环内均值，再在每组内对有值的中心等权平均。
  center <- cbind(fractions, Center_MDK = q$MDK_logCP10K, Center_NCL = q$NCL_logCP10K)
  colnames(center)[seq_len(ncol(fractions))] <- paste0("Fraction.", colnames(fractions))
  mm <- c(list("0" = center), lapply(graph, ring_metrics, spots = q, fractions = fractions))
  mm <- lapply(mm, function(m) m[is_center, , drop = FALSE])
  for (j in names(mm)) {
    ring <- as.integer(j)
    metrics <- mm[[j]]
    out[[length(out) + 1L]] <- center_summary(metrics, centers$state, ring, s)
    composition <- metrics[, startsWith(colnames(metrics), "Fraction."), drop = FALSE]
    non_epithelial <- composition[, !colnames(composition) %in% paste0("Fraction.", epithelial_types), drop = FALSE]
    # 先在同一中心的环内平均原比例，再以该中心的非上皮总占比作分母。
    share <- rowSums(non_epithelial)
    share[!complete.cases(composition)] <- NA_real_
    denominator <- cbind(`Non-epithelial share` = share)
    out[[length(out) + 1L]] <- center_summary(denominator, centers$state, ring, s, "Non-epithelial denominator")
    if (ncol(non_epithelial)) {
      conditional <- non_epithelial / ifelse(share > 0, share, NA_real_)
      out[[length(out) + 1L]] <- center_summary(conditional, centers$state, ring, s, "Non-epithelial composition")
    }
  }
  
  # 同一邻点可进入多个中心的环，记录两组使用的唯一邻点和交集。
  for (j in names(graph)) {
    a <- unique(unlist(graph[[j]][q$mapping_label %in% "Infected"], use.names = FALSE))
    b <- unique(unlist(graph[[j]][q$mapping_label %in% "Bystander"], use.names = FALSE))
    out[[length(out) + 1L]] <- data.frame(analysis = "Neighborhood overlap", s, ring = as.integer(j),
                                          metric = c("Shared unique neighbors", "Unique Infected neighbors", "Unique Bystander neighbors"),
                                          value = c(length(intersect(a, b)), length(a), length(b)))
  }
  plot_centers[[id]] <- bind_rows(lapply(names(mm), function(j) {
    data.frame(centers, mm[[j]], check.names = FALSE) %>%
      pivot_longer(-c(spot, state), names_to = "metric", values_to = "value") %>%
      mutate(ring = as.integer(j))
  }))
}


##### 精简输出 #####
summary_table <- bind_rows(out)
outdir <- path.expand(outdir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE, na = "")
info <- modifyList(info, list(run_time = Sys.time(), reuse_details_rds = reuse_details_rds,
                              mapping_files = mapping_files, mapping_object = mapping_object, mapping_label_col = mapping_label_col,
                              rings = rings, min_neighbors = min_neighbors, epithelial_types = epithelial_types,
                              plot_sc_violin = plot_sc_violin, sc_plot_file = if (is.null(sc_plot)) "" else sc_file,
                              session = capture.output(sessionInfo())))
info[c("thresholds", "main_threshold", "min_pairs", "tumor_min", "calipers", "tumor_labels", "mapping_method",
       "distribution_min_n", "violin_min_distinct", "violin_min_iqr_fraction", "plot_seed")] <- NULL
if (is.data.frame(info$sample_table)) info$sample_table$region_csv <- NULL
saveRDS(list(spatial = details, info = info), file.path(outdir, "details.rds"), compress = "gzip")

##### 出图：邻域组成优先，MDK与NCL单独展示 #####
status_colors <- c(Infected = "#C35B3F", Bystander = "#426B91")
cell_colors <- c(Epithelial = "#426B91", Fibroblast = "#91C482", Myeloid = "#C35B3F",
                 Tcell = "#785E96", Mast = "#A17418", Endothelial = "#397C73",
                 Bcell = "#739ABD", NK = "#6D6F75")
cell_names <- c(Tcell = "T cell", Bcell = "B cell", Mast = "Mast cell", NK = "NK cell")
cell_label <- function(x) {
  x <- as.character(x)
  use <- x %in% names(cell_names)
  x[use] <- cell_names[x[use]]
  x
}
type_order <- function(x) unique(x[order(match(x, names(cell_colors)), x)])
type_palette <- function(x) {
  extra <- setdiff(x, names(cell_colors))
  c(cell_colors, setNames(grDevices::hcl.colors(length(extra), "Dark 3"), extra))
}
ring_label <- function(x) factor(x, levels = 0:2, labels = c("Center", "Ring 1", "Ring 2"))
chart_theme <- theme_classic(base_size = 10.5, base_family = "Helvetica") +
  theme(axis.line = element_blank(),
        panel.border = element_rect(color = "#526170", fill = NA, linewidth = 0.35),
        strip.background = element_rect(fill = "#EDF2F6", color = "#526170", linewidth = 0.35),
        strip.text = element_text(face = "bold", color = "#22313F", margin = margin(4, 6, 4, 6)),
        panel.spacing = grid::unit(3, "mm"), axis.title = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 11, margin = margin(b = 9)),
        plot.margin = margin(6, 15, 12, 17))
report_page <- function(body, title) {
  header <- grid::textGrob(title, x = grid::unit(17, "pt"), y = 0.54, just = "left",
                           gp = grid::gpar(fontfamily = "Helvetica", fontface = "bold",
                                           fontsize = 18, col = "#22313F"))
  wrap_plots(wrap_elements(full = header), body, ncol = 1, heights = c(0.65, 6.85))
}

plot_composition <- function(comp, y_title) {
  types <- type_order(comp$cell_type)
  means <- pivot_longer(comp, c(Infected, Bystander), names_to = "state", values_to = "fraction") %>%
    mutate(state = factor(state, levels = names(status_colors)),
           cell_type = factor(cell_type, levels = types))
  bar_position <- position_dodge(width = 0.67, preserve = "single")
  p1 <- ggplot(means, aes(ring_label, 100 * fraction, fill = state)) +
    geom_col(width = 0.62, position = bar_position, alpha = 0.78,
             colour = "#526170", linewidth = 0.35, na.rm = TRUE) +
    facet_grid(cell_type ~ ., scales = "free_y", switch = "y", labeller = labeller(cell_type = cell_label)) +
    scale_fill_manual(values = status_colors, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = "Cell fractions", x = NULL, y = y_title) +
    chart_theme + theme(legend.position = "bottom", strip.placement = "outside")
  if (any(means$fraction == 0, na.rm = TRUE)) {
    p1 <- p1 + geom_text(aes(label = ifelse(fraction == 0, "0", "")),
                         position = bar_position, vjust = -0.4, size = 2.6,
                         colour = "#22313F", show.legend = FALSE, na.rm = TRUE)
  }
  comp$cell_type <- factor(comp$cell_type, levels = rev(types))
  delta <- 100 * comp$difference
  high <- max(1, abs(delta), na.rm = TRUE)
  delta_labels <- filter(comp, is.finite(difference))
  p2 <- ggplot(comp, aes(ring_label, cell_type, fill = 100 * difference)) +
    geom_tile(colour = "white", linewidth = 1.1) +
    scale_fill_gradient2(low = "#426B91", mid = "#F7F7F5", high = "#C35B3F", midpoint = 0, na.value = "#B8BDC3",
                         limits = c(-high, high), name = "Infected -\nBystander (pp)",
                         guide = if (nrow(delta_labels)) {
                           guide_colorbar(barheight = grid::unit(32, "mm"), barwidth = grid::unit(3, "mm"))
                         } else "none") +
    scale_y_discrete(labels = cell_label, expand = expansion(add = 0)) +
    scale_x_discrete(expand = expansion(add = 0)) +
    labs(title = "Composition differences", x = NULL, y = NULL) + chart_theme +
    theme(axis.ticks = element_blank(), legend.title = element_text(size = 9))
  if (nrow(delta_labels)) {
    p2 <- p2 + geom_text(data = delta_labels,
                         aes(label = sprintf("%+.1f", 100 * difference),
                             colour = abs(100 * difference) > 0.62 * high),
                         size = 3.45, show.legend = FALSE) +
      scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "#22313F"))
  }
  wrap_plots(p1, p2, widths = c(1.15, 1))
}

plot_center_values <- function(z, facets, colors = status_colors) {
  z <- z %>%
    group_by(across(all_of(c(facets, "state")))) %>%
    mutate(mark = case_when(
      n() == 1L ~ "bar",
      n_distinct(value) == 1L ~ "box",
      TRUE ~ "violin"
    )) %>% ungroup()
  groups <- split(z, z$mark)
  p <- ggplot(z, aes(state, value)) +
    scale_fill_manual(values = colors, guide = "none") +
    scale_x_discrete(drop = FALSE) + chart_theme + labs(x = NULL)
  if (!is.null(groups$bar)) {
    p <- p + geom_col(data = groups$bar, aes(fill = state), width = 0.56, alpha = 0.78,
                      colour = "#526170", linewidth = 0.35)
    if (any(groups$bar$value == 0)) {
      p <- p + geom_text(data = filter(groups$bar, value == 0), label = "0",
                         vjust = -0.4, size = 2.6, colour = "#22313F")
    }
  }
  if (!is.null(groups$violin)) {
    p <- p + geom_violin(data = groups$violin, aes(fill = state), width = 0.72, trim = TRUE,
                         scale = "width", alpha = 0.76, colour = "#526170", linewidth = 0.35) +
      stat_summary(data = groups$violin, fun = median, geom = "point",
                   shape = 95, size = 7, colour = "#22313F")
  }
  if (!is.null(groups$box)) {
    p <- p + geom_boxplot(data = groups$box, aes(fill = state), width = 0.45,
                          alpha = 0.78, colour = "#526170", linewidth = 0.35,
                          outlier.size = 0.9, outlier.alpha = 0.6)
  }
  p
}

plot_cell_centers <- function(z) {
  z$state <- factor(z$state, levels = names(status_colors))
  z$ring_label <- ring_label(z$ring)
  types <- type_order(z$cell_type)
  z$cell_type <- factor(z$cell_type, levels = types, labels = cell_label(types))
  z$value <- 100 * z$value
  plot_center_values(z, c("cell_type", "ring_label")) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.10))) +
    facet_grid(cell_type ~ ring_label, scales = "free_y") +
    labs(x = NULL, y = "Estimated proportion (%)") + chart_theme +
    theme(axis.text.x = element_text(size = 8.5))
}

plot_mdk_centers <- function(m) {
  m <- filter(m, is.finite(value))
  m$state <- factor(m$state, levels = names(status_colors))
  expr_order <- c("Center MDK", "Ring 1 NCL", "Ring 2 NCL")
  m$panel <- ifelse(m$metric == "Center_MDK", "Center MDK", paste("Ring", m$ring, "NCL"))
  expr <- filter(m, metric %in% c("Center_MDK", "Neighbor_NCL"))
  expr$panel <- factor(expr$panel, levels = expr_order)
  codetect <- filter(m, metric %in% c("MDK_center_NCL_neighbor", "NCL_center_MDK_neighbor"))
  codetect$panel <- ring_label(codetect$ring)
  codetect$direction <- factor(codetect$metric, levels = c("MDK_center_NCL_neighbor", "NCL_center_MDK_neighbor"),
                               labels = c("MDK center / NCL neighbors", "NCL center / MDK neighbors"))
  p1 <- plot_center_values(expr, "panel") + facet_wrap(~panel, nrow = 1, scales = "free_y") +
    labs(title = "Expression", y = "log1p(CP10K)")
  codetect$value <- 100 * codetect$value
  p2 <- plot_center_values(codetect, c("panel", "direction")) + facet_grid(panel ~ direction) +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100)) +
    labs(title = "Spatial co-detection", y = "Co-detected neighbors (%)")
  wrap_plots(p1, p2, ncol = 1, heights = c(0.88, 1.2))
}

plot_single_cell <- function(d) {
  z <- d %>%
    pivot_longer(c(MDK_logCP10K, NCL_logCP10K), names_to = "gene", values_to = "value") %>%
    filter(is.finite(value)) %>%
    mutate(gene = sub("_logCP10K$", "", gene))
  if (!nrow(z)) return(NULL)
  colors <- status_colors
  if (!all(z$state %in% names(status_colors))) {
    treatments <- unique(unname(sc_groups))
    colors <- setNames(rep("#426B91", length(treatments)), treatments)
    colors[virus_group] <- "#C35B3F"
  }
  z$state <- factor(z$state, levels = names(colors))
  z$sample_id <- factor(z$sample_id, levels = names(sc_groups))
  plot_center_values(z, c("gene", "sample_id"), colors) +
    facet_grid(gene ~ sample_id, scales = "free_y", switch = "y") +
    labs(y = "log1p(CP10K)") + theme(strip.placement = "outside")
}

##### 空间定位与表达图 #####
spatial_map_theme <- function() {
  theme_void(base_size = 10, base_family = "Helvetica") +
    theme(
      plot.title = element_text(face = "bold", size = 10, margin = margin(0, 0, 8, 0)),
      plot.margin = margin(4, 9, 4, 9),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      strip.background = element_rect(fill = "#E6EDF3", colour = NA),
      strip.text = element_text(face = "bold", size = 10, margin = margin(4, 4, 4, 4))
    )
}

spatial_pitch <- function(detail) {
  s <- detail$spots
  edge <- detail$rings[["1"]]
  dx <- s$x[match(unlist(edge, use.names = FALSE), s$spot)] -
    s$x[match(rep(names(edge), lengths(edge)), s$spot)]
  dy <- s$y[match(unlist(edge, use.names = FALSE), s$spot)] -
    s$y[match(rep(names(edge), lengths(edge)), s$spot)]
  distances <- sqrt(dx^2 + dy^2)
  if (any(is.finite(distances) & distances > 0)) {
    return(median(distances[is.finite(distances) & distances > 0]))
  }
  max(diff(range(s$x)), diff(range(s$y)), 1) / max(1, sqrt(nrow(s)))
}

plot_anchor_locations <- function(detail) {
  s <- detail$spots
  if (!nrow(s)) return(NULL)
  s$y <- -s$y
  pitch <- spatial_pitch(detail)
  groups <- c("Infected", "Bystander")
  s$anchor <- ifelse(s$QC_pass %in% TRUE & s$mapping_label %in% groups,
                     as.character(s$mapping_label), NA_character_)
  
  # 每组独立选择具有有效一环组成且最接近本组坐标中位数的示例中心。
  complete_spots <- rownames(detail$fractions)[complete.cases(detail$fractions)]
  selected <- vapply(groups, function(state) {
    centers <- s[which(s$anchor == state), , drop = FALSE]
    usable <- vapply(detail$rings[["1"]][centers$spot], function(ids) {
      sum(ids %in% complete_spots) >= min_neighbors
    }, logical(1))
    centers <- centers[usable, , drop = FALSE]
    if (!nrow(centers)) return(NA_character_)
    distance <- (centers$x - median(centers$x))^2 + (centers$y - median(centers$y))^2
    centers$spot[order(distance, centers$spot)[1L]]
  }, character(1))
  selected <- selected[!is.na(selected)]
  
  # 各局部图独立按窗口比例缩放，横纵使用同一个比例，分面统一对齐。
  panels <- c(unique(s$treatment)[1L], names(selected))
  window <- function(x, y) {
    span <- max(diff(x), diff(y), pitch)
    list(x = mean(x), y = mean(y), span = span)
  }
  orient <- function(d, w, panel) {
    d$plot_x <- (d$x - w$x) / w$span
    d$plot_y <- (d$y - w$y) / w$span
    d$panel <- rep(panel, nrow(d))
    d
  }
  mark_rings <- function(d, ids) {
    r1 <- unique(unlist(detail$rings[["1"]][ids], use.names = FALSE))
    r2 <- unique(unlist(detail$rings[["2"]][ids], use.names = FALSE))
    d$display <- ifelse(d$spot %in% r1, "Ring 1", ifelse(d$spot %in% r2, "Ring 2", NA_character_))
    d$display[d$spot %in% ids] <- NA_character_
    d
  }
  overview <- window(range(s$x) + c(-0.2, 0.2) * pitch, range(s$y) + c(-0.2, 0.2) * pitch)
  windows <- lapply(selected, function(id) {
    include <- unique(c(id, detail$rings[["1"]][[id]], detail$rings[["2"]][[id]]))
    d <- s[s$spot %in% include, ]
    list(x = range(d$x) + c(-0.65, 0.65) * pitch, y = range(d$y) + c(-1.2, 0.65) * pitch)
  })
  full <- orient(mark_rings(s, unname(selected)), overview, panels[1L])
  full$extra <- NA_character_
  crops <- bind_rows(lapply(names(selected), function(panel) {
    limits <- windows[[panel]]
    id <- selected[[panel]]
    # 局部图下方留出比例尺位置，邻域各点仍保留在绘图区内。
    d <- s[s$x >= limits$x[1L] & s$x <= limits$x[2L] &
             s$y >= limits$y[1L] + 0.55 * pitch & s$y <= limits$y[2L], ]
    d <- orient(mark_rings(d, id), window(limits$x, limits$y), panel)
    d$extra <- ifelse(d$spot != id, d$anchor, NA_character_)
    d$anchor[d$spot != id] <- NA_character_
    d
  }))
  maps <- bind_rows(full, crops)
  maps$panel <- factor(maps$panel, levels = panels)
  maps$zoom <- maps$panel != panels[1L]
  boxes <- bind_rows(lapply(windows, function(w) {
    data.frame(panel = panels[1L], xmin = (w$x[1L] - overview$x) / overview$span,
               xmax = (w$x[2L] - overview$x) / overview$span,
               ymin = (w$y[1L] - overview$y) / overview$span, ymax = (w$y[2L] - overview$y) / overview$span)
  }))
  bars <- bind_rows(lapply(names(windows), function(panel) {
    lim <- windows[[panel]]
    w <- window(lim$x, lim$y)
    # 传统Visium相邻spot中心相距100微米。
    orient(data.frame(x = lim$x[1L] + 0.35 * pitch, y = lim$y[1L] + 0.30 * pitch), w, panel) %>%
      mutate(end_x = plot_x + 2 * pitch / w$span, text_y = plot_y + 0.12 * pitch / w$span)
  }))
  cols <- c(Infected = "#C35B3F", Bystander = "#FFFFFF", "Ring 1" = "#92B7CC", "Ring 2" = "#DCE9F2")
  key <- which(names(cols) %in% c(maps$anchor, maps$display))
  map_legend <- if (length(key)) guide_legend(override.aes = list(shape = 21,
                                                                  colour = c("#78382B", "#426B91", "#92B7CC", "#C5D8E6")[key], size = 3.1, stroke = 0.65)) else "none"
  layout <- ggplot(maps, aes(plot_x, plot_y)) +
    geom_point(data = filter(maps, !zoom), size = 0.82, colour = "#D7DBDF") +
    geom_point(data = filter(maps, zoom), shape = 21, size = 2.3, colour = "#D1D5DA", fill = "#F1F3F4", stroke = 0.22) +
    geom_point(data = filter(maps, !zoom, !is.na(display)), aes(fill = display), shape = 21, size = 1.25,
               colour = "#7C98A9", stroke = 0.15) +
    geom_point(data = filter(maps, zoom, !is.na(display)), aes(fill = display), shape = 21, size = 3.25,
               colour = "#7291A5", stroke = 0.28) +
    geom_point(data = filter(maps, !is.na(extra)), aes(colour = extra), shape = 16, size = 1.0) +
    geom_point(data = filter(maps, !is.na(anchor), !zoom),
               aes(fill = anchor, colour = anchor), shape = 21, size = 1.9, stroke = 0.45) +
    geom_point(data = filter(maps, !is.na(anchor), zoom),
               aes(fill = anchor, colour = anchor), shape = 21, size = 4.0, stroke = 0.75) +
    scale_colour_manual(values = c(Infected = "#78382B", Bystander = "#426B91"), guide = "none") +
    scale_fill_manual(values = cols, breaks = names(cols)[key], name = NULL, guide = map_legend) +
    facet_wrap(~panel, nrow = 1) + coord_equal(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5), expand = FALSE) +
    spatial_map_theme() + theme(panel.spacing.x = grid::unit(5, "mm"))
  if (nrow(boxes)) {
    boxes$panel <- factor(boxes$panel, levels = panels)
    bars$panel <- factor(bars$panel, levels = panels)
    layout <- layout +
      geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = NA, colour = "#394A57", linetype = "dashed", linewidth = 0.65, inherit.aes = FALSE) +
      geom_segment(data = bars, aes(xend = end_x, yend = plot_y), colour = "#394A57", linewidth = 0.48) +
      geom_text(data = bars, aes(x = (plot_x + end_x) / 2, y = text_y), label = "200 \u00b5m",
                colour = "#394A57", size = 2.5, vjust = 0, family = "Helvetica")
  }
  
  examples <- bind_rows(lapply(names(selected), function(state) {
    ids <- detail$rings[["1"]][[selected[[state]]]]
    values <- detail$fractions[intersect(ids, complete_spots), , drop = FALSE]
    data.frame(state = state, cell_type = colnames(values), fraction = colMeans(values))
  }))
  if (nrow(examples)) {
    examples$state <- factor(examples$state, levels = rev(groups))
    examples$cell_type <- factor(examples$cell_type, levels = type_order(examples$cell_type))
    example_plot <- ggplot(examples, aes(100 * fraction, state, fill = cell_type)) +
      geom_col(width = 0.63, colour = "white", linewidth = 0.25, orientation = "y") +
      geom_text(data = distinct(examples, state), aes(x = -2.5, y = state, label = state, colour = state),
                inherit.aes = FALSE, hjust = 1, fontface = "bold", size = 3) +
      coord_cartesian(xlim = c(0, 100), clip = "off") +
      scale_x_continuous(breaks = seq(0, 100, 25), expand = expansion(mult = 0)) +
      scale_y_discrete(labels = NULL, expand = expansion(add = 0.55)) +
      scale_colour_manual(values = status_colors, guide = "none") +
      scale_fill_manual(values = type_palette(levels(examples$cell_type)), labels = cell_label, name = NULL) +
      labs(title = "Local first-ring composition", x = "Estimated proportion (%)", y = NULL) +
      chart_theme + theme(axis.ticks.y = element_blank(), plot.margin = margin(6, 15, 12, 75),
                          legend.key.size = grid::unit(3.5, "mm"), legend.text = element_text(size = 8.5))
    layout <- patchwork::wrap_plots(layout, example_plot, ncol = 1, heights = c(2.6, 1))
  }
  report_page(layout, "Spatial neighborhoods")
}


plot_expression_maps <- function(details) {
  d <- lapply(details, function(z) z$spots)
  if (!length(d) || !any(lengths(lapply(d, function(s) s$spot)) > 0L)) return(NULL)
  d <- bind_rows(d)
  id <- unique(d[, c("sample_id", "treatment"), drop = FALSE])
  id$panel <- ifelse(duplicated(id$treatment) | duplicated(id$treatment, fromLast = TRUE),
                     id$sample_id, id$treatment)
  d$panel <- factor(id$panel[match(d$sample_id, id$sample_id)], levels = id$panel)
  # Recenter the tissue footprints while retaining the original orientation and pixel aspect.
  d$plot_x <- d$x - ave(d$x, d$sample_id, FUN = function(x) mean(range(x)))
  d$plot_y <- -d$y + ave(d$y, d$sample_id, FUN = function(x) mean(range(x)))
  map <- pivot_longer(d, c(MDK_logCP10K, NCL_logCP10K), names_to = "gene", values_to = "expression") %>%
    mutate(gene = sub("_logCP10K$", "", gene), expression = ifelse(QC_pass, expression, NA_real_))
  high <- max(c(1, map$expression[is.finite(map$expression)]))
  p <- ggplot(map, aes(plot_x, plot_y, colour = expression)) + geom_point(shape = 16, size = 1.18) +
    scale_colour_gradientn(
      colours = c("#F3F7FB", "#91B3CB", "#2D5375"), na.value = "#B8BDC3",
      limits = c(0, high), name = "log1p(CP10K)",
      guide = guide_colorbar(barheight = grid::unit(32, "mm"), barwidth = grid::unit(3, "mm"))
    ) +
    facet_grid(gene ~ panel, switch = "y") + coord_equal(expand = TRUE) + spatial_map_theme() +
    theme(panel.spacing = grid::unit(5, "mm"), strip.text.y.left = element_text(angle = 0))
  report_page(p, "MDK and NCL expression")
}
pdf_file <- file.path(outdir, "plots.pdf")
if (capabilities("cairo")) {
  grDevices::cairo_pdf(pdf_file, width = 11, height = 7.5, family = "sans")
} else {
  pdf(pdf_file, width = 11, height = 7.5, useDingbats = FALSE)
}
tryCatch({
  ids <- names(plot_centers)
  main <- summary_table %>% filter(analysis %in% c("Composition", "Non-epithelial composition")) %>%
    select(sample_id, analysis, ring, metric, Infected, Bystander, difference) %>%
    mutate(ring_label = ring_label(ring), cell_type = sub("^Fraction[.]", "", metric)) %>%
    filter(!cell_type %in% epithelial_types)
  for (id in ids) {
    d <- details[[id]]
    title_suffix <- if (length(ids) > 1L) paste(" | ", id, sep = "") else ""
    for (section in c("Composition", "Non-epithelial composition")) {
      comp <- filter(main, sample_id == id, analysis == section)
      title <- if (section == "Composition") "Neighborhood composition" else section
      y_title <- if (section == "Composition") "Mean % of total" else "Mean % of non-epithelial"
      if (any(is.finite(comp$Infected) | is.finite(comp$Bystander))) {
        panel <- plot_composition(comp, y_title)
      } else {
        panel <- ggplot() + annotate("text", x = 0, y = 0, label = "No evaluable composition", colour = "#526170") +
          theme_void(base_family = "Helvetica")
      }
      print(report_page(panel, paste0(title, title_suffix)))
    }
    location <- plot_anchor_locations(d)
    if (!is.null(location)) print(location)
    dist <- plot_centers[[id]] %>%
      filter(ring != 0, is.finite(value), startsWith(metric, "Fraction.")) %>%
      mutate(cell_type = sub("^Fraction[.]", "", metric)) %>%
      filter(!cell_type %in% epithelial_types)
    if (nrow(dist)) {
      print(report_page(plot_cell_centers(dist), paste0("Cell-type proportions", title_suffix)))
    }
  }
  expression <- plot_expression_maps(details)
  if (!is.null(expression)) print(expression)
  for (id in ids) {
    m <- plot_centers[[id]]
    if (any(m$metric == "MDK_center_NCL_neighbor" & is.finite(m$value))) {
      title <- if (length(ids) > 1L) paste("MDK-NCL neighborhoods |", id) else "MDK-NCL neighborhoods"
      print(report_page(plot_mdk_centers(m), title))
    }
  }
  if (!is.null(sc_plot)) {
    single_cell <- plot_single_cell(sc_plot)
    if (!is.null(single_cell)) print(report_page(single_cell, "Single-cell MDK and NCL"))
  }
}, finally = dev.off())
message("完成。输出目录：", outdir)

