##### 清理运行环境 #####
rm(list = ls(all.names = TRUE))
graphics.off()
invisible(gc())

##### 加载依赖 #####
library(Seurat)
library(Matrix)
library(UCell)
library(BiocParallel)
library(edgeR)
library(statmod)
library(ggplot2)
library(patchwork)
library(curl)
# nichenetr在第二阶段加载；下载/安装失败不会抹掉第一阶段结果。

##### 路径、样本与分析参数 #####
input_file <- "~/rawdata/scRNA_virus/virus_3D/virus_3D_anno.RData"
input_object <- "sce"
analysis_label <- "D3"
outdir <- file.path(path.expand("~/OV"), "scRNA_virus", paste0("immune_NicheNet_", analysis_label))
sample_col <- "orig.ident"
celltype_col <- "celltype"
condition_levels <- c("Vehicle", "VG161") # 全部变化方向：VG161 - Vehicle。
sample_table <- data.frame(
  sample_id = c("Vehicle_3D_1", "Vehicle_3D_2", "VG161_3D_1", "VG161_3D_2"),
  replicate_id = c("Vehicle_3D_1", "Vehicle_3D_2", "VG161_3D_1", "VG161_3D_2"),
  condition = c("Vehicle", "Vehicle", "VG161", "VG161"), stringsAsFactors = FALSE)
# 默认是两组独立样本，不把编号1/2当作配对。技术拆库才填写同一replicate_id。
host_prefix <- ""
viral_pattern <- "^(VG161|HSV[-_]?1)([-_.]|$)"
viral_feature_ids <- character() # 无前缀病毒/载荷feature按真实参考注释补入。

# 可选：只读取已有细分注释，不重新聚类、不按单个marker猜CD8身份。
tcell_annotation_file <- "~/rawdata/scRNA_virus/virus_3D/virus_3D_tcell.RData"
tcell_annotation_object <- "sce"
receiver_definitions <- list(
  Tcell = list(parent = "Tcell", column = celltype_col, labels = "Tcell", lineage = "T"),
  CD8_T = list(parent = "Tcell", column = "t_type", labels = "CD8+ Tcell", lineage = "T"),
  Myeloid = list(parent = "Myeloid", column = celltype_col, labels = "Myeloid", lineage = "Myeloid"))
tumor_labels <- "Epithelial"
min_cells <- 20L
min_samples_per_condition <- 2L
min_signature_genes <- 5L
min_signature_coverage <- 0.70
ucell_max_rank <- 1500L
expression_fraction <- 0.10
filter_min_count <- 10
filter_min_total_count <- 15
de_fdr <- 0.05
target_min_abs_logFC <- 0.25
top_target_n <- 200L
min_target_genes <- 15L
min_background_genes <- 100L
min_candidate_ligands <- 5L
# 统一从响应目标与其背景排除待检轴端点；不从counts、表达筛选或候选配体中删除。
target_exclusion_genes <- c("MDK", "NCL", "LRP1", "SDC4", "PTPRZ1", "ALK")
target_gene <- "MDK"
run_nichenet <- TRUE
model_dir <- "~/rawdata/scRNA_virus/nichenet_v2"
model_source <- "official" # "local"仅用同版本已有文件，仍校验MD5；不自动替换模型。
download_timeout <- 1800L
download_retries <- 3L
seed <- 123L
top_ligands_plot <- 15L
top_ligands_links <- 10L
top_target_links <- 20L
condition_colors <- c(Vehicle = "#426B91", VG161 = "#C35B3F")

##### 探索性功能基因集：定义与解释见README #####
# 作者整理的小型starter集合，不冒充原文完整signature；不以PDCD1单基因定义耗竭。
signatures <- list(
  Cytotoxicity = c("NKG7", "PRF1", "GNLY", "GZMB", "GZMH", "CTSW", "FASLG"),
  Naive_memory = c("CCR7", "SELL", "TCF7", "LEF1", "IL7R", "LTB", "MAL"),
  Coinhibitory_activation = c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4", "TOX", "ENTPD1"),
  Proliferation = c("MKI67", "TOP2A", "TYMS", "PCNA", "STMN1", "UBE2C", "MCM5"),
  IFN_response = c("STAT1", "IRF1", "ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1"),
  MHC_I_processing = c("B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9"),
  MHC_II_presentation = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB", "CD74", "CIITA"),
  Inflammatory_response = c("IL1B", "TNF", "NFKBIA", "CXCL8", "CCL3", "CCL4", "PTGS2", "IL6"),
  Phagolysosomal_complement = c("C1QA", "C1QB", "C1QC", "CTSB", "CTSD", "LAMP1", "LAMP2", "FCER1G"))
signature_lineages <- list(T = c("Cytotoxicity", "Naive_memory", "Coinhibitory_activation", "Proliferation", "IFN_response", "MHC_I_processing"),
                           Myeloid = c("IFN_response", "MHC_I_processing", "MHC_II_presentation", "Inflammatory_response", "Phagolysosomal_complement"))

##### 通用函数与输入校验 #####
read_sce <- function(path, object_name) {
  path <- path.expand(path)
  if (!file.exists(path)) stop("找不到对象文件：", path)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  env <- new.env()
  nms <- load(path, envir = env)
  if (!object_name %in% nms) stop("文件中缺少对象：", object_name)
  get(object_name, envir = env)
}
bind_rows_base <- function(x) {
  x <- Filter(function(z) is.data.frame(z) && nrow(z) > 0L, x)
  if (length(x)) do.call(rbind, x) else data.frame()
}
is_viral_feature <- function(f) {
  symbol <- f
  if (nzchar(host_prefix)) symbol[startsWith(f, host_prefix)] <- substring(f[startsWith(f, host_prefix)], nchar(host_prefix) + 1L)
  hit <- if (nzchar(viral_pattern)) grepl(viral_pattern, f, ignore.case = TRUE) | grepl(viral_pattern, symbol, ignore.case = TRUE) else rep(FALSE, length(f))
  hit | f %in% viral_feature_ids | symbol %in% viral_feature_ids
}
required <- c("sample_id", "replicate_id", "condition")
if (!all(required %in% names(sample_table)) || !nrow(sample_table) || anyNA(sample_table[, required]) ||
    any(!nzchar(as.matrix(sample_table[, required]))) || anyDuplicated(sample_table$sample_id)) stop("样本表缺失、空值或sample_id重复。")
if (!identical(condition_levels, c("Vehicle", "VG161")) || !all(sample_table$condition %in% condition_levels)) stop("请保留Vehicle、VG161方向，并明确映射样本组。")
if (any(vapply(split(sample_table$condition, sample_table$replicate_id), function(z) length(unique(z)) != 1L, logical(1)))) stop("同一replicate_id跨条件：当前脚本是独立组设计，不支持重复测量。")
if (min_samples_per_condition < 2L || min_cells < 1L || min_target_genes < 2L || top_target_n < min_target_genes ||
    expression_fraction <= 0 || expression_fraction > 1 || min_signature_coverage <= 0 || min_signature_coverage > 1 ||
    de_fdr <= 0 || de_fdr >= 1 || target_min_abs_logFC < 0) stop("请核对细胞数、样本数、覆盖、目标大小和阈值。")
if (!model_source %in% c("official", "local")) stop("model_source应为official或local。")
set.seed(seed)
started_at <- Sys.time()
outdir <- path.expand(outdir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
sce <- read_sce(input_file, input_object)
if (!inherits(sce, "Seurat") || !"RNA" %in% names(sce@assays)) stop("输入须为含RNA原始counts的Seurat对象。")
meta <- sce[[]]
if (!all(c(sample_col, celltype_col) %in% names(meta)) || anyDuplicated(rownames(meta))) stop("metadata必要字段缺失或条码重复。")
if (!all(sample_table$sample_id %in% as.character(meta[[sample_col]]))) stop("样本表与对象不符；对象实际样本：", paste(unique(meta[[sample_col]]), collapse = ", "))
meta <- meta[!is.na(meta[[sample_col]]) & meta[[sample_col]] %in% sample_table$sample_id, , drop = FALSE]
meta$sample_id <- as.character(meta[[sample_col]])
meta$replicate_id <- sample_table$replicate_id[match(meta$sample_id, sample_table$sample_id)]
meta$condition <- factor(sample_table$condition[match(meta$sample_id, sample_table$sample_id)], levels = condition_levels)
annotation_audit <- list(status = "not_used", file = path.expand(tcell_annotation_file))
if (nzchar(tcell_annotation_file) && file.exists(path.expand(tcell_annotation_file))) {
  t_sce <- read_sce(tcell_annotation_file, tcell_annotation_object)
  if (!inherits(t_sce, "Seurat")) stop("T细分对象不是Seurat。")
  tm <- t_sce[[]]
  if (anyDuplicated(rownames(tm))) stop("T细分对象条码重复。")
  matched <- intersect(rownames(tm), rownames(meta)[meta[[celltype_col]] %in% "Tcell"])
  if (!length(matched)) stop("现有T细分对象与本次Tcell没有共同条码；不能跨D3/D5或改猜条码。")
  if (!sample_col %in% names(tm) || any(as.character(tm[matched, sample_col]) != meta[matched, "sample_id"])) stop("T细分对象的样本身份与主对象不一致。")
  if ("t_type" %in% names(tm)) {
    new_labels <- as.character(tm[matched, "t_type"])
    if ("t_type" %in% names(meta)) {
      old <- as.character(meta[matched, "t_type"])
      if (any(!is.na(old) & nzchar(old) & !is.na(new_labels) & old != new_labels)) stop("主对象和T细分对象的t_type冲突。")
    } else meta$t_type <- NA_character_
    meta[matched, "t_type"] <- new_labels
    annotation_audit$status <- "joined_existing_t_type"
    annotation_audit$n_matched <- length(matched)
    annotation_audit$labels <- table(new_labels, useNA = "ifany")
  } else annotation_audit$status <- "file_has_no_t_type_broad_T_only"
  rm(t_sce, tm)
} else annotation_audit$status <- "file_absent_broad_T_only"
message("T细分注释：", annotation_audit$status, "；不重新判定CD8或DC身份。")
print(table(meta$sample_id, meta[[celltype_col]], useNA = "ifany"))

##### 提取原始宿主counts：病毒排除与共同测量背景 #####
rna <- sce[["RNA"]]
if (inherits(rna, "Assay5")) {
  layer_names <- Layers(rna, search = "^counts($|\\.)")
  if (!length(layer_names)) stop("RNA缺少counts层，不能用data/SCT/integrated替代。")
  cells_by_layer <- lapply(layer_names, function(k) colnames(LayerData(rna, layer = k, fast = FALSE)))
  if (anyDuplicated(unlist(cells_by_layer, use.names = FALSE))) stop("counts层条码重复。")
  if (!all(rownames(meta) %in% unlist(cells_by_layer, use.names = FALSE))) stop("纳入细胞缺少原始counts。")
  layer_names <- layer_names[vapply(cells_by_layer, function(z) any(z %in% rownames(meta)), logical(1))]
  count_layers <- lapply(layer_names, function(k) {
    z <- LayerData(rna, layer = k, fast = FALSE)
    z[, intersect(colnames(z), rownames(meta)), drop = FALSE]
  })
} else if (inherits(rna, "Assay")) {
  z <- if (packageVersion("SeuratObject") >= "5.0.0") GetAssayData(sce, assay = "RNA", layer = "counts") else GetAssayData(sce, assay = "RNA", slot = "counts")
  if (!all(rownames(meta) %in% colnames(z))) stop("纳入细胞缺少原始counts。")
  layer_names <- "counts"
  count_layers <- list(z[, rownames(meta), drop = FALSE])
} else stop("不支持的RNA assay。")
names(count_layers) <- layer_names
for (k in layer_names) {
  z <- count_layers[[k]]
  if (anyDuplicated(rownames(z)) || anyDuplicated(colnames(z)) || is.null(rownames(z)) || is.null(colnames(z))) stop(k, " feature/条码不唯一。")
  if (!inherits(z, "sparseMatrix")) z <- Matrix(z, sparse = TRUE)
  z <- as(z, "dgCMatrix")
  if (any(!is.finite(z@x) | z@x < 0 | abs(z@x - round(z@x)) > 1e-6)) stop("需要非负有限整数counts。")
  count_layers[[k]] <- z
}
viral_exclusion_audit <- bind_rows_base(lapply(layer_names, function(k) {
  z <- count_layers[[k]]; ids <- rownames(z)[is_viral_feature(rownames(z))]
  data.frame(layer = rep(k, length(ids)), feature = ids, total_UMI = as.numeric(Matrix::rowSums(z[ids, , drop = FALSE])))
}))
host_features <- lapply(count_layers, function(z) {
  f <- rownames(z)
  if (nzchar(host_prefix)) f <- f[startsWith(f, host_prefix)]
  f[!is_viral_feature(f)]
})
common_features <- Reduce(intersect, host_features)
symbols <- if (nzchar(host_prefix)) substring(common_features, nchar(host_prefix) + 1L) else common_features
if (length(symbols) < min_background_genes || anyDuplicated(symbols) || anyNA(symbols) || any(!nzchar(symbols)) ||
    mean(grepl("^ENSG[0-9]", symbols)) > 0.5) stop("共同宿主symbol过少/重复，或输入主要为Ensembl ID；请先完成确定映射与物种分离。")
layer_audit <- data.frame(layer = layer_names, n_cells = vapply(count_layers, ncol, integer(1)),
                          n_host_features = lengths(host_features), n_common = length(symbols))
if (any(lengths(host_features) != length(symbols))) warning("不同counts层使用宿主feature交集；未测到的feature不补零。")
host_counts <- do.call(cbind, lapply(count_layers, function(z) z[common_features, , drop = FALSE]))
rownames(host_counts) <- symbols
meta <- meta[colnames(host_counts), , drop = FALSE]
host_umi <- Matrix::colSums(host_counts)
cell_qc <- aggregate(cbind(n_selected = rep(1L, nrow(meta)), n_valid_host = as.integer(host_umi > 0)),
                     by = meta[, c("sample_id", "replicate_id", "condition"), drop = FALSE], FUN = sum)
meta <- meta[host_umi > 0, , drop = FALSE]
host_counts <- host_counts[, rownames(meta), drop = FALSE]
rm(sce, rna, count_layers, host_umi, z)
invisible(gc())
replicate_table <- unique(sample_table[, c("replicate_id", "condition")])
replicate_table$condition <- factor(replicate_table$condition, levels = condition_levels)
row.names(replicate_table) <- replicate_table$replicate_id
receiver_cells <- list()
receiver_audit <- list()
for (receiver in names(receiver_definitions)) {
  spec <- receiver_definitions[[receiver]]
  present <- spec$column %in% names(meta)
  cells <- if (present) rownames(meta)[!is.na(meta[[spec$column]]) & meta[[spec$column]] %in% spec$labels & meta[[celltype_col]] %in% spec$parent] else character()
  receiver_cells[[receiver]] <- cells
  n <- vapply(replicate_table$replicate_id, function(id) sum(meta[cells, "replicate_id"] == id), integer(1))
  receiver_audit[[receiver]] <- data.frame(receiver = receiver, replicate_table, n_cells = n,
                                           passes_qc = n >= min_cells, annotation_status = if (present) "existing_label" else "annotation_column_absent", row.names = NULL)
}
receiver_audit <- do.call(rbind, receiver_audit)

##### 第一阶段：逐细胞功能评分，按独立样本汇总 #####
score_rows <- list(); coverage_rows <- list(); signature_expression_rows <- list()
for (receiver in names(receiver_cells)) {
  cells <- receiver_cells[[receiver]]
  spec <- receiver_definitions[[receiver]]
  nms <- signature_lineages[[spec$lineage]]
  sg <- lapply(signatures[nms], function(g) setdiff(unique(g), target_exclusion_genes))
  measured <- lapply(sg, intersect, y = rownames(host_counts))
  coverage <- data.frame(receiver = receiver, signature = names(sg), n_defined = lengths(sg),
                         n_measured = lengths(measured), coverage = lengths(measured) / lengths(sg), stringsAsFactors = FALSE)
  coverage$assessable <- coverage$n_measured >= min_signature_genes & coverage$coverage >= min_signature_coverage
  coverage$missing_genes <- vapply(names(sg), function(nm) paste(setdiff(sg[[nm]], measured[[nm]]), collapse = ";"), character(1))
  coverage_rows[[receiver]] <- coverage
  ok <- coverage$signature[coverage$assessable]
  scores <- NULL
  if (length(cells) && length(ok)) {
    scores <- UCell::ScoreSignatures_UCell(host_counts[, cells, drop = FALSE], features = measured[ok],
                                           maxRank = min(ucell_max_rank, nrow(host_counts)), ncores = 1, BPPARAM = SerialParam())
    if (!all(cells %in% rownames(scores))) stop("UCell输出细胞条码不完整。")
  }
  for (id in replicate_table$replicate_id) for (nm in nms) {
    cc <- cells[meta[cells, "replicate_id"] == id]
    enough <- length(cc) >= min_cells
    assessable <- nm %in% ok
    val <- if (enough && assessable) scores[cc, paste0(nm, "_UCell")] else NA_real_
    score_rows[[length(score_rows) + 1L]] <- data.frame(receiver = receiver, replicate_id = id,
                                                        condition = as.character(replicate_table[id, "condition"]), signature = nm, n_cells = length(cc),
                                                        mean_score = if (enough && assessable) mean(val) else NA_real_,
                                                        median_score = if (enough && assessable) median(val) else NA_real_,
                                                        status = if (!enough) "insufficient_cells" else if (!assessable) "insufficient_signature_coverage" else "assessed")
    for (g in sg[[nm]]) {
      present <- g %in% rownames(host_counts)
      signature_expression_rows[[length(signature_expression_rows) + 1L]] <- data.frame(receiver = receiver,
                                                                                        replicate_id = id, signature = nm, gene = g, n_cells = length(cc), feature_measured = present,
                                                                                        total_UMI = if (present && length(cc)) sum(host_counts[g, cc]) else NA_real_,
                                                                                        detected_fraction = if (present && length(cc)) mean(host_counts[g, cc] > 0) else NA_real_)
    }
  }
}
scores_summary <- do.call(rbind, score_rows)
signature_coverage <- do.call(rbind, coverage_rows)
signature_expression <- bind_rows_base(signature_expression_rows)
score_contrasts <- bind_rows_base(lapply(split(scores_summary, interaction(scores_summary$receiver, scores_summary$signature, drop = TRUE)), function(z) {
  a <- z$mean_score[z$condition == "Vehicle" & is.finite(z$mean_score)]
  b <- z$mean_score[z$condition == "VG161" & is.finite(z$mean_score)]
  enough <- length(a) >= min_samples_per_condition && length(b) >= min_samples_per_condition
  data.frame(receiver = z$receiver[1], signature = z$signature[1], n_Vehicle = length(a), n_VG161 = length(b),
             mean_Vehicle = if (length(a)) mean(a) else NA_real_, mean_VG161 = if (length(b)) mean(b) else NA_real_,
             difference_VG161_minus_Vehicle = if (enough) mean(b) - mean(a) else NA_real_,
             PValue = NA_real_, status = if (enough) "descriptive_no_test" else "insufficient_samples")
}))

##### 接收细胞pseudobulk与独立组edgeR：建立完整响应背景 #####
de_results <- list(); de_models <- list(); de_audit <- list()
for (receiver in names(receiver_cells)) {
  cells <- receiver_cells[[receiver]]
  ids <- replicate_table$replicate_id
  ra <- receiver_audit[receiver_audit$receiver == receiver, ]
  selected_ids <- ids[ra$passes_qc[match(ids, ra$replicate_id)]]
  mm <- replicate_table[selected_ids, , drop = FALSE]
  status <- "assessed"
  if (any(table(factor(mm$condition, levels = condition_levels)) < min_samples_per_condition)) status <- "insufficient_biological_samples"
  if (status != "assessed") {
    de_audit[[receiver]] <- data.frame(receiver = receiver, status = status, n_Vehicle = sum(mm$condition == "Vehicle"), n_VG161 = sum(mm$condition == "VG161"))
    next
  }
  cc <- cells[meta[cells, "replicate_id"] %in% selected_ids]
  membership <- sparseMatrix(i = seq_along(cc), j = match(meta[cc, "replicate_id"], selected_ids), x = 1,
                             dims = c(length(cc), length(selected_ids)))
  pb <- as.matrix(host_counts[, cc, drop = FALSE] %*% membership)
  colnames(pb) <- selected_ids
  design <- model.matrix(~ condition, data = mm)
  if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) stop("接收细胞独立组设计不可估计。")
  y <- DGEList(pb)
  keep <- filterByExpr(y, group = mm$condition, min.count = filter_min_count, min.total.count = filter_min_total_count)
  if (sum(keep) < min_background_genes) {
    de_audit[[receiver]] <- data.frame(receiver = receiver, status = "insufficient_tested_genes", n_Vehicle = sum(mm$condition == "Vehicle"), n_VG161 = sum(mm$condition == "VG161"))
    next
  }
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  tt <- glmQLFTest(fit, coef = "conditionVG161")
  de <- as.data.frame(topTags(tt, n = Inf, sort.by = "none")$table)
  de$gene <- rownames(de); de$receiver <- receiver
  de$rank_stat <- sign(de$logFC) * sqrt(de$F)
  de <- de[is.finite(de$rank_stat) & is.finite(de$logFC), , drop = FALSE]
  de_results[[receiver]] <- de
  de_models[[receiver]] <- list(pseudobulk = pb, samples = mm, design = design,
                                normalization = y$samples, logCPM = cpm(y, log = TRUE, prior.count = 2), tested_genes = de$gene)
  de_audit[[receiver]] <- data.frame(receiver = receiver, status = status, n_Vehicle = sum(mm$condition == "Vehicle"), n_VG161 = sum(mm$condition == "VG161"))
}
de_audit <- do.call(rbind, de_audit)

##### 绘图与阶段保存函数：保留样本点，不添加方法性caption #####
theme_research <- theme_bw(base_size = 10, base_family = "sans") + theme(
  panel.grid = element_blank(), panel.border = element_rect(colour = "#454545", linewidth = 0.4),
  axis.ticks = element_line(colour = "#454545", linewidth = 0.35),
  strip.background = element_rect(fill = "#E7EEF4", colour = "#454545", linewidth = 0.35),
  strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", size = 12),
  legend.title = element_blank())
empty_plot <- function(title, status) ggplot() + annotate("text", x = 0, y = 0, label = status, size = 4) +
  labs(title = title) + theme_void(base_family = "sans")
ligand_results <- data.frame(receiver = character(), source = character(), unit = character(), replicate_id = character(),
                             target_rule = character(), direction = character(), ligand = character(), auroc = double(), aupr = double(), aupr_corrected = double(), pearson = double(),
                             rank = double(), rank_percentile = double(), n_candidates = integer(), n_scored = integer(), n_targets = integer(),
                             n_background = integer(), n_FDR_supported = integer(), score_status = character())
model_audit <- data.frame(); candidate_availability <- data.frame(); target_tables <- list(); target_links <- data.frame()
prior_info <- list(status = "not_loaded"); model_coverage <- data.frame()
draw_plots <- function() {
  pdf(file.path(outdir, "plots.pdf"), width = 12, height = 9, onefile = TRUE, useDingbats = FALSE)
  on.exit(dev.off())
  pp <- receiver_audit
  pp$replicate_id <- factor(pp$replicate_id, levels = replicate_table$replicate_id)
  print(ggplot(pp, aes(replicate_id, n_cells, fill = condition)) + geom_col(width = 0.65, colour = "#454545", linewidth = 0.35) +
          facet_wrap(~ receiver, scales = "free_y") + scale_fill_manual(values = condition_colors) +
          labs(title = "Immune cells", x = NULL, y = "Cells with host counts") + theme_research +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top"))
  for (receiver in names(receiver_definitions)) {
    ss <- scores_summary[scores_summary$receiver == receiver, ]
    ss$replicate_id <- factor(ss$replicate_id, levels = replicate_table$replicate_id)
    if (!any(is.finite(ss$mean_score))) {
      print(empty_plot(receiver, "Not assessable: check annotation, cell counts and signature coverage.")); next
    }
    print(ggplot(ss, aes(replicate_id, mean_score, colour = condition)) + geom_point(size = 2.5, na.rm = TRUE) +
            facet_wrap(~ signature, ncol = 3, scales = "free_y") + scale_colour_manual(values = condition_colors, limits = condition_levels) +
            labs(title = receiver, x = NULL, y = "Mean UCell score") + theme_research +
            theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top"))
  }
  if (nrow(ligand_results)) draw_ligand_plots()
}
save_stage <- function(stage) {
  info <- list(stage = stage, started_at = as.character(started_at), finished_at = as.character(Sys.time()),
               input_file = path.expand(input_file), input_object = input_object, analysis_label = analysis_label,
               sample_table = sample_table, contrast = "VG161 - Vehicle", design = "independent groups; ~ condition",
               receiver_definitions = receiver_definitions, tumor_labels = tumor_labels, annotation_audit = annotation_audit,
               host_prefix = host_prefix, viral_pattern = viral_pattern, viral_feature_ids = viral_feature_ids,
               min_cells = min_cells, min_samples_per_condition = min_samples_per_condition,
               min_signature_genes = min_signature_genes, min_signature_coverage = min_signature_coverage,
               ucell_max_rank = min(ucell_max_rank, nrow(host_counts)), expression_fraction = expression_fraction,
               filter_min_count = filter_min_count, filter_min_total_count = filter_min_total_count,
               de_fdr = de_fdr, target_min_abs_logFC = target_min_abs_logFC, top_target_n = top_target_n,
               min_target_genes = min_target_genes, min_background_genes = min_background_genes,
               target_exclusion_genes = target_exclusion_genes, seed = seed, prior = prior_info,
               limits = c("Exploratory with two biological samples per condition.", "Broad Tcell/Myeloid states can reflect subtype composition.",
                          "NicheNet ranking has no ligand P/FDR and no activating/inhibitory sign.", "NCL specificity and surface localization are not established.",
                          "VG161 engineered payloads and shared host/payload feature names are not fully represented by this host prior."),
               package_versions = vapply(c("Seurat", "SeuratObject", "Matrix", "UCell", "edgeR", "ggplot2", "curl"),
                                         function(p) as.character(packageVersion(p)), character(1)), session = capture.output(sessionInfo()))
  if (requireNamespace("nichenetr", quietly = TRUE)) info$package_versions <- c(info$package_versions, nichenetr = as.character(packageVersion("nichenetr")))
  write.csv(scores_summary, file.path(outdir, "scores.csv"), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  write.csv(ligand_results, file.path(outdir, "ligands.csv"), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  saveRDS(list(info = info, scores = scores_summary, signatures = signatures, signature_coverage = signature_coverage,
               score_contrasts = score_contrasts, signature_expression = signature_expression,
               receiver_audit = receiver_audit, cell_qc = cell_qc, layer_audit = layer_audit, viral_exclusion_audit = viral_exclusion_audit,
               receiver_DE = de_results, receiver_models = de_models, DE_audit = de_audit,
               ligands = ligand_results, model_audit = model_audit, candidate_availability = candidate_availability,
               model_coverage = model_coverage, target_tables = target_tables, target_links = target_links),
          file.path(outdir, "results.rds"), compress = TRUE)
  draw_plots()
  message("已保存阶段：", stage, "；", normalizePath(outdir))
}
save_stage("immune_function_complete_NicheNet_pending")

##### 第二阶段：NicheNet官方先验下载与校验 #####
# 从本节续跑依赖前面完整对象；不要重新执行清理环境。该模型与前次GMT不同，LT文件约250MiB。
model_manifest <- data.frame(
  file = c("ligand_target_matrix_nsga2r_final.rds", "lr_network_human_21122021.rds"),
  bytes = c(262116605, 20724),
  md5 = c("b09606b04b2d4490418d9028c0e58b9f", "2a155f81e9ffffd5d5e709fe66bcc465"),
  stringsAsFactors = FALSE)
model_manifest$url <- paste0("https://zenodo.org/records/7074291/files/", model_manifest$file, "?download=1")
download_model <- function(url, destination, expected_bytes, expected_md5) {
  valid_file <- function(path) file.exists(path) && isTRUE(file.info(path)$size == expected_bytes) &&
    identical(unname(tools::md5sum(path)), expected_md5)
  if (valid_file(destination)) return(invisible(destination))
  if (model_source == "local") stop("本地先验缺失或校验失败：", destination)
  partial <- paste0(destination, ".part")
  last_error <- ""
  for (attempt in seq_len(download_retries)) {
    if (file.exists(partial) && file.info(partial)$size >= expected_bytes && !valid_file(partial)) unlink(partial)
    success <- tryCatch({
      if (!valid_file(partial)) {
        args <- list(urls = url, destfiles = partial, resume = TRUE, progress = TRUE,
                     connecttimeout = 30, low_speed_limit = 1024, low_speed_time = 120)
        timeout_argument <- if ("multi_timeout" %in% names(formals(curl::multi_download))) "multi_timeout" else "timeout"
        args[[timeout_argument]] <- download_timeout
        downloaded <- do.call(curl::multi_download, args)
        if (!isTRUE(downloaded$success[1])) stop("网络传输未完成；已保留.part。")
      }
      if (!valid_file(partial)) stop("先验大小或MD5不符，不接受部分下载。")
      if (!file.copy(partial, destination, overwrite = TRUE)) stop("不能写入先验缓存。")
      unlink(partial) # 只删除刚完成且已验证的中转文件。
      TRUE
    }, error = function(e) { last_error <<- conditionMessage(e); FALSE })
    if (success) return(invisible(destination))
    message("先验第", attempt, "次获取未完成：", last_error)
  }
  stop(last_error, "。第一阶段已保存，可从第二阶段重试，勿清理当前环境。", call. = FALSE)
}

##### 全候选排名函数与结果图 #####
draw_ligand_plots <- function() {
  for (receiver in unique(ligand_results$receiver)) for (source in unique(ligand_results$source)) {
    for (rule in c("FDR", "Top_QLF")) {
      zz <- ligand_results[ligand_results$receiver == receiver & ligand_results$source == source &
                             ligand_results$unit == "group" & ligand_results$target_rule == rule & is.finite(ligand_results$rank), ]
      if (!nrow(zz)) next
      z <- bind_rows_base(lapply(split(zz, zz$direction), function(a) {
        a <- a[order(a$rank), ]; a[a$rank <= top_ligands_plot | a$ligand == target_gene, ]
      }))
      z$highlight <- ifelse(z$ligand == target_gene, target_gene, "Other ligands")
      z$display <- factor(paste(z$direction, z$ligand, sep = "::"),
                          levels = unique(paste(z$direction, z$ligand, sep = "::"))[rev(seq_len(nrow(z)))])
      print(ggplot(z, aes(aupr_corrected, display, fill = highlight)) + geom_col(width = 0.7) +
              geom_vline(xintercept = 0, colour = "#656A70", linewidth = 0.25) +
              facet_wrap(~ direction, ncol = 2, scales = "free_y") +
              scale_y_discrete(labels = function(x) sub("^.*::", "", x)) +
              scale_fill_manual(values = setNames(c("#C35B3F", "#426B91"), c(target_gene, "Other ligands"))) +
              labs(title = paste(receiver, source, rule, sep = " | "), x = "AUPR corrected", y = NULL) +
              theme_research + theme(legend.position = "top"))
    }
    z <- ligand_results[ligand_results$receiver == receiver & ligand_results$source == source & ligand_results$unit == "sample", ]
    if (!nrow(z)) next
    finite_ranks <- z[is.finite(z$rank_percentile), ]
    if (!nrow(finite_ranks)) next
    by_ligand <- aggregate(rank_percentile ~ ligand, data = finite_ranks, FUN = mean)
    selected <- unique(c(head(by_ligand$ligand[order(by_ligand$rank_percentile)], top_ligands_plot), target_gene))
    grid <- expand.grid(ligand = selected, replicate_id = unique(z$replicate_id), direction = c("Up", "Down"), stringsAsFactors = FALSE)
    z <- merge(grid, z, by = c("ligand", "replicate_id", "direction"), all.x = TRUE, sort = FALSE)
    z$ligand <- factor(z$ligand, levels = rev(selected))
    z$replicate_id <- factor(z$replicate_id, levels = replicate_table$replicate_id)
    z$rank_label <- ifelse(is.finite(z$rank), as.character(z$rank), "NA")
    z$text_colour <- ifelse(is.finite(z$rank_percentile) & z$rank_percentile < 0.35, "white", "#454545")
    p <- ggplot(z, aes(replicate_id, ligand, fill = rank_percentile)) + geom_tile(colour = "white", linewidth = 0.4) +
      geom_text(aes(label = rank_label, colour = text_colour), size = 3) + scale_colour_identity() +
      facet_wrap(~ direction, ncol = 2) +
      scale_fill_gradient(low = "#2D5375", high = "#F3F7FB", limits = c(0, 1), na.value = "#B8BDC3", name = "Rank percentile") +
      labs(title = paste(receiver, source, "Sample response", sep = " | "), x = NULL, y = NULL) +
      theme_research + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right", legend.title = element_text())
    unavailable <- z[!is.na(z$available_this_sample) & !z$available_this_sample & is.finite(z$rank), ]
    if (nrow(unavailable)) p <- p + geom_point(data = unavailable, aes(shape = "Unavailable in sample"),
                                               colour = "#C35B3F", size = 2, stroke = 0.7, position = position_nudge(x = 0.24, y = 0.23)) +
      scale_shape_manual(values = c("Unavailable in sample" = 4), name = NULL)
    print(p)
  }
}

if (run_nichenet && length(de_models)) {
  if (!requireNamespace("nichenetr", quietly = TRUE)) stop("请安装nichenetr后从第二阶段续跑；免疫功能和接收端DE已保存。")
  library(nichenetr)
  model_dir <- path.expand(model_dir)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  model_manifest$path <- file.path(model_dir, model_manifest$file)
  for (i in seq_len(nrow(model_manifest))) download_model(model_manifest$url[i], model_manifest$path[i], model_manifest$bytes[i], model_manifest$md5[i])
  ligand_target_matrix <- readRDS(model_manifest$path[1])
  lr_network <- as.data.frame(readRDS(model_manifest$path[2]))
  if (!is.matrix(ligand_target_matrix) || !is.numeric(ligand_target_matrix) ||
      is.null(rownames(ligand_target_matrix)) || is.null(colnames(ligand_target_matrix)) ||
      anyDuplicated(rownames(ligand_target_matrix)) || anyDuplicated(colnames(ligand_target_matrix)) ||
      any(!is.finite(ligand_target_matrix)) || !all(c("from", "to") %in% names(lr_network))) stop("先验结构不符合human ligand-target/LR模型。")
  lr_network <- unique(lr_network[, c("from", "to")])
  lr_network$from <- as.character(lr_network$from); lr_network$to <- as.character(lr_network$to)
  lr_network <- lr_network[!is.na(lr_network$from) & !is.na(lr_network$to) & nzchar(lr_network$from) & nzchar(lr_network$to), ]
  prior_info <- list(status = "loaded", model = "NicheNet human final v2; Zenodo 7074291", files = model_manifest,
                     matrix_dimensions = dim(ligand_target_matrix), n_LR = nrow(lr_network))
  model_coverage <- data.frame(ligand = target_gene, receptor = c("NCL", "LRP1", "SDC4", "PTPRZ1", "ALK"))
  model_coverage$ligand_in_LT <- target_gene %in% colnames(ligand_target_matrix)
  model_coverage$pair_in_LR <- vapply(model_coverage$receptor, function(r) any(lr_network$from == target_gene & lr_network$to == r), logical(1))
  message("固定模型MDK–NCL覆盖：", model_coverage$pair_in_LR[model_coverage$receptor == "NCL"], "；不补造网络未覆盖的边。")
  
  ##### 各样本发送端和接收端表达：不按发送端DEG筛候选 #####
  expressed <- function(cells) {
    if (length(cells) < min_cells) return(character())
    rownames(host_counts)[Matrix::rowSums(host_counts[, cells, drop = FALSE] > 0) / length(cells) >= expression_fraction]
  }
  sender_expression <- list()
  for (id in replicate_table$replicate_id) {
    cc <- rownames(meta)[meta$replicate_id == id]
    by_type <- split(cc, as.character(meta[cc, celltype_col]))
    source_genes <- lapply(by_type, expressed)
    sender_expression[[id]] <- list(All_sources = unique(unlist(source_genes, use.names = FALSE)),
                                    Tumor = expressed(cc[meta[cc, celltype_col] %in% tumor_labels]))
  }
  candidate_rows <- list(); audit_rows <- list(); ligand_rows <- list(); link_rows <- list()
  model_index <- 0L
  for (receiver in names(de_models)) {
    de <- de_results[[receiver]]
    dm <- de_models[[receiver]]
    background <- setdiff(intersect(de$gene, rownames(ligand_target_matrix)), target_exclusion_genes)
    if (length(background) < min_background_genes) {
      audit_rows[[length(audit_rows) + 1L]] <- data.frame(receiver = receiver, source = "All", unit = "group", replicate_id = NA_character_,
                                                          target_rule = "all", direction = "Both", n_targets = NA_integer_, n_background = length(background), n_candidates = NA_integer_,
                                                          n_FDR_supported = NA_integer_, status = "insufficient_prior_background", model_id = NA_character_)
      next
    }
    if (length(background) < 5000L) warning(receiver, " 的receiver先验背景少于5000基因；100只是计算门槛，请检查物种、ID和覆盖。")
    rownames(de) <- de$gene
    lt <- ligand_target_matrix[background, , drop = FALSE]
    mm <- dm$samples
    receptor_expression <- lapply(rownames(mm), function(id) expressed(receiver_cells[[receiver]][meta[receiver_cells[[receiver]], "replicate_id"] == id]))
    names(receptor_expression) <- rownames(mm)
    sample_effect <- sweep(dm$logCPM[background, , drop = FALSE], 1L,
                           rowMeans(dm$logCPM[background, mm$condition == "Vehicle", drop = FALSE]), "-")
    for (source in c("All_sources", "Tumor")) {
      available <- lapply(rownames(mm), function(id) {
        compatible <- unique(lr_network$from[lr_network$to %in% receptor_expression[[id]]])
        intersect(colnames(lt), intersect(compatible, sender_expression[[id]][[source]]))
      })
      names(available) <- rownames(mm)
      # 两组各自求样本交集，再取并集：保留处理诱导与处理后消失的候选。
      common_by_condition <- lapply(condition_levels, function(cond) Reduce(intersect, available[rownames(mm)[mm$condition == cond]]))
      names(common_by_condition) <- condition_levels
      candidates <- sort(unique(unlist(common_by_condition, use.names = FALSE)))
      for (id in rownames(mm)) {
        all_prior_ligands <- sort(intersect(colnames(lt), unique(lr_network$from)))
        candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(receiver = receiver, source = source,
                                                                    replicate_id = id, condition = as.character(mm[id, "condition"]), ligand = all_prior_ligands,
                                                                    sender_expressed = all_prior_ligands %in% sender_expression[[id]][[source]],
                                                                    receiver_compatible = all_prior_ligands %in% lr_network$from[lr_network$to %in% receptor_expression[[id]]],
                                                                    available_this_sample = all_prior_ligands %in% available[[id]],
                                                                    supported_in_Vehicle = all_prior_ligands %in% common_by_condition$Vehicle,
                                                                    supported_in_VG161 = all_prior_ligands %in% common_by_condition$VG161,
                                                                    in_fixed_candidate_universe = all_prior_ligands %in% candidates)
      }
      
      ##### 分开评估FDR目标、预先声明的探索目标与逐样本响应 #####
      run_specs <- list()
      for (direction in c("Up", "Down")) {
        sign_use <- if (direction == "Up") 1 else -1
        ok_genes <- background[sign_use * de[background, "logFC"] >= target_min_abs_logFC &
                                 sign_use * de[background, "rank_stat"] > 0]
        ranked <- ok_genes[order(-sign_use * de[ok_genes, "rank_stat"], ok_genes)]
        run_specs[[length(run_specs) + 1L]] <- list(unit = "group", id = NA_character_, rule = "FDR", direction = direction,
                                                    genes = ranked[de[ranked, "FDR"] < de_fdr])
        run_specs[[length(run_specs) + 1L]] <- list(unit = "group", id = NA_character_, rule = "Top_QLF", direction = direction,
                                                    genes = head(ranked, top_target_n))
        for (id in rownames(mm)[mm$condition == "VG161"]) {
          lfc <- sample_effect[, id]
          keep <- names(lfc)[is.finite(lfc) & sign_use * lfc >= target_min_abs_logFC & sign_use * lfc > 0]
          keep <- keep[order(-sign_use * lfc[keep], keep)]
          run_specs[[length(run_specs) + 1L]] <- list(unit = "sample", id = id, rule = "Sample_top_effect", direction = direction,
                                                      genes = head(keep, top_target_n))
        }
      }
      for (spec in run_specs) {
        model_index <- model_index + 1L
        model_id <- paste0("M", model_index)
        targets <- spec$genes
        sign_use <- if (spec$direction == "Up") 1 else -1
        fdr_supported <- sum(de[targets, "FDR"] < de_fdr & sign_use * de[targets, "logFC"] > 0)
        status <- if (length(targets) < min_target_genes) "insufficient_targets" else if (length(background) - length(targets) < min_target_genes) "insufficient_non_targets" else
          if (length(candidates) < min_candidate_ligands) "insufficient_candidates" else "assessed"
        audit_rows[[length(audit_rows) + 1L]] <- data.frame(receiver = receiver, source = source, unit = spec$unit,
                                                            replicate_id = spec$id, target_rule = spec$rule, direction = spec$direction, n_targets = length(targets),
                                                            n_background = length(background), n_candidates = length(candidates), n_FDR_supported = fdr_supported, status = status, model_id = model_id)
        target_tables[[model_id]] <- data.frame(gene = targets, group_logFC = de[targets, "logFC"],
                                                group_P = de[targets, "PValue"], group_FDR = de[targets, "FDR"],
                                                sample_logCPM_difference = if (spec$unit == "sample") sample_effect[targets, spec$id] else rep(NA_real_, length(targets)))
        if (status != "assessed") next
        informative <- candidates[vapply(candidates, function(lig) stats::sd(lt[, lig]) > 0, logical(1))]
        activity <- if (length(informative)) as.data.frame(nichenetr::predict_ligand_activities(
          geneset = targets, background_expressed_genes = background, ligand_target_matrix = lt,
          potential_ligands = informative, single = TRUE)) else data.frame(test_ligand = character(), auroc = double(), aupr = double(), aupr_corrected = double(), pearson = double())
        out <- data.frame(ligand = candidates, stringsAsFactors = FALSE)
        ix <- match(out$ligand, activity$test_ligand)
        out$auroc <- activity$auroc[ix]; out$aupr <- activity$aupr[ix]
        out$aupr_corrected <- activity$aupr_corrected[ix]; out$pearson <- activity$pearson[ix]
        out$score_status <- ifelse(!out$ligand %in% informative, "constant_prior", ifelse(is.finite(out$aupr_corrected), "scored", "non_finite_score"))
        out <- out[order(-out$aupr_corrected, -out$pearson, out$ligand, na.last = TRUE), ]
        out$rank <- NA_real_
        scored <- which(out$score_status == "scored")
        out$rank[scored] <- seq_along(scored)
        out$rank_percentile <- if (length(scored) > 1L) (out$rank - 1) / (length(scored) - 1) else NA_real_
        out$receiver <- receiver; out$source <- source; out$unit <- spec$unit; out$replicate_id <- spec$id
        out$target_rule <- spec$rule; out$direction <- spec$direction; out$model_id <- model_id
        out$n_candidates <- length(candidates); out$n_scored <- length(scored); out$n_targets <- length(targets)
        out$n_background <- length(background); out$n_FDR_supported <- fdr_supported
        out$available_this_sample <- if (spec$unit == "sample") out$ligand %in% available[[spec$id]] else NA
        ligand_rows[[length(ligand_rows) + 1L]] <- out
        link_ligands <- intersect(unique(c(head(out$ligand[scored], top_ligands_links), target_gene)), candidates)
        for (lig in link_ligands) {
          w <- lt[targets, lig]; w <- sort(w[is.finite(w) & w > 0], decreasing = TRUE)
          w <- head(w, top_target_links)
          if (length(w)) link_rows[[length(link_rows) + 1L]] <- data.frame(model_id = model_id, receiver = receiver,
                                                                           source = source, ligand = lig, target = names(w), regulatory_potential = as.numeric(w))
        }
        message(receiver, "/", source, "/", spec$unit, "/", spec$rule, "/", spec$direction, ": ", length(candidates), "个候选。")
      }
    }
  }
  ranked_results <- bind_rows_base(ligand_rows)
  if (nrow(ranked_results)) ligand_results <- ranked_results
  model_audit <- bind_rows_base(audit_rows)
  candidate_availability <- bind_rows_base(candidate_rows)
  target_links <- bind_rows_base(link_rows)
  rm(ligand_target_matrix)
  invisible(gc())
  save_stage("complete")
} else {
  prior_info$status <- if (!run_nichenet) "disabled_by_configuration" else "no_receiver_with_sufficient_samples"
  save_stage(if (!run_nichenet) "immune_function_only" else "NicheNet_not_assessable")
}
message("完成：scores.csv、ligands.csv、plots.pdf、results.rds。免疫功能分数按样本展示；NicheNet分数没有配体P/FDR。")
