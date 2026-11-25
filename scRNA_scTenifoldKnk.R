##### 清理运行环境 #####
rm(list = ls(all.names = TRUE))
graphics.off()
invisible(gc())

##### 加载依赖 #####
if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))
library(Seurat)
library(Matrix)
library(scTenifoldKnk)
library(fgsea)
library(BiocParallel)
library(ggplot2)

##### 路径、样本与参数 #####
input_file <- "~/rawdata/scRNA_virus/virus_3D/virus_3D_anno.RData"
input_object <- "sce"
outdir <- file.path(path.expand("~/OV"), "scRNA_virus", "MDK_knockout_D3")
# 优先读取本地的一对GMT；任一缺失时改用Broad官方固定版本。
script_args <- grep("^--file=", commandArgs(), value = TRUE)
source_files <- unlist(lapply(sys.frames(), function(z) z$ofile), use.names = FALSE)
script_file <- if (length(source_files)) tail(source_files, 1L) else if (length(script_args)) sub("^--file=", "", script_args[1L]) else ""
code_dir <- if (nzchar(script_file)) dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)) else if (file.exists("code/scRNA_MDK_knockout.R")) normalizePath("code", winslash = "/") else getwd()
hallmark_gmt <- file.path(code_dir, "gene_sets", "hallmark.gmt")
gobp_gmt <- file.path(code_dir, "gene_sets", "gobp.gmt")
gmt_version <- "MSigDB 2026.1.Hs / human symbols"
sample_col <- "orig.ident"
celltype_col <- "celltype"
tumor_labels <- "Epithelial" # 沿用既有注释；Epithelial标签本身不验证恶性身份。
condition_levels <- c("Vehicle", "VG161")
sample_table <- data.frame(
  sample_id = c("Vehicle_3D_1", "Vehicle_3D_2", "VG161_3D_1", "VG161_3D_2"),
  replicate_id = c("Vehicle_3D_1", "Vehicle_3D_2", "VG161_3D_1", "VG161_3D_2"),
  condition = c("Vehicle", "Vehicle", "VG161", "VG161"), stringsAsFactors = FALSE)
# 两组独立动物；同一肿瘤的技术拆库才填写相同replicate_id。此处不使用infected字段。
host_prefix <- "" # 输入须已完成物种分离；例如human_MDK才配置human_。
viral_pattern <- "^(VG161|HSV[-_]?1)([-_.]|$)"
viral_feature_ids <- character() # 按真实参考补无前缀病毒/外源载荷ID；不删除宿主ISG。
target_gene <- "MDK"
seed <- 123L
min_cells_per_replicate <- 100L
max_cells_per_replicate <- 500L
min_detection_fraction <- 0.05 # 共同背景：每组均有至少5%平衡抽样细胞检出。
max_network_genes <- 3000L    # 合并平衡细胞vst排序，最多3000基因（含MDK）。
min_network_genes <- 200L
network_nNet <- 10L
network_nCells <- 500L       # 包内每网有放回bootstrap的细胞数，不是独立重复数。
network_nComp <- 3L
network_q <- 0.90
tensor_K <- 3L
manifold_nDim <- 2L
n_cores <- 2L
gene_fdr <- 0.05
pathway_fdr <- 0.05
pathway_min_size <- 15L
pathway_max_size <- 500L
gsea_nperm_simple <- 10000L
gsea_eps <- 1e-10
max_extra_immune_pathways <- 8L
display_min_NES <- 1.0 # 仅过滤非预设项的整行展示，不改变检验/FDR。
focus_pathways <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_COMPLEMENT",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
  "GOBP_RESPONSE_TO_VIRUS", "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_CHEMOKINE_MEDIATED_SIGNALING_PATHWAY", "GOBP_REGULATION_OF_T_CELL_ACTIVATION")
immune_regex <- paste(c("IMMUN", "LYMPHOCYTE", "LEUKOCYTE", "T_CELL", "B_CELL", "MYELOID",
                        "MACROPHAGE", "MONOCYTE", "DENDRITIC", "NEUTROPHIL", "NATURAL_KILLER", "NK_CELL",
                        "ANTIGEN", "CYTOKINE", "CHEMOKINE", "CHEMOTAXIS", "INFLAMMAT", "INTERFERON",
                        "COMPLEMENT", "PHAGOCYT", "RESPONSE_TO_VIRUS", "ANTIMICROBIAL", "DEFENSE_RESPONSE",
                        "INTERSPECIES", "MODULATION_OF_PROCESS_OF_ANOTHER_ORGANISM"), collapse = "|")
condition_colors <- c(Vehicle = "#426B91", VG161 = "#C35B3F")

##### 辅助函数与配置校验 #####
read_sce <- function(path, object_name) {
  path <- path.expand(path)
  if (!file.exists(path)) stop("找不到输入：", path)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  if (!object_name %in% loaded) stop("RData中没有对象：", object_name)
  get(object_name, envir = env)
}
is_viral_feature <- function(ids) {
  symbols <- ids
  if (nzchar(host_prefix)) {
    hit <- startsWith(ids, host_prefix)
    symbols[hit] <- substring(ids[hit], nchar(host_prefix) + 1L)
  }
  hit <- if (nzchar(viral_pattern)) grepl(viral_pattern, ids, ignore.case = TRUE) |
    grepl(viral_pattern, symbols, ignore.case = TRUE) else rep(FALSE, length(ids))
  hit | ids %in% viral_feature_ids | symbols %in% viral_feature_ids
}
read_gmt <- function(path) {
  if (!file.exists(path)) stop("找不到GMT：", path)
  rows <- strsplit(readLines(path, warn = FALSE, encoding = "UTF-8"), "\t", fixed = TRUE)
  if (!length(rows) || any(lengths(rows) < 3L)) stop("GMT格式不完整：", path)
  ids <- vapply(rows, `[`, character(1), 1L)
  if (anyDuplicated(ids) || any(!nzchar(ids))) stop("GMT集合名重复或为空。")
  out <- lapply(rows, function(z) unique(z[-c(1L, 2L)][nzchar(z[-c(1L, 2L)])]))
  names(out) <- ids
  out
}
bind_tables <- function(x) {
  x <- Filter(function(z) is.data.frame(z) && nrow(z), x)
  if (length(x)) do.call(rbind, x) else data.frame()
}
bh_valid <- function(p, n_total) {
  ans <- rep(NA_real_, length(p)); ok <- is.finite(p)
  if (any(ok)) ans[ok] <- p.adjust(p[ok], method = "BH", n = n_total)
  ans
}
required <- c("sample_id", "replicate_id", "condition")
if (!all(required %in% names(sample_table)) || !nrow(sample_table) ||
    anyNA(sample_table[, required]) || any(!nzchar(as.matrix(sample_table[, required]))) ||
    anyDuplicated(sample_table$sample_id)) stop("样本表缺字段、空值或sample_id重复。")
if (!identical(condition_levels, c("Vehicle", "VG161")) || !all(sample_table$condition %in% condition_levels)) stop("请明确Vehicle与VG161样本映射。")
if (any(vapply(split(sample_table$condition, sample_table$replicate_id), function(z) length(unique(z)) != 1L, logical(1)))) stop("同一动物跨组：当前是独立组分析。")
replicate_table <- unique(sample_table[, c("replicate_id", "condition")])
n_animals <- table(factor(replicate_table$condition, levels = condition_levels))
if (any(n_animals < 2L) || length(unique(as.integer(n_animals))) != 1L) stop("本流程要求两组相同数量、每组至少2个真实独立动物。")
numeric_parameters <- c(seed, min_cells_per_replicate, max_cells_per_replicate, max_network_genes,
                        min_network_genes, network_nNet, network_nCells, network_nComp, tensor_K, manifold_nDim,
                        n_cores, pathway_min_size, pathway_max_size, gsea_nperm_simple)
if (any(!is.finite(numeric_parameters) | numeric_parameters < 1 | numeric_parameters != round(numeric_parameters)) ||
    min_cells_per_replicate > max_cells_per_replicate || max_network_genes < min_network_genes ||
    min_network_genes <= network_nComp + 2L || network_nComp < 2L || pathway_min_size > pathway_max_size ||
    any(!is.finite(c(min_detection_fraction, network_q, gene_fdr, pathway_fdr, gsea_eps, display_min_NES))) ||
    min_detection_fraction <= 0 || min_detection_fraction > 1 || network_q <= 0 || network_q >= 1 ||
    gene_fdr <= 0 || gene_fdr >= 1 || pathway_fdr <= 0 || pathway_fdr >= 1 || gsea_eps <= 0 || display_min_NES < 0) stop("数值参数范围不正确。")
if (packageVersion("scTenifoldKnk") < "2.0.0") stop("本脚本要求scTenifoldKnk >=2.0.0，已验证2.0.0。请更新包。")
if (!all(c("dr_direction", "ma_method", "seed") %in% names(formals(scTenifoldKnk::scTenifoldKnk)))) stop("scTenifoldKnk接口与本脚本不兼容。")
set.seed(seed)
started_at <- Sys.time()
run_warnings <- character()
outdir <- path.expand(outdir)

##### 读取肿瘤细胞与原始RNA层 #####
sce <- read_sce(input_file, input_object)
if (!inherits(sce, "Seurat") || !"RNA" %in% names(sce@assays)) stop("输入必须为带RNA counts的Seurat对象。")
meta <- sce[[]]
if (!all(c(sample_col, celltype_col) %in% names(meta)) || anyDuplicated(rownames(meta))) stop("metadata字段缺失或条码重复。")
if (!all(sample_table$sample_id %in% as.character(meta[[sample_col]]))) stop("样本表与对象不一致。对象样本：", paste(unique(meta[[sample_col]]), collapse = ", "))
keep <- !is.na(meta[[sample_col]]) & meta[[sample_col]] %in% sample_table$sample_id &
  !is.na(meta[[celltype_col]]) & meta[[celltype_col]] %in% tumor_labels
meta <- meta[keep, , drop = FALSE]
if (!nrow(meta)) stop("所选样本中没有配置的肿瘤标签。")
meta$sample_id <- as.character(meta[[sample_col]])
meta$replicate_id <- sample_table$replicate_id[match(meta$sample_id, sample_table$sample_id)]
meta$condition <- sample_table$condition[match(meta$sample_id, sample_table$sample_id)]
rna <- sce[["RNA"]]
if (inherits(rna, "Assay5")) {
  layer_names <- Layers(rna, search = "^counts($|\\.)")
  if (!length(layer_names)) stop("RNA缺原始counts层，不能用data/SCT/integrated替代。")
  count_layers <- lapply(layer_names, function(k) {
    z <- LayerData(rna, layer = k, fast = FALSE)
    z[, intersect(colnames(z), rownames(meta)), drop = FALSE]
  })
  names(count_layers) <- layer_names
  count_layers <- Filter(function(z) ncol(z) > 0L, count_layers)
} else if (inherits(rna, "Assay")) {
  z <- if (packageVersion("SeuratObject") >= "5.0.0") GetAssayData(sce, assay = "RNA", layer = "counts") else GetAssayData(sce, assay = "RNA", slot = "counts")
  count_layers <- list(counts = z[, intersect(colnames(z), rownames(meta)), drop = FALSE])
} else stop("不支持的RNA assay。")
all_cells <- unlist(lapply(count_layers, colnames), use.names = FALSE)
if (anyDuplicated(all_cells) || !setequal(all_cells, rownames(meta))) stop("counts层条码重复或所选细胞缺counts。")
for (k in names(count_layers)) {
  z <- count_layers[[k]]
  if (is.null(rownames(z)) || is.null(colnames(z)) || anyNA(rownames(z)) || any(!nzchar(rownames(z))) ||
      anyDuplicated(rownames(z)) || anyDuplicated(colnames(z))) stop(k, " feature/条码缺失或重复。")
  z <- as(Matrix(z, sparse = TRUE), "dgCMatrix")
  if (any(!is.finite(z@x) | z@x < 0 | abs(z@x - round(z@x)) > 1e-6)) stop("RNA counts须为非负有限整数。")
  count_layers[[k]] <- z
}
host_features <- lapply(count_layers, function(z) {
  ids <- rownames(z)
  if (nzchar(host_prefix)) ids <- ids[startsWith(ids, host_prefix)]
  ids[!is_viral_feature(ids)]
})
common_features <- Reduce(intersect, host_features)
symbols <- if (nzchar(host_prefix)) substring(common_features, nchar(host_prefix) + 1L) else common_features
if (length(symbols) < min_network_genes || anyDuplicated(symbols) || mean(grepl("^ENSG[0-9]", symbols)) > 0.5) stop("共同宿主symbol过少/重复或主要为Ensembl ID；先核对物种与映射。")
if (!target_gene %in% symbols) stop("MDK没有在所有纳入counts层中共同测量；不能补零模拟敲除。")
layer_audit <- data.frame(layer = names(count_layers), n_cells = vapply(count_layers, ncol, integer(1)),
                          n_host_features = lengths(host_features), n_common_features = length(symbols), stringsAsFactors = FALSE)
viral_exclusion <- bind_tables(lapply(names(count_layers), function(k) {
  ids <- rownames(count_layers[[k]])[is_viral_feature(rownames(count_layers[[k]]))]
  data.frame(layer = rep(k, length(ids)), feature = ids, stringsAsFactors = FALSE)
}))
host_counts <- do.call(cbind, lapply(count_layers, function(z) z[common_features, , drop = FALSE]))
rownames(host_counts) <- symbols
meta <- meta[colnames(host_counts), , drop = FALSE]
host_umi <- Matrix::colSums(host_counts)
cell_audit <- data.frame(cell = rownames(meta), meta[, c("sample_id", "replicate_id", "condition")],
                         host_UMI = as.numeric(host_umi), valid_host = host_umi > 0, selected_network = FALSE,
                         stringsAsFactors = FALSE, row.names = NULL)
if (any(lengths(host_features) != length(symbols))) warning("counts层feature不同：使用共同宿主交集，未测项不补零。")
meta <- meta[host_umi > 0, , drop = FALSE]
host_counts <- host_counts[, rownames(meta), drop = FALSE]
rm(sce, rna, count_layers, z, host_umi)
invisible(gc())

##### 每动物等量抽样与共同基因背景 #####
n_valid <- table(factor(meta$replicate_id, levels = replicate_table$replicate_id))
if (any(n_valid < min_cells_per_replicate)) stop("有动物的有效肿瘤细胞不足", min_cells_per_replicate, "：", paste(names(n_valid), n_valid, sep = "=", collapse = ", "))
cells_per_replicate <- min(max_cells_per_replicate, min(n_valid))
set.seed(seed)
selected_cells <- unlist(lapply(replicate_table$replicate_id, function(id) {
  sample(rownames(meta)[meta$replicate_id == id], cells_per_replicate, replace = FALSE)
}), use.names = FALSE)
cell_audit$selected_network <- cell_audit$cell %in% selected_cells
balanced_counts <- host_counts[, selected_cells, drop = FALSE]
balanced_meta <- meta[selected_cells, , drop = FALSE]
detection <- sapply(condition_levels, function(g) Matrix::rowMeans(balanced_counts[, balanced_meta$condition == g, drop = FALSE] > 0))
eligible <- rownames(balanced_counts)[apply(detection >= min_detection_fraction, 1L, all)]
if (!target_gene %in% eligible) stop("MDK在至少一组的检出率不足，不能稳定建模。两组检出率：", paste(round(detection[target_gene, ], 4), collapse = ", "))
if (length(eligible) < min_network_genes) stop("两组共同检出的可建网基因不足。")
# 合并两组平衡细胞仅用于固定相同基因背景；不依据KO结果选择基因。
vst_object <- CreateSeuratObject(counts = balanced_counts[eligible, , drop = FALSE], min.cells = 0, min.features = 0)
vst_object <- FindVariableFeatures(vst_object, selection.method = "vst",
                                   nfeatures = min(max_network_genes, length(eligible)), verbose = FALSE)
vst_order <- VariableFeatures(vst_object)
network_genes <- unique(c(target_gene, head(setdiff(vst_order, target_gene), max_network_genes - 1L)))
network_genes <- sort(intersect(network_genes, eligible))
if (length(network_genes) < min_network_genes) stop("vst筛选后的共同基因不足。")
gene_selection <- data.frame(gene = rownames(balanced_counts),
                             detection_Vehicle = detection[, "Vehicle"], detection_VG161 = detection[, "VG161"],
                             eligible = rownames(balanced_counts) %in% eligible,
                             in_network = rownames(balanced_counts) %in% network_genes, stringsAsFactors = FALSE)
network_counts <- balanced_counts[network_genes, , drop = FALSE]
if (any(Matrix::colSums(network_counts) <= 0)) stop("有细胞在选定网络节点中总UMI为0，无法进行CPM。请核对节点覆盖，不自动补值。")
summary <- bind_tables(lapply(seq_len(nrow(replicate_table)), function(i) {
  id <- replicate_table$replicate_id[i]; cells <- rownames(meta)[meta$replicate_id == id]
  values <- as.numeric(host_counts[target_gene, cells])
  umi <- Matrix::colSums(host_counts[, cells, drop = FALSE])
  data.frame(replicate_id = id, condition = replicate_table$condition[i],
             sample_ids = paste(sample_table$sample_id[sample_table$replicate_id == id], collapse = ";"),
             n_tumor = sum(cell_audit$replicate_id == id), n_valid_host = length(cells),
             n_network_cells = cells_per_replicate, MDK_detection = mean(values > 0),
             MDK_mean_logCP10K = mean(log1p(1e4 * values / umi)), stringsAsFactors = FALSE)
}))
rm(vst_object, host_counts, balanced_counts)
invisible(gc())
message("每动物纳入", cells_per_replicate, "细胞；共同网络基因", length(network_genes), "。")

##### 获取并读取完整通路库 #####
# 可从此节续跑。依赖前面已定义的 read_gmt()、hallmark_gmt、gobp_gmt 和 gmt_version。
gmt_paths <- c(Hallmark = path.expand(hallmark_gmt), `GO:BP` = path.expand(gobp_gmt))
gmt_source <- "local"
resolved_gmt_version <- gmt_version
if (!all(file.exists(gmt_paths))) {
  message("本地GMT缺失：", paste(gmt_paths[!file.exists(gmt_paths)], collapse = ", "),
          "。改用Broad官方MSigDB 2026.1.Hs。")
  broad_base <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/"
  broad_manifest <- c(
    "h.all.v2026.1.Hs.symbols.gmt" = "367eec875967c2cfbf664a1a065b7b8d",
    "c5.go.bp.v2026.1.Hs.symbols.gmt" = "c1903b564ead3b31260418bdd0d22b10")
  gmt_cache_dir <- file.path(tempdir(), "msigdb_2026_1_Hs")
  dir.create(gmt_cache_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(gmt_cache_dir)) stop("无法创建GMT临时目录：", gmt_cache_dir, call. = FALSE)
  download_broad_gmt <- function(file_name, expected_md5) {
    path <- file.path(gmt_cache_dir, file_name)
    if (file.exists(path) && identical(unname(tools::md5sum(path)), expected_md5)) return(path)
    tmp <- tempfile(tmpdir = gmt_cache_dir, fileext = ".gmt")
    old_options <- options(timeout = max(180, getOption("timeout", 60)))
    on.exit({ unlink(tmp); options(old_options) }, add = TRUE)
    status <- utils::download.file(paste0(broad_base, file_name), tmp,
                                   method = "libcurl", mode = "wb", quiet = TRUE)
    if (status != 0L || !identical(unname(tools::md5sum(tmp)), expected_md5)) {
      stop("下载文件未通过MD5校验：", file_name)
    }
    if (!file.copy(tmp, path, overwrite = TRUE)) stop("无法写入GMT临时文件：", path)
    path
  }
  gmt_paths[] <- tryCatch(
    vapply(names(broad_manifest), function(file_name) download_broad_gmt(file_name, broad_manifest[[file_name]]), character(1)),
    error = function(e) stop("Broad GMT获取失败：", conditionMessage(e),
                             "。如服务器不能联网，上传配套的hallmark.gmt和gobp.gmt，设置两个本地路径后从本节继续。", call. = FALSE))
  gmt_source <- "Broad"
  resolved_gmt_version <- "MSigDB 2026.1.Hs / human symbols"
}
pathways <- c(read_gmt(gmt_paths[["Hallmark"]]), read_gmt(gmt_paths[["GO:BP"]]))
if (anyDuplicated(names(pathways))) stop("两个GMT的集合名重复。")
universe <- setdiff(network_genes, target_gene) # 排除敲除目标的自效应。
measured_pathways <- lapply(pathways, intersect, y = universe)
pathway_catalog <- data.frame(pathway = names(pathways),
                              collection = ifelse(startsWith(names(pathways), "HALLMARK_"), "H", "GO:BP"),
                              full_size = lengths(pathways), measured_size = lengths(measured_pathways),
                              coverage = lengths(measured_pathways) / lengths(pathways),
                              immune = grepl(immune_regex, names(pathways)) | names(pathways) %in% focus_pathways,
                              focus = names(pathways) %in% focus_pathways,
                              stringsAsFactors = FALSE)
pathway_catalog$eligible <- with(pathway_catalog, measured_size >= pathway_min_size &
                                   measured_size <= pathway_max_size & measured_size < length(universe))
test_pathways <- measured_pathways[pathway_catalog$eligible]
if (!length(test_pathways)) stop("共同网络背景中无可测试通路；核对基因ID和网络大小。")
if (any(!focus_pathways %in% names(pathways))) warning("部分预设通路不在当前GMT中，将保留NA状态：", paste(setdiff(focus_pathways, names(pathways)), collapse = ", "))

##### 分组构建WT网络并虚拟敲除MDK #####
ko_results <- setNames(vector("list", length(condition_levels)), condition_levels)
network_edges_by_group <- setNames(vector("list", length(condition_levels)), condition_levels)
for (j in seq_along(condition_levels)) {
  group <- condition_levels[j]
  message("开始 ", group, " MDK虚拟敲除……")
  x <- network_counts[, balanced_meta$condition == group, drop = FALSE]
  t0 <- Sys.time(); warnings_group <- character()
  ans <- tryCatch(withCallingHandlers(scTenifoldKnk::scTenifoldKnk(
    countMatrix = x, gKO = target_gene, qc = FALSE,
    nc_nNet = network_nNet, nc_nCells = network_nCells, nc_nComp = network_nComp,
    nc_q = network_q, td_K = tensor_K, ma_nDim = manifold_nDim,
    ma_method = "manifold", dr_direction = FALSE, nCores = n_cores, seed = seed),
    warning = function(w) { warnings_group <<- c(warnings_group, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) e)
  status <- "ok"; note <- ""; n_outgoing <- NA_integer_; outgoing_abs_weight <- NA_real_
  dr <- data.frame(gene = network_genes, distance = NA_real_, Z = NA_real_, FC = NA_real_,
                   p.value = NA_real_, p.adj = NA_real_, stringsAsFactors = FALSE)
  if (inherits(ans, "error")) {
    status <- "model_failed"; note <- conditionMessage(ans)
  } else {
    wt <- ans$tensorNetworks$WT
    needed <- c("gene", "distance", "Z", "FC", "p.value", "p.adj")
    if (is.null(wt) || !target_gene %in% rownames(wt) || !all(needed %in% names(ans$diffRegulation))) stop("包输出结构改变，请核对scTenifoldKnk版本。")
    outgoing <- as.numeric(wt[target_gene, setdiff(colnames(wt), target_gene), drop = FALSE])
    n_outgoing <- sum(is.finite(outgoing) & outgoing != 0)
    outgoing_abs_weight <- sum(abs(outgoing))
    if (any(!is.finite(outgoing))) {
      status <- "invalid_network"; note <- "MDK的WT出边含非有限值。"
    } else if (!n_outgoing) {
      status <- "no_MDK_outgoing_edges"; note <- "WT中MDK没有保留出边，当前网络无法评估MDK敲除。"
    } else {
      dr <- as.data.frame(ans$diffRegulation[, needed], stringsAsFactors = FALSE)
      if (anyDuplicated(dr$gene) || !setequal(dr$gene, network_genes)) stop("diffRegulation基因背景与输入不一致。")
      dr <- dr[match(network_genes, dr$gene), , drop = FALSE]
      if (any(!is.finite(dr$distance)) || any(dr$distance < 0) ||
          any(!is.finite(dr$p.adj)) || any(dr$p.adj < 0 | dr$p.adj > 1)) {
        status <- "invalid_perturbation"; note <- "距离或模型FDR无效，跳过富集。"
      } else if (!any(dr$distance[dr$gene != target_gene] > 0)) {
        status <- "no_perturbation"; note <- "非MDK基因的受扰距离全为0。"
      }
    }
  }
  if (status != "ok") {
    dr[, setdiff(names(dr), "gene")] <- NA_real_
    warning(group, "：", status, "；", note)
  }
  dr$condition <- group
  dr$target <- dr$gene == target_gene
  dr$rank_percentile <- NA_real_
  if (status == "ok") {
    # 仅保存实际WT网络中MDK的非零出边，精简结果无需保留完整邻接矩阵。
    neighbors <- setdiff(colnames(wt), target_gene)
    nonzero <- outgoing != 0
    network_edges_by_group[[group]] <- data.frame(condition = group, source = target_gene,
                                                  target = neighbors[nonzero], weight = outgoing[nonzero], stringsAsFactors = FALSE)
    k <- which(!dr$target)
    dr$rank_percentile[k] <- (rank(dr$distance[k], ties.method = "average") - 1) / (length(k) - 1)
  }
  dr$immune_member <- dr$gene %in% unique(unlist(measured_pathways[pathway_catalog$immune], use.names = FALSE))
  dr$status <- status
  ko_results[[group]] <- list(genes = dr, status = status, note = note,
                              n_MDK_outgoing = n_outgoing, MDK_outgoing_abs_weight = outgoing_abs_weight,
                              n_cells = ncol(x), n_animals = unname(n_animals[group]),
                              warnings = unique(warnings_group), elapsed_seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  run_warnings <- c(run_warnings, paste0(group, ": ", warnings_group))
  rm(ans, x)
  if (exists("wt", inherits = FALSE)) rm(wt)
  invisible(gc())
}

##### 全通路检验：正向受扰富集与辅助ORA #####
pathway_results <- list()
for (j in seq_along(condition_levels)) {
  group <- condition_levels[j]; model <- ko_results[[group]]
  tab <- pathway_catalog; tab$condition <- group
  tab$ES <- tab$NES <- tab$gsea_p <- tab$gsea_FDR <- tab$gsea_log2err <- NA_real_
  tab$leading_edge <- NA_character_
  tab$ORA_overlap <- tab$ORA_query_size <- NA_integer_
  tab$ORA_p <- tab$ORA_FDR <- NA_real_
  tab$ORA_genes <- NA_character_
  tab$status <- ifelse(tab$eligible, model$status, "outside_size_limits")
  if (model$status == "ok") {
    dr <- model$genes
    stats <- setNames(dr$distance[!dr$target], dr$gene[!dr$target])
    # 同距离时按基因名稳定排序；不为打破并列而添加伪随机噪声。
    stats <- stats[order(-stats, names(stats))]
    rank_ties <- sum(duplicated(stats)) / length(stats)
    set.seed(seed)
    fg <- withCallingHandlers(fgsea::fgseaMultilevel(pathways = test_pathways, stats = stats,
                                                     minSize = pathway_min_size, maxSize = pathway_max_size, scoreType = "pos",
                                                     nPermSimple = gsea_nperm_simple, eps = gsea_eps, BPPARAM = BiocParallel::SerialParam()),
                              warning = function(w) { run_warnings <<- c(run_warnings, paste0(group, ": ", conditionMessage(w))); invokeRestart("muffleWarning") })
    fg <- as.data.frame(fg)
    ii <- match(fg$pathway, tab$pathway)
    tab$ES[ii] <- fg$ES; tab$NES[ii] <- fg$NES; tab$gsea_p[ii] <- fg$pval
    if ("log2err" %in% names(fg)) tab$gsea_log2err[ii] <- fg$log2err
    tab$leading_edge[ii] <- vapply(fg$leadingEdge, paste, character(1), collapse = ";")
    # 全部可测试H+GO:BP一起BH；先检验，后取免疫主题，不能对子集重新校正。
    tab$gsea_FDR <- bh_valid(tab$gsea_p, n_total = sum(tab$eligible))
    tab$status[tab$eligible] <- ifelse(is.finite(tab$NES[tab$eligible]) & is.finite(tab$gsea_p[tab$eligible]), "tested", "gsea_unavailable")
    query <- dr$gene[!dr$target & is.finite(dr$p.adj) & dr$p.adj < gene_fdr]
    kk <- which(tab$eligible)
    tab$ORA_query_size[kk] <- length(query)
    tab$ORA_overlap[kk] <- vapply(measured_pathways[kk], function(gs) length(intersect(gs, query)), integer(1))
    tab$ORA_genes[kk] <- vapply(measured_pathways[kk], function(gs) paste(intersect(gs, query), collapse = ";"), character(1))
    # 空命中集是已观测无显著受扰基因，p=1；网络不可评估时上方保留NA。
    tab$ORA_p[kk] <- if (!length(query)) 1 else phyper(tab$ORA_overlap[kk] - 1L,
                                                       tab$measured_size[kk], length(universe) - tab$measured_size[kk], length(query), lower.tail = FALSE)
    tab$ORA_FDR <- bh_valid(tab$ORA_p, n_total = sum(tab$eligible))
    ko_results[[group]]$rank_tie_fraction <- rank_ties
    ko_results[[group]]$n_DR_non_target <- length(query)
  }
  tab$gsea_significant <- ifelse(is.finite(tab$gsea_FDR), tab$gsea_FDR < pathway_fdr, NA)
  tab$ORA_significant <- ifelse(is.finite(tab$ORA_FDR), tab$ORA_FDR < pathway_fdr, NA)
  pathway_results[[group]] <- tab
}
genes <- bind_tables(lapply(ko_results, `[[`, "genes"))
pathway_table <- bind_tables(pathway_results)

##### 固定关注项与整行展示选择 #####
immune_candidates <- pathway_catalog$pathway[pathway_catalog$immune & !pathway_catalog$focus]
candidate_summary <- bind_tables(lapply(immune_candidates, function(id) {
  z <- pathway_table[pathway_table$pathway == id, , drop = FALSE]
  valid <- is.finite(z$NES)
  if (!any(valid)) return(NULL)
  data.frame(pathway = id, max_NES = max(z$NES[valid]),
             min_FDR = if (any(is.finite(z$gsea_FDR))) min(z$gsea_FDR, na.rm = TRUE) else Inf)
}))
extra <- character()
if (nrow(candidate_summary)) {
  candidate_summary <- candidate_summary[candidate_summary$max_NES >= display_min_NES, , drop = FALSE]
  candidate_summary <- candidate_summary[order(candidate_summary$min_FDR, -candidate_summary$max_NES, candidate_summary$pathway), , drop = FALSE]
  extra <- head(candidate_summary$pathway, max_extra_immune_pathways)
}
display_ids <- unique(c(focus_pathways, extra))
display <- merge(expand.grid(pathway = display_ids, condition = condition_levels, stringsAsFactors = FALSE),
                 pathway_table, by = c("pathway", "condition"), all.x = TRUE, sort = FALSE)
missing <- is.na(display$status)
display$status[missing] <- "not_in_GMT"
display$focus <- display$pathway %in% focus_pathways
if (any(missing)) pathway_table <- rbind(pathway_table, display[missing, names(pathway_table), drop = FALSE])
display$display_order <- match(display$pathway, display_ids)
display <- display[order(display$display_order, match(display$condition, condition_levels)), , drop = FALSE]
display$label <- gsub("_", " ", sub("^(HALLMARK_|GOBP_)", "", display$pathway))
display$label <- paste0(ifelse(startsWith(display$pathway, "HALLMARK_"), "H: ", "GO: "), display$label)
display$label <- vapply(display$label, function(z) paste(strwrap(z, width = 42L), collapse = "\n"), character(1))
display$label <- factor(display$label, levels = rev(unique(display$label)))
display$condition <- factor(display$condition, levels = condition_levels)
display$FDR_label <- ifelse(is.finite(display$gsea_FDR), ifelse(display$gsea_FDR < .001, "<0.001", sprintf("%.3f", display$gsea_FDR)), "NA")
delta <- data.frame(pathway = display_ids, stringsAsFactors = FALSE)
for (group in condition_levels) {
  z <- display[display$condition == group, , drop = FALSE]
  delta[[group]] <- z$NES[match(delta$pathway, z$pathway)]
}
delta$delta_NES <- delta$VG161 - delta$Vehicle
delta$label <- display$label[match(delta$pathway, display$pathway)]

##### 绘图与精简保存 #####
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
summary$model_status <- vapply(summary$condition, function(g) ko_results[[g]]$status, character(1))
summary$n_network_genes <- length(network_genes)
summary$n_MDK_outgoing <- vapply(summary$condition, function(g) ko_results[[g]]$n_MDK_outgoing, integer(1))
pathway_table$displayed <- pathway_table$pathway %in% display_ids
for (name in c("summary", "genes")) write.csv(get(name), file.path(outdir, paste0(name, ".csv")), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
write.csv(pathway_table, file.path(outdir, "pathways.csv"), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
finished_at <- Sys.time()
info <- list(started_at = started_at, finished_at = finished_at,
             input_file = path.expand(input_file), input_object = input_object, sample_table = sample_table,
             annotation = list(sample_col = sample_col, celltype_col = celltype_col, tumor_labels = tumor_labels,
                               infection_state_filter = FALSE, species_assumption = "Human host symbols; verify upstream species assignment"),
             parameters = list(target = target_gene, seed = seed, min_cells_per_replicate = min_cells_per_replicate,
                               max_cells_per_replicate = max_cells_per_replicate, cells_per_replicate = cells_per_replicate,
                               detection_fraction = min_detection_fraction, max_network_genes = max_network_genes, min_network_genes = min_network_genes,
                               gene_selection = "Pooled balanced cells, shared detection filter, Seurat vst; MDK forced, then excluded from enrichment",
                               nc_nNet = network_nNet, nc_nCells = network_nCells, nc_nComp = network_nComp, nc_q = network_q,
                               td_K = tensor_K, ma_nDim = manifold_nDim, ma_method = "manifold", dr_direction = FALSE, nCores = n_cores,
                               gene_fdr = gene_fdr, pathway_fdr = pathway_fdr, minSize = pathway_min_size, maxSize = pathway_max_size,
                               scoreType = "pos", gsea_nPermSimple = gsea_nperm_simple, gsea_eps = gsea_eps,
                               focus_pathways = focus_pathways, immune_regex = immune_regex, display_min_NES = display_min_NES,
                               max_extra_immune_pathways = max_extra_immune_pathways, condition_colors = condition_colors, host_prefix = host_prefix,
                               viral_pattern = viral_pattern, viral_feature_ids = viral_feature_ids),
             gene_sets = list(source = gmt_source, version = resolved_gmt_version, paths = gmt_paths, md5 = tools::md5sum(gmt_paths)),
             figures = list(top_genes_per_group = 20L, rank_metric = "distance", volcano_axes = "Z and -log10(model FDR)",
                            max_WT_neighbors_per_group = 8L, network_selection = "Nonzero MDK WT edges to genes below model FDR, ranked by absolute weight",
                            functional_display = "Up to 6 terms per group per method, union shown with both groups; ORA hit fraction > 0, GSEA NES >= 1",
                            order = "Top genes, volcano, WT network, functional ORA and GSEA, immune delta NES, GSEA and ORA"),
             interpretation = c("Pooled group models: no animal-level KO treatment test",
                                "Unsigned distance: enrichment indicates perturbation concentration, not activation/suppression",
                                "Delta NES is descriptive model sensitivity; not an interaction test or expression difference",
                                "scTenifoldKnk FC is not expression fold change; gene model FDR is not biological-replicate DE FDR",
                                "GSEA and ORA each use BH over all testable H+GO:BP terms, separately per group",
                                "Network-gene selection limits coverage and conditions the enrichment null",
                                "Tumor-cell results do not simulate immune receiver-cell responses"),
             warnings = unique(run_warnings), session = capture.output(sessionInfo()))
model_audit <- lapply(ko_results, function(z) z[setdiff(names(z), "genes")])
results <- list(info = info, summary = summary, genes = genes, pathways = pathway_table,
                display = display, delta = delta, model_audit = model_audit, cell_audit = cell_audit,
                layer_audit = layer_audit, viral_exclusion = viral_exclusion, gene_selection = gene_selection,
                network_genes = network_genes, enrichment_universe = universe,
                network_edges = bind_tables(network_edges_by_group))
saveRDS(results, file.path(outdir, "results.rds"), compress = "xz")
writeLines(c("D3 tumor-cell MDK virtual knockout", capture.output(str(info, max.level = 3)),
             "", "MODEL AUDIT", capture.output(str(model_audit)), "", "SESSION", info$session),
           file.path(outdir, "info.txt"), useBytes = TRUE)
# 表格和RDS先写入，绘图可直接从results.rds续跑。
source(file.path(code_dir, "scRNA_MDK_knockout_figures.R"), local = TRUE, encoding = "UTF-8")
write_mdk_plots(results, file.path(outdir, "plots.pdf"))
message("已保存：", normalizePath(outdir, winslash = "/"))
if (any(vapply(ko_results, function(z) z$status != "ok", logical(1)))) warning("部分模型不可评估，详见results.rds$model_audit与info.txt；NA不解释为无效应。")

