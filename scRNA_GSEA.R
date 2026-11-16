##### 清理运行环境 #####
rm(list = ls(all.names = TRUE))
graphics.off()
invisible(gc())

##### 加载依赖 #####
library(Seurat)
library(Matrix)
library(edgeR)
library(statmod)
library(fgsea)
library(BiocParallel)
library(ggplot2)
library(patchwork)

##### 输入路径与参数 #####
# 第一步：肿瘤细胞 Infected - Bystander 的功能变化及 MDK 排序贡献。
# 不运行细胞通讯，不重判感染，不修改或另存大型 Seurat 对象。
input_file <- "~/rawdata/scRNA_virus/virus_3D/virus_3D_anno.RData"
input_object <- "sce"
analysis_label <- "D3"
outdir <- file.path(path.expand("~/OV"), "scRNA_virus", paste0("GSEA_", analysis_label))

sample_col <- "orig.ident"
celltype_col <- "celltype"
state_col <- "infected"
tumor_labels <- c("Epithelial")
state_levels <- c("Bystander", "Infected") # 正 logFC / 正 NES 对应 Infected。

# replicate_id 是真实独立生物样本的阻断单位。下面暂用文库名，运行前核对。
# 技术拆库可填同一 replicate_id；同鼠两侧/不同日期应先分开运行。
# D5 注射侧示例 sample_id: VG161_5D_L1、VG161_5D_L2，勿直接混入 R 侧。
sample_table <- data.frame(
  sample_id = c("VG161_3D_1", "VG161_3D_2"),
  replicate_id = c("VG161_3D_1", "VG161_3D_2"),
  stringsAsFactors = FALSE
)
replicate_unit <- "biological sample; confirm sample_table against animal records"
host_prefix <- ""                  # 已分出人源宿主且无前缀时留空；如 human_ 则填写。
viral_pattern <- "^(VG161|HSV[-_]?1)([-_.]|$)" # 从宿主分析排除带病毒前缀的骨架及载荷。
viral_feature_ids <- character()   # 无前缀/特殊命名时补充精确feature名，如c("UL29", "ICP27")。
target_gene <- "MDK"
min_cells_per_state <- 20L         # 每个独立样本的两状态均须满足；不是生物学阈值。
min_pairs <- 2L                    # 不能用细胞数代替样本重复；两对仍属探索。
filter_min_count <- 10
filter_min_total_count <- 15
seed <- 123L

# 默认从Broad官方只获取两份GMT（合计约4.69 MiB），不调用Zenodo整包下载。
# 文件仅留在当前R临时目录，同一会话重跑时复用通过校验的文件。
# 可选模式："GMT"读取自己指定的两份文件；"msigdbr"使用原在线入口。
gene_set_source_mode <- "Broad"
gene_set_dir <- file.path(tempdir(), "gsea_broad_2026_1")
hallmark_gmt <- ""                 # 仅GMT模式需填写本地Hallmark文件路径。
gobp_gmt <- ""                     # 仅GMT模式需填写本地GO:BP文件路径。
gmt_version <- "MSigDB 2026.1.Hs"
download_timeout <- 180L
download_retries <- 3L
min_set_size <- 15L                # 与本次排序背景取交集后的基因数。
max_set_size <- 500L
gsea_eps <- 1e-10
gsea_simple_n <- 1000L             # fgsea 内部抽样，不是生物样本置换。
fdr_cutoff <- 0.05
remove_target_sensitivity <- TRUE  # 固定其他基因排序，仅删除 MDK 后重跑完整家族。

# 完整检验后才用这些条目组织展示，不据此缩小检验背景或重算更小家族 FDR。
focus_pathways <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "GOBP_RESPONSE_TO_VIRUS", "GOBP_RESPONSE_TO_CYTOKINE",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION", "GOBP_LEUKOCYTE_MIGRATION",
  "GOBP_REGULATION_OF_T_CELL_ACTIVATION", "GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS",
  "GOBP_CELL_ADHESION"
)
overview_per_direction <- 6L        # 每集合、每方向最多展示几项。
max_enrichment_curves <- 8L
state_colors <- c(Bystander = "#426B91", Infected = "#C35B3F")

##### 配置与输入校验 #####
if (!identical(state_levels, c("Bystander", "Infected"))) {
  stop("请保留 state_levels 的顺序；若原标签不同，请先在 metadata 中明确映射。")
}
if (!all(c("sample_id", "replicate_id") %in% names(sample_table)) ||
    !nrow(sample_table) || anyNA(sample_table[, c("sample_id", "replicate_id")]) ||
    any(!nzchar(sample_table$sample_id)) || any(!nzchar(sample_table$replicate_id)) ||
    anyDuplicated(sample_table$sample_id)) stop("sample_table 的样本/独立重复 ID 缺失或重复。")
if (min_pairs < 2L || min_cells_per_state < 1L || min_set_size < 2L ||
    max_set_size < min_set_size || fdr_cutoff <= 0 || fdr_cutoff >= 1 ||
    max_enrichment_curves < 1L) stop("请核对细胞数、重复数、集合大小、FDR 或作图参数。")
if (!gene_set_source_mode %in% c("Broad", "GMT", "msigdbr")) stop("gene_set_source_mode应为Broad、GMT或msigdbr。")
if (gene_set_source_mode == "GMT") {
  if (!nzchar(hallmark_gmt) || !nzchar(gobp_gmt)) stop("GMT模式请同时填写 Hallmark 和 GO:BP 文件路径。")
  if (any(!file.exists(path.expand(c(hallmark_gmt, gobp_gmt))))) stop("找不到离线基因集，请核对hallmark_gmt与gobp_gmt。")
}
if (gene_set_source_mode == "Broad" &&
    (!nzchar(gene_set_dir) || !is.finite(download_timeout) || download_timeout < 1 ||
     !is.finite(download_retries) || download_retries < 1 || download_retries != as.integer(download_retries))) {
  stop("请核对Broad临时目录、下载超时及重试次数。")
}
input_file <- path.expand(input_file)
if (!file.exists(input_file)) stop("找不到输入文件：", input_file)
if (normalizePath(input_file, mustWork = TRUE) %in%
    normalizePath(file.path(path.expand(outdir), c("summary.csv", "genes.csv", "results.rds", "plots.pdf")),
                  mustWork = FALSE)) stop("输入文件不能与输出同名。")
started_at <- Sys.time()
set.seed(seed)
read_env <- new.env()
loaded <- tryCatch(suppressWarnings(load(input_file, envir = read_env)), error = function(e) NULL)
if (is.null(loaded)) {
  sce <- readRDS(input_file)
} else {
  if (!input_object %in% loaded) stop("输入文件中没有对象：", input_object)
  sce <- get(input_object, envir = read_env)
}
rm(read_env)
if (!inherits(sce, "Seurat") || !"RNA" %in% names(sce@assays)) stop("需要含 RNA counts 的 Seurat 对象。")
meta <- sce[[]]
required_cols <- c(sample_col, celltype_col, state_col)
if (!all(required_cols %in% names(meta))) stop("metadata 缺少：", paste(setdiff(required_cols, names(meta)), collapse = ", "))
if (is.null(rownames(meta)) || anyNA(rownames(meta)) || anyDuplicated(rownames(meta))) stop("metadata 条码不唯一。")
sample_ids <- as.character(meta[[sample_col]])
if (!all(sample_table$sample_id %in% sample_ids)) {
  stop("sample_table 与真实 metadata 不符，缺少：", paste(setdiff(sample_table$sample_id, sample_ids), collapse = ", "),
       "。请先查看 table(sce[[sample_col]][, 1])；不能以条码前缀代替核对。")
}
selected <- !is.na(meta[[celltype_col]]) & meta[[celltype_col]] %in% tumor_labels &
  !is.na(sample_ids) & sample_ids %in% sample_table$sample_id &
  !is.na(meta[[state_col]]) & meta[[state_col]] %in% state_levels
meta_use <- meta[selected, , drop = FALSE]
if (!nrow(meta_use)) stop("所选样本中没有带 Infected/Bystander 标签的目标肿瘤细胞。")
meta_use$sample_id <- as.character(meta_use[[sample_col]])
meta_use$replicate_id <- sample_table$replicate_id[match(meta_use$sample_id, sample_table$sample_id)]
meta_use$state <- factor(as.character(meta_use[[state_col]]), levels = state_levels)
message("仅分析配置样本的肿瘤细胞；沿用现有感染标签。请确认 replicate_id 是独立生物样本。")
print(table(meta_use$sample_id, meta_use$state))

##### 显式提取原始 RNA counts，保留稀疏结构 #####
# 只保留纳入细胞所在的层。不同层缺失的 feature 不是观测到的零：
# DE 使用这些层共有的宿主 feature，避免 JoinLayers 填零造成假差异。
rna <- sce[["RNA"]]
if (inherits(rna, "Assay5")) {
  layer_names <- Layers(rna, search = "^counts($|\\.)")
  if (!length(layer_names)) stop("RNA 没有原始 counts 层；不能使用 data/SCT/integrated 代替。")
  layer_cell_names <- lapply(layer_names, function(k) colnames(LayerData(rna, layer = k, fast = FALSE)))
  if (anyDuplicated(unlist(layer_cell_names, use.names = FALSE))) stop("RNA counts 层之间条码重复，请先明确原始计数层。")
  if (!all(rownames(meta_use) %in% unlist(layer_cell_names, use.names = FALSE))) stop("所选细胞缺少 RNA counts，不能补零。")
  take_layers <- vapply(layer_cell_names, function(z) any(z %in% rownames(meta_use)), logical(1))
  layer_names <- layer_names[take_layers]
  count_layers <- lapply(layer_names, function(k) {
    z <- LayerData(rna, layer = k, fast = FALSE)
    z[, intersect(colnames(z), rownames(meta_use)), drop = FALSE]
  })
} else if (inherits(rna, "Assay")) {
  z <- if (packageVersion("SeuratObject") >= "5.0.0") GetAssayData(sce, assay = "RNA", layer = "counts") else
    GetAssayData(sce, assay = "RNA", slot = "counts")
  if (!all(rownames(meta_use) %in% colnames(z))) stop("所选细胞缺少 RNA counts，不能补零。")
  layer_names <- "counts"
  count_layers <- list(z[, rownames(meta_use), drop = FALSE])
  rm(z)
} else stop("不支持的 RNA assay 类型：", paste(class(rna), collapse = ", "))
names(count_layers) <- layer_names
for (k in layer_names) {
  z <- count_layers[[k]]
  if (is.null(rownames(z)) || is.null(colnames(z)) || anyNA(rownames(z)) ||
      anyNA(colnames(z)) || anyDuplicated(rownames(z)) || anyDuplicated(colnames(z))) stop(k, "：feature/条码缺失或重复。")
  if (!inherits(z, "sparseMatrix")) z <- Matrix(z, sparse = TRUE)
  z <- as(z, "dgCMatrix")
  if (any(!is.finite(z@x) | z@x < 0 | abs(z@x - round(z@x)) > 1e-6)) stop(k, "：需要有限、非负整数 RNA counts。")
  count_layers[[k]] <- z
}
# 排除病毒feature发生在DE、TMM和GSEA排序之前；不删除VG161处理样本或宿主抗病毒应答基因。
is_viral_feature <- function(f) {
  symbol <- f
  if (nzchar(host_prefix)) {
    prefixed <- startsWith(f, host_prefix)
    symbol[prefixed] <- substring(f[prefixed], nchar(host_prefix) + 1L)
  }
  by_pattern <- if (nzchar(viral_pattern)) {
    grepl(viral_pattern, f, ignore.case = TRUE) | grepl(viral_pattern, symbol, ignore.case = TRUE)
  } else rep(FALSE, length(f))
  by_pattern | f %in% viral_feature_ids | symbol %in% viral_feature_ids
}
viral_exclusion_audit <- do.call(rbind, lapply(layer_names, function(k) {
  z <- count_layers[[k]]
  ids <- rownames(z)[is_viral_feature(rownames(z))]
  data.frame(layer = rep(k, length(ids)), feature = ids,
             total_UMI = as.numeric(Matrix::rowSums(z[ids, , drop = FALSE])), stringsAsFactors = FALSE)
}))
message("排除已识别病毒feature：", length(unique(viral_exclusion_audit$feature)), " 个；保留所选样本的宿主基因。")
host_features <- lapply(count_layers, function(z) {
  f <- rownames(z)
  if (nzchar(host_prefix)) f <- f[startsWith(f, host_prefix)]
  f[!is_viral_feature(f)]
})
common_features <- Reduce(intersect, host_features)
if (length(common_features) < 100L) stop("共同宿主 feature 少于100，请核对物种前缀/计数层。")
symbols <- if (nzchar(host_prefix)) substring(common_features, nchar(host_prefix) + 1L) else common_features
if (anyNA(symbols) || any(!nzchar(symbols)) || anyDuplicated(symbols)) stop("去前缀后 gene symbol 不唯一；请先提供确定的 feature 到 symbol 映射。")
if (mean(grepl("^ENSG[0-9]", symbols)) > 0.5) stop("当前主要是 Ensembl ID；请先按固定注释映射为人 gene symbol，再运行。")
layer_audit <- data.frame(layer = layer_names, n_cells = vapply(count_layers, ncol, integer(1)),
                          n_host_features = lengths(host_features), n_shared_features = length(common_features),
                          n_viral_excluded = vapply(count_layers, function(z) sum(is_viral_feature(rownames(z))), integer(1)),
                          target_present = vapply(host_features, function(f) {
                            (if (nzchar(host_prefix)) paste0(host_prefix, target_gene) else target_gene) %in% f
                          }, logical(1)), stringsAsFactors = FALSE)
if (any(lengths(host_features) != length(common_features))) warning("计数层 feature 集不同：DE 仅用共有宿主 feature，详情见 results.rds 的 layer_audit。")

##### 按独立样本与感染状态汇总 pseudobulk #####
replicate_ids <- unique(sample_table$replicate_id)
pb_meta <- expand.grid(state = state_levels, replicate_id = replicate_ids, KEEP.OUT.ATTRS = FALSE,
                       stringsAsFactors = FALSE)
pb_meta$state <- factor(pb_meta$state, levels = state_levels)
pb_meta$pb_id <- paste0("PB", seq_len(nrow(pb_meta)))
pb_meta$n_cells <- 0L
pb_meta$host_UMI <- 0
pb_meta$target_UMI <- 0
pb_meta$target_detected_n <- 0L
pb_meta$target_log_sum <- 0
pb_meta$target_assessable <- TRUE
pb_counts <- matrix(0, nrow = length(common_features), ncol = nrow(pb_meta),
                    dimnames = list(symbols, pb_meta$pb_id))
cell_audit <- vector("list", length(count_layers))
for (i in seq_along(count_layers)) {
  z <- count_layers[[i]]
  m <- meta_use[colnames(z), , drop = FALSE]
  # 使用所有纳入层共有的宿主基因作为统一分母，避免 feature 缺失改变跨层 CP10K。
  host_umi <- Matrix::colSums(z[common_features, , drop = FALSE])
  valid <- is.finite(host_umi) & host_umi > 0
  cell_audit[[i]] <- data.frame(sample_id = m$sample_id, replicate_id = m$replicate_id,
                                state = m$state, valid_host = valid, stringsAsFactors = FALSE)
  if (!any(valid)) next
  z <- z[, valid, drop = FALSE]
  m <- m[valid, , drop = FALSE]
  host_umi <- host_umi[valid]
  pb_index <- match(paste(m$replicate_id, m$state, sep = "\r"),
                    paste(pb_meta$replicate_id, pb_meta$state, sep = "\r"))
  membership <- sparseMatrix(i = seq_along(pb_index), j = pb_index, x = 1,
                             dims = c(length(pb_index), nrow(pb_meta)))
  # 稀疏基因×细胞乘稀疏细胞×样本，只有小型 pseudobulk 转为普通矩阵。
  pb_counts <- pb_counts + as.matrix(z[common_features, , drop = FALSE] %*% membership)
  tab <- tabulate(pb_index, nbins = nrow(pb_meta))
  pb_meta$n_cells <- pb_meta$n_cells + tab
  pb_meta$host_UMI <- pb_meta$host_UMI + as.numeric(crossprod(membership, host_umi))
  target_feature <- if (nzchar(host_prefix)) paste0(host_prefix, target_gene) else target_gene
  if (target_feature %in% rownames(z)) {
    target_counts <- as.numeric(z[target_feature, ])
    pb_meta$target_UMI <- pb_meta$target_UMI + as.numeric(crossprod(membership, target_counts))
    pb_meta$target_detected_n <- pb_meta$target_detected_n + as.numeric(crossprod(membership, target_counts > 0))
    pb_meta$target_log_sum <- pb_meta$target_log_sum + as.numeric(crossprod(membership, log1p(1e4 * target_counts / host_umi)))
  } else pb_meta$target_assessable[tab > 0] <- FALSE
}
cell_audit <- do.call(rbind, cell_audit)
sample_audit <- aggregate(cbind(n_selected = rep(1L, nrow(cell_audit)), n_valid_host = as.integer(cell_audit$valid_host)),
                          by = cell_audit[, c("sample_id", "replicate_id", "state")], FUN = sum)
rm(count_layers, z, membership, cell_audit, rna, sce, meta, meta_use)
invisible(gc()) # 后续只需小型 pseudobulk，及时释放大型对象与细胞矩阵。
pb_meta$library_size <- colSums(pb_counts)
pb_meta$target_assessable <- pb_meta$target_assessable & pb_meta$n_cells > 0
for (nm in c("target_UMI", "target_detected_n", "target_log_sum")) {
  pb_meta[[nm]][!pb_meta$target_assessable] <- NA_real_
}
pb_meta$target_mean_logCP10K <- ifelse(pb_meta$n_cells > 0 & pb_meta$target_assessable,
                                       pb_meta$target_log_sum / pb_meta$n_cells, NA_real_)
pb_meta$target_detected_fraction <- ifelse(pb_meta$n_cells > 0 & pb_meta$target_assessable,
                                           pb_meta$target_detected_n / pb_meta$n_cells, NA_real_)
pb_meta$passes_qc <- pb_meta$n_cells >= min_cells_per_state & pb_meta$library_size > 0
pair_ok <- vapply(replicate_ids, function(id) all(pb_meta$passes_qc[pb_meta$replicate_id == id]), logical(1))
included_pairs <- replicate_ids[pair_ok]
pb_meta$included <- pb_meta$replicate_id %in% included_pairs
if (length(included_pairs) < min_pairs) {
  print(pb_meta[, c("replicate_id", "state", "n_cells", "library_size", "passes_qc")])
  stop("同时满足两状态条件的独立样本仅 ", length(included_pairs), " 个；不能改用逐细胞 P 值替代。请核对样本和阈值。")
}
if (length(included_pairs) < 3L) warning("只有两对独立样本：样本级效应及GSEA仍属探索，需要独立验证。")
if (any(!pair_ok)) warning("部分样本因一侧细胞数/计数不足而整对排除；详情见 pb_meta。")
pb_use <- droplevels(pb_meta[pb_meta$included, , drop = FALSE])
pb_use$block <- factor(pb_use$replicate_id, levels = included_pairs)
pb_use$state <- factor(pb_use$state, levels = state_levels)
counts_use <- pb_counts[, pb_use$pb_id, drop = FALSE]

##### 配对 edgeR 差异分析与全基因排序 #####
design <- model.matrix(~ block + state, data = pb_use)
rownames(design) <- pb_use$pb_id
if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) stop("配对设计不满秩或没有剩余自由度。")
y <- DGEList(counts = counts_use)
# 对配对设计，表达过滤按待比较的 state，而非阻断样本项。
keep_gene <- filterByExpr(y, group = pb_use$state,
                          min.count = filter_min_count, min.total.count = filter_min_total_count)
if (sum(keep_gene) < 100L) stop("通过表达过滤的宿主基因少于100，无法进行可靠的排序富集。")
y <- y[keep_gene, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
test <- glmQLFTest(fit, coef = "stateInfected")
de <- topTags(test, n = Inf, sort.by = "none")$table
de$gene <- rownames(de)
de$rank_stat <- sign(de$logFC) * sqrt(pmax(de$F, 0))
finite_rank <- is.finite(de$rank_stat)
if (sum(finite_rank) < 100L || all(de$rank_stat[finite_rank] == 0)) stop("有限、非全零的排序统计量不足。")
# 并列值按 gene symbol 固定顺序，不加随机噪声制造精度。
rank_order <- order(-de$rank_stat[finite_rank], de$gene[finite_rank], method = "radix")
ranks <- setNames(de$rank_stat[finite_rank][rank_order], de$gene[finite_rank][rank_order])
if (any(is_viral_feature(names(ranks)))) stop("GSEA排序仍含已识别病毒feature，请核对排除规则。")
tie_fraction <- mean(duplicated(unname(ranks)))
if (tie_fraction > 0.1) warning("排序并列比例超过10%，需结合记录谨慎解释基因位置和 leading edge。")
de_all <- data.frame(gene = rownames(counts_use), total_count = rowSums(counts_use),
                     passed_filter = keep_gene, stringsAsFactors = FALSE)
for (nm in c("logFC", "logCPM", "F", "PValue", "FDR", "rank_stat")) de_all[[nm]] <- de[[nm]][match(de_all$gene, de$gene)]
de_all$rank_position <- match(de_all$gene, names(ranks))
de_all$rank_percentile <- (de_all$rank_position - 1) / (length(ranks) - 1)
de_all$is_target <- de_all$gene == target_gene
target_status <- if (!target_gene %in% symbols) "not_in_shared_features" else if (!target_gene %in% de$gene) {
  "expression_filtered"
} else if (!target_gene %in% names(ranks)) "non_finite_statistic" else "ranked"
target_de <- de_all[de_all$gene == target_gene, , drop = FALSE]
message(target_gene, " status: ", target_status)
print(target_de)

##### 完整 Hallmark 与 GO:BP 基因集 #####
download_broad_gmt <- function(url, destination, expected_md5, expected_sets) {
  if (file.exists(destination) && identical(unname(tools::md5sum(destination)), expected_md5)) {
    message("复用本次R会话已校验的基因集：", basename(destination))
    return(invisible(destination))
  }
  old_options <- options(timeout = max(download_timeout, getOption("timeout", 60)))
  on.exit(options(old_options), add = TRUE)
  last_error <- ""
  for (attempt in seq_len(download_retries)) {
    temporary <- tempfile(pattern = "download_", tmpdir = dirname(destination), fileext = ".gmt")
    success <- tryCatch({
      code <- utils::download.file(url, destfile = temporary, method = "libcurl", mode = "wb", quiet = FALSE)
      if (!identical(as.integer(code), 0L)) stop("下载返回非零状态。")
      if (!identical(unname(tools::md5sum(temporary)), expected_md5)) stop("文件校验失败，可能是截断下载。")
      fields <- strsplit(readLines(temporary, warn = FALSE, encoding = "UTF-8"), "\t", fixed = TRUE)
      if (length(fields) != expected_sets || any(lengths(fields) < 3L)) stop("GMT集合数或格式不符。")
      if (anyDuplicated(vapply(fields, `[`, character(1), 1L))) stop("GMT集合名称重复。")
      if (!file.copy(temporary, destination, overwrite = TRUE)) stop("无法写入基因集临时目录。")
      TRUE
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      FALSE
    })
    unlink(temporary) # 仅删除本函数刚创建的中转文件，不接受半个下载作为GSEA输入。
    if (success) return(invisible(destination))
    message("第", attempt, "次下载未完成：", last_error)
    if (attempt < download_retries) Sys.sleep(1)
  }
  stop("Broad基因集下载失败：", last_error,
       "。已生成的DE/ranks仍在当前会话，可在连接恢复后从本节重试；不要重跑清理环境。", call. = FALSE)
}
if (gene_set_source_mode == "Broad") {
  broad_base <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/"
  broad_files <- data.frame(
    file = c("h.all.v2026.1.Hs.symbols.gmt", "c5.go.bp.v2026.1.Hs.symbols.gmt"),
    expected_sets = c(50L, 7538L),
    expected_md5 = c("367eec875967c2cfbf664a1a065b7b8d", "c1903b564ead3b31260418bdd0d22b10"),
    stringsAsFactors = FALSE)
  gene_set_dir <- path.expand(gene_set_dir)
  dir.create(gene_set_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(gene_set_dir)) stop("无法创建基因集临时目录：", gene_set_dir)
  broad_files$url <- paste0(broad_base, broad_files$file)
  broad_files$local <- file.path(gene_set_dir, broad_files$file)
  for (i in seq_len(nrow(broad_files))) {
    download_broad_gmt(broad_files$url[i], broad_files$local[i], broad_files$expected_md5[i], broad_files$expected_sets[i])
  }
  hallmark_gmt <- broad_files$local[1L]
  gobp_gmt <- broad_files$local[2L]
  gmt_version <- "MSigDB 2026.1.Hs" # 与固定URL及校验值共同锁定版本。
  message("Broad基因集已就绪：Hallmark 50项、GO:BP 7538项，约4.69 MiB临时文件。")
}
read_gmt <- function(path, collection) {
  path <- path.expand(path)
  if (!file.exists(path)) stop("找不到 GMT：", path)
  fields <- strsplit(readLines(path, warn = FALSE, encoding = "UTF-8"), "\t", fixed = TRUE)
  if (!length(fields) || any(lengths(fields) < 3L)) stop("GMT 应逐行包含名称、说明、至少一个gene symbol。")
  nms <- vapply(fields, `[`, character(1), 1L)
  if (any(!nzchar(nms)) || anyDuplicated(nms)) stop("GMT 基因集名称缺失或重复。")
  sets <- lapply(fields, function(z) unique(z[-c(1L, 2L)][nzchar(z[-c(1L, 2L)])]))
  names(sets) <- nms
  list(sets = sets, meta = data.frame(pathway = nms, collection = collection,
                                      description = vapply(fields, `[`, character(1), 2L),
                                      db_version = gmt_version, stringsAsFactors = FALSE))
}
if (gene_set_source_mode %in% c("Broad", "GMT")) {
  if (!nzchar(gmt_version)) stop("离线 GMT 请填写 gmt_version，并确认两份文件来自同一版本。")
  h <- read_gmt(hallmark_gmt, "H")
  go <- read_gmt(gobp_gmt, "GO:BP")
  pathways <- c(h$sets, go$sets)
  pathway_meta <- rbind(h$meta, go$meta)
  gene_set_source <- list(mode = gene_set_source_mode, paths = path.expand(c(hallmark_gmt, gobp_gmt)),
                          md5 = tools::md5sum(path.expand(c(hallmark_gmt, gobp_gmt))), version = gmt_version)
  if (gene_set_source_mode == "Broad") {
    gene_set_source$urls <- broad_files$url
    gene_set_source$bytes <- file.info(broad_files$local)$size
    gene_set_source$download_timeout <- download_timeout
    gene_set_source$download_retries <- download_retries
  }
} else {
  library(msigdbr)
  if (packageVersion("msigdbr") < "10.0.0") stop("请更新 msigdbr，使用 collection/subcollection 的现代接口，或改用离线 GMT。")
  catalog <- tryCatch(msigdbr_collections(db_species = "HS"), error = function(e) {
    stop("在线基因集读取失败（", conditionMessage(e), "）。请使用配套离线GMT，并填写gene_set_dir；这一步不需要重算前面的DE。", call. = FALSE)
  })
  if (!all(c("gs_collection", "gs_subcollection") %in% names(catalog)) ||
      !any(catalog$gs_collection == "H") ||
      !any(catalog$gs_collection == "C5" & catalog$gs_subcollection == "GO:BP")) stop("当前数据库缺少人 Hallmark 或 GO:BP。")
  online_sets <- tryCatch(list(
    h = msigdbr(db_species = "HS", species = "Homo sapiens", collection = "H"),
    go = msigdbr(db_species = "HS", species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
  ), error = function(e) {
    stop("在线基因集读取失败（", conditionMessage(e), "）。请改用配套离线GMT并从基因集分节继续。", call. = FALSE)
  })
  h <- online_sets$h
  go <- online_sets$go
  sets_df <- rbind(as.data.frame(h), as.data.frame(go))
  sets_df$analysis_collection <- ifelse(sets_df$gs_collection == "H", "H", "GO:BP")
  pathways <- lapply(split(sets_df$gene_symbol, sets_df$gs_name), function(z) unique(z[!is.na(z) & nzchar(z)]))
  pathway_meta <- unique(sets_df[, c("gs_name", "analysis_collection", "gs_description", "db_version")])
  names(pathway_meta) <- c("pathway", "collection", "description", "db_version")
  gene_set_source <- list(mode = "msigdbr", version = unique(sets_df$db_version), catalog = catalog)
  rm(sets_df)
}
if (anyDuplicated(names(pathways)) || anyDuplicated(pathway_meta$pathway)) stop("跨集合存在重名，请核对基因集来源。")
pathway_meta <- pathway_meta[match(names(pathways), pathway_meta$pathway), , drop = FALSE]
if (anyNA(pathway_meta$collection)) stop("部分基因集缺少集合注释。")
matched_genes <- intersect(names(ranks), unique(unlist(pathways, use.names = FALSE)))
mapping_fraction <- length(matched_genes) / length(ranks)
if (length(matched_genes) < 100L || mapping_fraction < 0.1) stop("人基因集匹配过低，请核对gene symbol、物种及基因集版本。")
if (mapping_fraction < 0.5) warning("少于一半排序基因能匹配所用基因集，详情已记录；勿通过转大写强制混合物种。")
rm(h, go)

##### 双向 GSEA 与集合完整性记录 #####
run_gsea <- function(stats, seed_use) {
  overlap <- lapply(pathways, intersect, y = names(stats))
  sizes <- lengths(overlap)
  hit_weight <- vapply(overlap, function(g) sum(abs(stats[g])), numeric(1))
  status <- ifelse(sizes >= length(stats), "covers_rank_universe", ifelse(sizes < min_set_size, "below_min_size",
                                                                          ifelse(sizes > max_set_size, "above_max_size",
                                                                                 ifelse(hit_weight == 0, "zero_hit_weight", "tested"))))
  result <- pathway_meta
  result$database_size <- lengths(pathways)
  result$size <- sizes
  result$test_status <- status
  for (nm in c("ES", "NES", "pval", "log2err", "padj_family", "padj_global")) result[[nm]] <- NA_real_
  result$leadingEdge <- I(rep(list(character()), nrow(result)))
  eligible <- names(pathways)[status == "tested"]
  if (length(eligible)) {
    set.seed(seed_use)
    res <- as.data.frame(fgseaMultilevel(pathways = pathways[eligible], stats = stats,
                                         minSize = min_set_size, maxSize = max_set_size,
                                         eps = gsea_eps, scoreType = "std", gseaParam = 1,
                                         nPermSimple = gsea_simple_n, BPPARAM = SerialParam()))
    ix <- match(res$pathway, result$pathway)
    for (nm in c("ES", "NES", "pval", "log2err")) result[[nm]][ix] <- res[[nm]]
    result$leadingEdge[ix] <- res$leadingEdge
    # 校正家族是完整、满足大小和权重条件的集合；保留计算失败的NA，不当成阴性。
    for (family in unique(result$collection)) {
      j <- which(result$collection == family & result$test_status == "tested")
      result$padj_family[j] <- p.adjust(result$pval[j], method = "BH", n = length(j))
    }
    j <- which(result$test_status == "tested")
    result$padj_global[j] <- p.adjust(result$pval[j], method = "BH", n = length(j))
    result$test_status[j[!is.finite(result$pval[j]) | !is.finite(result$NES[j])]] <- "fgsea_unresolved"
  }
  result
}
gsea <- run_gsea(ranks, seed)
gsea$target_in_database_set <- vapply(pathways, function(z) target_gene %in% z, logical(1))
gsea$target_in_rank <- target_gene %in% names(ranks)
gsea$target_in_leading_edge <- vapply(seq_len(nrow(gsea)), function(i) {
  if (target_status != "ranked" || gsea$test_status[i] != "tested") return(NA)
  target_gene %in% gsea$leadingEdge[[i]]
}, logical(1))
# hit_weight 是加权游走中该排序项占集合命中权重的比例，不是生物学效应份额。
gsea$target_hit_weight <- vapply(seq_along(pathways), function(i) {
  if (target_status != "ranked") return(NA_real_)
  gs <- intersect(pathways[[i]], names(ranks))
  if (!target_gene %in% gs) return(0)
  den <- sum(abs(ranks[gs]))
  if (den > 0) abs(ranks[target_gene]) / den else NA_real_
}, numeric(1))
gsea$focused <- gsea$pathway %in% focus_pathways
gsea$significant_family <- !is.na(gsea$padj_family) & gsea$padj_family < fdr_cutoff
if (!any(gsea$test_status == "tested")) warning("没有可评估的GSEA集合；仍保存完整状态和MDK记录。")

##### MDK 单个排序项去除敏感性 #####
# 不是虚拟敲除：不重拟合DE、不改变其他基因rank。完整H+GO:BP重跑，防止缩小FDR家族。
gsea_without_target <- NULL
if (remove_target_sensitivity && target_status == "ranked") {
  ranks_without_target <- ranks[names(ranks) != target_gene]
  gsea_without_target <- run_gsea(ranks_without_target, seed)
}
summary_table <- gsea
summary_table$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
summary_table$target_gene <- target_gene
summary_table$target_status <- target_status
for (nm in c("size", "ES", "NES", "pval", "padj_family", "padj_global")) {
  summary_table[[paste0(nm, "_without_target")]] <- if (is.null(gsea_without_target)) NA_real_ else gsea_without_target[[nm]]
}
summary_table$status_without_target <- if (is.null(gsea_without_target)) {
  if (target_status != "ranked") "not_assessed_target_not_ranked" else "disabled"
} else gsea_without_target$test_status
summary_table$leadingEdge_without_target <- if (is.null(gsea_without_target)) NA_character_ else
  vapply(gsea_without_target$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
summary_table$direction_same_without_target <- ifelse(is.finite(summary_table$NES) & is.finite(summary_table$NES_without_target),
                                                      sign(summary_table$NES) == sign(summary_table$NES_without_target), NA)
summary_table$significant_without_target <- ifelse(is.finite(summary_table$padj_family_without_target),
                                                   summary_table$padj_family_without_target < fdr_cutoff, NA)
target_pathways <- summary_table[summary_table$target_in_database_set, , drop = FALSE]
missing_focus <- setdiff(focus_pathways, summary_table$pathway)

##### 选择图中展示条目，不改变全结果和多重校正 #####
eligible <- gsea[gsea$test_status == "tested" & is.finite(gsea$NES), , drop = FALSE]
nes_limit <- if (nrow(eligible)) max(1, abs(eligible$NES)) else 1
overview_names <- character()
for (family in c("H", "GO:BP")) for (direction in c(-1, 1)) {
  z <- eligible[eligible$collection == family & sign(eligible$NES) == direction, , drop = FALSE]
  z <- z[order(z$padj_family, -abs(z$NES), z$pathway, na.last = TRUE), , drop = FALSE]
  overview_names <- c(overview_names, head(z$pathway, overview_per_direction))
}
mdk_display <- eligible[eligible$target_in_database_set, , drop = FALSE]
mdk_display <- mdk_display[order(mdk_display$padj_family, -abs(mdk_display$NES), mdk_display$pathway), , drop = FALSE]
curve_names <- unique(c(head(mdk_display$pathway, 3L), focus_pathways[focus_pathways %in% eligible$pathway], overview_names))
curve_names <- head(curve_names, max_enrichment_curves)

##### 绘图：样本、通路概览、MDK 核心位置与敏感性 #####
theme_research <- theme_bw(base_size = 10, base_family = "sans") +
  theme(panel.grid = element_blank(), panel.border = element_rect(colour = "#454545", linewidth = 0.4),
        axis.line = element_blank(), axis.ticks = element_line(colour = "#454545", linewidth = 0.35),
        strip.background = element_rect(fill = "#E7EEF4", colour = "#454545", linewidth = 0.35),
        strip.text = element_text(face = "bold", size = 10), legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        legend.position = "right")
short_name <- function(z, width = 65L) {
  z <- gsub("^(HALLMARK_|GOBP_)", "", z)
  vapply(gsub("_", " ", z), function(x) paste(strwrap(x, width), collapse = "\n"), character(1))
}
empty_plot <- function(title, text) {
  ggplot() + annotate("text", x = 0, y = 0, label = text, size = 4, hjust = 0.5) +
    xlim(-1, 1) + ylim(-1, 1) + labs(title = title) + theme_void(base_family = "sans")
}
plot_pathways <- function(nms, title) {
  z <- eligible[match(intersect(nms, eligible$pathway), eligible$pathway), , drop = FALSE]
  if (!nrow(z)) return(empty_plot(title, "No assessable pathways in this display selection."))
  z <- z[order(z$collection, z$NES), , drop = FALSE]
  z$label <- factor(z$pathway, levels = z$pathway)
  ggplot(z, aes(NES, label, fill = NES)) + geom_col(width = 0.7, colour = "#454545", linewidth = 0.2,
                                                    aes(alpha = significant_family)) + geom_vline(xintercept = 0, colour = "#454545", linewidth = 0.3) +
    facet_grid(collection ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(x) short_name(x)) + scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1),
                                                                              labels = c(`FALSE` = paste0("FDR >= ", fdr_cutoff), `TRUE` = paste0("FDR < ", fdr_cutoff))) +
    scale_fill_gradient2(low = "#426B91", mid = "#F7F7F5", high = "#C35B3F", midpoint = 0,
                         limits = c(-nes_limit, nes_limit)) +
    labs(title = title, x = "NES", y = NULL) +
    theme_research
}

# 图内统计图例：左列为指标，右列为本通路的结果；MDK线型与曲线中的位置线一致。
# 方法、指标定义和多重校正口径写在README，不作为PDF底部caption。
enrichment_key <- function(nes, fdr, target_label, member, corner) {
  fmt <- function(x) if (is.finite(x)) formatC(x, digits = 3, format = "g") else "NA"
  right <- grepl("right", corner)
  top <- grepl("top", corner)
  key_colour <- if (member) "#C35B3F" else "#6D6F75"
  children <- list(
    grid::rectGrob(gp = grid::gpar(fill = "white", col = "#D3D8DD", lwd = 0.6)),
    grid::textGrob(c("NES", "Family FDR", target_gene), x = 0.06,
                   y = c(0.81, 0.55, 0.23), just = "left",
                   gp = grid::gpar(fontfamily = "sans", fontsize = 9, col = "#555B61")),
    grid::textGrob(c(fmt(nes), fmt(fdr)), x = 0.94, y = c(0.81, 0.55), just = "right",
                   gp = grid::gpar(fontfamily = "sans", fontsize = 10, fontface = "bold", col = "#252A30")),
    grid::textGrob(target_label, x = 0.94, y = 0.23, just = "right",
                   gp = grid::gpar(fontfamily = "sans", fontsize = 8.5, lineheight = 1.05, col = key_colour))
  )
  if (member) children <- c(children, list(
    grid::segmentsGrob(x0 = 0.22, x1 = 0.32, y0 = 0.23, y1 = 0.23,
                       gp = grid::gpar(col = key_colour, lwd = 1, lty = "dashed"))))
  # annotation_custom会设置最外层viewport；额外包一层，保留图例自身的毫米尺寸和四角定位。
  grid::grobTree(grid::grobTree(children = do.call(grid::gList, children),
                                vp = grid::viewport(x = if (right) 0.98 else 0.02, y = if (top) 0.96 else 0.04,
                                                    width = grid::unit(66, "mm"), height = grid::unit(23, "mm"),
                                                    just = c(if (right) "right" else "left", if (top) "top" else "bottom"))))
}

enrichment_corner <- function(walk, target_position = NA_integer_) {
  # 在四角中优先选曲线经过最少的位置，并避开MDK竖线；白底图例保持数值清晰。
  x <- (walk$position - min(walk$position)) / diff(range(walk$position))
  yr <- range(c(0, walk$ES))
  y <- (walk$ES - yr[1]) / max(diff(yr), .Machine$double.eps)
  tx <- (target_position - min(walk$position)) / diff(range(walk$position))
  corners <- c("top right", "top left", "bottom right", "bottom left")
  score <- vapply(corners, function(corner) {
    x_inside <- if (grepl("right", corner)) x > 0.64 else x < 0.36
    y_inside <- if (grepl("top", corner)) y > 0.50 else y < 0.50
    target_inside <- is.finite(tx) && if (grepl("right", corner)) tx > 0.64 else tx < 0.36
    mean(x_inside & y_inside) + 0.5 * target_inside
  }, numeric(1))
  if (min(score) == 0) return(list(corner = corners[which.min(score)], reserve = FALSE))
  # 四角均有曲线时，额外保留上方空间，而不是把最少遮挡当作没有遮挡。
  list(corner = if (is.finite(tx) && tx > 0.64) "top left" else "top right", reserve = TRUE)
}

enrichment_plot <- function(nm) {
  gs <- intersect(pathways[[nm]], names(ranks))
  positions <- which(names(ranks) %in% gs)
  n <- length(ranks)
  weight <- abs(ranks[positions])
  step <- rep(-1 / (n - length(positions)), n)
  step[positions] <- weight / sum(weight)
  walk <- data.frame(position = seq_len(n), ES = cumsum(step), rank_stat = unname(ranks))
  row <- gsea[match(nm, gsea$pathway), , drop = FALSE]
  mdkr <- match(target_gene, names(ranks))
  member <- target_gene %in% gs
  target_label <- if (!target_gene %in% names(ranks)) "Not ranked" else if (!member) "Not a member" else
    if (is.na(row$target_in_leading_edge)) "Member\nLeading edge: NA" else
      if (isTRUE(row$target_in_leading_edge)) "Member\nIn leading edge" else "Member\nOutside leading edge"
  key_position <- enrichment_corner(walk, if (member) mdkr else NA_integer_)
  p1 <- ggplot(walk, aes(position, ES)) + geom_line(colour = "#426B91", linewidth = 0.55) +
    geom_hline(yintercept = 0, colour = "#9E9E9E", linewidth = 0.25) +
    labs(title = short_name(nm, 80L), x = NULL, y = "Running ES") + theme_research +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  hits <- data.frame(position = positions)
  p2 <- ggplot(hits, aes(position)) + geom_segment(aes(xend = position, y = 0, yend = 1), linewidth = 0.25, colour = "#656A70") +
    xlim(1, n) + labs(x = NULL, y = NULL) + theme_void()
  if (member) {
    p1 <- p1 + geom_vline(xintercept = mdkr, colour = "#C35B3F", linewidth = 0.4, linetype = "dashed")
    p2 <- p2 + geom_segment(data = data.frame(position = mdkr), aes(xend = position, y = 0, yend = 1),
                            colour = "#C35B3F", linewidth = 1)
  }
  if (key_position$reserve) p1 <- p1 + scale_y_continuous(expand = expansion(mult = c(0.05, 0.75)))
  p1 <- p1 + annotation_custom(enrichment_key(row$NES, row$padj_family, target_label, member, key_position$corner))
  p3 <- ggplot(walk, aes(position, rank_stat)) + geom_area(fill = "#B8BDC3", colour = NA) +
    labs(x = "Gene rank (left: Infected; right: Bystander)", y = "Rank statistic") + theme_research +
    theme(axis.text.y = element_text(size = 7), axis.title.y = element_text(size = 8))
  p1 / p2 / p3 + plot_layout(heights = c(3, 0.35, 1))
}
dir.create(path.expand(outdir), recursive = TRUE, showWarnings = FALSE)
outdir <- path.expand(outdir)
write.csv(de_all, file.path(outdir, "genes.csv"), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
pdf(file.path(outdir, "plots.pdf"), width = 12, height = 10, onefile = TRUE, useDingbats = FALSE)
tryCatch({
  p_counts <- ggplot(pb_meta, aes(state, n_cells, fill = state)) + geom_col(width = 0.65, colour = "#454545", linewidth = 0.35) +
    geom_hline(aes(yintercept = min_cells_per_state, linetype = "Minimum"), colour = "#656A70", linewidth = 0.3) +
    facet_wrap(~ replicate_id, ncol = 3) + scale_fill_manual(values = state_colors) +
    scale_linetype_manual(values = c(Minimum = "dashed"), labels = paste0("Minimum: ", min_cells_per_state, " cells")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    guides(fill = "none") + labs(title = "Tumor cells", x = NULL, y = "Cells with host counts") +
    theme_research +
    (if (packageVersion("ggplot2") >= "3.5.0") theme(legend.position = "inside", legend.position.inside = c(0.98, 0.98)) else
      theme(legend.position = c(0.98, 0.98))) + theme(legend.justification = c(1, 1),
                                                      legend.background = element_rect(fill = "white", colour = NA), legend.margin = margin(2, 4, 2, 4),
                                                      legend.key.height = grid::unit(3, "mm"))
  target_plot_data <- rbind(data.frame(pb_use[, c("replicate_id", "state")], metric = "MDK expression", value = pb_use$target_mean_logCP10K),
                            data.frame(pb_use[, c("replicate_id", "state")], metric = "MDK detected (%)", value = 100 * pb_use$target_detected_fraction))
  target_plot_data$metric <- sub("MDK", target_gene, target_plot_data$metric, fixed = TRUE)
  line_key <- interaction(target_plot_data$replicate_id, target_plot_data$metric, drop = TRUE)
  n_valid_points <- ave(as.integer(is.finite(target_plot_data$value)), line_key, FUN = sum)
  paired_plot_data <- target_plot_data[n_valid_points == length(state_levels) & is.finite(target_plot_data$value), , drop = FALSE]
  p_target <- ggplot(target_plot_data, aes(state, value, group = replicate_id)) +
    geom_line(data = paired_plot_data, colour = "#C4C7CA", linewidth = 0.25, na.rm = TRUE) +
    geom_point(aes(colour = state), size = 2, na.rm = TRUE) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) + scale_colour_manual(values = state_colors) +
    labs(title = target_gene, x = NULL, y = NULL) +
    theme_research + theme(legend.position = "none")
  print(p_counts / p_target + plot_layout(heights = c(1.2, 1)))
  print(plot_pathways(overview_names, "Pathway overview"))
  print(plot_pathways(focus_pathways, "Immune and tissue programs"))
  if (nrow(mdk_display)) print(plot_pathways(head(mdk_display$pathway, 14L), paste0(target_gene, "-annotated programs"))) else
    print(empty_plot(paste0(target_gene, "-annotated programs"), "No assessable target-containing sets.\nCheck target status and the complete results."))
  if (nrow(mdk_display) && !is.null(gsea_without_target)) {
    z <- summary_table[match(head(mdk_display$pathway, 12L), summary_table$pathway), , drop = FALSE]
    z$label <- factor(z$pathway, levels = rev(z$pathway))
    z$unassessed_label <- ifelse(z$status_without_target == "below_min_size", "Below size cutoff", "Not assessable")
    p_sensitivity <- ggplot(z, aes(y = label)) +
      geom_segment(aes(x = NES, xend = NES_without_target, yend = label), colour = "#C4C7CA", linewidth = 0.4, na.rm = TRUE) +
      geom_point(aes(x = NES, shape = "Full rank"), colour = "#426B91", size = 2.3) +
      geom_point(aes(x = NES_without_target, shape = "Without target"), colour = "#C35B3F", size = 2.3, na.rm = TRUE) +
      geom_text(data = z[!is.finite(z$NES_without_target), , drop = FALSE],
                aes(x = NES, label = unassessed_label), nudge_y = 0.2, hjust = 1, size = 2.6, colour = "#6D6F75") +
      scale_shape_manual(values = c("Full rank" = 16, "Without target" = 1)) +
      scale_y_discrete(labels = function(x) short_name(x)) +
      labs(title = paste0(target_gene, " rank-item sensitivity"), x = "NES", y = NULL) +
      theme_research
    print(p_sensitivity)
  }
  if (length(curve_names)) for (i in seq(1L, length(curve_names), by = 2L)) {
    nms <- curve_names[i:min(i + 1L, length(curve_names))]
    ps <- lapply(nms, enrichment_plot)
    print(wrap_plots(ps, ncol = 1))
  }
}, finally = dev.off())

##### 保存精简结果与运行记录 #####
packages <- c("Seurat", "SeuratObject", "Matrix", "edgeR", "limma", "statmod", "fgsea", "BiocParallel", "ggplot2", "patchwork")
if (gene_set_source_mode == "msigdbr") packages <- c(packages, "msigdbr")
info <- list(
  started_at = as.character(started_at), finished_at = as.character(Sys.time()),
  input_file = normalizePath(input_file), input_object = input_object, analysis_label = analysis_label,
  input_size = file.info(input_file)$size, contrast = "Infected - Bystander", replicate_unit = replicate_unit,
  sample_table = sample_table, tumor_labels = tumor_labels,
  metadata_columns = c(sample = sample_col, celltype = celltype_col, state = state_col),
  host_prefix = host_prefix, viral_pattern = viral_pattern, viral_feature_ids = viral_feature_ids,
  n_viral_features_excluded = length(unique(viral_exclusion_audit$feature)),
  n_viral_features_in_rank = sum(is_viral_feature(names(ranks))), min_cells_per_state = min_cells_per_state,
  min_pairs = min_pairs, included_pairs = included_pairs, n_pairs = length(included_pairs),
  filter_min_count = filter_min_count, filter_min_total_count = filter_min_total_count,
  rank_method = "sign(logFC) * sqrt(single-coefficient edgeR QL F); ties ordered by symbol", tie_fraction = tie_fraction,
  min_set_size = min_set_size, max_set_size = max_set_size, seed = seed, eps = gsea_eps, nPermSimple = gsea_simple_n,
  fdr_cutoff = fdr_cutoff, gene_set_source = gene_set_source, mapping_fraction = mapping_fraction,
  n_ranked = length(ranks), n_matched = length(matched_genes), target_gene = target_gene, target_status = target_status,
  focus_pathways = focus_pathways, missing_focus = missing_focus, overview_names = overview_names, curve_names = curve_names,
  remove_target_sensitivity = remove_target_sensitivity,
  package_versions = setNames(vapply(packages, function(p) as.character(packageVersion(p)), character(1)), packages),
  session = capture.output(sessionInfo()),
  limits = c("Infection labels are inherited, not revalidated here.",
             "Epithelial annotation must be validated as malignant tumor cells upstream.",
             "GO T-cell terms in tumor cells are not measurements of T-cell function.",
             "Leading edge and rank-item removal do not establish MDK causality.",
             "Small biological n and residual within-tumor state heterogeneity limit inference.")
)
results <- list(info = info, sample_audit = sample_audit, layer_audit = layer_audit,
                viral_exclusion_audit = viral_exclusion_audit, pb_meta = pb_meta,
                design = design, pseudobulk_counts = counts_use, normalization = y$samples,
                genes = de_all, ranks = ranks, pathways = pathways, pathway_meta = pathway_meta,
                gsea = gsea, gsea_without_target = gsea_without_target, target_de = target_de,
                target_pathways = target_pathways, summary = summary_table)
saveRDS(results, file.path(outdir, "results.rds"), compress = TRUE)
message("完成：", normalizePath(outdir), "；输出 genes.csv、summary.csv、plots.pdf、results.rds。")
