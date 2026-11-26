##### 清理运行环境 #####
rm(list = ls(all.names = TRUE))
graphics.off()
invisible(gc())

##### 加载依赖 #####
if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))
library(fgsea)
library(BiocParallel)
library(ggplot2)

##### 路径与分析参数 #####
script_args <- grep("^--file=", commandArgs(), value = TRUE)
source_files <- unlist(lapply(sys.frames(), function(z) z$ofile), use.names = FALSE)
script_file <- if (length(source_files)) tail(source_files, 1L) else if (length(script_args)) sub("^--file=", "", script_args[1L]) else ""
code_dir <- if (nzchar(script_file)) dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)) else if (file.exists("code/scRNA_MDK_perturbation.R")) normalizePath("code", winslash = "/") else getwd()
project_dir <- dirname(code_dir)
reference_root <- if (.Platform$OS.type == "windows") file.path(project_dir, "_analysis") else
  path.expand("~/rawdata/scRNA_virus")
source_dir <- file.path(reference_root, "mdk_perturbation_sources")
datasets <- c("GSE227286", "GSE227287")
reference_paths <- setNames(file.path(source_dir, paste0(datasets, ".tsv")), datasets)
gmt_paths <- c(Hallmark = file.path(source_dir, "hallmark_mouse.gmt"),
               GO_BP = file.path(source_dir, "gobp_mouse.gmt"))
ortholog_file <- file.path(source_dir, "orthologs.tsv")
string_file <- file.path(source_dir, "string_associations.tsv")
d3_root <- if (.Platform$OS.type == "windows") "D:/XLH/溶瘤病毒/MDK/scRNA" else file.path(path.expand("~/OV"), "scRNA_virus")
d3_file <- file.path(d3_root, "immune_NicheNet_D3", "results.rds")
function_file <- if (.Platform$OS.type == "windows") file.path(code_dir, "new_function.R") else path.expand("~/code/new_function.R")
outdir <- file.path(project_dir, "OV","scRNA_virus", "MDK_perturbation_D3")
receivers <- c("Myeloid", "Tcell")
receiver_contrast <- "VG161 - Vehicle"
directions <- c("Up", "Down")
gene_fdr <- 0.05
min_abs_log2FC <- 0.25
top_signature <- 100L
min_size <- 15L
max_size <- 500L
min_projection_overlap <- 15L
seed <- 123L
eps <- 1e-10
n_perm_simple <- 1000L
display_top_genes <- 10L
display_top_pathways <- 6L
immune_pathways <- c(
  "Interferon alpha response" = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "Interferon gamma response" = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "IL6 JAK STAT3 signaling" = "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "Inflammatory response" = "HALLMARK_INFLAMMATORY_RESPONSE",
  "Myeloid leukocyte activation" = "GOBP_MYELOID_LEUKOCYTE_ACTIVATION",
  "Dendritic cell differentiation" = "GOBP_DENDRITIC_CELL_DIFFERENTIATION",
  "Myeloid dendritic cell activation" = "GOBP_MYELOID_DENDRITIC_CELL_ACTIVATION",
  "Antigen processing and presentation" = "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "Presentation of exogenous peptides via MHC II" = "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II",
  "Phagocytosis" = "GOBP_PHAGOCYTOSIS",
  "Complement activation" = "GOBP_COMPLEMENT_ACTIVATION",
  "T cell activation" = "GOBP_T_CELL_ACTIVATION",
  "T cell mediated cytotoxicity" = "GOBP_T_CELL_MEDIATED_CYTOTOXICITY")

##### 读入已有统计和检查数据格式 #####
started_at <- Sys.time()
paths <- c(reference_paths, gmt_paths, orthologs = ortholog_file, STRING = string_file, D3 = d3_file, functions = function_file)
if (any(!file.exists(paths))) stop("缺少文件：", paste(paths[!file.exists(paths)], collapse = ", "), call. = FALSE)
mdk_functions <- new.env(parent = environment())
source(function_file, local = mdk_functions, encoding = "UTF-8")
if (!exists("write_mdk_perturbation_plots", envir = mdk_functions, mode = "function", inherits = FALSE)) {
  stop("请先将MDK公开扰动绘图函数合并到：", function_file, call. = FALSE)
}
d3 <- readRDS(d3_file)
if (!identical(d3$info$analysis_label, "D3") || !identical(d3$info$contrast, receiver_contrast) ||
    !all(receivers %in% names(d3$receiver_DE))) {
  stop("D3输入缺少指定接收端或VG161−Vehicle对比。", call. = FALSE)
}
reference <- lapply(reference_paths, read.delim, stringsAsFactors = FALSE, check.names = FALSE,
                    na.strings = c("", "NA", "NaN"), fileEncoding = "UTF-8")
required <- c("gene_id", "mouse_symbol", "log2FC", "stat", "pvalue", "padj")
if (!all(vapply(reference, function(x) all(required %in% names(x)) && nrow(x) > 0L, logical(1)))) {
  stop("公开统计表需包含：", paste(required, collapse = ", "), call. = FALSE)
}
orthologs <- read.delim(ortholog_file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
if (!all(c("mouse_symbol", "human_symbol") %in% names(orthologs))) stop("同源表缺少鼠和人基因符号字段。", call. = FALSE)
orthologs <- unique(orthologs[c("mouse_symbol", "human_symbol")])
orthologs <- orthologs[!is.na(orthologs$mouse_symbol) & !is.na(orthologs$human_symbol) &
                         nzchar(orthologs$mouse_symbol) & nzchar(orthologs$human_symbol), , drop = FALSE]
ambiguous_mouse <- orthologs$mouse_symbol[duplicated(orthologs$mouse_symbol)]
ambiguous_human <- orthologs$human_symbol[duplicated(orthologs$human_symbol)]
orthologs_1to1 <- orthologs[!orthologs$mouse_symbol %in% ambiguous_mouse &
                              !orthologs$human_symbol %in% ambiguous_human, , drop = FALSE]
if (!nrow(orthologs_1to1)) stop("同源表未留下严格一对一符号映射。", call. = FALSE)
network <- read.delim(string_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
network_fields <- c("accession", "gene_a", "gene_b", "string_id_a", "string_id_b", "score", "experimental_score", "database_score")
if (!all(network_fields %in% names(network)) || anyNA(network[network_fields]) ||
    any(!network$accession %in% datasets | network$gene_a == network$gene_b | network$score < .4 |
        network$score > 1 | (network$experimental_score <= 0 & network$database_score <= 0))) {
  stop("STRING表的基因对、证据字段或分值无效。", call. = FALSE)
}

##### 构建公开鼠源排序与原始通路库 #####
gene_tables <- vector("list", length(datasets))
names(gene_tables) <- datasets
mouse_ranks <- gene_tables
rank_methods <- setNames(character(length(datasets)), datasets)
for (dataset in datasets) {
  x <- reference[[dataset]]
  for (col in c("log2FC", "stat", "pvalue", "padj")) x[[col]] <- as.numeric(x[[col]])
  if (any(x$pvalue < 0 | x$pvalue > 1, na.rm = TRUE) || any(x$padj < 0 | x$padj > 1, na.rm = TRUE)) {
    stop(dataset, "的P值或调整P值超出[0,1]。", call. = FALSE)
  }
  rank_methods[[dataset]] <- if (any(is.finite(x$stat))) "reported_stat" else "signed_neg_log10_pvalue"
  x$rank_stat <- if (rank_methods[[dataset]] == "reported_stat") x$stat else
    sign(x$log2FC) * -log10(pmax(x$pvalue, .Machine$double.xmin))
  ranked <- which(!is.na(x$mouse_symbol) & nzchar(x$mouse_symbol) & is.finite(x$rank_stat) & is.finite(x$log2FC))
  ranked <- ranked[order(-abs(x$rank_stat[ranked]), x$gene_id[ranked], na.last = TRUE)]
  ranked <- ranked[!duplicated(x$mouse_symbol[ranked])]
  if (!length(ranked)) stop(dataset, "未找到具可用统计量的鼠基因符号。", call. = FALSE)
  x$used_in_rank <- seq_len(nrow(x)) %in% ranked
  x$dataset <- dataset
  gene_tables[[dataset]] <- x
  mouse_ranks[[dataset]] <- setNames(x$rank_stat[ranked], x$mouse_symbol[ranked])
}
genes <- do.call(rbind, lapply(gene_tables, function(x) {
  columns <- c("dataset", required, "rank_stat", "used_in_rank", "baseMean")
  if (!"baseMean" %in% names(x)) x$baseMean <- NA_real_
  x[columns]
}))
rownames(genes) <- NULL
for (dataset in datasets) {
  edges <- network[network$accession == dataset, , drop = FALSE]
  if (!all(c(edges$gene_a, edges$gene_b) %in% names(mouse_ranks[[dataset]]))) {
    stop(dataset, "的 STRING 边包含无法关联排序表的鼠基因。", call. = FALSE)
  }
}

gmt <- lapply(gmt_paths, function(path) {
  lines <- strsplit(readLines(path, encoding = "UTF-8", warn = FALSE), "\t", fixed = TRUE)
  if (!length(lines) || any(lengths(lines) < 3L)) stop("GMT格式无效：", path, call. = FALSE)
  gene_sets <- lapply(lines, function(fields) unique(fields[-c(1L, 2L)][nzchar(fields[-c(1L, 2L)])]))
  setNames(gene_sets, vapply(lines, `[`, character(1), 1L))
})
collections <- rep(names(gmt), lengths(gmt))
gene_sets <- do.call(c, unname(gmt))
if (is.null(names(gene_sets)) || anyDuplicated(names(gene_sets))) stop("鼠GMT中通路名缺失或重复。", call. = FALSE)

##### 公开鼠源签名的全库GSEA #####
external_paths <- vector("list", length(datasets))
names(external_paths) <- datasets
for (dataset in datasets) {
  ranks <- mouse_ranks[[dataset]]
  ranks <- ranks[order(-ranks, names(ranks))]
  overlap <- vapply(gene_sets, function(x) sum(x %in% names(ranks)), integer(1))
  eligible <- overlap >= min_size & overlap <= max_size & overlap < length(ranks)
  message(dataset, "：", length(ranks), " 个排序基因，", sum(eligible), " 条可检验鼠源通路。")
  result <- data.frame(dataset = dataset, collection = collections, pathway = names(gene_sets),
                       size = lengths(gene_sets), n_overlap = unname(overlap), eligible = unname(eligible),
                       ES = NA_real_, NES = NA_real_, P = NA_real_, FDR = NA_real_, log2err = NA_real_,
                       leading_edge = NA_character_, status = ifelse(eligible, "not_estimated", "outside_size_filter"), stringsAsFactors = FALSE)
  if (any(eligible)) {
    set.seed(seed)
    fit <- fgsea::fgseaMultilevel(pathways = gene_sets[eligible], stats = ranks, minSize = min_size,
                                  maxSize = max_size, scoreType = "std", eps = eps, nPermSimple = n_perm_simple, BPPARAM = BiocParallel::SerialParam())
    idx <- match(fit$pathway, result$pathway)
    result$ES[idx] <- fit$ES
    result$NES[idx] <- fit$NES
    result$P[idx] <- fit$pval
    result$log2err[idx] <- fit$log2err
    result$leading_edge[idx] <- vapply(fit$leadingEdge, paste, character(1), collapse = ";")
    ok <- which(is.finite(result$P))
    result$FDR[ok] <- p.adjust(result$P[ok], method = "BH", n = sum(eligible))
    result$status[ok] <- "tested"
  }
  external_paths[[dataset]] <- result
}
pathways <- do.call(rbind, external_paths)
rownames(pathways) <- NULL
pathways$immune_focus <- pathways$pathway %in% unname(immune_pathways)

##### 从公开扰动独立选取上下调签名并映射 #####
signature_specs <- expand.grid(direction = directions, dataset = datasets, stringsAsFactors = FALSE)
signature_tables <- vector("list", nrow(signature_specs))
signature_stats <- signature_specs
signature_stats$n_mouse_significant <- signature_stats$n_mouse_signature <- signature_stats$n_mapped <- 0L
for (i in seq_len(nrow(signature_specs))) {
  dataset <- signature_specs$dataset[i]
  direction <- signature_specs$direction[i]
  direction_sign <- if (direction == "Up") 1 else -1
  x <- gene_tables[[dataset]]
  eligible_genes <- which(x$used_in_rank & is.finite(x$padj) & x$padj < gene_fdr &
                            direction_sign * x$log2FC >= min_abs_log2FC)
  order_genes <- eligible_genes[order(-abs(x$rank_stat[eligible_genes]), x$mouse_symbol[eligible_genes])]
  x <- x[head(order_genes, top_signature), c("dataset", "gene_id", "mouse_symbol", "log2FC", "stat", "rank_stat", "pvalue", "padj"), drop = FALSE]
  x$direction <- rep(direction, nrow(x))
  x$signature_rank <- seq_len(nrow(x))
  x$human_symbol <- orthologs_1to1$human_symbol[match(x$mouse_symbol, orthologs_1to1$mouse_symbol)]
  x$mapping_status <- ifelse(is.na(x$human_symbol), "no_unique_1to1_symbol", "mapped_1to1")
  signature_tables[[i]] <- x
  signature_stats$n_mouse_significant[i] <- length(eligible_genes)
  signature_stats$n_mouse_signature[i] <- nrow(x)
  signature_stats$n_mapped[i] <- sum(!is.na(x$human_symbol))
}
signature_genes <- do.call(rbind, signature_tables)
rownames(signature_genes) <- NULL

##### 映射到真实D3接收端排序 #####
projection <- merge(expand.grid(receiver = receivers, i = seq_len(nrow(signature_specs)), stringsAsFactors = FALSE),
                    transform(signature_stats, i = seq_len(nrow(signature_stats))), by = "i", sort = FALSE)
projection <- projection[order(match(projection$receiver, receivers), projection$i), , drop = FALSE]
projection$n_overlap <- 0L
projection$n_receiver_ranked <- NA_integer_
projection$eligible <- FALSE
projection$ES <- projection$NES <- projection$P <- projection$FDR <- projection$log2err <- NA_real_
projection$overlap_genes <- projection$leading_edge <- NA_character_
projection$status <- "insufficient_overlap"
for (receiver in receivers) {
  de <- d3$receiver_DE[[receiver]]
  if (!all(c("gene", "rank_stat") %in% names(de)) || anyDuplicated(de$gene)) stop("D3接收端缺少唯一基因排序：", receiver, call. = FALSE)
  valid <- !is.na(de$gene) & nzchar(de$gene) & is.finite(de$rank_stat)
  ranks <- setNames(de$rank_stat[valid], de$gene[valid])
  ranks <- ranks[order(-ranks, names(ranks))]
  receiver_rows <- which(projection$receiver == receiver)
  projection$n_receiver_ranked[receiver_rows] <- length(ranks)
  sets <- lapply(projection$i[receiver_rows], function(i) {
    ids <- signature_tables[[i]]$human_symbol
    intersect(ids[!is.na(ids)], names(ranks))
  })
  projection$n_overlap[receiver_rows] <- lengths(sets)
  projection$overlap_genes[receiver_rows] <- vapply(sets, paste, character(1), collapse = ";")
  eligibility <- lengths(sets) >= min_projection_overlap & lengths(sets) < length(ranks)
  message(receiver, "：", length(ranks), " 个D3排序基因，", sum(eligibility), " 个可检验外部签名。")
  projection$eligible[receiver_rows] <- eligibility
  projection$status[receiver_rows[eligibility]] <- "not_estimated"
  if (!any(eligibility)) next
  names(sets) <- as.character(receiver_rows)
  set.seed(seed)
  fit <- fgsea::fgseaMultilevel(pathways = sets[eligibility], stats = ranks,
                                minSize = min_projection_overlap, maxSize = top_signature, scoreType = "std", eps = eps,
                                nPermSimple = n_perm_simple, BPPARAM = BiocParallel::SerialParam())
  idx <- as.integer(fit$pathway)
  projection$ES[idx] <- fit$ES
  projection$NES[idx] <- fit$NES
  projection$P[idx] <- fit$pval
  projection$log2err[idx] <- fit$log2err
  projection$leading_edge[idx] <- vapply(fit$leadingEdge, paste, character(1), collapse = ";")
}
# 2个扰动时间/体系 × 2个签名方向 × 2个接收种类为预指定投影检验家族。
ok <- which(is.finite(projection$P))
projection$FDR[ok] <- p.adjust(projection$P[ok], method = "BH", n = sum(projection$eligible))
projection$status[ok] <- "tested"
projection$orientation <- ifelse(is.finite(projection$NES) & projection$NES != 0,
                                 ifelse(sign(projection$NES) == ifelse(projection$direction == "Up", 1, -1), "concordant", "opposite"), NA_character_)
projection$i <- NULL
rownames(projection) <- NULL

##### 输出与绘图 #####
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
info <- list(started_at = started_at, finished_at = Sys.time(),
             inputs = data.frame(input = names(paths), path = unname(paths), md5 = unname(tools::md5sum(paths))),
             dataset_description = c(GSE227286 = "Day 3 FLT3LG bone marrow cultures, rMDK on days 0-3",
                                     GSE227287 = "Day 9 FLT3LG bone-marrow-derived dendritic cells, rMDK on days 7-9"),
             external_contrast = "rMDK - control", external_species = "Mus musculus", rank_methods = rank_methods,
             receiver_contrast = receiver_contrast, receiver_primary = "Myeloid", receiver_secondary = "Tcell",
             receiver_species = "Homo sapiens symbol assumption consistent with D3 NicheNet input; annotations require verification",
             sample_table = d3$info$sample_table, receiver_rank = "Original receiver edgeR signed rank_stat",
             gene_fdr = gene_fdr, min_abs_log2FC = min_abs_log2FC, top_signature = top_signature,
             min_size = min_size, max_size = max_size, min_projection_overlap = min_projection_overlap,
             external_family = "BH across all size-eligible Hallmark and GO BP terms per mouse dataset",
             projection_family = "Joint BH across size-eligible members of eight prespecified receiver by dataset by direction tests",
             n_projection_family_eligible = sum(projection$eligible),
             mapping = "Exact symbol match, strict one-to-one only, selected mouse signatures fixed before mapping",
             seed = seed, eps = eps, n_perm_simple = n_perm_simple,
             display_top_genes = display_top_genes, display_top_pathways = display_top_pathways,
             immune_pathways = immune_pathways, network = "STRING v12.0 mouse functional associations among original top 25 up and top 25 down; combined score >= 0.4 and experimental or database score > 0",
             interpretation = "Receiver agreement is observational association with a mouse MDK signature; contextual T-cell projection does not test direct MDK effects",
             versions = list(R = R.version.string, fgsea = as.character(packageVersion("fgsea")), ggplot2 = as.character(packageVersion("ggplot2"))))
results <- list(info = info, genes = genes, pathways = pathways, projection = projection,
                network = network, signature_genes = signature_genes, signature_summary = signature_stats,
                score_contrasts = d3$score_contrasts)
utils::write.csv(genes, file.path(outdir, "genes.csv"), row.names = FALSE, na = "NA")
utils::write.csv(pathways, file.path(outdir, "pathways.csv"), row.names = FALSE, na = "NA")
utils::write.csv(projection, file.path(outdir, "summary.csv"), row.names = FALSE, na = "NA")
saveRDS(results, file.path(outdir, "results.rds"), compress = "xz")
mdk_functions$write_mdk_perturbation_plots(results, file.path(outdir, "plots.pdf"))
message("已保存MDK公开扰动参考与D3映射：", normalizePath(outdir, winslash = "/"))
