##### 清理运行环境 #####
rm(list = ls(all.names = TRUE))
graphics.off()
invisible(gc())

##### 加载依赖 #####
library(Seurat)
library(Matrix)
library(CellChat)
library(future)
library(ggplot2)
library(patchwork)
if (packageVersion("ggplot2") >= "4.0.0" && packageVersion("patchwork") < "1.3.2") {
  stop("当前ggplot2为4.x，patchwork需至少1.3.2；请先更新patchwork，再从本脚本开头运行。")
}

##### 路径、样本与分析参数 #####
# 旧代码为D5；这里默认衔接当前正在分析的D3完整主注释对象。
analysis_label <- "D3"
if (!grepl("^D[0-9]+$", analysis_label)) stop("analysis_label请使用D3、D5等格式。")
data_tag <- paste0(sub("^D", "", analysis_label), "D")
input_file <- file.path("~/rawdata", "scRNA_virus", paste0("virus_", data_tag), paste0("virus_", data_tag, "_anno.RData"))
input_object <- "sce"
outdir <- file.path(path.expand("~/OV"), "scRNA_virus", paste0("CellChat_", analysis_label))
sample_col <- "orig.ident"
celltype_col <- "celltype"
state_col <- "infected"
study_samples <- paste0(rep(c("Vehicle", "VG161"), each = 2L), "_", data_tag, "_", rep(1:2, 2L))
sample_table <- data.frame(
  sample_id = study_samples,
  replicate_id = study_samples,
  condition = c("Vehicle", "Vehicle", "VG161", "VG161"), stringsAsFactors = FALSE)
condition_levels <- c("Vehicle", "VG161")
tumor_labels <- "Epithelial"
immune_receivers <- c("Tcell", "Myeloid")
state_levels <- c("Bystander", "Infected")
state_groups <- c("Tumor_Bystander", "Tumor_Infected")
db_scope <- "All"                       # All=全部数据库类别；Secreted仅分泌通讯。
database_override <- NULL               # 默认使用安装版human DB；仅在明确验证的自定义DB时替换。
expression_mode <- "host_logCP10K"       # 默认从原始counts排除病毒后按全部共同宿主UMI归一化。
probability_scale <- "shared_max"        # shared_max用于可比网络；native_per_object仅复现官方逐对象尺度。
run_overexpression_diagnostics <- FALSE # 候选固定，过表达预筛选不参与本轮不同网络的候选选择。
min_cells <- 20L                        # 每动物、每细胞大类/肿瘤状态。
min_animals <- 2L
communication_p <- 0.05                 # CellChat细胞标签置换P，非动物差异P/FDR。
nboot <- 100L
average_type <- "triMean"
population_size <- FALSE                # 本轮比较表达相关通讯，不乘捕获细胞比例。
seed <- 123L
host_prefix <- ""
viral_pattern <- "^(VG161|HSV[-_]?1)([-_.]|$)"
viral_feature_ids <- character()
focus_pathway <- "MK"
focus_ligand <- "MDK"
focus_receptor <- "NCL"
expression_genes <- c("MDK", "NCL")
save_cellchat_objects <- FALSE          # 默认只保留紧凑网络/表，不重复保存表达矩阵。
condition_colors <- c(Vehicle = "#426B91", VG161 = "#C35B3F")
state_colors <- c(Tumor_Bystander = "#426B91", Tumor_Infected = "#C35B3F")

##### 通用函数与配置校验 #####
bind_rows_base <- function(x) {
  x <- Filter(function(z) is.data.frame(z) && nrow(z), x)
  if (length(x)) do.call(rbind, x) else data.frame()
}
read_seurat <- function(path, object_name) {
  path <- path.expand(path)
  if (!file.exists(path)) stop("找不到输入文件：", path)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  e <- new.env(parent = emptyenv()); nms <- load(path, envir = e)
  if (!object_name %in% nms) stop("输入中缺少对象：", object_name)
  get(object_name, envir = e)
}
required <- c("sample_id", "replicate_id", "condition")
if (!all(required %in% names(sample_table)) || anyNA(sample_table[, required]) ||
    any(!nzchar(as.matrix(sample_table[, required]))) || anyDuplicated(sample_table$sample_id)) stop("样本表缺失/重复。")
if (!all(sample_table$condition %in% condition_levels) ||
    any(vapply(split(sample_table$condition, sample_table$replicate_id), function(z) length(unique(z)) != 1L, logical(1)))) {
  stop("动物与条件映射不明确；Vehicle/VG161为独立动物设计。")
}
if (!db_scope %in% c("All", "Secreted") || min_cells < 3L || min_animals < 2L || nboot < 20L ||
    communication_p <= 0 || communication_p >= 1) stop("请核对DB范围、细胞/动物数与置换参数。")
if (!identical(state_levels, c("Bystander", "Infected")) || !identical(state_groups, c("Tumor_Bystander", "Tumor_Infected"))) stop("状态名称配置不符。")
if (!expression_mode %in% c("host_logCP10K", "RNA_data") || !probability_scale %in% c("shared_max", "native_per_object")) stop("表达归一化或概率缩放配置有误。")
future::plan("sequential")
set.seed(seed)
started_at <- Sys.time()
outdir <- path.expand(outdir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

##### 读取完整对象和共同宿主表达背景 #####
sce <- read_seurat(input_file, input_object)
if (!inherits(sce, "Seurat") || !"RNA" %in% names(sce@assays)) stop("需要含RNA assay的完整主Seurat对象。")
meta <- sce[[]]
if (!all(c(sample_col, celltype_col) %in% names(meta)) || anyDuplicated(rownames(meta))) stop("缺少样本/大类列或条码重复。")
if (!all(sample_table$sample_id %in% as.character(meta[[sample_col]]))) stop("样本表与对象不匹配。")
meta <- meta[!is.na(meta[[sample_col]]) & meta[[sample_col]] %in% sample_table$sample_id, , drop = FALSE]
if (anyNA(meta[[celltype_col]]) || any(!nzchar(trimws(as.character(meta[[celltype_col]]))))) stop("纳入细胞有缺失大类注释。")
meta$sample_id <- as.character(meta[[sample_col]])
meta$replicate_id <- sample_table$replicate_id[match(meta$sample_id, sample_table$sample_id)]
meta$condition <- sample_table$condition[match(meta$sample_id, sample_table$sample_id)]
meta$celltype <- as.character(meta[[celltype_col]])
state_column_present <- state_col %in% names(meta)
meta$state <- if (state_column_present) as.character(meta[[state_col]]) else NA_character_
# 不读取T细分对象，不使用t_type/t_function/cluster重命名。
rna <- sce[["RNA"]]
expression_layer <- if (expression_mode == "host_logCP10K") "counts" else "data"
if (inherits(rna, "Assay5")) {
  layer_names <- Layers(rna, search = paste0("^", expression_layer, "($|\\.)"))
  if (!length(layer_names)) stop("RNA缺少所选表达层：", expression_layer)
  layers <- lapply(layer_names, function(k) {
    z <- LayerData(rna, layer = k, fast = FALSE)
    z[, intersect(colnames(z), rownames(meta)), drop = FALSE]
  })
  keep <- vapply(layers, ncol, integer(1)) > 0L
  layers <- layers[keep]; layer_names <- layer_names[keep]
} else {
  z <- if (packageVersion("SeuratObject") >= "5.0.0") GetAssayData(sce, assay = "RNA", layer = expression_layer) else GetAssayData(sce, assay = "RNA", slot = expression_layer)
  layers <- list(z[, intersect(colnames(z), rownames(meta)), drop = FALSE]); layer_names <- expression_layer
}
if (!length(layers) || anyDuplicated(unlist(lapply(layers, colnames), use.names = FALSE)) ||
    !setequal(unlist(lapply(layers, colnames), use.names = FALSE), rownames(meta))) stop("RNA表达层条码重叠或覆盖不完整。")
layers <- lapply(layers, function(z) {
  if (is.null(rownames(z)) || anyDuplicated(rownames(z))) stop("RNA feature缺失/重复。")
  z <- as(Matrix(z, sparse = TRUE), "dgCMatrix")
  if (any(!is.finite(z@x) | z@x < 0)) stop("RNA data必须是非负有限log-normalized表达；不能用scale.data/integrated。")
  if (expression_mode == "host_logCP10K" && any(abs(z@x - round(z@x)) > 1e-6)) stop("host_logCP10K需要RNA原始非负整数counts。")
  z
})
features_by_layer <- lapply(layers, function(z) {
  f <- rownames(z)
  if (nzchar(host_prefix)) f <- f[startsWith(f, host_prefix)]
  symbol <- if (nzchar(host_prefix)) substring(f, nchar(host_prefix) + 1L) else f
  f[!grepl(viral_pattern, f, ignore.case = TRUE) & !grepl(viral_pattern, symbol, ignore.case = TRUE) &
      !f %in% viral_feature_ids & !symbol %in% viral_feature_ids]
})
common_features <- Reduce(intersect, features_by_layer)
symbols <- if (nzchar(host_prefix)) substring(common_features, nchar(host_prefix) + 1L) else common_features
if (length(symbols) < 100L || anyDuplicated(symbols) || mean(grepl("^ENSG[0-9]", symbols)) > 0.5) stop("共同宿主gene symbol不足/重复；先核对人鼠来源与ID。")
layer_audit <- data.frame(layer = layer_names, n_cells = vapply(layers, ncol, integer(1)),
                          n_features = vapply(layers, nrow, integer(1)), n_host = lengths(features_by_layer), n_common = length(symbols))
data_input <- do.call(cbind, lapply(layers, function(z) z[common_features, , drop = FALSE]))
rownames(data_input) <- symbols
meta <- meta[colnames(data_input), , drop = FALSE]
valid <- Matrix::colSums(data_input) > 0
zero_expression_cells <- rownames(meta)[!valid]
data_input <- data_input[, valid, drop = FALSE]; meta <- meta[valid, , drop = FALSE]
if (expression_mode == "host_logCP10K") {
  # 分母为完整共同宿主feature，而非只取信号基因；病毒增加本身不机械稀释宿主表达。
  host_library_size <- Matrix::colSums(data_input)
  expression_dimnames <- dimnames(data_input)
  data_input <- as(data_input %*% Diagonal(x = 10000 / host_library_size), "dgCMatrix")
  dimnames(data_input) <- expression_dimnames
  data_input@x <- log1p(data_input@x)
}
if (!identical(colnames(data_input), rownames(meta))) stop("归一化后表达条码与metadata顺序不一致。")
rm(sce, rna, layers)
invisible(gc())

##### 固定共同数据库和LR候选背景 #####
DB <- if (is.null(database_override)) CellChatDB.human else database_override
if (db_scope == "Secreted") DB <- subsetDB(DB, search = "Secreted Signaling")
db_interaction <- as.data.frame(DB$interaction, stringsAsFactors = FALSE)
if (!all(c("ligand", "receptor", "pathway_name") %in% names(db_interaction))) stop("CellChatDB缺少必要字段。")
if (!"interaction_name" %in% names(db_interaction)) db_interaction$interaction_name <- rownames(db_interaction)
for (field in intersect(c("interaction_name", "ligand", "receptor", "pathway_name", "annotation"), names(db_interaction))) {
  db_interaction[[field]] <- as.character(db_interaction[[field]])
}
if (anyDuplicated(db_interaction$interaction_name)) stop("数据库interaction_name不唯一。")
rownames(db_interaction) <- db_interaction$interaction_name
expand_complex <- function(x) {
  if (x %in% rownames(DB$complex)) {
    a <- as.character(unlist(DB$complex[x, , drop = FALSE], use.names = FALSE))
    a[!is.na(a) & nzchar(a)]
  } else x
}
required_genes <- lapply(seq_len(nrow(db_interaction)), function(i) unique(c(
  expand_complex(db_interaction$ligand[i]), expand_complex(db_interaction$receptor[i]))))
measured_lr <- vapply(required_genes, function(g) all(g %in% rownames(data_input)), logical(1))
db_audit <- data.frame(interaction_name = db_interaction$interaction_name,
                       pathway = db_interaction$pathway_name, ligand = db_interaction$ligand, receptor = db_interaction$receptor,
                       measured = measured_lr, missing_genes = vapply(required_genes, function(g) paste(setdiff(g, rownames(data_input)), collapse = ";"), character(1)))
fixed_LR <- db_interaction[measured_lr, , drop = FALSE]
if (!nrow(fixed_LR)) stop("没有共同测到的配体-受体；请核对数据库、物种及feature。")
# 不按各对象自己的overexpressed结果改变候选宇宙；所有对象计算同一套完整可测LR。
# 普通computeCommunProb通过LRsig读取该固定表；不补造数据库没有的MDK-NCL边。
reference_meta <- data.frame(labels = meta$celltype, samples = factor(meta$replicate_id), row.names = rownames(meta))
reference_cellchat <- createCellChat(object = data_input, meta = reference_meta, group.by = "labels")
reference_cellchat@DB <- DB
reference_cellchat <- subsetData(reference_cellchat)
signaling_features <- rownames(reference_cellchat@data.signaling)
global_signaling_max <- max(reference_cellchat@data.signaling)
if (!is.finite(global_signaling_max) || global_signaling_max <= 0) stop("数据库信号基因没有有效表达。")
rm(reference_cellchat)

##### 统一CellChat内部表达缩放，保留其余官方公式和置换 #####
# 官方函数每对象除以自己的max；直接预除共同常数会被再次抵消。
# 复制当前安装版函数，仅参数化这一处分母；不修改包namespace，不添加虚构表达锚点。
compute_probability_fixed <- CellChat::computeCommunProb
original_probability_body <- body(compute_probability_fixed)
replacement_count <- 0L
replace_scale <- function(x, reverse = FALSE) {
  target <- if (reverse) quote(data.use <- data/global_signaling_max) else quote(data.use <- data/max(data))
  replacement <- if (reverse) quote(data.use <- data/max(data)) else quote(data.use <- data/global_signaling_max)
  if (is.call(x) && identical(x, target)) {
    if (!reverse) replacement_count <<- replacement_count + 1L
    return(replacement)
  }
  if (is.call(x) || is.expression(x) || is.pairlist(x)) for (i in seq_along(x)) {
    if (!is.null(x[[i]]) && !identical(x[[i]], quote(expr = ))) x[i] <- list(replace_scale(x[[i]], reverse))
  }
  x
}
if (probability_scale == "shared_max") {
  body(compute_probability_fixed) <- replace_scale(original_probability_body)
  if (replacement_count != 1L) stop("当前CellChat内部缩放实现已变化；需核对computeCommunProb，不能静默套用固定尺度。")
  if (!identical(replace_scale(body(compute_probability_fixed), reverse = TRUE), original_probability_body)) {
    stop("固定尺度实现的逆向核对失败；不得改变官方函数的其他公式或参数。")
  }
  scale_environment <- new.env(parent = environment(CellChat::computeCommunProb))
  scale_environment$global_signaling_max <- global_signaling_max
  environment(compute_probability_fixed) <- scale_environment
}
source_file <- tempfile(fileext = ".R")
writeLines(deparse(original_probability_body), source_file)
probability_function_md5 <- unname(tools::md5sum(source_file))
unlink(source_file)

##### 各动物细胞资格与肿瘤状态配对 #####
replicate_table <- unique(sample_table[, c("replicate_id", "condition")])
main_groups <- unique(c(tumor_labels[tumor_labels %in% meta$celltype], immune_receivers[immune_receivers %in% meta$celltype], sort(unique(meta$celltype))))
main_counts <- as.data.frame(table(replicate_id = factor(meta$replicate_id, levels = replicate_table$replicate_id),
                                   group = factor(meta$celltype, levels = main_groups)), stringsAsFactors = FALSE)
names(main_counts)[3] <- "n_cells"
main_counts$replicate_id <- as.character(main_counts$replicate_id); main_counts$group <- as.character(main_counts$group)
main_counts$condition <- replicate_table$condition[match(main_counts$replicate_id, replicate_table$replicate_id)]
main_counts$eligible <- main_counts$n_cells >= min_cells
main_key <- paste(main_counts$replicate_id, main_counts$group, sep = "\r")
meta$main_eligible <- main_counts$eligible[match(paste(meta$replicate_id, meta$celltype, sep = "\r"), main_key)]
vg_ids <- replicate_table$replicate_id[replicate_table$condition == "VG161"]
state_counts <- expand.grid(replicate_id = vg_ids, state = state_levels, stringsAsFactors = FALSE)
state_counts$n_cells <- vapply(seq_len(nrow(state_counts)), function(i) sum(meta$replicate_id == state_counts$replicate_id[i] &
                                                                              meta$celltype %in% tumor_labels & !is.na(meta$state) & meta$state == state_counts$state[i]), integer(1))
complete_pairs <- vg_ids[vapply(vg_ids, function(id) all(state_counts$n_cells[state_counts$replicate_id == id] >= min_cells), logical(1))]
state_counts$pair_eligible <- state_counts$replicate_id %in% complete_pairs
meta$state_group <- meta$celltype
meta$state_group[meta$celltype %in% tumor_labels] <- NA_character_
for (i in seq_along(state_levels)) {
  ii <- meta$celltype %in% tumor_labels & !is.na(meta$state) & meta$state == state_levels[i]
  meta$state_group[ii] <- state_groups[i]
}
state_group_levels <- unique(c(state_groups, setdiff(main_groups, tumor_labels)))
meta$state_eligible <- meta$replicate_id %in% complete_pairs & !is.na(meta$state_group) &
  (meta$celltype %in% tumor_labels | meta$main_eligible)
state_status <- if (!state_column_present) "state_column_absent" else if (length(complete_pairs) >= min_animals) "assessed" else "insufficient_complete_pairs"

##### 单网络推断与完整边表：缺失为NA，实测无支持为0 #####
all_edges <- list(); all_pathways <- list(); object_audit <- list(); compact_networks <- list(); cellchat_objects <- list()
make_edge_grid <- function(object_id, scope, condition, replicate_id, groups, active_groups, object_status) {
  grid <- expand.grid(source = groups, target = groups, interaction_name = db_interaction$interaction_name, stringsAsFactors = FALSE)
  ix <- match(grid$interaction_name, db_interaction$interaction_name)
  grid$object_id <- object_id; grid$scope <- scope; grid$condition <- condition; grid$replicate_id <- replicate_id
  grid$pathway <- db_interaction$pathway_name[ix]; grid$ligand <- db_interaction$ligand[ix]; grid$receptor <- db_interaction$receptor[ix]
  grid$assessed <- measured_lr[ix] & grid$source %in% active_groups & grid$target %in% active_groups & object_status %in% c("assessed", "zero_signal")
  grid$status <- ifelse(!measured_lr[ix], "required_gene_not_measured", ifelse(!grid$source %in% active_groups | !grid$target %in% active_groups,
                                                                               "group_below_min_cells_or_absent", object_status))
  grid$prob_raw <- grid$prob <- grid$pval <- NA_real_
  grid$supported <- NA
  grid
}
infer_network <- function(object_id, scope, condition, replicate_id, cells, labels, groups) {
  labels <- setNames(as.character(labels), cells)
  counts <- table(factor(labels, levels = groups))
  active <- groups[counts >= min_cells]
  use <- cells[labels %in% active]
  status <- if (length(active) < 2L) "insufficient_groups" else "assessed"
  object <- NULL
  if (status == "assessed") {
    cc_meta <- data.frame(labels = factor(labels[use], levels = active), samples = factor(meta[use, "replicate_id"]), row.names = use)
    object <- createCellChat(object = data_input[, use, drop = FALSE], meta = cc_meta, group.by = "labels")
    object@DB <- DB
    object <- subsetData(object)
    # 记录过表达信息供诊断，实际候选仍固定为全部共同可测LR，避免条件依赖的预筛选。
    if (run_overexpression_diagnostics) {
      object <- identifyOverExpressedGenes(object)
      object <- identifyOverExpressedInteractions(object)
    }
    object@LR$LRsig <- fixed_LR
    if (max(object@data.signaling) <= 0) status <- "zero_signal" else {
      object <- compute_probability_fixed(object, type = average_type,
                                          raw.use = TRUE, population.size = population_size, nboot = nboot, seed.use = seed)
      object <- filterCommunication(object, min.cells = min_cells)
      object <- computeCommunProbPathway(object, thresh = communication_p)
      object <- aggregateNet(object, thresh = communication_p)
    }
  }
  grid <- make_edge_grid(object_id, scope, condition, replicate_id, groups, active, status)
  if (status == "zero_signal") {
    grid$prob_raw[grid$assessed] <- 0; grid$pval[grid$assessed] <- 1
  } else if (status == "assessed") {
    net <- object@net
    inferred <- expand.grid(source = dimnames(net$prob)[[1]], target = dimnames(net$prob)[[2]],
                            interaction_name = dimnames(net$prob)[[3]], stringsAsFactors = FALSE)
    # expand.grid第一维变化最快，与R数组展开顺序一致。
    inferred$prob_raw <- as.numeric(net$prob); inferred$pval <- as.numeric(net$pval)
    key <- function(z) paste(z$source, z$target, z$interaction_name, sep = "\r")
    j <- match(key(grid), key(inferred))
    ii <- which(grid$assessed)
    if (anyNA(j[ii])) stop("固定LR未完整返回；不能把未计算的边填0。")
    grid$prob_raw[ii] <- inferred$prob_raw[j[ii]]; grid$pval[ii] <- inferred$pval[j[ii]]
  }
  grid$supported[grid$assessed] <- grid$prob_raw[grid$assessed] > 0 & grid$pval[grid$assessed] < communication_p
  grid$prob[grid$assessed] <- ifelse(grid$supported[grid$assessed], grid$prob_raw[grid$assessed], 0)
  all_edges[[object_id]] <<- grid
  # 每列为一个LR、每行为source-target；按pathway做矩阵rowsum，避免数万小data.frame。
  n_pair <- length(groups)^2L
  path_levels <- unique(db_interaction$pathway_name)
  path_sum <- function(value) {
    value[!is.finite(value)] <- 0
    t(rowsum(t(matrix(value, nrow = n_pair)), group = db_interaction$pathway_name, reorder = FALSE))
  }
  n_measured <- as.numeric(path_sum(as.integer(grid$assessed)))
  pathway <- expand.grid(source = groups, target = groups, pathway = path_levels, stringsAsFactors = FALSE)
  pathway$object_id <- object_id; pathway$scope <- scope; pathway$condition <- condition; pathway$replicate_id <- replicate_id
  pathway$n_LR_defined <- rep(as.integer(table(factor(db_interaction$pathway_name, levels = path_levels))), each = n_pair)
  pathway$n_LR_measured <- n_measured
  pathway$assessed <- n_measured > 0L
  pathway$complete_coverage <- n_measured == pathway$n_LR_defined
  pathway$prob <- ifelse(pathway$assessed, as.numeric(path_sum(grid$prob)), NA_real_)
  pathway$prob_raw <- ifelse(pathway$assessed, as.numeric(path_sum(grid$prob_raw)), NA_real_)
  pathway$n_links <- ifelse(pathway$assessed, as.numeric(path_sum(as.numeric(grid$supported))), NA_integer_)
  pathway$status <- ifelse(!pathway$assessed, "not_assessable", ifelse(pathway$complete_coverage, "assessed", "partial_LR_coverage"))
  all_pathways[[object_id]] <<- pathway
  mat <- function(column) {
    x <- if (column == "count") as.numeric(grid$supported) else grid$prob
    x[!is.finite(x)] <- 0
    value <- rowSums(matrix(x, nrow = n_pair))
    value[rowSums(matrix(grid$assessed, nrow = n_pair)) == 0L] <- NA_real_
    matrix(value, length(groups), length(groups), dimnames = list(groups, groups))
  }
  compact_networks[[object_id]] <<- list(object_id = object_id, scope = scope, condition = condition,
                                         replicate_id = replicate_id, status = status, groups = groups, active_groups = active,
                                         n_cells = setNames(as.integer(counts), groups), count = mat("count"), weight = mat("weight"))
  object_audit[[object_id]] <<- data.frame(object_id = object_id, scope = scope, condition = condition,
                                           replicate_id = replicate_id, n_cells = length(use), n_groups = length(active), n_LR = nrow(fixed_LR),
                                           status = status, stringsAsFactors = FALSE)
  if (save_cellchat_objects && !is.null(object)) cellchat_objects[[object_id]] <<- object
  message("完成网络：", object_id, "；", length(use), "细胞，", length(active), "个大类，", status)
  invisible(NULL)
}

##### Vehicle/VG161：合并展示网络与每动物网络 #####
for (condition in condition_levels) {
  cc <- rownames(meta)[meta$condition == condition & meta$main_eligible]
  infer_network(paste0("pooled__", condition), "treatment_pooled", condition, NA_character_,
                cc, meta[cc, "celltype"], main_groups)
}
for (id in replicate_table$replicate_id) {
  cc <- rownames(meta)[meta$replicate_id == id & meta$main_eligible]
  infer_network(paste0("animal__", id), "treatment_animal", replicate_table$condition[match(id, replicate_table$replicate_id)], id,
                cc, meta[cc, "celltype"], main_groups)
}

##### VG161肿瘤状态：相同免疫池及四种肿瘤内部方向 #####
cc <- rownames(meta)[meta$state_eligible]
infer_network("state_pooled__VG161", "state_pooled", "VG161", NA_character_,
              cc, meta[cc, "state_group"], state_group_levels)
for (id in vg_ids) {
  cc <- rownames(meta)[meta$replicate_id == id & meta$state_eligible]
  infer_network(paste0("state_animal__", id), "state_animal", "VG161", id,
                cc, meta[cc, "state_group"], state_group_levels)
}
edges <- bind_rows_base(all_edges)
pathways <- bind_rows_base(all_pathways)
object_audit <- bind_rows_base(object_audit)
rm(all_edges, all_pathways)

##### 动物层描述性比较：不把细胞置换P当处理差异P #####
summarize_treatment <- function(data, keys) {
  data <- data[data$scope == "treatment_animal", , drop = FALSE]
  if (!nrow(data)) return(data.frame())
  by_animal <- split(data, data$replicate_id)
  template <- by_animal[[1]][, keys, drop = FALSE]
  key <- function(z) do.call(paste, c(z[keys], sep = "\r"))
  template_key <- key(template)
  values <- vapply(by_animal, function(z) {
    j <- match(template_key, key(z)); x <- z$prob[j]
    x[is.na(j) | !z$assessed[j]] <- NA_real_; x
  }, numeric(nrow(template)))
  values <- matrix(values, nrow = nrow(template), dimnames = list(NULL, names(by_animal)))
  cond <- replicate_table$condition[match(colnames(values), replicate_table$replicate_id)]
  get_stats <- function(x) {
    valid <- is.finite(x); n <- rowSums(valid); x[!valid] <- 0
    total <- rowSums(x); average <- total / pmax(n, 1L); average[n == 0L] <- NA_real_
    variance <- (rowSums(x^2) - total^2 / pmax(n, 1L)) / pmax(n - 1L, 1L)
    spread <- sqrt(pmax(variance, 0)); spread[n < 2L] <- NA_real_
    lo <- do.call(pmin, c(lapply(seq_len(ncol(x)), function(j) ifelse(valid[, j], x[, j], Inf)), na.rm = TRUE))
    hi <- do.call(pmax, c(lapply(seq_len(ncol(x)), function(j) ifelse(valid[, j], x[, j], -Inf)), na.rm = TRUE))
    lo[n == 0L] <- hi[n == 0L] <- NA_real_
    list(n = n, mean = average, SD = spread, min = lo, max = hi)
  }
  a <- get_stats(values[, cond == "Vehicle", drop = FALSE])
  b <- get_stats(values[, cond == "VG161", drop = FALSE])
  enough <- a$n >= min_animals & b$n >= min_animals
  cbind(template, data.frame(n_Vehicle = a$n, n_VG161 = b$n,
                             mean_Vehicle = a$mean, mean_VG161 = b$mean, SD_Vehicle = a$SD, SD_VG161 = b$SD,
                             difference_VG161_minus_Vehicle = ifelse(enough, b$mean - a$mean, NA_real_),
                             range_relation = ifelse(!enough, "insufficient_animals", ifelse(b$min > a$max, "VG161_above_all_Vehicle",
                                                                                             ifelse(b$max < a$min, "VG161_below_all_Vehicle", "overlap_or_equal"))),
                             PValue = NA_real_, status = ifelse(enough, "descriptive_no_test", "insufficient_animals")))
}
paired_sender_differences <- function(data, keys) {
  z <- data[data$scope == "state_animal" & data$source %in% state_groups, , drop = FALSE]
  a <- z[z$source == state_groups[1], , drop = FALSE]
  b <- z[z$source == state_groups[2], , drop = FALSE]
  key_cols <- c("replicate_id", "target", keys)
  key <- function(x) do.call(paste, c(x[key_cols], sep = "\r"))
  j <- match(key(a), key(b))
  ok <- a$assessed & !is.na(j) & b$assessed[j]
  out <- a[, key_cols, drop = FALSE]
  out$Bystander <- a$prob; out$Infected <- b$prob[j]
  out$assessed <- ok
  out$difference_Infected_minus_Bystander <- ifelse(ok, out$Infected - out$Bystander, NA_real_)
  out$target_scope <- ifelse(out$target %in% immune_receivers, "immune", ifelse(out$target %in% state_groups, "tumor_internal", "other"))
  out$status <- ifelse(ok, "assessed", "not_assessable")
  out
}
summarize_sender <- function(data, keys) {
  if (!nrow(data)) return(data.frame())
  key_cols <- c("target", "target_scope", keys)
  by_animal <- split(data, data$replicate_id)
  template <- by_animal[[1]][, key_cols, drop = FALSE]
  key <- function(z) do.call(paste, c(z[key_cols], sep = "\r"))
  template_key <- key(template)
  values <- vapply(by_animal, function(z) {
    j <- match(template_key, key(z)); value <- z$difference_Infected_minus_Bystander[j]
    value[is.na(j) | !z$assessed[j]] <- NA_real_; value
  }, numeric(nrow(template)))
  values <- matrix(values, nrow = nrow(template))
  valid <- is.finite(values); n <- rowSums(valid)
  values[!valid] <- 0
  total <- rowSums(values); average <- total / pmax(n, 1L); average[n == 0L] <- NA_real_
  spread <- sqrt(pmax((rowSums(values^2) - total^2 / pmax(n, 1L)) / pmax(n - 1L, 1L), 0))
  spread[n < 2L] <- NA_real_
  enough <- n >= min_animals
  direction <- ifelse(!enough, "insufficient_pairs", ifelse(rowSums(values > 0 & valid) == n, "all_pairs_up",
                                                            ifelse(rowSums(values < 0 & valid) == n, "all_pairs_down", ifelse(rowSums(values == 0 & valid) == n, "all_pairs_zero", "mixed_or_zero"))))
  cbind(template, data.frame(n_pairs = n, mean_difference_Infected_minus_Bystander = average,
                             SD_difference = spread, direction_consistency = direction, PValue = NA_real_,
                             status = ifelse(enough, "descriptive_no_test", "insufficient_pairs")))
}
contrasts <- list(
  pathways = summarize_treatment(pathways, c("source", "target", "pathway")),
  interactions = summarize_treatment(edges, c("source", "target", "interaction_name", "pathway", "ligand", "receptor")))
sender_differences <- list(
  pathways = paired_sender_differences(pathways, "pathway"),
  interactions = paired_sender_differences(edges, c("interaction_name", "pathway", "ligand", "receptor")))
sender_contrasts <- list(
  pathways = summarize_sender(sender_differences$pathways, "pathway"),
  interactions = summarize_sender(sender_differences$interactions, c("interaction_name", "pathway", "ligand", "receptor")))

##### MDK/NCL原表达与紧凑结果：保留零值和未测状态 #####
expression_cells <- bind_rows_base(lapply(unique(expression_genes), function(g) {
  measured <- g %in% rownames(data_input)
  data.frame(cell = rownames(meta), sample_id = meta$sample_id, replicate_id = meta$replicate_id,
             condition = meta$condition, celltype = meta$celltype, state_group = meta$state_group,
             main_eligible = meta$main_eligible, state_eligible = meta$state_eligible, gene = g,
             measured = measured, value = if (measured) as.numeric(data_input[g, rownames(meta)]) else NA_real_, stringsAsFactors = FALSE)
}))
expression_summary <- bind_rows_base(lapply(split(expression_cells, paste(expression_cells$replicate_id, expression_cells$celltype,
                                                                          expression_cells$gene, sep = "\r")), function(z) {
                                                                            x <- z$value[is.finite(z$value)]
                                                                            eligible <- all(z$main_eligible)
                                                                            data.frame(replicate_id = z$replicate_id[1], condition = z$condition[1], celltype = z$celltype[1], gene = z$gene[1],
                                                                                       n_cells = nrow(z), eligible = eligible, feature_measured = all(z$measured),
                                                                                       mean_log_normalized = if (eligible && length(x)) mean(x) else NA_real_,
                                                                                       detected_fraction = if (eligible && length(x)) mean(x > 0) else NA_real_,
                                                                                       status = if (!eligible) "insufficient_cells" else if (!all(z$measured)) "gene_not_measured" else "assessed")
                                                                          }))
merged_ids <- paste0("pooled__", condition_levels)
individual_ids <- paste0("animal__", replicate_table$replicate_id)
state_individual_ids <- paste0("state_animal__", vg_ids)
info <- list(started_at = as.character(started_at), analysis_label = analysis_label,
             input_file = path.expand(input_file), input_object = input_object, outdir = outdir, sample_table = sample_table,
             scope = "broad_celltypes_only", celltype_column = celltype_col, sample_column = sample_col,
             state_column = state_col, tumor_labels = tumor_labels, immune_receivers = immune_receivers,
             celltype_order = main_groups, tumor_groups = state_groups, state_groups = state_groups,
             state_status = state_status, complete_pairs = complete_pairs, database_scope = db_scope,
             database_interactions = nrow(db_interaction), fixed_candidates = nrow(fixed_LR),
             global_signaling_max = global_signaling_max, scale_patch_replacements = replacement_count,
             probability_scale = probability_scale, probability_function_md5 = probability_function_md5,
             probability_function_original = deparse(original_probability_body),
             scale_method = if (probability_scale == "shared_max") "Installed computeCommunProb copied; only data/max(data) denominator replaced with global signaling maximum" else "Native per-object maximum; comparison uses object-specific scales",
             candidate_method = "All DB ligand-receptor pairs with required ligand/receptor subunits measured in the common RNA feature background",
             min_cells = min_cells, min_animals = min_animals, communication_p = communication_p,
             nboot = nboot, average_type = average_type, population_size = population_size, raw_use = TRUE, seed = seed,
             focus_pathway = focus_pathway, focus_ligand = focus_ligand, focus_receptor = focus_receptor,
             host_prefix = host_prefix, viral_pattern = viral_pattern, viral_feature_ids = viral_feature_ids,
             expression_mode = expression_mode,
             expression_unit = if (expression_mode == "host_logCP10K") "log1p(10000 * RNA count / all-common-host count total); viral features excluded from numerator/background and denominator" else "Existing RNA LogNormalize data; viral features removed but normalization denominator not recomputed",
             inference = "CellChat permutation P is not animal-level treatment differential P; animal comparisons are descriptive with biological repeats",
             save_cellchat_objects = save_cellchat_objects,
             versions = vapply(c("Seurat", "SeuratObject", "CellChat", "Matrix", "ggplot2"), function(p) as.character(packageVersion(p)), character(1)),
             session = capture.output(sessionInfo()))
results <- list(info = info,
                audit = list(objects = object_audit, database = db_audit, layers = layer_audit,
                             main_counts = main_counts, state_counts = state_counts, zero_expression_cells = zero_expression_cells),
                merged = setNames(compact_networks[merged_ids], condition_levels),
                individual = setNames(compact_networks[individual_ids], replicate_table$replicate_id),
                state_merged = compact_networks[["state_pooled__VG161"]],
                state_individual = setNames(compact_networks[state_individual_ids], vg_ids),
                edges = edges, pathways = pathways,
                mk_edges = edges[edges$pathway == focus_pathway, , drop = FALSE],
                focus_edges = edges[edges$ligand == focus_ligand & edges$receptor == focus_receptor, , drop = FALSE],
                intratumor = list(edges = edges[edges$scope %in% c("state_pooled", "state_animal") & edges$source %in% state_groups & edges$target %in% state_groups, , drop = FALSE],
                                  pathways = pathways[pathways$scope %in% c("state_pooled", "state_animal") & pathways$source %in% state_groups & pathways$target %in% state_groups, , drop = FALSE]),
                contrasts = contrasts, sender_differences = sender_differences, sender_contrasts = sender_contrasts,
                expression = expression_summary, expression_cells = expression_cells)
if (save_cellchat_objects) results$cellchat_objects <- cellchat_objects
# 后续由绘图和保存段导出；不写回输入sce，不保存重复的大型Seurat对象。

##### 绘图参数：标题色条、组色和细胞谱系色分别管理 #####
header_colors <- c(Overview = "#276C9F", Pathways = "#F0B174", Immune = "#38923C",
                   Tumor = "#785E96", MK = "#A7BDDB", MDK_NCL = "#397C73", MDK = "#785E96", NCL = "#397C73")
celltype_colors <- c(Epithelial = "#426B91", Fibroblast = "#F0B174", Endothelial = "#38923C",
                     Tcell = "#C35B3F", Myeloid = "#785E96", Mast = "#A17418", Other = "#6D6F75",
                     Tumor_Bystander = "#426B91", Tumor_Infected = "#C35B3F")
pathways_per_page <- 28L
plots_width <- 12
plots_height <- 9
na_color <- "#B8BDC3"
plot_theme <- theme_bw(base_size = 10, base_family = "sans") + theme(
  panel.grid = element_blank(), panel.border = element_rect(colour = "#454545", linewidth = 0.4),
  axis.ticks = element_line(colour = "#454545", linewidth = 0.35),
  strip.text = element_text(face = "bold", size = 10), strip.background = element_rect(linewidth = 0.35),
  plot.title = element_text(face = "bold", size = 13), legend.position = "right")
cc_strip <- function(p, label, module) {
  fill <- if (grepl("^#", module)) module else unname(header_colors[module])
  rgb <- grDevices::col2rgb(fill)[, 1]
  foreground <- if (sum(rgb * c(0.299, 0.587, 0.114)) > 155) "#252525" else "white"
  p$data$.cc_panel <- label
  p + facet_wrap(~ .cc_panel) + theme(strip.background = element_rect(fill = fill, colour = "#454545", linewidth = 0.35),
                                      strip.text = element_text(colour = foreground, face = "bold"))
}
cc_empty <- function(label, status, module = "Overview") {
  p <- ggplot(data.frame(x = 0, y = 0, label = status), aes(x, y, label = label)) +
    geom_text(size = 4, colour = "#6D6F75") + theme_void(base_family = "sans")
  cc_strip(p, label, module)
}
cc_limit <- function(x, symmetric = FALSE) {
  x <- x[is.finite(x)]
  upper <- if (length(x)) max(if (symmetric) abs(x) else x) else 0
  if (!is.finite(upper) || upper <= 0) upper <- 1
  if (symmetric) c(-upper, upper) else c(0, upper)
}
cc_label <- function(x) gsub("^Tumor_", "", x)
cc_sum <- function(x) if (any(is.finite(x))) sum(x[is.finite(x)]) else NA_real_
cc_aggregate <- function(data, keys, value = "prob") {
  if (!nrow(data)) return(data.frame())
  rows <- split(data, do.call(paste, c(data[keys], sep = "\r")))
  bind_rows_base(lapply(rows, function(z) cbind(z[1, keys, drop = FALSE], data.frame(value = cc_sum(z[[value]])))))
}
cc_network_table <- function(data, value = "prob") {
  cc_aggregate(data, c("object_id", "scope", "condition", "replicate_id", "source", "target"), value)
}
cc_heatmap <- function(data, label, module, limits, value = "value", diverging = FALSE, legend_name = "Inferred score") {
  if (!nrow(data)) return(cc_empty(label, "Not assessable", module))
  levels_use <- unique(c(results$info$celltype_order, as.character(data$source), as.character(data$target)))
  levels_use <- levels_use[levels_use %in% c(as.character(data$source), as.character(data$target))]
  data$source <- factor(data$source, levels = rev(levels_use))
  data$target <- factor(data$target, levels = levels_use)
  data$value_plot <- data[[value]]
  p <- ggplot(data, aes(target, source, fill = value_plot)) + geom_tile(colour = "white", linewidth = 0.25) +
    scale_x_discrete(labels = cc_label, drop = FALSE) + scale_y_discrete(labels = cc_label, drop = FALSE) +
    labs(x = "Receiver", y = "Sender") + plot_theme +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.title = element_text(size = 9))
  if (diverging) p <- p + scale_fill_gradient2(low = "#426B91", mid = "#F7F7F5", high = "#C35B3F",
                                               midpoint = 0, limits = limits, na.value = na_color, name = "Difference") else
                                                 p <- p + scale_fill_gradientn(colours = c("#F3F7FB", "#91B3CB", "#2D5375"),
                                                                               limits = limits, na.value = na_color, name = legend_name)
  missing <- data[!is.finite(data$value_plot), , drop = FALSE]
  if (nrow(missing)) p <- p + geom_point(data = missing, aes(shape = "Not assessable"), colour = "#6D6F75", size = 1.5) +
    scale_shape_manual(values = c("Not assessable" = 4), name = NULL)
  cc_strip(p, label, module)
}
cc_circle <- function(data, label, module, limits) {
  if (!nrow(data)) return(cc_empty(label, "Not assessable", module))
  types <- unique(c(as.character(data$source), as.character(data$target)))
  types <- unique(c(results$info$celltype_order[results$info$celltype_order %in% types], types))
  angle <- seq(pi / 2, pi / 2 + 2 * pi, length.out = length(types) + 1L)[seq_along(types)]
  vertices <- data.frame(celltype = types, x = cos(angle), y = sin(angle))
  palette <- celltype_colors
  extra <- setdiff(types, names(palette))
  if (length(extra)) palette <- c(palette, setNames(grDevices::hcl.colors(length(extra), "Dark 3"), extra))
  lines <- data[is.finite(data$value) & data$value > 0 & data$source != data$target, , drop = FALSE]
  a <- match(lines$source, vertices$celltype); b <- match(lines$target, vertices$celltype)
  lines$x <- vertices$x[a]; lines$y <- vertices$y[a]; lines$xend <- vertices$x[b]; lines$yend <- vertices$y[b]
  distance <- sqrt((lines$xend - lines$x)^2 + (lines$yend - lines$y)^2)
  ux <- (lines$xend - lines$x) / distance; uy <- (lines$yend - lines$y) / distance
  lines$x <- lines$x + 0.07 * ux; lines$y <- lines$y + 0.07 * uy
  lines$xend <- lines$xend - 0.09 * ux; lines$yend <- lines$yend - 0.09 * uy
  self <- data[is.finite(data$value) & data$value > 0 & data$source == data$target, , drop = FALSE]
  loops <- bind_rows_base(lapply(seq_len(nrow(self)), function(i) {
    j <- match(self$source[i], vertices$celltype)
    theta <- seq(angle[j] + pi * 0.9, angle[j] + pi * 2.7, length.out = 50)
    data.frame(x = 1.15 * vertices$x[j] + 0.14 * cos(theta), y = 1.15 * vertices$y[j] + 0.14 * sin(theta),
               value = self$value[i], source = self$source[i], id = i)
  }))
  p <- ggplot(vertices, aes(x, y))
  if (nrow(lines)) p <- p + geom_curve(data = lines, aes(x = x, y = y, xend = xend, yend = yend, linewidth = value, colour = source),
                                       inherit.aes = FALSE, alpha = 0.72, curvature = 0.14,
                                       arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed"))
  if (nrow(loops)) p <- p + geom_path(data = loops, aes(x, y, group = id, linewidth = value, colour = source), inherit.aes = FALSE,
                                      alpha = 0.72, arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed"))
  vertices$available <- vapply(vertices$celltype, function(x) any(is.finite(data$value[data$source == x | data$target == x])), logical(1))
  p <- p + geom_point(aes(fill = celltype), shape = 21, size = 4.8, colour = "#454545", stroke = 0.4) +
    geom_text(aes(x = 1.38 * x, y = 1.38 * y, label = cc_label(celltype)), size = 3.3) +
    scale_fill_manual(values = palette, name = "Cell type") +
    scale_colour_manual(values = palette, limits = types, guide = "none") +
    scale_linewidth_continuous(limits = limits, range = c(0.25, 2.2), name = "Inferred score") +
    coord_equal(xlim = c(-1.8, 1.8), ylim = c(-1.65, 1.65), clip = "off") +
    theme_void(base_family = "sans") + theme(legend.position = "right", plot.margin = margin(8, 18, 8, 18))
  unavailable <- vertices[!vertices$available, , drop = FALSE]
  if (nrow(unavailable)) p <- p + geom_point(data = unavailable, aes(shape = "Not assessable"),
                                             fill = na_color, colour = "#454545", size = 4.8) + scale_shape_manual(values = c("Not assessable" = 4), name = NULL)
  if (!nrow(lines) && !nrow(self)) p <- p + annotate("text", x = 0, y = 0,
                                                     label = if (any(is.finite(data$value))) "No inferred edges" else "Not assessable", colour = "#6D6F75", size = 3.5)
  cc_strip(p, label, module)
}
cc_difference_matrix <- function(data) {
  a <- data[data$condition == condition_levels[1], , drop = FALSE]
  b <- data[data$condition == condition_levels[2], , drop = FALSE]
  j <- match(paste(a$source, a$target), paste(b$source, b$target))
  a$value <- b$value[j] - a$value
  a
}
cc_flow <- function(data, pooled_scope) {
  z <- data[data$scope == pooled_scope, , drop = FALSE]
  if (!nrow(z)) return(data.frame())
  # 全通路流量只比较两组均可评估的source-target对，不把缺失群当作零通讯。
  if (pooled_scope == "treatment_pooled") {
    a <- z[z$condition == condition_levels[1], , drop = FALSE]
    b <- z[z$condition == condition_levels[2], , drop = FALSE]
    key <- function(x) paste(x$source, x$target, x$pathway, sep = "\r")
    j <- match(key(a), key(b)); ok <- is.finite(a$prob) & !is.na(j) & is.finite(b$prob[j])
    z <- rbind(a[ok, , drop = FALSE], b[j[ok], , drop = FALSE])
    if (!nrow(z)) return(data.frame())
    cc_aggregate(z, c("condition", "pathway"))
  } else cc_aggregate(z, c("source", "target", "pathway"))
}
cc_flow_pages <- function(data, label, module = "Pathways") {
  if (!nrow(data)) { print(cc_empty(label, "No jointly assessable pathways", module)); return(invisible(NULL)) }
  by_path <- split(data$value, data$pathway)
  order_use <- names(sort(vapply(by_path, function(z) if (all(is.finite(z))) max(z) - min(z) else -Inf, numeric(1)), decreasing = TRUE))
  for (page in split(order_use, ceiling(seq_along(order_use) / pathways_per_page))) {
    z <- data[data$pathway %in% page, , drop = FALSE]
    z$pathway <- factor(z$pathway, levels = rev(page))
    z$condition <- factor(z$condition, levels = condition_levels)
    p <- ggplot(z, aes(value, pathway, colour = condition, shape = condition)) +
      geom_point(size = 2.2, position = position_dodge(width = 0.35), na.rm = TRUE) +
      scale_colour_manual(values = condition_colors, drop = FALSE) + scale_shape_manual(values = c(16, 17), drop = FALSE) +
      scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.08))) +
      labs(x = "Inferred information flow", y = NULL, colour = NULL, shape = NULL) + plot_theme
    print(cc_strip(p, label, module))
  }
}
cc_state_pathway_pages <- function(data, targets, label, module) {
  z <- data[data$scope == "state_pooled" & data$source %in% state_groups & data$target %in% targets, , drop = FALSE]
  if (!nrow(z) || !any(is.finite(z$prob))) {
    print(cc_empty(label, "No assessable pathways", module)); return(invisible(NULL))
  }
  strength <- vapply(split(z$prob, z$pathway), function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_, numeric(1))
  order_use <- names(sort(strength[is.finite(strength)], decreasing = TRUE))
  columns <- as.vector(outer(cc_label(state_groups), cc_label(targets), paste, sep = " -> "))
  z$pair <- factor(paste(cc_label(z$source), cc_label(z$target), sep = " -> "), levels = columns)
  limits <- cc_limit(z$prob)
  for (page in split(order_use, ceiling(seq_along(order_use) / pathways_per_page))) {
    a <- z[z$pathway %in% page, , drop = FALSE]
    a$pathway <- factor(a$pathway, levels = rev(page))
    p <- ggplot(a, aes(pair, pathway, fill = prob)) + geom_tile(colour = "white", linewidth = 0.2) +
      scale_fill_gradientn(colours = c("#F3F7FB", "#91B3CB", "#2D5375"), limits = limits, na.value = na_color,
                           name = "Inferred score") + scale_x_discrete(drop = FALSE) +
      labs(x = NULL, y = NULL) + plot_theme + theme(axis.text.x = element_text(angle = 35, hjust = 1))
    missing <- a[!is.finite(a$prob), , drop = FALSE]
    if (nrow(missing)) p <- p + geom_point(data = missing, aes(shape = "Not assessable"), size = 1.2, colour = "#6D6F75") +
      scale_shape_manual(values = c("Not assessable" = 4), name = NULL)
    print(cc_strip(p, label, module))
  }
}
cc_bubbles <- function(data, label, module, xlevels, colors, limits) {
  if (!nrow(data) || !any(is.finite(data$prob))) return(cc_empty(label, "Not assessable", module))
  data$x <- factor(data$display_group, levels = xlevels)
  data$interaction <- factor(gsub("_", "-", data$interaction_name), levels = rev(unique(gsub("_", "-", data$interaction_name))))
  p <- ggplot(data, aes(x, interaction)) +
    geom_point(data = data[is.finite(data$prob) & data$prob > 0, , drop = FALSE], aes(size = prob, colour = display_group)) +
    geom_text(data = data[is.finite(data$prob) & data$prob == 0, , drop = FALSE], label = "0", size = 2.7, colour = "#7C8388") +
    scale_size_area(max_size = 6, limits = limits, name = "Inferred score") + scale_colour_manual(values = colors, limits = xlevels) +
    scale_x_discrete(labels = cc_label, drop = FALSE) +
    labs(x = NULL, y = NULL) + plot_theme + theme(legend.position = "right")
  missing <- data[!is.finite(data$prob), , drop = FALSE]
  if (nrow(missing)) p <- p + geom_point(data = missing, aes(shape = "Not assessable"), colour = na_color, size = 2) +
    scale_shape_manual(values = c("Not assessable" = 4), name = NULL)
  cc_strip(p, label, module)
}
cc_animal_bars <- function(data, label, module, colors, groups, paired = FALSE) {
  if (!nrow(data) || !any(is.finite(data$value))) return(cc_empty(label, "Not assessable", module))
  data$group <- factor(data$display_group, levels = groups)
  summaries <- bind_rows_base(lapply(groups, function(g) {
    x <- data$value[data$display_group == g & is.finite(data$value)]
    data.frame(group = factor(g, levels = groups), n = length(x), mean = if (length(x)) mean(x) else NA_real_,
               SD = if (length(x) >= 2L) sd(x) else NA_real_)
  }))
  summaries$lower <- summaries$mean - summaries$SD; summaries$upper <- summaries$mean + summaries$SD
  p <- ggplot(summaries, aes(group, mean, fill = group)) +
    geom_col(width = 0.62, alpha = 0.78, colour = "#454545", linewidth = 0.35, na.rm = TRUE) +
    geom_errorbar(data = summaries[is.finite(summaries$SD), , drop = FALSE], aes(ymin = lower, ymax = upper), width = 0.12, linewidth = 0.35) +
    scale_fill_manual(values = colors, limits = groups) + scale_x_discrete(labels = cc_label, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    labs(x = NULL, y = "Animal mean inferred score") + plot_theme + theme(legend.position = "none")
  points <- data[is.finite(data$value), , drop = FALSE]
  if (paired && nrow(points)) {
    complete_ids <- names(which(table(points$replicate_id) == length(groups)))
    paired_points <- points[points$replicate_id %in% complete_ids, , drop = FALSE]
    p <- p + geom_line(data = paired_points, aes(x = group, y = value, group = replicate_id), inherit.aes = FALSE,
                       colour = "#C4C7CA", linewidth = 0.25)
  }
  p <- p + geom_point(data = points, aes(group, value, colour = group), inherit.aes = FALSE,
                      position = if (paired) position_identity() else position_jitter(width = 0.05, height = 0, seed = seed), size = 2) +
    scale_colour_manual(values = colors, limits = groups)
  bad <- summaries[summaries$n < min_animals, , drop = FALSE]
  if (nrow(bad)) p <- p + geom_text(data = bad, aes(x = group, y = Inf, label = ifelse(n == 0, "NA", "n = 1 animal")),
                                    inherit.aes = FALSE, vjust = 1.2, colour = "#6D6F75", size = 3)
  cc_strip(p, label, module)
}

cc_expression_violin <- function(data, gene, label, colors, groups, module, xfield, y_limits = NULL) {
  data <- data[is.finite(data$value), , drop = FALSE]
  if (!nrow(data)) return(cc_empty(label, "Not measured or no eligible cells", module))
  data$group <- factor(data[[xfield]], levels = groups)
  densities <- bind_rows_base(Filter(function(z) nrow(z) >= 2L && length(unique(z$value)) >= 2L,
                                     split(data, data$group, drop = TRUE)))
  p <- ggplot(data, aes(group, value, fill = group))
  if (nrow(densities)) p <- p + geom_violin(data = densities, scale = "width", trim = TRUE,
                                            alpha = 0.78, colour = "#454545", linewidth = 0.35)
  p <- p + geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.35, alpha = 0.9) +
    scale_fill_manual(values = colors, limits = groups) + scale_x_discrete(drop = FALSE, labels = cc_label) +
    labs(x = NULL, y = if (results$info$expression_mode == "host_logCP10K") "log1p(host CP10K)" else "RNA log-normalized expression") +
    plot_theme + theme(legend.position = "none")
  n <- vapply(groups, function(g) length(unique(data$replicate_id[as.character(data$group) == g])), integer(1))
  missing <- data.frame(group = factor(groups[n < 2L], levels = groups),
                        label = ifelse(n[n < 2L] == 0L, "NA", "n = 1 animal"))
  if (nrow(missing)) p <- p + geom_text(data = missing, aes(group, Inf, label = label),
                                        inherit.aes = FALSE, vjust = 1.2, size = 3, colour = "#6D6F75")
  if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
  cc_strip(p, label, module)
}
cc_internal_bars <- function(data, label, module) {
  a <- data[data$scope == "state_animal" & data$source %in% state_groups & data$target %in% state_groups, , drop = FALSE]
  a$display_group <- paste(cc_label(a$source), cc_label(a$target), sep = " -> ")
  groups <- as.vector(outer(cc_label(state_groups), cc_label(state_groups), paste, sep = " -> "))
  colors <- setNames(rep(unname(state_colors), 2L), groups)
  a$value <- a$prob
  cc_animal_bars(a, label, module, colors, groups)
}

##### 全通路、MK和MDK-NCL：分开输出，沿用同组合并网络 #####
draw_cellchat_plots <- function() {
  treatment <- results$pathways[results$pathways$scope == "treatment_pooled", , drop = FALSE]
  overall <- cc_network_table(treatment)
  pdf(file.path(outdir, "plots.pdf"), width = plots_width, height = plots_height, onefile = TRUE, useDingbats = FALSE)
  tryCatch({
    network_panels <- lapply(condition_levels, function(g) cc_circle(overall[overall$condition == g, , drop = FALSE],
                                                                     g, "Overview", cc_limit(overall$value)))
    print(wrap_plots(network_panels, ncol = 2, guides = "collect") + plot_annotation(title = "All pathways"))
    heat_panels <- lapply(condition_levels, function(g) cc_heatmap(overall[overall$condition == g, , drop = FALSE],
                                                                   g, "Overview", cc_limit(overall$value)))
    print(wrap_plots(heat_panels, ncol = 2, guides = "collect") + plot_annotation(title = "All pathways"))
    counts <- cc_network_table(treatment, "n_links")
    print(wrap_plots(lapply(condition_levels, function(g) cc_heatmap(counts[counts$condition == g, , drop = FALSE],
                                                                     g, "Overview", cc_limit(counts$value), legend_name = "Supported LR entries")), ncol = 2, guides = "collect") +
            plot_annotation(title = "Number of inferred interactions"))
    difference <- cc_difference_matrix(overall)
    print(cc_heatmap(difference, "VG161 - Vehicle", "Overview", cc_limit(difference$value, TRUE), diverging = TRUE))
    cc_flow_pages(cc_flow(results$pathways, "treatment_pooled"), "All pathways")
    cc_state_pathway_pages(results$pathways, immune_receivers, "All pathways | Tumor states -> immune", "Immune")
    cc_state_pathway_pages(results$pathways, state_groups, "All pathways | Tumor internal communication", "Tumor")
    state <- results$intratumor$pathways[results$intratumor$pathways$scope == "state_pooled", , drop = FALSE]
    if (nrow(state)) {
      tumor_net <- cc_network_table(state)
      print(cc_heatmap(tumor_net, "Tumor internal communication", "Tumor", cc_limit(tumor_net$value)))
    } else print(cc_empty("Tumor internal communication", results$info$state_status, "Tumor"))
  }, finally = dev.off())
  
  pdf(file.path(outdir, "MK.pdf"), width = plots_width, height = plots_height, onefile = TRUE, useDingbats = FALSE)
  tryCatch({
    mk <- results$pathways[results$pathways$pathway == focus_pathway, , drop = FALSE]
    mk_limits <- cc_limit(mk$prob[mk$scope %in% c("treatment_pooled", "state_pooled")])
    mk_lr_limits <- cc_limit(results$mk_edges$prob[results$mk_edges$scope %in% c("treatment_pooled", "state_pooled")])
    pooled <- mk[mk$scope == "treatment_pooled", , drop = FALSE]
    pooled$value <- pooled$prob
    print(wrap_plots(lapply(condition_levels, function(g) cc_heatmap(pooled[pooled$condition == g, , drop = FALSE],
                                                                     g, "MK", mk_limits)), ncol = 2, guides = "collect") + plot_annotation(title = focus_pathway))
    print(wrap_plots(lapply(condition_levels, function(g) cc_circle(pooled[pooled$condition == g, , drop = FALSE],
                                                                    g, "MK", mk_limits)), ncol = 2, guides = "collect") + plot_annotation(title = focus_pathway))
    state <- mk[mk$scope == "state_pooled", , drop = FALSE]; state$value <- state$prob
    print(cc_heatmap(state, "MK | VG161 tumor states", "MK", mk_limits))
    for (target in unique(c(immune_receivers, tumor_labels))) {
      z <- results$mk_edges[results$mk_edges$scope == "treatment_pooled" & results$mk_edges$source %in% tumor_labels &
                              results$mk_edges$target == target, , drop = FALSE]
      z$display_group <- z$condition
      print(cc_bubbles(z, paste("MK | Tumor ->", target), "MK", condition_levels, condition_colors, mk_lr_limits))
      a <- mk[mk$scope == "treatment_animal" & mk$source %in% tumor_labels & mk$target == target, , drop = FALSE]
      a$display_group <- a$condition; a$value <- a$prob
      print(cc_animal_bars(a, paste("MK | Tumor ->", target), "MK", condition_colors, condition_levels))
    }
    for (target in immune_receivers) {
      z <- results$mk_edges[results$mk_edges$scope == "state_pooled" & results$mk_edges$source %in% state_groups &
                              results$mk_edges$target == target, , drop = FALSE]
      z$display_group <- z$source
      print(cc_bubbles(z, paste("MK | Tumor states ->", target), "Immune", state_groups, state_colors, mk_lr_limits))
      a <- mk[mk$scope == "state_animal" & mk$source %in% state_groups & mk$target == target, , drop = FALSE]
      a$display_group <- a$source; a$value <- a$prob
      print(cc_animal_bars(a, paste("MK | Tumor states ->", target), "Immune", state_colors, state_groups, paired = TRUE))
    }
    z <- mk[mk$scope == "state_pooled" & mk$source %in% state_groups & mk$target %in% state_groups, , drop = FALSE]
    z$value <- z$prob
    print(cc_heatmap(z, "MK | Tumor internal communication", "Tumor", mk_limits))
    print(cc_internal_bars(mk, "MK | Tumor internal animal summary", "Tumor"))
  }, finally = dev.off())
  
  pdf(file.path(outdir, "MDK_NCL.pdf"), width = plots_width, height = plots_height, onefile = TRUE, useDingbats = FALSE)
  tryCatch({
    pair <- results$focus_edges
    pair_limits <- cc_limit(pair$prob[pair$scope %in% c("treatment_pooled", "state_pooled")])
    p <- pair[pair$scope == "treatment_pooled", , drop = FALSE]; p$value <- p$prob
    print(wrap_plots(lapply(condition_levels, function(g) cc_heatmap(p[p$condition == g, , drop = FALSE],
                                                                     g, "MDK_NCL", pair_limits)), ncol = 2, guides = "collect") + plot_annotation(title = "MDK-NCL"))
    z <- pair[pair$scope == "state_pooled", , drop = FALSE]; z$value <- z$prob
    print(cc_heatmap(z, "MDK-NCL | VG161 tumor states", "MDK_NCL", pair_limits))
    for (target in unique(c(immune_receivers, tumor_labels))) {
      a <- pair[pair$scope == "treatment_animal" & pair$source %in% tumor_labels & pair$target == target, , drop = FALSE]
      a$display_group <- a$condition; a$value <- a$prob
      print(cc_animal_bars(a, paste("MDK-NCL | Tumor ->", target), "MDK_NCL", condition_colors, condition_levels))
    }
    for (target in immune_receivers) {
      a <- pair[pair$scope == "state_animal" & pair$source %in% state_groups & pair$target == target, , drop = FALSE]
      a$display_group <- a$source; a$value <- a$prob
      print(cc_animal_bars(a, paste("MDK-NCL | Tumor states ->", target), "Immune", state_colors, state_groups, paired = TRUE))
    }
    z <- pair[pair$scope == "state_pooled" & pair$source %in% state_groups & pair$target %in% state_groups, , drop = FALSE]
    z$value <- z$prob
    print(cc_heatmap(z, "MDK-NCL | Tumor internal communication", "Tumor", pair_limits))
    print(cc_internal_bars(pair, "MDK-NCL | Tumor internal animal summary", "Tumor"))
    for (gene in expression_genes) {
      z <- results$expression_cells[results$expression_cells$gene == gene & results$expression_cells$main_eligible, , drop = FALSE]
      gene_limits <- cc_limit(z$value)
      types <- unique(z$celltype)
      plots <- lapply(types, function(type) {
        a <- z[z$celltype == type, , drop = FALSE]
        cc_expression_violin(a, gene, paste(gene, type, sep = " | "), condition_colors, condition_levels,
                             if (type %in% names(celltype_colors)) unname(celltype_colors[type]) else gene, "condition", gene_limits)
      })
      if (length(plots)) print(wrap_plots(plots, ncol = 3) + plot_annotation(title = gene)) else
        print(cc_empty(gene, "No eligible cells", gene))
      a <- results$expression_cells[results$expression_cells$gene == gene & results$expression_cells$state_eligible &
                                      results$expression_cells$celltype %in% tumor_labels & results$expression_cells$state_group %in% state_groups, , drop = FALSE]
      print(cc_expression_violin(a, gene, paste(gene, "| Tumor states"), state_colors, state_groups, gene, "state_group", gene_limits))
    }
  }, finally = dev.off())
}

##### 保存紧凑结果与比较表 #####
results$info$header_colors <- header_colors
results$info$celltype_colors <- celltype_colors
results$info$plot_scope <- "Whole-network totals sum assessable fixed-universe entries. All-pathway pages retain zero scores; fully unassessable pathways remain in results tables but are not ranked in the figures. State overview compares both tumor senders to the same immune pool and all four tumor-internal directions."
results$info$plot_references <- c(
  "Shen et al., Gut, doi:10.1136/gutjnl-2025-335904, Figures 3/5/6/7: narrow colored headers, separate condition colors",
  "CellChat, Nat Commun 2021, doi:10.1038/s41467-021-21246-9, Figure 1: network, pathway and LR views",
  "Tang et al., Cancer Cell 2026, doi:10.1016/j.ccell.2026.05.012, Figure 3: source-target matrix and summary encodings")
results$info$finished_at <- as.character(Sys.time())
saveRDS(results, file.path(outdir, "results.rds"), compress = TRUE)
summary_parts <- list()
for (kind in names(results$contrasts)) {
  z <- results$contrasts[[kind]]
  if (nrow(z)) { z$comparison <- "VG161_minus_Vehicle"; z$level <- kind; summary_parts[[length(summary_parts) + 1L]] <- z }
}
for (kind in names(results$sender_contrasts)) {
  z <- results$sender_contrasts[[kind]]
  if (nrow(z)) { z$comparison <- "Infected_minus_Bystander"; z$level <- kind; summary_parts[[length(summary_parts) + 1L]] <- z }
}
all_fields <- unique(unlist(lapply(summary_parts, names)))
summary_table <- bind_rows_base(lapply(summary_parts, function(z) {
  for (field in setdiff(all_fields, names(z))) z[[field]] <- NA
  z[, all_fields, drop = FALSE]
}))
write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
draw_cellchat_plots()
message("完成：", normalizePath(outdir, winslash = "/"), "；plots.pdf、MK.pdf、MDK_NCL.pdf、summary.csv、results.rds。")
