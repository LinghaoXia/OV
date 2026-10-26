library(readxl)
library(DESeq2)
library(ggplot2)
library(ggrepel)


rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "TCT"
input_dir <- "~/rawdata/TCT/CRISPR"
out_dir = paste0("~/OV/",gene)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# -------------------- 1. 文件和重点基因 --------------------
targets <- c(
  "IRF1", "ZSCAN2", "ANXA7", "ACYP1", "USP44",
  "RHBDF1", "HDAC11", "CPSF4L", "KCNC2", "GATSL3"
)

files <- c(
  Control = "对照组.xlsx",
  D2 = "治疗组-D2.xlsx",
  D3 = "治疗组-D3.xlsx",
  D7 = "治疗组-D7.xlsx"
)
input_files <- file.path(input_dir, files)

sample_cols <- paste0("count", 1:10)

# 统计及配色阈值，可在这里统一修改
fdr_cutoff <- 0.05
enriched_z_cutoff <- 1.5
depleted_z_cutoff <- -1.5

# 低计数过滤：至少3个样本中计数达到10
min_count <- 10
min_samples <- 3

# -------------------- 2. 读取并合并40个样本 --------------------
tabs <- lapply(
  input_files,
  function(f) as.data.frame(read_excel(f))
)
names(tabs) <- names(files)

genes <- as.character(tabs[[1]]$gene_id)


matrices <- lapply(names(files), function(g) {
  
  x <- tabs[[g]]
  
  required_columns <- c("gene_id", sample_cols)
  
  if (!all(required_columns %in% names(x))) {
    stop("文件缺少gene_id或count1-count10：", files[[g]])
  }
  
  if (
    anyDuplicated(x$gene_id) ||
    !setequal(as.character(x$gene_id), genes)
  ) {
    stop("基因集合不一致或存在重复gene_id：", files[[g]])
  }
  
  # 按基因ID对齐，不依赖原表行顺序
  a <- as.matrix(
    x[
      match(genes, as.character(x$gene_id)),
      sample_cols
    ]
  )
  
  if (
    !is.numeric(a) ||
    any(!is.finite(a)) ||
    any(a < 0) ||
    any(abs(a - round(a)) > 1e-8)
  ) {
    stop(
      files[[g]],
      "：count列必须是非负、无缺失的原始整数计数。"
    )
  }
  
  if (any(a > .Machine$integer.max)) {
    stop("计数超出R整数范围：", files[[g]])
  }
  
  colnames(a) <- paste(g, sample_cols, sep = "_")
  a
})

cts <- do.call(cbind, matrices)
rownames(cts) <- genes

# 已验证为整数；这里仅转换存储类型
storage.mode(cts) <- "integer"

if (any(colSums(cts) == 0)) {
  stop("存在总计数为0的样本，请检查输入。")
}

meta <- data.frame(
  group = factor(
    rep(names(files), each = 10),
    levels = c("Control", "D2", "D3", "D7")
  ),
  row.names = colnames(cts)
)

stopifnot(identical(colnames(cts), rownames(meta)))

write.csv(
  data.frame(sample = rownames(meta), meta),
  file.path(out_dir, "sample_metadata.csv"),
  row.names = FALSE
)

# -------------------- 3. DESeq2分析 --------------------
# 不先平均10个样本，也不要求处理组全部非零
keep <- rowSums(cts >= min_count) >= min_samples

if (sum(keep) < 3) {
  stop("通过低计数过滤的基因不足。")
}

message(
  "输入基因数：", nrow(cts),
  "；过滤后基因数：", sum(keep)
)

dds <- DESeqDataSetFromMatrix(
  countData = cts[keep, , drop = FALSE],
  colData = meta,
  design = ~ group
)

# 针对零值较多的计数矩阵使用poscounts
dds <- DESeq(dds, sfType = "poscounts")

saveRDS(
  dds,
  file.path(out_dir, "DESeq2_model.rds")
)

write.csv(
  data.frame(
    sample = names(sizeFactors(dds)),
    size_factor = sizeFactors(dds)
  ),
  file.path(out_dir, "size_factors.csv"),
  row.names = FALSE
)

# -------------------- 4. 三个时间点分别与对照比较 --------------------

days <- c("D2", "D3", "D7")
contrast_results <- setNames(vector("list", length(days)), days)

dat <- data.frame(
  gene = rownames(dds),
  stringsAsFactors = FALSE
)

for (day in days) {
  
  r <- as.data.frame(
    results(
      dds,
      contrast = c("group", day, "Control"),
      alpha = fdr_cutoff
    )
  )
  
  r <- r[dat$gene, , drop = FALSE]
  contrast_results[[day]] <- r
  
  dat[[paste0("LFC_", day)]] <- r$log2FoldChange
  dat[[paste0("SE_", day)]] <- r$lfcSE
  dat[[paste0("P_", day)]] <- r$pvalue
  dat[[paste0("FDR_", day)]] <- r$padj
}

# 保留默认离群值处理和独立过滤行为
# FDR为NA时，不能当作显著，也不能当作0

# -------------------- 5. 计算LFC Z-score --------------------

# 三个比较使用相同的有效基因集合
finite <- apply(
  dat[paste0("LFC_", days)],
  1,
  function(x) all(is.finite(x))
)

for (day in days) {
  
  x <- dat[[paste0("LFC_", day)]][finite]
  
  if (length(x) < 3 || sd(x) == 0) {
    stop("有效LFC数量或变异不足，无法计算Z-score。")
  }
  
  dat[[paste0("Z_", day)]] <- NA_real_
  
  dat[finite, paste0("Z_", day)] <- as.numeric(
    scale(x)
  )
}

# -------------------- 6. 匹配10个重点基因 --------------------

# 忽略大小写，例如Irf1匹配IRF1
key <- toupper(trimws(dat$gene))

if (anyDuplicated(key)) {
  stop("忽略大小写后出现重复基因名，请核查注释。")
}

dat$target <- targets[match(key, targets)]
dat$is_target <- !is.na(dat$target)

status <- data.frame(
  target = targets,
  status = "Ready for plotting"
)

status$status[
  !targets %in% toupper(trimws(genes))
] <- "Not found in input"

status$status[
  targets %in% toupper(trimws(genes[!keep]))
] <- "Removed by low-count filter"

status$status[
  targets %in% dat$target[!finite & dat$is_target]
] <- "Non-finite model effect"

write.csv(
  status,
  file.path(out_dir, "target_status.csv"),
  row.names = FALSE
)

print(status)

if (any(status$status != "Ready for plotting")) {
  warning(
    "部分重点基因无法绘图，",
    "请查看target_status.csv；不会填入虚构数值。"
  )
}

# -------------------- 7. 三个时间点的共同/特异分类 --------------------

for (day in days) {
  
  f <- dat[[paste0("FDR_", day)]]
  z <- dat[[paste0("Z_", day)]]
  lfc <- dat[[paste0("LFC_", day)]]
  
  known <- is.finite(f) &
    is.finite(z) &
    is.finite(lfc)
  
  state <- rep("Not selected", nrow(dat))
  
  state[
    known &
      f < fdr_cutoff &
      z > enriched_z_cutoff &
      lfc > 0
  ] <- "Enriched"
  
  state[
    known &
      f < fdr_cutoff &
      z < depleted_z_cutoff &
      lfc < 0
  ] <- "Depleted"
  
  state[!known] <- "Unavailable"
  
  dat[[paste0("Call_", day)]] <- state
}

dat$three_timepoint_pattern <- apply(
  dat[paste0("Call_", days)],
  1,
  function(v) {
    
    if (any(v == "Unavailable")) {
      return("Unclassified: unavailable result")
    }
    
    up <- days[v == "Enriched"]
    down <- days[v == "Depleted"]
    
    if (length(up) > 0 && length(down) > 0) {
      return("Opposite directions across timepoints")
    }
    
    if (length(up) > 0) {
      return(
        paste0(
          "Enriched: ",
          paste(up, collapse = "+")
        )
      )
    }
    
    if (length(down) > 0) {
      return(
        paste0(
          "Depleted: ",
          paste(down, collapse = "+")
        )
      )
    }
    
    "Not selected at any timepoint"
  }
)

write.csv(
  dat,
  file.path(out_dir, "all_gene_results.csv"),
  row.names = FALSE,
  na = "NA"
)

write.csv(
  dat[dat$is_target, ],
  file.path(out_dir, "target_gene_results.csv"),
  row.names = FALSE,
  na = "NA"
)

write.csv(
  dat[
    c(
      "gene",
      paste0("Call_", days),
      "three_timepoint_pattern"
    )
  ],
  file.path(out_dir, "three_timepoint_overlap.csv"),
  row.names = FALSE
)

# -------------------- 8. 每个两两比较的配色分类 --------------------

# “X only/Y only”仅针对所在子图的两个时间点

make_pair <- function(xday, yday) {
  
  d <- data.frame(
    gene = dat$gene,
    target = dat$target,
    x = dat[[paste0("Z_", xday)]],
    y = dat[[paste0("Z_", yday)]],
    LFCx = dat[[paste0("LFC_", xday)]],
    LFCy = dat[[paste0("LFC_", yday)]],
    FDRx = dat[[paste0("FDR_", xday)]],
    FDRy = dat[[paste0("FDR_", yday)]]
  )
  
  d <- d[
    is.finite(d$x) &
      is.finite(d$y) &
      is.finite(d$LFCx) &
      is.finite(d$LFCy),
  ]
  
  known <- is.finite(d$FDRx) &
    is.finite(d$FDRy)
  
  px <- known &
    d$FDRx < fdr_cutoff &
    d$x > enriched_z_cutoff &
    d$LFCx > 0
  
  py <- known &
    d$FDRy < fdr_cutoff &
    d$y > enriched_z_cutoff &
    d$LFCy > 0
  
  nx <- known &
    d$FDRx < fdr_cutoff &
    d$x < depleted_z_cutoff &
    d$LFCx < 0
  
  ny <- known &
    d$FDRy < fdr_cutoff &
    d$y < depleted_z_cutoff &
    d$LFCy < 0
  
  d$category <- "Not selected"
  
  d$category[px & !py] <- "Enriched: X only"
  d$category[!px & py] <- "Enriched: Y only"
  d$category[px & py] <- "Enriched: both"
  
  d$category[nx & !ny] <- "Depleted: X only"
  d$category[!nx & ny] <- "Depleted: Y only"
  d$category[nx & ny] <- "Depleted: both"
  
  d$category[
    (px & ny) | (nx & py)
  ] <- "Opposite directions"
  
  # FDR缺失时，不判为另一组“无命中”
  d$category[!known] <- "FDR unavailable"
  
  d$comparison <- paste0(
    xday, " (X) / ", yday, " (Y)"
  )
  
  d
}

d <- rbind(
  make_pair("D2", "D3"),
  make_pair("D2", "D7"),
  make_pair("D3", "D7")
)

if (nrow(d) == 0) {
  stop("没有可绘制的有效效应。")
}

# -------------------- 9. 类似参考图的分类配色 --------------------

palette <- c(
  "Enriched: both"      = "#E95448",
  "Enriched: X only"    = "#E3AD3D",
  "Enriched: Y only"    = "#D793CB",
  "Depleted: both"      = "#1877CF",
  "Depleted: X only"    = "#73A6B3",
  "Depleted: Y only"    = "#AAA0D7",
  "Opposite directions" = "#71548B",
  "Not selected"       = "#C9C9C9",
  "FDR unavailable"    = "#999999"
)

d$category <- factor(
  d$category,
  levels = names(palette)
)

hit <- d[
  !d$category %in% c(
    "Not selected",
    "FDR unavailable"
  ),
]

lab <- d[!is.na(d$target), ]

lim <- max(abs(c(d$x, d$y))) * 1.18

# -------------------- 10. 绘制三个子图 --------------------

p <- ggplot(d, aes(x, y)) +
  
  geom_hline(
    yintercept = 0,
    linetype = 3,
    colour = "grey55"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 3,
    colour = "grey55"
  ) +
  
  geom_point(
    aes(colour = category),
    size = 0.6,
    alpha = 0.25,
    show.legend = TRUE
  ) +
  
  geom_point(
    data = hit,
    aes(colour = category),
    size = 1.6,
    alpha = 0.9
  ) +
  
  # 重点基因：保留分类颜色，外加黑圈
  geom_point(
    data = lab,
    shape = 21,
    fill = NA,
    colour = "black",
    size = 3,
    stroke = 0.7
  ) +
  
  # 无论是否显著，均标注可用的重点基因
  geom_text_repel(
    data = lab,
    aes(label = target),
    colour = "black",
    fontface = "bold",
    size = 3.1,
    max.overlaps = Inf,
    seed = 42,
    min.segment.length = 0,
    box.padding = 0.45,
    point.padding = 0.3,
    max.time = 3,
    max.iter = 20000
  ) +
  
  scale_colour_manual(
    values = palette,
    drop = FALSE,
    name = NULL
  ) +
  
  facet_wrap(
    ~ comparison,
    nrow = 1
  ) +
  
  coord_equal(
    xlim = c(-lim, lim),
    ylim = c(-lim, lim)
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  
  guides(
    colour = guide_legend(
      nrow = 3,
      override.aes = list(
        alpha = 1,
        size = 3
      )
    )
  ) +
  
  labs(
    title = paste0(
      "CRISPR SCREEN"
    ),
    x = paste0(
      "X timepoint vs untreated endpoint control: ",
      "LFC Z-score"
    ),
    y = paste0(
      "Y timepoint vs untreated endpoint control: ",
      "LFC Z-score"
    ),
    caption = paste(
      paste0(
        "Enriched: FDR < ", fdr_cutoff,
        ", Z > ", enriched_z_cutoff,
        ", LFC > 0. Depleted: FDR < ", fdr_cutoff,
        ", Z < ", depleted_z_cutoff,
        ", LFC < 0."
      ),
      sep = "\n"
    )
  )

if (interactive()) print(p)

ggsave(
  file.path(out_dir, "shared_hits_three_panels.png"),
  p,
  width = 18,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(out_dir, "shared_hits_three_panels.pdf"),
  p,
  width = 18,
  height = 7
)

write.csv(
  d,
  file.path(out_dir, "shared_hits_classification.csv"),
  row.names = FALSE,
  na = "NA"
)

# -------------------- 11. 重点基因三时间点热图 --------------------

# 热图使用原始模型log2FC，三个时间点使用相同色标
heat <- do.call(
  rbind,
  lapply(days, function(day) {
    
    x <- dat[dat$is_target, ]
    
    data.frame(
      gene = x$target,
      day = rep(day, nrow(x)),
      LFC = x[[paste0("LFC_", day)]],
      FDR = x[[paste0("FDR_", day)]]
    )
  })
)

heat$day <- factor(
  heat$day,
  levels = days
)

heat$gene <- factor(
  heat$gene,
  levels = rev(targets)
)

heat$mark <- ifelse(
  !is.na(heat$FDR) &
    heat$FDR < fdr_cutoff,
  "*",
  ""
)

if (any(is.finite(heat$LFC))) {
  
  heat_limit <- max(
    abs(heat$LFC[is.finite(heat$LFC)])
  )
  
  if (heat_limit == 0) {
    heat_limit <- 1
  }
  
  p2 <- ggplot(
    heat,
    aes(day, gene, fill = LFC)
  ) +
    
    geom_tile(
      colour = "white",
      linewidth = 0.7
    ) +
    
    geom_text(
      aes(label = mark),
      colour = "black",
      size = 6
    ) +
    
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-heat_limit, heat_limit),
      na.value = "grey85",
      name = "log2 FC"
    ) +
    
    scale_y_discrete(drop = FALSE) +
    
    theme_minimal(base_size = 12) +
    
    theme(
      panel.grid = element_blank()
    ) +
    
    labs(
      title = paste0(
        "Target Genes"
      ),
      x = NULL,
      y = NULL,
      caption = paste0(
        "* BH FDR < ", fdr_cutoff,
        " versus untreated endpoint control"
      )
    )
  
  if (interactive()) print(p2)
  
  ggsave(
    file.path(out_dir, "target_genes_heatmap.png"),
    p2,
    width = 6,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    file.path(out_dir, "target_genes_heatmap.pdf"),
    p2,
    width = 6,
    height = 6
  )
}


# -------------------- 12. 通路参数及富集 --------------------
organism <- "mouse"             # 本批输入为小鼠符号；人源文库改为"human"
gene_fdr <- 0.05
min_abs_lfc <- 0                # 新增火山图/ORA：FDR + LFC方向；不使用Z门槛
pathway_fdr <- 0.05
top_n <- 10                     # 每张通路图最多显示的显著通路数
scatter_z_enriched <- 1.5       # 保留原三联散点图的判定规则
scatter_z_depleted <- -1.5
run_reactome_ora <- TRUE
run_kegg_ora <- TRUE
common_essential_file <- NULL  # 可选：同物种共同必需基因Entrez ID，每行一个、无表头

#extension_dir <- file.path(out_dir, "volcano_reactome")
extension_dir <- file.path(out_dir, "volcano_kegg")
dir.create(extension_dir, recursive=TRUE, showWarnings=FALSE)

classify <- function(lfc, fdr, cut=0.05, min_lfc=0) {
  ok <- is.finite(lfc)&is.finite(fdr)
  ans <- rep("Not selected", length(lfc))
  ans[ok&fdr<cut&lfc>min_lfc] <- "Enriched"
  ans[ok&fdr<cut&lfc< -min_lfc] <- "Depleted"
  ans[!ok] <- "Unavailable"
  ans
}

comparison_tables <- setNames(vector("list", length(days)), days)

for (day in days) {
  folder <- file.path(extension_dir, paste0(day, "_vs_Control"))
  dir.create(folder, recursive=TRUE, showWarnings=FALSE)
  
  r <- contrast_results[[day]]
  r$gene <- rownames(r)
  r$target <- targets[match(toupper(trimws(r$gene)), targets)]
  r$category <- classify(
    r$log2FoldChange, r$padj, gene_fdr, min_abs_lfc
  )
  
  r <- r[order(r$padj, -abs(r$log2FoldChange), na.last=TRUE),]
  comparison_tables[[day]] <- r
  
  for (nm in c("All", "Significant", "Enriched", "Depleted")) {
    sub <- switch(
      nm,
      All=r,
      Significant=r[r$category %in% c("Enriched", "Depleted"),],
      Enriched=r[r$category=="Enriched",],
      Depleted=r[r$category=="Depleted",]
    )
    write.csv(
      sub, file.path(folder, paste0(nm, "_genes.csv")),
      row.names=FALSE, na="NA"
    )
  }
  
  target_table <- r[!is.na(r$target),]
  target_table$volcano_available <-
    is.finite(target_table$log2FoldChange)&is.finite(target_table$padj)
  
  write.csv(target_table, file.path(folder, "Target_genes.csv"),
            row.names=FALSE, na="NA")
  
  missing_targets <- setdiff(
    targets, target_table$target[target_table$volcano_available]
  )
  
  writeLines(missing_targets, file.path(folder, "Targets_not_plottable.txt"))
  
  if (length(missing_targets))
    message(day, ": absent/unavailable targets: ",
            paste(missing_targets, collapse=", "))
  
  v <- r[is.finite(r$log2FoldChange)&is.finite(r$padj),]
  v$neglogFDR <- -log10(pmax(v$padj, .Machine$double.xmin))
  v$category <- factor(
    v$category, levels=c("Enriched", "Depleted", "Not selected")
  )
  vlab <- v[!is.na(v$target),]
  
  p_volcano <- ggplot(v, aes(log2FoldChange, neglogFDR)) +
    geom_point(aes(colour=category),
               size=1, alpha=0.5, show.legend=TRUE) +
    geom_hline(
      yintercept=-log10(gene_fdr), linetype=3, colour="grey50"
    ) +
    geom_vline(
      xintercept=unique(c(-min_abs_lfc, min_abs_lfc)),
      linetype=3, colour="grey50"
    ) +
    geom_point(
      data=vlab, shape=21, fill=NA, colour="black",
      size=3, stroke=0.7, show.legend=FALSE
    ) +
    geom_text_repel(
      data=vlab, aes(label=target),
      colour="black", fontface="bold", size=3.3,
      max.overlaps=Inf, seed=42, min.segment.length=0,
      box.padding=0.5, max.time=3, show.legend=FALSE
    ) +
    scale_colour_manual(
      values=c(
        Enriched="#E95448",
        Depleted="#1877CF",
        "Not selected"="#BDBDBD"
      ),
      drop=FALSE, name=NULL
    ) +
    theme_classic(base_size=12) +
    guides(
      colour=guide_legend(override.aes=list(alpha=1, size=3))
    ) +
    labs(
      title=paste(day, "vs Control"),
      x="log2 fold change",
      y="-log10(BH FDR)",
    )
  
  ggsave(file.path(folder, "Volcano.png"),
         p_volcano, width=9, height=7, dpi=300)
  ggsave(file.path(folder, "Volcano.pdf"),
         p_volcano, width=9, height=7)
}


# -------------------- 13. 通路富集 --------------------
run_enrichment <- function() {
  org_pkg <- if (organism=="mouse") "org.Mm.eg.db" else "org.Hs.eg.db"
  kegg_organism <- if (organism=="mouse") "mmu" else "hsa"
  
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages(
      "BiocManager", repos="https://cloud.r-project.org", lib=user_lib
    )
  
  needed <- c("AnnotationDbi", "clusterProfiler", org_pkg)
  #needed <- c("AnnotationDbi", "ReactomePA", org_pkg)
  missing <- needed[
    !vapply(needed, requireNamespace, logical(1), quietly=TRUE)
  ]
  
  if (length(missing))
    BiocManager::install(
      missing, ask=FALSE, update=FALSE, lib=user_lib
    )
  
  missing <- needed[
    !vapply(needed, requireNamespace, logical(1), quietly=TRUE)
  ]
  
  if (length(missing))
    stop("通路依赖未成功安装：", paste(missing, collapse=", "))
  
  OrgDb <- getExportedValue(org_pkg, org_pkg)
  
  # 官方SYMBOL优先，唯一ALIAS回退；日期名和映射歧义均记录并排除。
  map_genes <- function(ids, db) {
    ids <- as.character(ids)
    
    date_like <- grepl(
      "^[0-9]{1,2}-(Mar|Sep)$|^[0-9]{4}-[0-9]{2}-[0-9]{2}",
      ids
    )
    
    m <- data.frame(
      gene=ids,
      ENTREZID=NA_character_,
      source=NA_character_,
      status=ifelse(date_like, "date_like_not_mapped", "unmapped")
    )
    
    for (kt in c("SYMBOL", "ALIAS")) {
      query <- intersect(
        m$gene[m$status=="unmapped"],
        AnnotationDbi::keys(db, keytype=kt)
      )
      if (!length(query)) next
      
      a <- suppressMessages(
        AnnotationDbi::select(
          db, keys=query, keytype=kt, columns="ENTREZID"
        )
      )
      
      a <- unique(data.frame(
        gene=as.character(a[[kt]]),
        ENTREZID=as.character(a$ENTREZID)
      ))
      
      a <- a[!is.na(a$ENTREZID)&nzchar(a$ENTREZID),]
      n <- table(a$gene)
      
      m$status[m$gene %in% names(n)[n>1]] <- "ambiguous_one_to_many"
      a <- a[a$gene %in% names(n)[n==1],]
      
      j <- match(a$gene, m$gene)
      m$ENTREZID[j] <- a$ENTREZID
      m$source[j] <- kt
      m$status[j] <- "mapped"
    }
    
    good <- which(m$status=="mapped")
    dup <- duplicated(m$ENTREZID[good])|
      duplicated(m$ENTREZID[good], fromLast=TRUE)
    
    m$status[good[dup]] <- "multiple_inputs_same_entrez_excluded"
    m
  }
  
  mapping <- map_genes(rownames(dds), OrgDb)
  
  write.csv(
    mapping, file.path(extension_dir, "gene_ID_mapping.csv"),
    row.names=FALSE
  )
  
  capture.output(
    AnnotationDbi::metadata(OrgDb),
    file=file.path(extension_dir, "annotation_metadata.txt")
  )
  
  message("唯一映射成功：",
          sum(mapping$status=="mapped"), " / ", nrow(mapping))
  
  essential <- character()
  
  if (!is.null(common_essential_file)) {
    essential <- trimws(
      readLines(path.expand(common_essential_file), warn=FALSE)
    )
    essential <- unique(
      essential[nzchar(essential)&!grepl("^#", essential)]
    )
    
    if (!all(grepl("^[0-9]+$", essential)))
      stop("必需基因文件应为Entrez ID，每行一个、无表头。")
  }
  
  run_reactome <- function(query, bg, folder, day, direction) {
    dir.create(folder, recursive=TRUE, showWarnings=FALSE)
    
    stale <- file.path(
      folder,
      c("ERROR.txt", "analysis_status.txt", "Reactome_result.rds",
        "effective_background_Entrez.txt")
    )
    invisible(file.remove(stale[file.exists(stale)]))
    
    bg <- unique(bg)
    query <- intersect(unique(query), bg)
    
    writeLines(query, file.path(folder, "input_Entrez.txt"))
    writeLines(bg, file.path(folder, "background_Entrez.txt"))
    
    plot_title <- paste(day, direction, "genes: Reactome")
    note <- NULL
    effective_n <- NA_integer_
    
    tab <- data.frame(
      ID=character(), Description=character(), pvalue=numeric(),
      p.adjust=numeric(), geneID=character(), Count=integer()
    )
    
    if (!length(query)) {
      note <- "No mapped input genes for this direction"
      
    } else if (!length(bg)) {
      note <- "No mapped background genes"
      
    } else {
      obj <- tryCatch(
        {
          
          ans <- clusterProfiler::enrichKEGG(
            gene=query,
            universe=bg,
            organism=kegg_organism,
            keyType="ncbi-geneid",
            pvalueCutoff=1,
            pAdjustMethod="BH",
            qvalueCutoff=1,
            minGSSize=10,
            maxGSSize=500,
            use_internal_data=FALSE
          )
          
          # 保留原来将结果中Entrez ID转换为基因名的功能
          if (!is.null(ans)) {
            ans <- clusterProfiler::setReadable(
              ans,
              OrgDb=OrgDb,
              keyType="ENTREZID"
            )
          }
          
          ans
          
        }, 
        # ReactomePA::enrichPathway(
        #   gene=query,
        #   universe=bg,
        #   organism=organism,
        #   pvalueCutoff=1,
        #   pAdjustMethod="BH",
        #   qvalueCutoff=1,
        #   minGSSize=10,
        #   maxGSSize=500,
        #   readable=TRUE
        # ),
        error=function(e)e
      )
      
      if (inherits(obj, "error")) {
        note <- paste("Analysis failed:", conditionMessage(obj))
        writeLines(note, file.path(folder, "ERROR.txt"))
        
      } else if (is.null(obj)) {
        note <- "No testable Reactome pathways after annotation"
        
      } else {
        saveRDS(obj, file.path(folder, "Reactome_result.rds"))
        
        if (isS4(obj) && "universe" %in% methods::slotNames(obj)) {
          effective <- as.character(methods::slot(obj, "universe"))
          effective_n <- length(effective)
          
          writeLines(
            effective,
            file.path(folder, "effective_background_Entrez.txt")
          )
        }
        
        candidate <- as.data.frame(obj)
        
        if (nrow(candidate)) {
          tab <- candidate
        } else {
          note <- "No testable Reactome pathways after annotation"
        }
      }
    }
    
    write.csv(
      tab, file.path(folder, "Reactome_all.csv"),
      row.names=FALSE, na="NA"
    )
    
    sig <- tab[
      is.finite(tab$p.adjust)&tab$p.adjust<pathway_fdr,
      ,drop=FALSE
    ]
    sig <- sig[order(sig$p.adjust),,drop=FALSE]
    
    write.csv(
      sig, file.path(folder, "Reactome_significant.csv"),
      row.names=FALSE, na="NA"
    )
    
    plot_subtitle <- paste0(
      "Mapped query: ", length(query),
      "; requested background: ", length(bg),
      "; BH FDR < ", pathway_fdr
    )
    
    if (!nrow(sig)) {
      if (is.null(note))
        note <- "No pathways pass the BH FDR threshold"
      
      p_pathway <- ggplot() +
        annotate(
          "text", x=0, y=0,
          label=paste(strwrap(note, 60), collapse="\n"), size=4
        ) +
        theme_void() +
        labs(title=plot_title, subtitle=plot_subtitle)
      
    } else {
      z <- head(sig, top_n)
      z$score <- -log10(pmax(z$p.adjust, .Machine$double.xmin))
      
      labels <- setNames(
        vapply(
          z$Description,
          function(s)paste(strwrap(s, 55), collapse="\n"),
          character(1)
        ),
        z$ID
      )
      
      z$ID <- factor(z$ID, levels=rev(z$ID))
      
      p_pathway <- ggplot(z, aes(score, ID)) +
        geom_col(
          fill=if(direction=="Enriched") "#DF827B" else "#126BB5",
          width=0.7
        ) +
        scale_y_discrete(labels=labels) +
        geom_vline(
          xintercept=-log10(pathway_fdr),
          linetype=3, colour="grey50"
        ) +
        theme_classic(base_size=11) +
        labs(
          x="-log10(pathway BH-adjusted P value)",
          y=NULL,
          title=plot_title,
          subtitle=plot_subtitle,
          caption=paste(
            "Over-representation among selected knockout targets;",
            "not pathway activation/inhibition."
          )
        )
    }
    
    height <- if(nrow(sig)) {
      max(4, 0.48*min(nrow(sig), top_n)+2)
    } else {
      4
    }
    
    ggsave(
      file.path(folder, "Reactome_barplot.png"),
      p_pathway, width=10, height=height, dpi=300
    )
    ggsave(
      file.path(folder, "Reactome_barplot.pdf"),
      p_pathway, width=10, height=height
    )
    
    writeLines(
      if(is.null(note)) "Completed" else note,
      file.path(folder, "analysis_status.txt")
    )
    
    message(
      day, " / ", direction, ": ",
      nrow(sig), " significant pathways; status: ",
      if (is.null(note)) "Completed" else note
    )
    
    data.frame(
      comparison=paste0(day, "_vs_Control"),
      direction=direction,
      mapped_input=length(query),
      mapped_background=length(bg),
      effective_background=effective_n,
      tested_pathways=nrow(tab),
      significant_pathways=nrow(sig),
      status=if(is.null(note)) "Completed" else note
    )
  }
  
  summary_rows <- list()
  
  for (day in days) {
    r <- comparison_tables[[day]]
    r$ENTREZID <- mapping$ENTREZID[match(r$gene, mapping$gene)]
    r$mapping_status <- mapping$status[match(r$gene, mapping$gene)]
    
    folder <- file.path(extension_dir, paste0(day, "_vs_Control"))
    
    write.csv(
      r, file.path(folder, "All_genes_with_mapping.csv"),
      row.names=FALSE, na="NA"
    )
    
    # 背景是该比较有可用FDR/LFC、唯一映射成功的基因，
    # 不是全基因组或仅显著基因。
    eligible <- is.finite(r$padj)&is.finite(r$log2FoldChange)&
      r$mapping_status=="mapped"
    
    bg <- unique(r$ENTREZID[eligible])
    
    for (direction in c("Enriched", "Depleted")) {
      query <- unique(r$ENTREZID[eligible&r$category==direction])
      b <- bg
      
      if (direction=="Depleted" && length(essential)) {
        query <- setdiff(query, essential)
        b <- setdiff(b, essential)
      }
      
      summary_rows[[paste(day, direction, sep="_")]] <- run_reactome(
        query, b,
        file.path(folder, paste0("Reactome_", direction)),
        day, direction
      )
    }
  }
  
  write.csv(
    do.call(rbind, summary_rows),
    file.path(extension_dir, "enrichment_run_summary.csv"),
    row.names=FALSE
  )
}

options(timeout = 600)
#run_reactome_ora
if (run_kegg_ora) {
  run_enrichment()
}

