library(tidyverse)
library(patchwork)
library(dplyr)
library(ggplot2)
library(MNet)
library(data.table)
library(do)

rm(list=ls())
gc()
getwd()
setwd('~/rawdata')
gene = "L2HGDH"
outdir = paste0("~/OV/",gene)

metab <- fread("~/rawdata/L-2h_rnaseq/L2HGDH-metabolism.csv")
##处理数据
metab <- as.data.frame(metab)
rownames(metab) <- metab$`Sample Name`
metab <- metab[,-1]
metab <- as.data.frame(t(metab))

metab <- metab %>% arrange(desc(Group))
group <- metab$Group
pair_id <- metab$Pair
metab <- metab[,-c(1:2)]
metab %>% mutate(across(where(is.character), as.numeric))  -> metab
metab <- as.data.frame(t(metab))

group <- ifelse(group == "T", "tumor", "normal")

# ##计算差异——不配对
# result_mlimma_all <- mlimma(metab,group)

##计算差异——配对
library(limma)
design <- model.matrix(~ 0+group + pair_id)
fit <- lmFit(metab,design)
fit <- eBayes(fit)
result_mlimma_all <- topTable(fit,
                           adjust = 'fdr',
                           coef = "groupnormal",
                           n = Inf)
result_mlimma_all$name <- rownames(result_mlimma_all)

write.table(result_mlimma_all,paste0(outdir,"/result_mlimma_all.txt"),quote=F,row.names=F,sep="\t")   

##绘制热图
result_mlimma_filter <- result_mlimma_all %>%
  dplyr::filter(abs(logFC) > 1) %>%
  dplyr::filter(`adj.P.Val` < 0.05)%>% 
  arrange(desc(logFC))%>%
  head(n=10)


if ("L-2-Hydroxyglutaric acid" %in% result_mlimma_filter$name) {
  dat_filter <- metab %>%
    tibble::rownames_to_column(var="label") %>%
    dplyr::filter(label %in% result_mlimma_filter$name) %>%
    tibble::column_to_rownames("label") 
}else{
  dat_filter <- metab %>%
    tibble::rownames_to_column(var="label") %>%
    dplyr::filter(label %in% c("L-2-Hydroxyglutaric acid",result_mlimma_filter$name)) %>%
    tibble::column_to_rownames("label") 
}

##不能出现0
dat_filter[apply(dat_filter, MARGIN = c(1, 2), FUN = function(x) x == 0)] <- 1
pdf(file = paste0(outdir,"/Metabolism-",gene,"-heatmap.pdf"),width = 20,height = 10)
p_heatmap <- MNet::pHeatmap(dat_filter,group,clustering_distance_cols ="manhattan",
                            clustering_method="ward.D",fontsize_row=12)
dev.off()



