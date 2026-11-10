# Visium 感染邻域与 MDK–NCL

本说明对应 [SCT_MDK_NCL.R](SCT_MDK_NCL.R)。默认分析 D3 注射侧。VG161 切片用于比较 Infected 和 Bystander 的邻域组成，Vehicle 和 VG161 均用于展示 MDK、NCL 的空间表达和检查原始病毒信号。

## 输入与运行

### 已跑过上一版

1. 把 `reuse_details_rds` 填为上一版的 `details.rds`。脚本复用该文件的 Visium 条码、阵列坐标、宿主质控、MDK/NCL 表达和 RCTD 比例，重新汇总两类中心的邻域。如果旁边有旧的 `summary.csv`，脚本保留其中适用的单细胞表达汇总。默认还从 `sc_file` 读取原始单细胞计数绘制最后一页，空间分析继续复用缓存。如果文件不可用，或将 `plot_sc_violin` 设为 `FALSE`，PDF 只输出空间分析部分。
2. 在 `mapping_files` 中按空间样本名填写 CMAP 的结果文件。D3 默认读取 `~/rawdata/SCT/analysis/virus_3D_1/virus_3D_anno.RData` 的 `sct`。用原 `SCT.R` **完成 CMAP 后最后保存的版本**。原脚本在更早的步骤也保存过同名文件。
3. 核对 `outdir`，在服务器进入脚本所在目录运行下面的命令。

```bash
Rscript SCT_MDK_NCL.R
```

### 从原始数据运行

把 `reuse_details_rds` 留空，仍然设置 `mapping_files`。再核对以下输入。

| 配置 | 内容 |
| --- | --- |
| `sc_file`、`sc_object` | 带原始 `RNA` counts 的 Seurat RData 或 RDS。D3 的 RData 对象名为 `sce`。Seurat v5 多个 counts 层按细胞条码合并，原层缺失的特征在诊断中保留为缺失。 |
| `sample_col`、`celltype_col`、`sc_groups` | 单细胞样本列、细胞类型列及样本到处理组的映射。默认是 `orig.ident`、`celltype` 和四个 D3 单细胞样本。 |
| `plot_sc_violin` | 默认 `TRUE`，在有可用原始单细胞计数时绘制上皮细胞 MDK 和 NCL 的逐细胞分布，放在同一份 PDF 的最后一页。 |
| `st_root`、`samples$data_dir` | Space Ranger 结果目录。每张切片需要 `filtered_feature_bc_matrix.h5`，以及 `spatial/` 中的 `tissue_positions.csv` 或旧版 `tissue_positions_list.csv`。 |
| `samples$sample_id`、`samples$treatment`、`virus_group` | 切片名、处理组及需要比较 Infected 与 Bystander 的组名。处理组默认是 VG161。不同切片各用自己的 CMAP 文件。 |
| `samples$composition_csv` | 留空运行 RCTD。已有比例时填 CSV，包含唯一的 `spot` 列和所有细胞类型的比例列。比例在 0–1，每行加和接近 1。导入前去除残差等非细胞类型列。Seurat 标签转移分数或 CMAP 的最高概率标签各有自己的用途，组成输入需要完整比例。 |
| `epithelial_types` | 上皮类别名称，默认是 `Epithelial`。用于非上皮分母、图上显示类别和可选的上皮单细胞汇总。若导入的比例表把上皮拆成多个类别，在此列出对应的全部名称。脚本在组成中找不到任何匹配名称时会提醒核对。 |
| `sc_host_prefix`、`st_host_prefix` | 人鼠混合参考中用于选取人源宿主基因的前缀。D3 默认使用无前缀基因名，需保留 `MDK`、`NCL`。 |

脚本加载 `Seurat`、`Matrix`、`dplyr`、`tidyr`、`ggplot2` 和 `patchwork`。读取 H5 还需要 `hdf5r`。重新拟合组成时需要 `spacexr`，支持 GitHub 和 Bioconductor 两种接口，安装命令在脚本开头。复用缓存或导入完整比例时可以跳过 RCTD 拟合。读取既有 CMAP 结果时，运行本脚本只需要相应的 R 输入文件。

### CMAP 文件格式与分组定义

- **Infected、Bystander** 沿用原 `SCT.R` 中单细胞的 `sce$infected`。原流程在 VG161 单细胞中用该列细分上皮类别，其他细胞类型保留原注释。原单细胞脚本在 VG161 组以骨架病毒特征信号大于 0 标记 Infected、信号为 0 标记 Bystander。这里按映射结果保存的原标签分析。
- **CMAP** 将单细胞转录组映射到空间位置。原 `SCT.R` 在每个 Visium spot 中取 `Probability` 最大的映射细胞，保存该代表细胞的类别到 `sct$celltype_new`。本脚本按切片和原始 spot 条码读入这列。两个状态的拼写和大小写原样保留。
- **分析单位**是 Visium spot。Infected 和 Bystander 指代表细胞映射为这两类的中心 spot。VG161 切片中通过质控且映射为这两类的 spot 全部纳入，各指标再按可评估性汇总。其他标签和缺少映射的 spot 通过质控后可以计入邻域。Vehicle 参与表达图和病毒质控。
- **支持格式**包括包含 `sct` 的 RData、直接保存 Seurat 对象的 RDS，以及含 `spot`、`celltype_new` 的 CSV。对象名和标签列分别由 `mapping_object`、`mapping_label_col` 指定。CSV 如带 `sample_id`，须与配置中的切片名一致。条码必须唯一且能在该切片找到。不同 Visium 切片可能复用条码，须从原流程确认文件和切片的对应关系。

如果希望用体积更小的 CSV，在原 `SCT.R` 的 CMAP 步骤完成后，可从当时的 `sct` 导出。

```r
mapping <- sct[[]][, "celltype_new", drop = FALSE]
mapping$spot <- rownames(mapping)
mapping$sample_id <- "VG161_3D_1"
write.csv(mapping[, c("spot", "sample_id", "celltype_new")], "VG161_3D_1_mapping.csv", row.names = FALSE)
```

## 邻域和组间比较是怎么计算的

1. **spot 质控和坐标。** `in_tissue = 1` 表示组织内 spot。UMI 是去重后的分子计数。脚本要求宿主 UMI 大于零且至少达到 `min_host_UMI`，默认 100。阵列坐标 `array_row`、`array_col` 用于建图，原图像像素坐标用于绘图，保留切片方向。复用缓存时沿用当时的质控和坐标，改变宿主 UMI 门槛需要从原始输入重跑。
2. **Center、Ring 1、Ring 2。** Center 是用于比较的一个 Infected 或 Bystander spot，单独汇总。Ring 1 是传统 Visium 六边形阵列上一步可到达的有效 spot，最多 6 个。相对坐标为 `(0, ±2)` 和 `(±1, ±1)`。Ring 2 是沿有效 Ring 1 边最短走两步的 spot，完整阵列内部最多 12 个。Ring 2 排除中心和所有 Ring 1 点，也去掉多条路径重复到达的点。例如从中心经一个一环点走到一个非一环点，该点属于 Ring 2。通过质控的 spot 才入图，组织边缘和缺失点会截断路径。没有有效边时不跨空洞补邻居。Ring 2 按图距离计算，局部图的 200 µm 比例尺按传统 Visium 的 100 µm spot 中心间距换算。
3. **RCTD 组成与环内均值。** RCTD 用有细胞类型注释的单细胞表达作为参考，估计每个 spot 在参考涵盖的细胞类型之间的相对比例。重新拟合时，两张切片共用同一份参考，至少需要 500 个共享宿主基因。`reference_min_cells`、`reference_max_cells` 默认每类至少 25 个、最多抽 1,500 个参考细胞。参考宿主深度至少 100 UMI，`seed` 默认 123，`rctd_cores` 默认 1。full 模式允许多个细胞类型同时存在。病毒特征和 MDK/NCL 排除在 RCTD 拟合基因之外。环内把组成可评估的邻点等权平均，每个中心在组内等权。缺失组成的邻点不计入组成均值，少于 `min_neighbors` 个可用邻点时，该中心的对应环组成记为缺失。默认需要 2 个邻点。
4. **同片分组比较。** 直接按映射标签比较两组，无需为 Infected 找对应的 Bystander，也无需用中心上皮比例设门槛。每张切片按环和指标分别汇总，组内各可评估中心等权。某个中心缺少某一环的可用数据时，其他环和可评估指标仍可使用。某组在一个指标上没有可评估中心时，保留另一组的估计，组间差记为缺失。
5. **两套分母。** 第一套以全部估计组成为分母，主图突出非上皮类别，汇总表也保留上皮。第二套以非上皮成分为分母。对每个中心和每一环，先计算环内各细胞类型的原比例均值，再把每种非上皮细胞的比例除以该环的非上皮比例总和。Center 直接用中心 spot 的比例。随后按组对各可评估中心等权平均。例如一环有 90% 上皮和 2% 髓系，髓系占全部组成的 2%，占非上皮组成的 20%。非上皮分母为 0 或缺失的中心仅在对应环的非上皮比例中记为缺失。分母大于 0 时，某个类别实际为 0 的估计值照常纳入。非上皮占比很小时，归一化后的比例应结合分母诊断一起看。
6. **邻域重叠。** 邻域可能重叠，按切片和环记录两组各自使用的唯一邻点数，以及两组共同使用的唯一邻点数。当前 D3 每个处理组一张切片，结果描述这张切片中两类映射中心周围的组成。两组在切片上分布的区域不同也可能带来组成差异。

## MDK–NCL 和病毒质控

- **表达量**采用原始宿主计数计算 `log1p(10,000 × 基因 UMI / 宿主 UMI)`。中心展示该 spot 的 MDK 或 NCL 表达，邻域表达取可用邻点的归一化值均值。每环至少有 `min_neighbors` 个通过空间质控的点才计算邻域表达和共检出。空间图统一色标，计数为 0 的 spot 保留在色标内，灰色表示不可评估。
- **共检出**以某个中心和指定环为单位。如果中心 MDK 有 UMI，计算该环中 NCL 有 UMI 的邻点比例。如果中心 MDK 为 0，该比例记为 0。反向指标交换中心和邻点的两个基因。两种方向分别报告，表示空间 RNA 共检出。NCL RNA 不提供细胞表面蛋白定位或受体结合证据。原 CMAP 分组沿用其上游基因选择，解读 MDK–NCL 结果时也要考虑映射与表达对比可能共享输入特征。
- **病毒质控**按 `viral_pattern` 查找 VG161 的 UL、ICP、US、IRL 骨架特征。`viral_genes_use` 可以指定精确基因名。`all_viral_pattern` 用于从宿主计数中排除 VG161 特征。缺少完整病毒特征时，`n_assessable` 反映可评估数量。空间病毒 UMI 仅用于质控。
- **单细胞表达**按样本和细胞大类汇总 MDK、NCL 的平均归一化表达与原始 UMI 检出比例。原始 `sce` 有 `infected` 列时，另汇总已有标签的 VG161 上皮细胞。单细胞分布图以 `epithelial_types` 选取上皮细胞，以每个细胞原始 RNA 计数计算的 `log1p(10,000 × 基因 UMI / 宿主 UMI)` 为纵轴。如果 VG161 上皮同时有明确的 Infected 和 Bystander 标签，则按这两个标签出图。其余情况下按原处理组标签 Vehicle 和 VG161 展示上皮细胞。各样本分别分面，同一基因跨样本使用同一纵轴尺度，表达为零的可评估细胞也进入分布。

## 输出与验证

- `summary.csv` 中，先看 `Mapping QC` 的各标签总数 `n` 和通过空间质控数 `value`。`Composition`、`Non-epithelial composition` 和 `MDK_NCL` 按切片、环和指标给出各组可评估中心数 `n_Infected`、`n_Bystander`，各组均值 `Infected`、`Bystander`，以及组间差 `difference`（Infected 减 Bystander）。`Non-epithelial denominator` 给出各组非上皮分母的均值和有效数。分母的有效数减去同组同环非上皮比例的有效数，即该环分母为 0 的中心数。组成在表中用 0–1 表示，图中乘以 100，组间组成差用百分点。
- `plots.pdf` 前两页分别使用全组成和非上皮分母。每组每环每类细胞只有一个均值，使用并列条形图，组间差用热图。随后是 VG161 空间定位、全组成分母下的非上皮细胞分布、Vehicle 和 VG161 的 MDK/NCL 空间表达，以及 MDK–NCL 表达和共检出。空间定位页为每个有可用一环组成的组选一个局部示例，选择最接近该组空间坐标中位数的中心。局部堆叠条保留完整的一环细胞组成。启用单细胞分布图且数据可用时，最后一页展示上皮细胞的 MDK 和 NCL 表达。
- **分布图的图型**。空间逐中心和单细胞逐细胞的分布均用小提琴图展示有不同取值的多观测组。小提琴限制在观测范围，横线标记中位数。一个观测值用条形，零值用数字标在基线。同组有多个观测且取值完全相同时，箱线图收缩为一条横线，标出共同取值。
- `details.rds` 保存每张切片的坐标和计数指标、组成、映射标签、邻接条码及 `centers` 表。中心表用 `spot`、`state` 记录中心条码和映射标签。`info` 记录映射来源、参数和运行环境。
- 已用合成数据验证原始 H5 与 RData、旧缓存与 CSV/RDS 两类入口。非上皮分母通过独立计算核对了中心与两环的组均值、差值和有效数，也核对了零分母和真实零比例的处理。已渲染检查并列均值条形图及空间和单细胞小提琴图，并测试了单细胞标签缺失、单细胞文件缺失的情况。真实 D3 的 CMAP 文件未在本机获得，真实标签下的两组中心数和邻域结果需在服务器运行确认。

## 参考文献

1. 项目原流程的 CMAP 标签生成及写回。 [原 SCT.R](https://github.com/LinghaoXia/OV/blob/scRNA_virus/SCT.R) · [单细胞注释脚本](https://github.com/LinghaoXia/OV/blob/scRNA_virus/scRNA_all.R)
2. Ke 等，*High-resolution mapping of single cells in spatial context*，Nature Communications，2025。 [论文](https://doi.org/10.1038/s41467-025-61667-4) · [CMAP 官方代码](https://github.com/SuoLab-GZLab/CMAP)
3. Cable 等，*Robust decomposition of cell type mixtures in spatial transcriptomics*，Nature Biotechnology，2022。 [论文](https://doi.org/10.1038/s41587-021-00830-w) · [GitHub 版 RCTD](https://github.com/dmcable/spacexr) · [Bioconductor 版 RCTD](https://bioconductor.org/packages/spacexr/)
4. 10x Genomics 的阵列坐标、文件结构与 spot 间距。 [Space Ranger 空间输出](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/spatial-outputs) · [Visium 术语](https://www.10xgenomics.com/support/software/space-ranger/latest/getting-started/space-ranger-glossary)
5. Pei、Min 等，*Spatial mapping of transcriptomic plasticity in metastatic pancreatic cancer*，Nature，2025。参考图 5 的定位、局部图、组成和分布。 [图 5](https://www.nature.com/articles/s41586-025-08927-x/figures/5)
6. Mo、Liu 等，Nature，2024。参考 Visium 细胞比例和表达分开呈现的图 3 和扩展图 9。 [论文](https://doi.org/10.1038/s41586-024-08087-4)
7. ggplot2 图层定义与参数。 [小提琴图](https://ggplot2.tidyverse.org/reference/geom_violin.html) · [箱线图](https://ggplot2.tidyverse.org/reference/geom_boxplot.html)
