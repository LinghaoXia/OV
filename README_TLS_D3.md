# D3注射侧：TLS分析步骤1—6

脚本：`SCT_TLS_D3.R`。默认处理Vehicle_3D_1、VG161_3D_1各一张切片，始终围绕TLS相关程序和TLS候选区域。独立样本数为每组1，因此本轮只做描述性比较，不计算治疗效应P值。

## 运行前修改

脚本最前面清理环境、关闭图形设备、回收内存，随后重新加载输入。路径及参数集中在开头，默认如下：

```r
data_root <- "~/rawdata/SCT/virus_3D"
outdir <- "~/OV/SCT/TLS_D3"
# samples$data_dir分别为上述目录下Vehicle_3D_1、VG161_3D_1。
```

默认直接读各自的`filtered_feature_bc_matrix.h5`，需要配套`spatial/`目录中的位置表、比例尺和PNG。也支持把样本表的`input_type`改为`object`并分别填写`object_file`；原始`data_dir`仍用于核对切片身份和空间坐标。保留输入中的原始全基因counts，不用SCT/predictions替代。

依赖Seurat、Matrix、dplyr、ggplot2、patchwork、AUCell；png、jsonlite、igraph、sp、BiocParallel等由对应依赖提供。若服务器没有AUCell，先单独安装一次：

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("AUCell", ask = FALSE, update = FALSE)
source("SCT_TLS_D3.R", encoding = "UTF-8")
```

## 做了哪些分析

1. **输入与QC**：更新旧Seurat对象；检查原始counts、原始位置表、图像是否对齐。基础QC默认与上一版一致。两遍顺序读取，避免同时保留两份大型Seurat对象。
2. **统一AUCell**：使用两片共同保留的宿主基因，固定基因顺序、随机种子和前5%的排名截断。检出基因不足统一截断值的spot标不可评价，不把零表达随机排名当信号。主评分为作者TLS22清单，辅助为Imprint29及B/T、组织化、GC相关面板。
3. **表达候选ROI**：两片共用199个按平均表达/检出率匹配的随机背景基因集。默认要求TLS22检出≥3基因、相对背景尾部比例≤0.01，并有≥3个直接相邻spot。用Visium阵列六邻域连通，保留孔洞，不做k近邻跨接、补洞或凸包填充。结果仅为待复核表达候选。
4. **H&E独立候选并集**：可先在R里圈现有H&E，或导入已有多边形。H&E候选即使转录支持弱也保留。表达/H&E区域按共享spot归并，记录`expression/HE_only/both`来源，全部标为`pending_review`。

筛查阈值是本项目的探索性参数，不是论文验证的TLS诊断标准。匹配背景尾部比例不是治疗P值、不是TLS概率；未做多重比较校正，也不代表ROI层面的FDR。随机背景不完全保留基因间相关性；表达筛查也允许单一淋巴细胞聚集入选，因此仍须在后续审核B/T结构和形态。`results.rds$background_matching`保存每个目标基因与抽中背景的表达/检出率匹配情况。`results.rds$sensitivity`保存尾部阈值0.005/0.01/0.025与最少相邻spot数2/3/4的检查结果，不据此挑选最符合预期的一套参数。

汇总同时报告通过QC后可计算AUCell的比例，PDF显示检出基因分布与统一K值。TLS可能只占小区域，全片中位数即使有局部差异也可能均为0，因此汇总保留均值、中位数、非零评分spot比例，最后图展示均值和表达候选数。若两片被低复杂度门控排除的比例不同，不能把剩余spot的评分差直接写成治疗效应。

主面板TLS22要求完整保留，缺失则停止；辅助Imprint29按所有样本共有的同一批基因评分，要求覆盖≥80%，覆盖明细随结果保存。HNSCC原文写`PTLP`，尚无核实过的`PLTP`更正，因此不静默替换。若矩阵无PTLP，明确记录实际28/29覆盖，不能称完整29基因复现。AUCell只是复用作者基因清单的本项目评分，不是Science的iSTAR模型。

## 如何完成H&E圈定

默认`draw_HE_ROI <- FALSE`，先生成匿名H&E页、统一评分图和表达候选轮廓。此时`HE_review_done=FALSE`，表达候选为零不能写成“本片无TLS”。匿名H&E页在表达图前，建议先独立查看；不将低分候选直接排除。

若在交互式R/RStudio中希望先圈图，把`draw_HE_ROI <- TRUE`后从头运行：先看只显示S01/S02的H&E，输入`y`后依次左键圈顶点，右键或Esc结束；输入`n`结束该切片。图面隐藏处理组，表达图尚未展示；控制台和脚本包含样本信息，因此这不等同正式的分组盲法评价。多边形坐标直接记录为原始fullres像素。

也可设置`HE_roi_file`读CSV，每个多边形多个顶点行：

```text
sample_id,roi_id,x,y
Vehicle_3D_1,H1,1000,2000
Vehicle_3D_1,H1,1600,2000
Vehicle_3D_1,H1,1600,2600
```

上述数字仅示意，不能用于真实样本；坐标必须是fullres、原点左上，不是PDF坐标。绘制过的多边形保存在结果中，下次可直接复用，不必另存CSV：

```r
HE_roi_file <- "~/OV/SCT/TLS_D3/results.rds"
draw_HE_ROI <- FALSE
```

PDF默认采用hires PNG，缺失时退回lowres；仅PDF预览下采样至最长边约1200像素，原图保持不变。每个模块内部跨样本共用色阶，不用不同模块的颜色比较细胞比例。现有PNG分辨率不足时，成熟度保留未知，不能凭预览图给出GC诊断。复用旧RDS中的人工ROI时核对坐标表与H&E文件指纹，避免把旧轮廓套到新切片。

## 前半段步骤1—3的三个输出

- `summary.csv`：每样本一行，包括QC、可评分数量、各面板覆盖、评分摘要、表达/H&E/合并候选数、H&E是否已审核。
- `plots.pdf`：匿名原始H&E、两片统一色阶的QC/面板图、候选ROI轮廓、最后的分组描述图。橙色为表达网格轮廓，蓝色为H&E多边形；均不是已确认TLS边界。
- `results.rds`：保存小型spot表、ROI成员和轮廓、人工多边形、少量marker counts、缺失基因审计、敏感性检查、参数、时间、包版本及来源。不保存Seurat对象、全基因counts、排名矩阵或图像。

```r
tls <- readRDS("~/OV/SCT/TLS_D3/results.rds")
tls$summary       # 样本汇总
tls$rois          # 待复核候选，尚不是最终TLS计数
tls$coverage      # 完整/部分面板以及缺失基因
tls$sensitivity   # 候选对预设阈值的稳定性
tls$info          # 参数、时间、解释边界、sessionInfo
```

前半段不输出TLS物理面积、密度或成熟度标签；网格轮廓只定位需要复核的捕获位置。后半段在证据及人工资料明确后才汇总对应指标。

## 已接在脚本后面的步骤4—6

从`##### 第4—6步：续接现有结果，不再次清理运行环境 #####`开始执行即可；不必重新加载大Seurat对象或重算AUCell。清理环境仍只在整个脚本最前面进行。

后半段默认沿用当前`outdir`，读取其中`results.rds`。若在新的R会话仅运行后半段，默认目录仍为`~/OV/SCT/TLS_D3`；需要时在后半段输入块改`review_dir`。不需要AUCell即可运行后半段。

### 第4步：局部证据及审核表

每个候选保留原始成员spot，并按原始六邻域外扩两圈。局部背景排除其他候选及QC不合格spot，不跨缺失位置或通过QC失败点连接。输出核心/邻域的marker原始UMI、检出数、QC后宿主归一化表达，以及B/T、组织化、FDC相关及GC相关证据。没有该基因、没有捕获spot或RNA质量不可评价，与实际测得零表达分别记录。

`review.pdf`包括候选局部H&E、marker空间图、原始signature构成及比较结果。原始UMI占比不是精确AUCell贡献率。图像文件仍在服务器的原始位置时可直接显示；搬家后可填写`review_image_files`命名向量。代码核对原图MD5及尺寸，不会拿不同图像按旧比例尺静默拉伸；图像不可用时保留坐标图并明确标记不能评价形态。

`review.csv`是唯一新增的人工回填表，包含两类记录：

| 记录/字段 | 填写方式 |
|---|---|
| `record_type=sample` | 每张切片一行；即使没有候选也保留 |
| `whole_section_HE_reviewed` | 整张现有H&E确已完成检查后才填TRUE；未看完或分辨率不足保留FALSE/未知 |
| `record_type=roi` | 每个已有候选一行，不自行添加或更改sample_id/roi_id/review_key |
| `HE_reviewed` | 该区域是否已查看原H&E，TRUE/FALSE |
| `decision` | pending、uncertain、rejected、TLS_supported |
| `morphology` | not_assessed、uncertain、no_aggregate、aggregate、follicle_like、GC_like |
| `BT_support` | 综合当前数据是否有B/T组成与结构支持，未知NA；不是仅有T信号就填TRUE |
| `n_TLS_structures` | 支持TLS后人工数独立结构；一个候选可以包含多个结构，未数清保留NA |
| `structures_nonoverlapping` | 样本行：多个支持ROI涉及的结构已经去重时填TRUE |
| `maturity` | 默认not_assessed；也可undetermined、early_like、follicle_like、GC_like。后三者须逐结构ROI且有相应依据 |
| `GC_structure_support` | GC样类别需明确结构支持，不能只根据BCL6/MKI67高表达填写 |
| `notes` | 记录依据与不确定性；B标记未检出/不可评价却支持TLS时必须解释额外证据 |

自动证据列`RNA_B_support`等只说明原始数据，不会替你把候选升级为TLS。审核结论与成熟相关标签都不自动填写。如果看图发现原候选之外的新结构，应回到第3步录入H&E多边形，再按新候选版本建立审核表。

### 第5步：审核后的结构数量、覆盖及成熟相关指标

整片H&E及全部候选没有解决时，全片TLS结构数保持NA；不是把未审核当0。先明确“支持的候选ROI数”，再根据人工填写的独立结构数量计算“TLS结构数”，两者不同。只有审核完成且未发现支持结构时才可报告0。

`supported_ROI_capture_fraction_QC`表示支持ROI覆盖的捕获spot比例，不是TLS实际面积。物理指标只在有真实测量时计算：样本行`tissue_area_mm2`为可评估组织面积；ROI行`TLS_area_mm2`为该ROI内已支持TLS结构的实际总面积，不能填候选网格面积。多个ROI的面积已去重后样本行`areas_nonoverlapping=TRUE`。未提供这些资料时，密度和面积比例保持NA，不影响其他步骤。

成熟相关类别要求逐结构ROI；同一个ROI含多个结构时先拆分，不能直接给整个区域套一个成熟标签。这里的`*_like`是当前证据支持的描述，不等同经过独立验证的病理成熟诊断。

### 第6步：两个样本的描述性比较

保留未调整评分，同时在两片共同的UMI/检出基因数分层内，用相同权重比较TLS22及Imprint评分。记录保留的spot比例；这只是技术敏感性检查，不能补回约40%未评分spot，也不是治疗效应的因果校正。每组只有一片，不计算spot或ROI充当动物重复的组间P值。

### 后半段只新增两个文件

- `review.csv`：首次生成后填写人工列；重新运行后半段会读取并保留合法审核内容。
- `review.pdf`：局部证据图及描述性比较。

原`plots.pdf`保留。`summary.csv`补充审核、结构数量及分层评分；详细证据与运行信息写入原`results.rds$followup`，不另存大型对象或更多中间文件。

```r
tls <- readRDS("~/OV/SCT/TLS_D3/results.rds")
tls$followup$roi_evidence
tls$followup$gene_evidence
tls$followup$structure_summary
tls$followup$depth_summary
tls$followup$descriptive_comparison
```

同版本可重复运行，审核内容按sample_id、roi_id及候选指纹精确匹配。若改变输入、轮廓或筛查版本，旧表不会自动套用；请保留旧表，将`review_file`指向一个新的明确文件名后再建立审核表。

## 验证范围

已用两份合成标准Visium数据验证步骤1—3的object/Assay5与默认10x HDF5/VisiumV2读取、统一AUCell、多块评分、连续及孤立信号、全零信号、空spot和低复杂度、组织孔洞、H&E-only候选、缺失必要基因报错、原始输入不变和前半段三个输出，并检查生成PDF。

后半段测试覆盖pending不记零、人工支持门控、一个ROI含多个结构、成熟度不自动赋值、缺基因/无捕获/全QC失败、审核后零结构、旧审核拒绝、幂等复跑及只新增两个文件；已在用户提供的真实D3结果副本上执行后半段并核对数值，未修改提供的原始文件。原始H&E仍在服务器的路径，局部原图应在服务器运行后查看。测试验证实现，不验证探索性阈值的生物学诊断性能。交互式鼠标圈图未通过自动化模拟实际鼠标操作。

## 参考

- [Science 2026 TLS图谱](https://doi.org/10.1126/science.adz2742)；[作者TLS22清单](https://raw.githubusercontent.com/Coolgenome/Pan-Cancer_Spatial_Atlas_TLS/main/data/TLS_signature/markers.tsv)。
- [Nature Communications 2025，HNSCC TLS-imprint/AUCell应用](https://doi.org/10.1038/s41467-025-59341-w)。
- [AUCell官方教程](https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html)。
