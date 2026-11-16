# 肿瘤细胞感染相关GSEA与MDK排序贡献分析

对应脚本：`scRNA_GSEA.R`。这是已有单细胞对象之后的第一步功能分析，默认比较D3肿瘤细胞中的 **Infected − Bystander**，同时记录MDK在基因排序、基因集注释及leading edge中的位置。

脚本独立运行，开头会清理当前R全局环境并关闭图形设备。它沿用对象已有的细胞类型和感染标签，不重新判定感染，不运行CellChat，不修改输入文件，不另存大型Seurat对象。`Epithelial`是否确为恶性肿瘤细胞、感染标签是否可靠，须由前置注释和原始实验记录支持。

本分析能说明感染状态相关的肿瘤转录程序及MDK这一排序项的关系，不能证明MDK引起这些变化，也不能证明MDK–NCL结合、T细胞耗竭或T细胞功能恢复。

## 1. 六个步骤及运行顺序

|步骤|用途与实现|需要先确认的事项|
|---|---|---|
|1 读取与核对输入|读取`.RData/.Rda`指定对象，或读取`.rds`；检查Seurat、RNA assay及三个metadata字段，按样本、肿瘤标签和两种感染状态选细胞|路径、对象名、真实样本表、恶性细胞身份和已有感染标签|
|2 提取原始宿主计数|读取RNA原始counts；兼容传统Assay和含`counts`/`counts.*`的Assay5；保留稀疏结构，在共同背景、pseudobulk及DE之前按正则和精确ID排除病毒特征并保存审计|宿主物种、feature前缀、病毒及外源载荷feature命名；不能用data、SCT或integrated替代|
|3 汇总配对pseudobulk|按`replicate_id × state`相加原始counts；两状态均达细胞数要求且有宿主计数才纳入整对|`replicate_id`是真实生物重复，技术拆库应归同一重复；同鼠两侧和不同时点不混成独立重复|
|4 配对edgeR及完整排序|表达过滤、TMM、稳健离散度/QL拟合；模型`~ block + state`，检验`stateInfected`；产生所有合格宿主基因的有符号排序|配对设计满秩、有剩余自由度；独立样本数不足不能用大量细胞弥补|
|5 完整GSEA及MDK记录|默认直接下载并校验Broad官方固定版本GMT；对完整Hallmark和GO:BP中符合可测大小/权重要求的集合运行双向fgsea，记录MDK成员/leading edge/命中权重；可选删除MDK排序项重跑|服务器能访问Broad官方数据地址；基因符号与基因集一致；展示清单不缩小检验家族|
|6 绘图与精简保存|保存样本审计、基因和通路结果、MDK表达/检出率、富集曲线及敏感性图；共四个输出文件|结果方向、NA和未测试状态均需结合审计表解释；没有重新运行真实对象前不报告真实结论|

建议按脚本中文分节顺序逐段运行，先读样本计数和MDK状态，再进行耗时的GSEA与作图。脚本也可整段运行。

## 2. 输入文件、字段及独立重复

### 2.1 Seurat对象

默认输入为`~/rawdata/scRNA_virus/virus_3D/virus_3D_anno.RData`，对象名为`sce`。`~`由当前R运行账户解释；在本地Windows和服务器上可能对应不同目录。换用`.rds`时仍填写`input_file`，`input_object`只在`.RData/.Rda`读取时使用。

输入必须包含：

- `RNA` assay的原始、非负、有限整数counts。不是归一化表达、整合值或SCT残差。
- 唯一且能与counts列名对应的细胞条码。
- 人宿主gene symbol；如有固定前缀，可配置`host_prefix`去掉。脚本不自动做Ensembl到symbol映射、不把鼠基因转大写、不混合人鼠同源基因。
- 已完成必要的细胞质量控制、物种分离、双细胞检查及细胞类型/感染状态注释。

metadata默认列如下；名称不同只修改开头对应配置。

|配置项|默认字段|含义与处理|
|---|---|---|
|`sample_col`|`orig.ident`|原始样本/文库标识；值必须与`sample_table$sample_id`逐字匹配|
|`celltype_col`|`celltype`|细胞类型；默认仅选`Epithelial`，不把该标签自动验证为恶性|
|`state_col`|`infected`|已有感染标签；仅纳入`Bystander`及`Infected`，其他标签和缺失值不纳入|

`Bystander`是当前对象中的操作标签，不保证细胞从未接触病毒；`Infected`也不自动等于生产性感染。脚本不复核或改变这些标签。这里的比较是同一肿瘤样本内两种状态的关联，不能直接解释成随机分配病毒处理的因果效应。

### 2.2 样本表与D3/D5运行

开头`sample_table`至少包含两列：

|字段|必须如何填写|
|---|---|
|`sample_id`|metadata中实际存在的样本/文库名，每行一个且不得重复|
|`replicate_id`|该文库所属的真实独立生物样本；同一肿瘤的技术拆库填写相同ID，汇总时合并|

默认暂列`VG161_3D_1`和`VG161_3D_2`，两个`replicate_id`也暂同文库名。**这不是对真实动物对应关系的确认，运行前必须对照实际记录。** `replicate_unit`同时改成清楚的单位说明，保存在运行记录中。

一个完整配对是同一个`replicate_id`同时具有合格的Bystander和Infected pseudobulk。两者并非同一细胞的治疗前后配对；它们是同一生物样本内两种细胞状态的比较。细胞、技术孔、文库、切片和ROI均不能替代独立动物/生物样本。

默认仅运行D3。D5应另改输入、`analysis_label`和样本表，并按L/R分别运行，例如仅纳入`VG161_5D_L1`、`VG161_5D_L2`的一次分析。先核对L/R是否分别为注射侧与对侧，不从字母推断。不要把同鼠左右瘤或不同日期当作互相独立的配对一起放入当前模型；当前脚本没有嵌套动物、多侧/多时间交互模型。

`sample_id`在对象中不存在会停止，不会猜条码前缀代替。样本真实存在但所选肿瘤状态不足时，后续按整对规则排除并记录。

### 2.3 原始层取交集与缺失

Assay5只提取含纳入细胞的原始counts层。层间条码重复、选定细胞没有counts、计数非整数或负值时停止。

不同层的feature集合可能不相同。脚本使用这些层**共有的宿主feature**进入DE，不使用`JoinLayers`把未测feature补为零。`layer_audit`记录每层宿主feature数、共同feature数及MDK是否存在。只有某层实际包含该feature的零计数才能解释为该层中的观测零。

细胞是否具有有效宿主计数、MDK细胞表达的CP10K分母和DE计数背景，均以这一共同宿主feature集合为准，使跨层分母一致。若交集明显缩小，应先追溯上游计数保存过程，而不是仅通过放宽阈值继续。

MDK缺失于任何纳入层时，它不进入共同feature背景；应结合`target_status`与`layer_audit`判断。用到缺失MDK层的样本/状态组，其MDK累计计数、检出细胞数、表达和检出率全部为NA，不能保留部分层总和冒充完整测量，也不能补零再画下降趋势。其他完整测到MDK的组可以保留描述值，但不能因此把MDK强行放回共同DE背景。若实际RNA行名主要是Ensembl ID，脚本会停止；应在前置步骤用固定注释完成确定映射，保留映射来源和重复处理规则。

### 2.4 病毒基因排除与宿主反应基因保留

病毒feature在原始counts层按名字识别，先从宿主feature集合排除，再形成共同宿主背景、计数分母、pseudobulk、DE及GSEA排序。不是仅在作图时隐藏，也不是仅删除不能匹配MSigDB的病毒名字。

默认`viral_pattern = "^(VG161|HSV[-_]?1)([-_.]|$)"`，忽略大小写匹配以`VG161`或`HSV1/HSV-1/HSV_1`开头，后接分隔符或名称结束的feature。`viral_feature_ids`为额外的精确ID列表，默认空向量。识别会同时检查原始feature名和去除已配置`host_prefix`后的名称；精确列表按名称大小写匹配。正则和精确列表按“满足任一条件即排除”补充使用。

- 若原始矩阵使用无统一前缀的病毒ID，例如单独的`UL...`、`US...`或自定义载荷名，默认规则不能保证覆盖；应对照参考注释补入`viral_feature_ids`或经核对的正则。不要凭猜测扩大到所有短字母开头的宿主基因。
- 若原始病毒feature另带前缀，脚本仅会尝试去除已配置的`host_prefix`再检查；其他未知前缀应按真实完整名称填写精确列表或修改正则。脚本不从gene symbol推断序列物种，也不重新比对。
- **IFIT、ISG15、MX、OAS等宿主干扰素/抗病毒反应基因应保留。** 它们是感染相关宿主程序的测量，不是病毒编码基因；“与病毒反应有关”不构成删除依据。
- 外源载荷若与宿主同名基因在上游计数时已合为同一行，仅靠名字不能区分来源。需要先追溯参考及计数规则，不能把混合行自动称为纯宿主表达。
- 查看`results.rds$viral_exclusion_audit`的实际排除清单及`layer_audit$n_viral_excluded`。零个排除项可能表示输入已分离，也可能表示命名规则不匹配，须结合原始feature核对。最终rank断言不含**规则已识别**的病毒feature；这个断言不能保证识别未知命名的病毒基因。

## 3. 开头需要修改的配置

优先检查以下几项，再运行其他代码：

|配置|默认值/用途|何时修改|
|---|---|---|
|`input_file`、`input_object`|输入文件与RData对象名|换文件或对象名；不把输出文件填作输入|
|`analysis_label`、`outdir`|默认D3；默认输出`~/OV/scRNA_virus/GSEA_D3`|D5_L、D5_R另设稳定目录，避免覆盖其他比较|
|三个metadata字段、`tumor_labels`|默认`orig.ident/celltype/infected`与`Epithelial`|实际注释列或标签不同|
|`sample_table`、`replicate_unit`|纳入文库及真实阻断单位|每次均与实验记录核对|
|`host_prefix`|默认空；若人特征为`human_MDK`等可填`human_`|输入不是已分离的人宿主无前缀数据|
|`viral_pattern`、`viral_feature_ids`|默认`^(VG161\|HSV[-_]?1)([-_.]\|$)`及空精确ID向量|按原始骨架/载荷feature名核对；无前缀的病毒ID必须补精确列表或规则；保留宿主ISGs|
|`target_gene`|默认`MDK`|确需追踪另一人gene symbol时修改|
|`gene_set_source_mode`|默认`"Broad"`；另允许`"GMT"`和`"msigdbr"`|默认直接读Broad官方固定版本；仅在明确改用离线或msigdbr入口时切换|
|`gene_set_dir`|默认`file.path(tempdir(), "gsea_broad_2026_1")`，保存Broad临时缓存|临时盘不可用时可改成可写目录；临时文件是否长期保留由运行环境决定|
|`hallmark_gmt/gobp_gmt/gmt_version`|两路径默认空；`Broad`模式在基因集分节赋实际缓存路径及固定版本|只有`GMT`模式需手填两文件路径与真实版本；路径留空本身不再切换来源|
|`download_timeout/download_retries`|默认180秒/3次|Broad下载超时与最多尝试次数；原R超时更长时保留更长值，下载成功会立即继续|
|`focus_pathways`|预先指定的免疫/组织展示清单|只改变展示，不改变完整检验；名称须精确匹配当前库|

`state_levels`固定为`c("Bystander", "Infected")`，决定比较方向。当前脚本会拒绝其他顺序；若原对象用别的标签，先在前置注释中明确映射，不直接调换因子顺序。

默认`gene_set_source_mode <- "Broad"`，直接在原脚本的“完整Hallmark与GO:BP基因集”分节下载Broad官方**MSigDB 2026.1.Hs**两份人gene symbol GMT；不需要另一个下载脚本、手工上传文件或msigdbr/Zenodo入口。下载后校验固定MD5和50个Hallmark、7,538个GO:BP集合数，再由同一个GMT读取流程处理。它们是公共注释，不是用户单细胞数据，也没有按MDK、免疫主题或显著性预筛。

`GMT`为可选本地入口，可继续使用随附固定版本文件（见`gene_sets/README.md`）。两份文件须为同一版本的人gene symbol格式，每行“集合名、说明、成员基因”，制表符分隔；填写版本字符串不能替代检查文件来源，脚本不会自动识别自换GMT是否删过集合。`msigdbr`为另一可选入口，使用人数据库`db_species="HS"`和`species="Homo sapiens"`，分别取`H`及`C5/GO:BP`。**现在按`gene_set_source_mode`选择来源，不再靠把两条GMT路径留空来切换。**

## 4. 阈值、统计口径及适用范围

### 4.1 纳入及表达过滤

|参数/检查|默认值|作用与依据|
|---|---:|---|
|`min_cells_per_state`|20|每个独立样本的每种状态至少20个具有有效宿主计数的细胞；任一状态不合格则整对排除。是最低可测性门槛，不是生物学分界，也不保证20个细胞就稳定|
|`min_pairs`|2|至少两对独立样本。两对会警告且仅作探索；增加细胞不能增加生物重复数。不满足时停止，不能改用逐细胞P值替代|
|宿主计数|大于0|无有效宿主计数的细胞不参与汇总；有细胞但共同宿主pseudobulk总量为0的状态也不合格|
|`filter_min_count`|10|`filterByExpr`依据库大小换算CPM阈值，并要求在足够数量的pseudobulk中达到；不是要求每个细胞10UMI|
|`filter_min_total_count`|15|另外要求基因在纳入pseudobulk中的合计原始计数达到最低水平，去除极稀疏基因|
|共同feature/保留rank最低数|100|脚本的输入防错门槛，低于100停止；不是“100个基因即可可靠GSEA”的充分条件|
|基因集匹配比例|低于0.1停止；低于0.5警告|排序基因中能匹配H或GO:BP成员的比例，另要求至少100个匹配；用于检查物种/ID/库错误，不是富集检验阈值|
|并列rank比例|超过0.1警告|并列值按gene symbol固定顺序，不添加随机噪声；多并列时位置与leading edge需谨慎解释|

配对模型是`~ block + state`，`block`为`replicate_id`；样本项吸收样本间基线差异。表达过滤按待比较的`state`分组，依照edgeR配对设计指导，不按阻断项过滤。

TMM校正计数库的组成差异，离散度和quasi-likelihood（QL，准似然）模型用来描述计数与生物重复之间的不确定性；`robust=TRUE`用于减小异常高变基因对估计的影响。它们不能修复错误的样本对应、细胞混入、病毒标签或不足的独立重复。当前模型不调整感染状态内的所有潜在肿瘤亚状态，也不能把观察性差异变成MDK因果效应。

### 4.2 全基因排序与GSEA

`rank_stat = sign(logFC) × sqrt(F)`，其中F来自edgeR单系数QL检验。F本身非负，乘以logFC的符号才能形成双向排序。这是构造的有符号QL排序统计量，不是edgeR原生输出的t统计量。

- `logFC > 0`、正rank：Infected中相对更高；负值：Bystander中相对更高。
- 所有通过表达过滤、统计量有限的共同宿主基因进入排序，包括不显著基因；不按P值/FDR或logFC阈值先筛DEG。
- 未进入某个基因集的排序基因仍保留在游走背景中，不先把排序截成Hallmark/GO成员的并集。
- `rank_position`从1开始，越靠前越偏Infected；`rank_percentile`从0（最前）到1（最后），不是表达百分比。

GSEA沿完整排序考察一个基因集的成员是否集中在某一端。ES为加权游走的极值，NES为经过对应随机集合零分布归一化的ES。**正/负NES仅表示这些基因的转录排序倾向，不直接等于过程被活化/抑制，更不是T细胞功能增强/受损。** 特别是在本脚本只分析的肿瘤上皮中，出现T细胞GO词不代表直接测到了T细胞反应。

|GSEA参数|默认值|含义|
|---|---:|---|
|`min_set_size/max_set_size`|15/500|集合与本次rank背景交集后的大小范围；过小集合易由少数排序项主导，过大集合缺少特异性。是预设分析范围，修改后应整体重跑并记录|
|`scoreType`、`gseaParam`|`std`、1|双向GSEA；命中步长按rank绝对值加权|
|`gsea_eps`|1e-10|P值估计下界；低于此精度的数值不应继续作精细大小比较，不代表样本间证据可无限精确|
|`gsea_simple_n`|1000|fgsea算法内部简单抽样阶段的次数，不是细胞数、动物数，也不是动物标签置换次数|
|`seed`|123|固定随机抽样；仍应保存包/数据版本。不同版本及输入变更可能改变结果|
|`fdr_cutoff`|0.05|用于标记BH调整后P值的展示阈值；未过阈值不能写成已证差异，接近阈值不宜二分夸大|

### 4.3 FDR家族及计算状态

`genes.csv`中的`FDR`是基因层面的edgeR BH校正；`summary.csv`中的`padj_family/padj_global`是通路层面的BH校正，不能互换。

GO条目之间常有成员重叠。同一组核心基因可对应多个显著条目，不能把这些条目数当作相互独立的生物学证据；阅读时应同时检查leading edge及条目定义。

- `padj_family`：分别在完整的H和完整的GO:BP中，对满足大小及非零命中权重条件的全部集合校正。用于默认图形的显著性标记。
- `padj_global`：将H与GO:BP中所有可测试集合合并校正，供需要跨两个集合家族统一判断时查看。
- `focus_pathways`、概览Top列表及MDK成员集合是在完整计算之后选展示，不对这些子集另算较小家族的FDR。
- `below_min_size/above_max_size`为交集大小未达范围；`covers_rank_universe`为集合覆盖整个rank背景，此时游走统计量不可定义；`zero_hit_weight`为集合命中rank绝对值和为0；`fgsea_unresolved`为算法未获得有限P/NES。它们不是“阴性通路”。
- 计算失败的NA保留；校正时仍按该次完整可测试家族数量计数，不把NA擅自填0或1。去MDK后集合资格可能变化，两个run的合格家族数不保证完全相同。

## 5. 怎样观察MDK，怎样避免循环解释

### 5.1 表达、检出与DE状态

MDK表达图每个点代表一个独立生物样本的一种状态，灰线仅连接同一真实样本。表达指标为逐细胞`log1p(10000 × MDK_count / 共同宿主feature总count)`后求均值，检出率为MDK原始count大于0的有效细胞比例；它们是描述性指标，不等于pseudobulk的edgeR logFC、logCPM或每细胞绝对分子数。共同背景分母不同于未经feature交集限制的整个人转录组总UMI，需按此口径报告。

`target_status`区分：

|状态|含义|
|---|---|
|`ranked`|MDK进入完整rank，可以检查其排序位置、成员及leading edge|
|`not_in_shared_features`|MDK不在所有纳入层共有的宿主feature中；检查`layer_audit`，不能解释为表达为0|
|`expression_filtered`|MDK在共同背景中，但未通过edgeR最低表达过滤；保留原始计数信息，统计量和rank为NA|
|`non_finite_statistic`|MDK统计量不有限，不能进入GSEA排序；不是无差异|

### 5.2 注释成员、leading edge与命中权重

|结果字段|定义|允许与不允许的解释|
|---|---|---|
|`target_in_database_set`|原始数据库成员表是否列出MDK|说明注释关系；不等于本批数据证实MDK作用于该过程|
|`target_in_rank`|MDK是否进入本次完整rank；对所有集合相同|说明可评估性；不是某个集合中特异表达|
|`target_in_leading_edge`|MDK是否落在该集合产生极端ES的核心排序区段|说明这一排序中的核心成员；不证明MDK是上游驱动或必要因子；未可评估时为NA|
|`target_hit_weight`|权重1时，`abs(rank_MDK) / sum(abs(rank_集合所有可测成员))`|只是MDK占集合总命中步长权重的比例；不是占ES、NES或生物效应的比例，也不是T细胞贡献比例|

MDK可为数据库成员但未进入rank，或进入rank但不在某集合leading edge。MDK在rank但不属于该集合时`target_hit_weight=0`；MDK未进入rank时该字段为NA。即使属于leading edge，也应同时查看其他核心基因和样本级表达，避免由单个目标基因反推整个机制。

### 5.3 删除MDK排序项的敏感性

默认`remove_target_sensitivity=TRUE`且MDK状态为`ranked`时，从rank中删去MDK，保留其他基因的统计量不变，重新运行完整H与GO:BP。**不重新拟合DE、不重归一化、不改动细胞，也不是虚拟敲除或MDK敲低实验。** 若MDK不在rank则不运行，并输出明确状态。

检查同一集合的方向、ES/NES、P/FDR、成员交集大小和leading edge是否稳定；`direction_same_without_target`只判断两次有限NES是否同号。删除一个排序项会改变总背景、集合大小和随机零分布，不能把ΔNES或“过/不过0.05”解释为MDK的因果贡献。如果集合由15个可测成员变成14个，它将成为`below_min_size`，应写“未再达到预设大小阈值”，不能写“去MDK后富集消失”。

不属于某集合的MDK被删去，也可能通过背景变化影响数值。只有两对样本或结果接近阈值时，这种数值变化尤其不宜被过度解释。真实MDK因果仍需独立扰动、回补及相应功能验证。

## 6. 依赖、安装与运行

显式依赖：`Seurat`、`Matrix`、`edgeR`、`statmod`、`fgsea`、`BiocParallel`、`ggplot2`、`patchwork`。默认`Broad`模式由R自身的HTTPS下载功能读取官方GMT，需要能访问`data.broadinstitute.org`，不加载msigdbr；`GMT`模式不下载基因集；只有`msigdbr`模式额外需要该包及其数据源连接。Seurat会使用其匹配的SeuratObject；edgeR使用limma等依赖。`ggplot2`需支持`linewidth`参数（3.4.0及以后）。`msigdbr`模式要求`msigdbr >= 10.0.0`的现代接口；实际数据库版本另记录，不能仅靠包名判断。

先在当前能正常读取原对象的R环境检查缺包，再只安装缺失项。不要为了本脚本盲目升级用户Seurat、SeuratObject或R。Bioconductor包须安装与当前R对应的发行版本；较旧R若无法安装当前包，可使用其兼容发行版并核查接口，或对基因集改用固定版本离线GMT。脚本本身不会自动升级软件或调用`UpdateSeuratObject`。

首次缺包时，可在R控制台参考以下安装方式；已有的包不用重复装：

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "fgsea", "BiocParallel"), ask = FALSE, update = FALSE)
install.packages(c("statmod", "ggplot2", "patchwork"))
# 只有主动改用gene_set_source_mode="msigdbr"时才需要：install.packages("msigdbr")
```

`Seurat/SeuratObject/Matrix`优先复用既有匹配环境。如果安装提示R版本或依赖不兼容，先检查`R.version.string`、`packageVersion()`及Bioconductor版本，不通过强制混装解决。基因集来源切换不要求升级这些包；首次安装R包仍需可用的软件源或服务器已有依赖。msigdbr 26.1.1曾修复同源映射缓存问题；本脚本采用人库到人符号，不进行跨物种转换。

### 6.1 默认Broad官方入口与缓存

保持开头`gene_set_source_mode <- "Broad"`即可。下载、校验和读取都在原脚本内完成；无需先运行外部下载工具。下载对象为固定版本的：

- `h.all.v2026.1.Hs.symbols.gmt`：完整Hallmark，50个集合。
- `c5.go.bp.v2026.1.Hs.symbols.gmt`：完整GO:BP，7,538个集合。

两份文件合计约**4.69 MiB**，缓存在`file.path(tempdir(), "gsea_broad_2026_1")`。这是当前R会话的临时目录，可能位于服务器临时盘；不是输出目录，也不是用户原始数据目录。相同会话已有且通过MD5的文件可复用；结束会话或系统清理后可能需要重新下载。下载过程还需为临时文件保留少量额外空间，分析的R内存需求另由输入对象与结果决定。

默认每份文件**最多3次尝试，每次下载超时设为180秒或原R更长的超时值**。成功即继续，不会强制等待180秒；失败、未下载完整或校验不符不能作为有效GMT继续使用。`download_timeout/download_retries`可在开头调整，慢速连接可增加超时。固定MD5用于核对文件内容，下载时再核对集合数和格式。它们不是通路显著性阈值。

通过校验后，本节设置实际`hallmark_gmt`、`gobp_gmt`与版本并沿用同一`read_gmt`流程。`results.rds$info$gene_set_source`记录`mode="Broad"`、官方URLs、缓存路径、MD5与版本。临时缓存将来可能不存在，但实际分析的`pathways`成员清单保存在结果中。

### 6.2 可选GMT和msigdbr入口

**已有本地GMT时**，设置`gene_set_source_mode <- "GMT"`，再填写：

```r
gene_set_dir <- file.path(path.expand("~/rawdata"), "gene_sets")
hallmark_gmt <- file.path(gene_set_dir, "hallmark.gmt")
gobp_gmt <- file.path(gene_set_dir, "gobp.gmt")
gmt_version <- "MSigDB 2026.1.Hs"
```

此模式下将随附两份GMT上传到该服务器目录，或改为自己的文件路径；不用上传`_analysis/gsea_gene_sets.rds`或用户单细胞对象作为基因集。`file.exists()`检查两文件可读，`tools::md5sum()`可与对应文件来源记录核对。随附导出版GMT与Broad原始GMT即使成员相同，也可能因排序或说明字段不同而具有不同MD5，不混用两种文件的散列值。

**明确需要msigdbr时**，设置`gene_set_source_mode <- "msigdbr"`并安装相应包。此模式可能仍需访问Zenodo等数据源；其包版本和数据库版本分别记录。切换来源只改变基因集入口，不改变既有统计与病毒排除逻辑；不同数据库版本或成员差异仍可能改变结果。

### 6.3 下载失败后续跑GSEA的条件

若原msigdbr/Zenodo出现`Could not resolve host`或连接超时，改为新版默认`gene_set_source_mode <- "Broad"`即可从官方固定版本读取；不是通过升级Seurat解决。Broad地址本身也需要服务器HTTPS连接，最多尝试后仍失败应检查报错及连接，或明确切换`GMT`，脚本不应把失败的部分文件当成完成结果。

1. 在当前会话执行下面四行Broad配置；若选`GMT`则改为该模式并配置真实路径和版本。

   ```r
   gene_set_source_mode <- "Broad"
   gene_set_dir <- file.path(tempdir(), "gsea_broad_2026_1")
   download_timeout <- 180L
   download_retries <- 3L
   ```

2. **只有当当前会话已经完成配对edgeR和完整rank，并有病毒排除函数、`viral_exclusion_audit`及其他保存所需的上游对象，且确认原始counts阶段已正确排除病毒feature时**，才能从`##### 完整 Hallmark 与 GO:BP 基因集 #####`开始，依顺序重跑其后全部代码到保存结果。不要只复制下载几行而跳过本节的读取和完整性校验。
3. 续跑时不要执行开头清理环境，也不要用整段`source()`当作续跑，因为它会清空现存rank。若要完整重算，则使用下面的正常运行命令。
4. 若更改了病毒排除、样本表、计数层等上游配置，不能确认旧rank背景，或旧版失败会话缺少新版审计对象，应从头重跑。不能只从最终rank删除病毒名字：库大小、TMM、DE和MDK分母都可能需要重新计算。
5. 新R会话没有上游对象时同样完整运行。当前脚本不提供从旧`results.rds`自动恢复DE拟合的入口，不从最终rank伪造原始层排除审计。

### 6.4 运行命令

从项目根目录运行：

```sh
Rscript --vanilla code/scRNA_GSEA.R
```

或在R/RStudio中将工作目录设为项目根目录，运行：

```r
source("code/scRNA_GSEA.R", encoding = "UTF-8")
```

也可打开脚本逐节执行。脚本默认输出路径不依赖当前工作目录，但上述脚本路径需要从项目根目录定位。每次运行会覆盖同一`outdir`中的四个同名结果；不同日期/侧别请先改`analysis_label/outdir`。

### 6.5 仅更新本轮图形格式

本轮无需新增或修改开头配置，沿用已确认的输入、样本表与输出目录。若当前R会话保留主脚本完成分析后的全部对象，可从`##### 选择图中展示条目，不改变全结果和多重校正 #####`开始顺序运行到末尾：只重新选择展示、绘图和保存，不重新拟合DE或运行GSEA。不要执行开头清理环境。它会重写原来的四个输出，其中主要变化为`plots.pdf`；真实统计表应保持一致。若已经清空会话，则按正常流程完整运行；当前脚本没有只靠旧`results.rds`自动恢复全部绘图对象的专用入口。

## 7. 输出四个文件及阅读顺序

|文件|内容|建议先看什么|
|---|---|---|
|`genes.csv`|所有共同宿主feature的原始合计、过滤状态；合格基因的edgeR logFC/logCPM/F/P/FDR、rank及MDK标记|找到MDK行，确认`passed_filter`、方向和rank；未通过过滤的统计列为NA|
|`summary.csv`|完整H/GO:BP成员集合，原始/可测大小、测试状态、ES/NES/P、两类FDR、leading edge、MDK记录及去MDK结果|先分清`collection/test_status`，再看family/global FDR；不只读图中Top条目|
|`plots.pdf`|样本细胞数及MDK描述图、概览、预设关注条目、MDK注释集合、排序项敏感性及富集曲线|先核对纳入样本与缺失，再读方向和显著性图例；最多曲线数由配置控制|
|`results.rds`|紧凑结果：`info`、`sample_audit`、`layer_audit`、`viral_exclusion_audit`、`pb_meta`、设计矩阵、小型pseudobulk、归一化记录、完整rank/集合/结果及MDK子表|优先核对真实n、病毒实际排除清单、QC排除、层feature覆盖、版本、参数及缺失关注集合|

结果文件不含大型Seurat对象，也不另外生成时间戳文件或重复中间文件。运行时间、参数、输入路径/大小、包版本、MSigDB/GMT来源及`sessionInfo()`集中保存在`results.rds$info`。

可以在完成后新会话查看：

```r
gsea_result <- readRDS("~/OV/scRNA_virus/GSEA_D3/results.rds")
gsea_result$info$n_pairs
gsea_result$pb_meta
gsea_result$layer_audit
gsea_result$viral_exclusion_audit
gsea_result$target_de
gsea_result$target_pathways
```

### 图形约定与缺失信息

- 脚本导出的各页不附方法、阈值定义或指标释义等底部长caption。这些解释集中在本README第4–5节和`results.rds$info`；图中继续保留坐标单位、分组、颜色/形状/线型图例及必要的实际统计值。
- Bystander为钢蓝`#426B91`，Infected为陶红`#C35B3F`；这两种颜色在此表示感染状态，不是Vehicle/VG161实验组。
- 配对线为浅灰，每个点只代表相应独立样本/状态的描述值，不画逐细胞P值。
- 通路条形以NES发散色标表示方向，各通路页共用以0为中心的对称颜色范围，实色/淡色由`padj_family`是否小于阈值区分。概览按每家族每方向的family FDR优先选取，默认最多6项，**可包含不显著结果**；关注页同样不保证所有条目显著。
- 默认最多8条富集曲线，优先若干MDK成员集合，再纳入关注及概览条目；曲线中灰色短线为集合成员，MDK为成员时用陶红短线/虚线突出，并保留相应线型图例。曲线名称写在标题；实际NES、family FDR及MDK状态放在图内白底紧凑统计块，分行对齐。脚本保守检查四角的曲线和MDK竖线占用，优先使用空白角落；四角均有曲线时增加上方留白，不覆盖曲线。统计块固定为66 × 23 mm，适用于默认12 × 10英寸PDF及每页最多两条曲线；修改页尺寸或布局后须重新检查。
- MDK缺失或算法不可评估不画成0；图页缺少点、条目或敏感性端点时，先查看CSV/RDS的NA与状态。关注清单中当前库没有的名称记录在`info$missing_focus`；不自动改用相似名字。
- 图形不加A/B/C面板编号，采用白底、无衬线字体及细框。最终用于论文前仍需按目标尺寸检查文本和长通路名的可读性。

图内紧凑标记及图形元素按以下口径阅读；这些定义不在图底重复排长文字：

|图内项目|本README中的解释|
|---|---|
|`NES`实际数值|对应当前集合的归一化富集分数；正值偏Infected、负值偏Bystander。不能直接称过程活化/抑制，见第4.2节|
|`Family FDR`实际数值|该集合所属完整H或GO:BP可测试家族内的BH调整P值；不是只对当前图中集合校正，也不是跨H/GO统一校正，见第4.3节|
|MDK状态|区分leading edge、集合成员但不在leading edge、非该可测集合成员，以及不可评估状态。leading edge只是排序核心身份，不是因果作用，见第5.2节|
|NA或不可评估状态|说明相应数据、rank或统计量不可评估；不把它改写为0、未表达或已证明无效，见第2.3、4.3及5.1节|
|细胞数页的虚线|每样本每状态的最低细胞数门槛，默认20；这是QC标线，不是统计显著性，见第4.1节|
|MDK表达和检出率|共同宿主背景下的逐细胞log1p(CP10K)均值及原始count大于0的比例，见第5.1节|
|概览实色/淡色|分别表示family FDR小于/不小于配置阈值；默认0.05。图中Top条目也可不显著，见第4.3节|
|敏感性页的实心/空心点|完整rank与删除MDK排序项后的结果；点间连线仅帮助对应同一集合，不表示生物学配对干预，见第5.3节|

## 8. 实际验证范围与使用限制

- 已按脚本逐节核对README中的输入字段、流程、默认参数和四项输出定义。
- **合成数据端到端验证已通过。** 在R 4.4.0、Seurat 5.5.0、SeuratObject 5.4.0、edgeR 4.4.2和fgsea 1.30.0环境下，使用真实Seurat对象测试传统Assay与拆分Assay5两种计数结构；两者的pseudobulk、DE和rank一致。没有另行安装Seurat 4版本运行，旧包环境仍需按实际依赖验证。
- 已核对植入的差异方向及GSEA方向、零宿主计数细胞的审计/排除、输入metadata保留、某层缺MDK时不补零、MDK整体缺失与表达过滤状态的区分、仅一对合格样本时明确停止，以及去MDK使集合15→14时记录`below_min_size`而非宣称富集消失。
- **在线基因集入口已验证。** 使用msigdbr 26.1.1读取MSigDB 2026.1.Hs的50个Hallmark及7,538个GO:BP集合，并用合成计数完成完整在线入口的DE、GSEA、MDK敏感性和四文件输出。这些合成数据仅用于验证程序，不构成本项目生物学结果。
- **离线GMT及扩展病毒排除已通过端到端合成验证。** 配套GMT与源缓存的完整成员逐通路一致；测试中禁止网络入口，加入高计数的VG161/HSV前缀病毒特征及精确列表指定的无前缀病毒特征，确认其在pseudobulk、DE和GSEA rank中均被排除、审计正确，宿主IFIT1保留；独立重算的宿主pseudobulk与脚本一致。还验证了GMT缺失时在读取真实对象和DE之前明确停止。原有六类边界验收在本轮修改后再次通过。服务器文件上传与路径可用性仍需在运行端核对。
- **本轮Broad官方默认入口已通过端到端合成验证。** 在禁用msigdbr调用的测试环境中，由原脚本直接下载两份官方GMT并完成四文件输出；共4,921,236字节、7,588个集合，逐通路成员与既有完整GMT参照一致，配对DE和rank完全相同，病毒排除继续生效。另验证同一会话复用已校验文件且不改写缓存，以及错误MD5被拒绝、失败文件不作为有效输入。默认最多3次尝试的逻辑已核对；没有模拟所有网络故障，也未在用户服务器实际运行，不能据此宣称该服务器连接已恢复。
- 此前已渲染检查在线入口合成运行的7页PDF，核对长通路名、图例、配对图、MDK标记及富集曲线；该检查针对本轮格式修改之前的版本。
- **本轮图形修改已重绘验证。** 使用已完成的合成分析结果，只执行新版展示和导出段，生成并渲染检查9页主PDF及3页补充显示检查；已移除各页方法caption，检查了正/负NES、MDK核心成员/非成员/未进入rank、MDK线型图例，以及四角均有曲线时的上方留白。另以专门的图例显示样例核对“Member / Outside leading edge”、科学计数法FDR和NA，未将这些指定显示值作为真实分析结果。新旧`genes.csv`和`summary.csv`内容一致；本次未重新拟合DE或重算GSEA，也未重跑用户真实对象。真实数据的曲线、标题及改过的页尺寸仍需按最终输出复核。
- **真实对象重算：尚未完成。** 当前没有重新计算用户原始单细胞对象，也没有获得本项目真实GSEA结果；本文不预判MDK方向或通路显著性。
- 样本表的独立性、Epithelial恶性身份、原感染标签、人鼠/病毒比对、实际层feature覆盖及最终图形可读性均需结合真实数据进一步验证。
- 若仅有两对样本，报告样本级效应与不确定性，结论保留为探索。GSEA的算法P值不会创造更多生物重复，也不能替代后续独立重复和功能实验。

## 参考文献

1. Chen Y, Lun ATL, Smyth GK. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. *F1000Research*. 2016;5:1438. [DOI:10.12688/f1000research.8987.2](https://doi.org/10.12688/f1000research.8987.2)。本脚本采用edgeR QL流程；配对设计与表达过滤参见[edgeR用户指南](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)。
2. Squair JW, et al. Confronting false discoveries in single-cell differential expression. *Nature Communications*. 2021;12:5692. [DOI:10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2)。用于生物重复和样本级分析口径。
3. Subramanian A, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*. 2005;102:15545–15550. [DOI:10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)；[GSEA官方用户指南](https://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/)。
4. Korotkevich G, et al. fgsea: Fast Gene Set Enrichment Analysis. [软件DOI:10.18129/B9.bioc.fgsea](https://doi.org/10.18129/B9.bioc.fgsea)；[官方教程](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html)。运行所用版本由`info$package_versions`记录。
5. Liberzon A, et al. The Molecular Signatures Database hallmark gene set collection. *Cell Systems*. 2015;1:417–425. [DOI:10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)。
6. MSigDB人基因集官方说明：[集合总览](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)、[GO:BP集合](https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=GO%3ABP)；Broad官方[2026.1.Hs发布目录](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/)、[Hallmark GMT](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt)、[GO:BP GMT](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/c5.go.bp.v2026.1.Hs.symbols.gmt)。实际成员和数据库版本以本次保存结果为准。
7. Dolgalev I. msigdbr官方文档：[读取接口](https://igordot.github.io/msigdbr/reference/msigdbr.html)、[集合查询](https://igordot.github.io/msigdbr/reference/msigdbr_collections.html)、[版本记录](https://igordot.github.io/msigdbr/news/index.html)。
8. Seurat/SeuratObject官方文档：[Assay5的Layers/LayerData](https://satijalab.github.io/seurat-object/reference/Layers.html)、[GetAssayData](https://satijalab.github.io/seurat-object/reference/AssayData.html)。用于提取原始RNA计数，不据此执行重新注释或整合。
