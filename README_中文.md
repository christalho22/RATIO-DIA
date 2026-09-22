# RATIO-DIA 使用说明

RATIO-DIA用于将DIA-MS/MS形成的过度连接分子网络拆分为拓扑聚焦模块，并结合
不同部位的化学丰度和生物活性梯度筛选优先候选节点。

## 输入文件

- MGF：每个谱图的`SCANS`必须与丰度表的`Alignment ID`一致。
- 丰度表：包含`Alignment ID`、`Average Rt(min)`、`Average Mz`和A-E列。
- 活性表：长格式，包含`Group`和`Response`；组别为Model、Control及A-E。

## 安装与运行

```bash
pip install -e .
ratio-dia --mgf examples/demo.mgf --abundance-table examples/abundance_table.csv --activity-table examples/activity_response.csv --output outputs/demo --min-cosine 0.70 --min-matched-peaks 6 --partition-method relation-aware --input-scope edge-connected --pair-motif-weight 0.25 --motif-k 15 --motif-minimum 0.10 --fragmentation-floor 0.10 --topology-floor 0.25 --motif-scale 0.80 --louvain-resolution 2.25 --random-seed 20260914
```

主分析使用`relation-aware`。在Modified cosine原始网络上分别提取ECV局部
拓扑支持、单碎片TF-IDF和碎片对共现TF-IDF证据，并形成非互惠的motif邻域
并集；原始谱图边接受连续重加权，motif关系以加性权重补充。随后通过固定参数
网格进行加权Louvain多分辨率划分。模块规模约束在候选方案筛选阶段应用，而不
是在Louvain优化内部强制执行。

在结构鉴定之前，将前体离子*m/z* 571.2803、518.2184、559.2807等相关特征
离子的共同聚类表现用于参数优化。筛选顺序依次为：单一模块中的目标离子数、
目标模块紧凑度、目标离子对共聚类比例、>10节点模块的覆盖量和加权模块度。
结构身份、部位丰度和活性标签不进入网络构建；活性信息仅在模块确定后映射。
`ecv-consensus`和`ecv-hierarchy`仅保留为旧版复现及敏感性分析。

运行后将得到完整相似度边表、模块节点表、模块内边表、模块汇总表、活性候选
节点表及活性候选聚焦网络。Cytoscape导入时，以`Source_Scan`和`Target_Scan`
作为边的两端，以`Alignment_ID`作为节点属性匹配键。

## 结果解释

`Module_ID=0`表示节点未进入满足最小规模条件的正式模块。Tier 1和Tier 2用于
候选优先排序，不代表结构确证或活性因果关系。DIA共碎裂、加合物、源内碎片、
聚集离子和共洗脱异构体仍需通过标准品、HRMS/MS、NMR及独立活性实验验证。

完整参数说明和输入格式见英文版[README](README.md)及`docs`目录。
