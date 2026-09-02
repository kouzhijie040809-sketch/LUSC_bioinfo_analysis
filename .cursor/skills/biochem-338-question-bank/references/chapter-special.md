# 章节特殊处理

## 对比题自动识别

如果两个知识点具有明显对立关系，自动建立比较题，并建立 `comparison_group`。

例如：

```text
竞争性抑制 vs 非竞争性抑制
DNA vs RNA
α-螺旋 vs β-折叠
糖酵解 vs 糖异生
转录 vs 翻译
```

生成比较维度（只保留实际适用维度）：

```text
定义 位置 底物 关键酶 机制 能量 调控 意义
```

## 代谢章节

必须额外建立：

```text
pathway
location
substrate
product
key_enzyme
rate_limiting_enzyme
energy
cofactor
regulation
hormonal_regulation
relationship_with_other_pathways
physiological_significance
```

每条代谢途径必须能够回答：

```text
在哪里发生？
从什么开始？
生成什么？
关键步骤是什么？
关键酶是什么？
消耗/产生什么？
如何调控？
与哪些途径联系？
有什么意义？
```

代谢论述优先逻辑：

```text
场所 → 底物 → 关键反应 → 关键酶 → 能量变化 → 调控 → 与其他代谢途径联系 → 生理意义
```

## 酶章节

必须特别记录：

```text
enzyme_name
substrate
product
reaction
specificity
kinetics
Km
Vmax
inhibition
regulation
cofactor
```

尤其关注：

```text
Km
Vmax
竞争性抑制
非竞争性抑制
反竞争性抑制
Lineweaver-Burk
```

## 蛋白质章节

重点建立：

```text
amino_acid
peptide_bond
primary_structure
secondary_structure
tertiary_structure
quaternary_structure
disulfide_bond
denaturation
renaturation
structure-function_relationship
```

重点建立逻辑链：`结构 → 性质 → 功能`

## 核酸章节

重点建立：

```text
DNA_structure
RNA_structure
replication
transcription
translation
mutation
DNA_RNA_comparison
```

尤其关注：结构差异、复制机制、转录机制、翻译机制、调控。
