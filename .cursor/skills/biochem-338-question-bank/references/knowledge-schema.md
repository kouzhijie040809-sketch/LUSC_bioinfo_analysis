# 知识点体系

## 三级知识结构

每章建立三级知识结构。最终题目必须尽可能绑定到 Level 3。

### Level 1：一级主题

例如：`蛋白质`

### Level 2：核心知识模块

例如：

```text
氨基酸
蛋白质一级结构
蛋白质高级结构
蛋白质结构与功能
```

### Level 3：可考知识点

例如：

```text
肽键
二硫键
蛋白质一级结构
α-螺旋
β-折叠
蛋白质变性
蛋白质复性
```

## 知识点属性

每个知识点建立以下字段：

```yaml
knowledge_id:
chapter:
topic:
subtopic:
knowledge_name:
definition:
mechanism:
structure:
function:
relationship:
comparison:
regulation:
clinical_or_application:
common_question_angle:
source:
importance:
frequency:
mastery:
```

### importance

```text
核心
重要
一般
边缘
```

### frequency

```text
高频
中频
低频
未知
```

### mastery

```text
未学习
初步掌握
基本掌握
熟练
```

默认新抽取的知识点 `mastery = 未学习`。

## 知识点 ID

按章编号，例如第 1 章：`K001`、`K002`……  
跨章引用时可用 `C01-K001` 形式，但章内文件仍用短 ID。

## 客观题知识点提取

即使题目是选择题、填空题、判断题，也必须提取其背后的知识点。

例如：

```text
题目：下列哪种氨基酸含有硫元素？
```

不能只保存题目。必须抽取：

```text
知识点：含硫氨基酸
关联知识：Met / Cys / 氨基酸结构 / 氨基酸分类
```

并将其加入知识点库。

## 重要性评分

可以建立 `importance_score`，建议：

```text
真题出现 +5
多个资料出现 +3
课后题出现 +2
可形成简答题 +2
可形成论述题 +3
属于章节核心知识 +2
仅选择题出现 +1
```

最终不要机械依赖分数。必须结合实际内容判断。

## 题目优先级

```text
P0 必须掌握：高频真题 / 多资料重复 / 核心基础 / 高频主观题
P1 重点掌握：多次出现 / 常见简答或论述角度 / 与多个知识点相关
P2 需要理解：中低频 / 较少出现 / 可能作为综合题组成部分
P3 了解即可：边缘知识 / 很少出现 / 缺乏主观题价值
```

## 跨章节关系

必须主动识别跨章节关系，建立 `cross_chapter_relationship`。

例如：

```text
蛋白质 → 酶 → 糖代谢 → 生物氧化 → 脂质代谢 → 氨基酸代谢
```

跨章节知识优先用于：论述题、综合题、高难度简答题。
