# 题目、真题、去重与题族

## 题目 ID 与字段

每道题建立唯一 ID：`Q000001`、`Q000002`……

```yaml
question_id:
chapter:
knowledge_points:
question_type:
source_type:
source_file:
source_page:
original_text:
normalized_text:
difficulty:
frequency:
question_family:
answer_status:
potential_subjective_type:
```

## 题型分类

统一使用：

```text
名词解释
简答题
论述题
选择题
填空题
判断题
其他
```

如果一道题可以同时属于多个类型：保留其原始类型，并增加 `potential_subjective_type`。

例如：原始为选择题，潜在主观题类型为简答题。

## 真题识别规则

只有满足资料明确证据时，才能标记 `source_type = 真题`。

证据包括：

- 明确年份
- 明确学校
- 明确考试名称
- 明确「真题」「考研真题」等标识
- 资料上下文明确说明题目来源

如果只有题目，没有来源信息：`source_type = 未知`。

禁止根据题目风格猜测年份。

## 真题信息

如果确定是真题，记录：

```yaml
exam_year:
exam_university:
exam_subject:
question_number:
source_file:
source_page:
original_text:
```

例如：`2022年 / 上海交通大学 / 338生物化学 / 名词解释第2题`

上交 338 真题与其他院校真题必须分开存放、分开标记。

## 题目标准化

去重之前建立：

- `original_text`：永远不能被覆盖
- `normalized_text`：用于匹配，例如「竞争性抑制的概念及动力学特点」

原题：

```text
什么是酶的竞争性抑制？其动力学特点是什么？
```

标准化：

```text
竞争性抑制的概念及动力学特点
```

## 题目去重（三级）

### Level 1：完全重复

文字基本完全相同。

处理：保留一个主记录，记录多个来源。

### Level 2：表述不同但本质相同

例如「解释竞争性抑制」与「什么是竞争性抑制？它有哪些特点？」

归为同一题族。

### Level 3：知识点相同但考查角度不同

例如「什么是竞争性抑制？」与「比较竞争性抑制与非竞争性抑制。」

不能直接删除。必须归入：同一知识点、不同题型/题族。

## 题族系统（核心）

不要把题目看成孤立的题。必须建立：

```text
知识点 → 考查角度 → 题族 → 具体题目
```

例如：

```text
竞争性抑制
├── 定义题
├── 特点题
├── 动力学题
├── 与非竞争性抑制比较
├── Lineweaver-Burk图分析
└── 综合论述题
```

## 题族命名规则

题族名称应该描述「考什么」，而不是简单复制题目。

推荐角度：

```text
定义类 结构类 功能类 机制类 特点类 比较类
过程类 调控类 意义类 实验原理类 应用类 综合分析类
```

例如：

```text
酶的竞争性抑制—定义与特点题族
酶的竞争性抑制—动力学变化题族
糖酵解—过程与关键步骤题族
糖酵解—调控机制题族
```

题族 ID：按章 `TF001`、`TF002`……

## 题目难度

```text
★        单一知识点直接记忆
★★       单一知识点理解
★★★      多个知识点结合
★★★★     机制分析或综合比较
★★★★★    跨章节综合、复杂机制或高难度论述
```

## AI 生成题规则

AI 只能在已有知识点和题族基础上生成题目。

必须标记：`source_type = AI生成题`

并记录：

```yaml
generated_from:
generation_reason:
related_original_questions:
```

例如：`generated_from: 竞争性抑制`；`related_original_questions: Q00121, Q00456, Q00981`

所有 AI 生成内容必须带 `[AI生成]`。原题不得添加该标记。

## AI 题目生成原则

优先模拟：已有真题考法 + 课后题考法 + 高频知识点 + 知识点之间的联系。

禁止无依据制造非常偏门的题目。

## 真题改编

如果根据真题生成相似题：

必须标记 `source_type = 真题改编`。

不得标记成 `真题`。
