# 上海交通大学 338 生物化学主观题库

这个目录是 **上交 338 生物化学** 的学习资料整理工作区，与仓库里的 LUSC 生信分析互不影响。

目标：把教材知识点总结、课后习题、考研真题里混杂的题型，按章节整理成可以背诵和默写的主观题：

- 名词解释
- 简答题
- 论述题

选择题、判断题、填空题不会作为最终练习题，但其中的知识点会被提取进知识库。

---

## 你怎么用（推荐流程）

### 1. 放入学习资料

把 PDF / Word / Markdown / 图片笔记放到：

```text
sjtu-338-biochem/sources/
```

建议文件名能看出内容，例如：

```text
sources/
├── 教材总结-蛋白质.pdf
├── 课后习题-全册.pdf
└── 上交338真题.pdf
```

来源不清楚也没关系，Agent 会登记为「未确定」，不会擅自猜成真题。

### 2. 启动题库构建 Agent

在 Cursor Agent 对话中输入：

```text
/biochem-338-question-bank
从第1章开始
```

或明确指定专用 Agent：

```text
/biochem-338-builder
从第1章开始
```

云端 Agent（Cloud Agent）同样适用：新建 Agent 时把上述两句话作为任务说明，并确保本仓库已包含 `.cursor/skills/` 与 `sjtu-338-biochem/`（先合并本 PR）。

Agent 只会处理 **一章**。第 1 章全部输出并质检后会停下来等你确认，不会自动进入第 2 章。

### 3. 确认后再下一章

第 1 章没问题后，对 Agent 说：

```text
开始第2章
```

### 4. 按章节背题

题库建成后，不要再让 Agent 重新“发明”一整章题。改用复习 skill：

```text
/biochem-338-review
给我第1章所有高频名词解释
```

其他常用指令：

```text
给我第1章上交338真题
给我第1章所有简答题
给我第1章最重要的10个题族
按照真题风格给我出5道论述题
检测我第1章哪些知识点没有掌握
把第1章所有题按照 P0/P1/P2 排序
给我第1章背诵版答案
随机抽我10道主观题
根据我做错的题生成变式题
```

---

## 目录说明

```text
sjtu-338-biochem/
├── README.md                 ← 本说明
├── source_manifest.md        ← 所有输入资料清单（Agent 维护）
├── sources/                  ← 你放入的原始学习资料
└── knowledge-base/
    ├── _index.md             ← 各章完成状态
    └── Chapter_01/           ← 每章完整题库（处理该章后生成）
        ├── 01_知识框架.md
        ├── 02_核心知识点.md
        ├── 03_原始题目库.md
        ├── 04_真题库.md
        ├── 05_题族.md
        ├── 06_名词解释.md
        ├── 07_简答题.md
        ├── 08_论述题.md
        ├── 09_易错点.md
        ├── 10_高频考点.md
        ├── 11_章节覆盖率.md
        └── processing_log.md
```

---

## Agent 不会做的事

- 编造真题年份、学校、题号
- 把 AI 生成题标成历年真题或课后原题
- 删除或改写资料里的原题原文
- 为了覆盖率 100% 硬编一堆偏题
- 一次做完全部章节
- 改动仓库里的 LUSC 生信分析文件

---

## 相关 Cursor 配置

| 类型 | 路径 | 作用 |
| --- | --- | --- |
| 建库 Skill | `.cursor/skills/biochem-338-question-bank/` | 读取→分类→去重→题族→主观题→答案 |
| 复习 Skill | `.cursor/skills/biochem-338-review/` | 从已有题库抽题背诵 |
| 专用 Agent | `.cursor/agents/biochem-338-builder.md` | 只负责建库的 Subagent |
| 规则 | `.cursor/rules/biochem-338.mdc` | 在本目录工作时自动约束 Agent |
