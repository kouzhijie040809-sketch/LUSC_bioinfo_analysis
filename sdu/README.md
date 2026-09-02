# sdu 工作区（上交 338 生物化学）

本目录对应用户电脑上的 `D:\sdu`。

云端 Agent **不能直接读取** 你电脑硬盘。它只能读取本 GitHub 仓库里的 `sdu/`。因此：

- 题库工作位置已经迁到 `sdu/`
- 学习资料仍需复制/上传到 `sdu/生化资料/338（生物化学）资料/`，路径才和本地一致

若要让 Agent **直接读 D:\sdu 里的原文件**：请用 Cursor 桌面版打开 `D:\sdu`，并把本仓库的 `.cursor/` 复制到 `D:\sdu\.cursor\`，再用本地 Agent。不要用当前这个绑定 LUSC 仓库的 Cloud Agent。

---

## 目录（与本地对齐）

```text
sdu/                                      ← 对应 D:\sdu
├── README.md
├── AGENT_PROMPT.md
├── source_manifest.md
├── 生化资料/
│   └── 338（生物化学）资料/              ← 对应 D:\sdu\生化资料\338（生物化学）资料
│       └── （把 PDF / 笔记放这里）
└── 338-题库/                             ← 生成的主观题库（不要和原资料混放）
    ├── _index.md
    └── Chapter_01/
```

---

## 放入资料后怎么开始第 1 章

```text
/biochem-338-question-bank
从第1章开始
```

Agent 只处理第 1 章，完成后停止，等你看效果。
