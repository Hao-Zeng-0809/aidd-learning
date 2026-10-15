# 第 2 周练习：用 RDKit 计算药物分子的理化性质

> 对应 [12周执行计划.md](../../12周执行计划.md) 第 2 周 D1–D5
> 难度：入门–中等　建议用时：4–5 小时（可拆成 3 天）

---

## 背景

你手上有 25 个**你药化课上一定学过的药**（`data/drugs.csv`）。
现在要用 RDKit 把它们从「结构」变成「数字」，并回答一个药学问题：

> **这些药符不符合 Lipinski 五规则？如果不符合，为什么它还能上市？**

这不是编程练习，是**用计算工具回答药学问题**。编程只是手段。

---

## 任务

### 第一步：计算理化性质

读入 `data/drugs.csv`（25 个药），为每个分子计算下面这些量：

| 列名 | 怎么算 | 物理意义 |
|---|---|---|
| `MolWt` | `Descriptors.MolWt` | 分子量 |
| `MolLogP` | `Descriptors.MolLogP` | 脂水分配系数（渗透性） |
| `TPSA` | `Descriptors.TPSA` | 拓扑极性表面积 |
| `NumHDonors` | `Descriptors.NumHDonors` | 氢键供体数 |
| `NumHAcceptors` | `Descriptors.NumHAcceptors` | 氢键受体数 |
| `NumRotatableBonds` | `Descriptors.NumRotatableBonds` | 可旋转键（分子柔性） |
| `RingCount` | `Descriptors.RingCount` | 环数 |
| `FractionCSP3` | `Descriptors.FractionCSP3` | sp3 碳比例（三维性） |
| `QED` | `Descriptors.qed` | 类药性综合评分（0–1，越大越好） |

### 第二步：判断 Lipinski 五规则

**五规则**（口服药物经验规律）：

| 条件 | 阈值 | 怎么判断 |
|---|---|---|
| 分子量 | ≤ 500 | `MolWt > 500` 算违反 |
| logP | ≤ 5 | `MolLogP > 5` 算违反 |
| 氢键供体 | ≤ 5 | `NumHDonors > 5` 算违反 |
| 氢键受体 | ≤ 10 | `NumHAcceptors > 10` 算违反 |

- 新增列 `violations`：违反了几条（0–4 的整数）
- 新增列 `passes_ro5`：`violations <= 1` 为 `True`（**注意：违反 1 条通常仍可接受**）

### 第三步：输出

**`drugs_properties.csv`** —— 保留原 4 列，加上上面 11 个新列，共 15 列。

**列顺序**（`check.py` 会检查）：

```
name_cn, name_en, smiles, drug_class,                                  ← 原 4 列
MolWt, MolLogP, TPSA, NumHDonors, NumHAcceptors,                       ← 第二步的表
NumRotatableBonds, RingCount, FractionCSP3, QED,
violations, passes_ro5                                                 ← 最后两列
```

> 💡 别依赖中间过程留下的顺序。最后统一写一句
> `df = df[列名列表]` 显式指定，就不会错。

**`drugs_grid.png`** —— 分子结构网格图，4 列 × 7 行（25 个分子），
图例用**英文名**（原因见下面的坑）。

### 第四步：思考题（写进你的仓库 README，这是重点）

1. **哪些药违反了 Lipinski 规则？各违反几条？** 运行 `check.py` 能看到这份名单。
2. **阿托伐他汀（Atorvastatin）违反了哪两条？它为什么还能口服上市？**
   提示：想想给药剂量、以及五规则本来是对"口服小分子"的**统计经验**，不是硬性门槛。
3. **`FractionCSP3` 高的分子有什么特点？** 想一想它和"三维性""成药性"的关系。
   提示：查一下 "flatland" 或 "escaping from flatland" 这个概念。
4. **`MolLogP` 超过 5 意味着什么实际问题？** 从吸收、分布的角度说。
5. 在这 25 个药里，**logP 最高的是哪个？最低的是哪个？这符合你的直觉吗？**

---

## 动手

```powershell
conda activate chemai
cd C:\Users\Aurora\Documents\deepseek-harness\default-workspace\exercises\week2-rdkit
python template.py      # 填空
python check.py         # 自检
```

---

## ⚠️ 三个必须知道的坑

### 坑 1：RDKit 画图**不支持中文**

如果你把图例写成 `"阿司匹林 Aspirin"`，输出的图上**只会显示 "Aspirin"**，
中文部分被静默丢弃（不报错，只是没有）。因为 RDKit 自带字体没有中文字形。

**所以图例请用英文名。** 报告里的图注用中文说明即可。

（如果你一定要在图上显示中文，需要用 `rdMolDraw2D` 指定一个中文 TTF 字体，
那是进阶内容，这周不做。）

### 坑 2：`MolFromSmiles` 可能返回 `None`

解析失败时返回 `None`，后面所有计算都会崩。
**必须**判断：`if mol is None: 跳过并计数`，最后打印跳过了几个。

### 坑 3：`Draw.MolsToGridImage` 还是 `Draw.MolToGridImage`

- 多个分子用 **`MolsToGridImage`**（有 s）
- 单个分子用 `MolToImage`

拼错会得到 `AttributeError`。另外网格图的参数是 `molsPerRow` 和 `subImgSize`。

---

## 常见坑

**通用报错查这份表**：[pandas 调试速查表](../../guides/pandas调试速查表.md)

| 现象 | 原因 |
|---|---|
| `ImportError: DLL load failed ... rdMolDraw2D` | cairo 依赖问题。**已经帮你修好了**，修复记录见 [项目 README](../../projects/qsar-egfr/README.md) 的坑 6 |
| 图上中文不见了 | 见坑 1，正常现象，用英文 |
| `AttributeError: module 'rdkit.Chem.Descriptors' has no attribute 'QED'` | 是**小写** `Descriptors.qed`。RDKit 的 API 在版本间会改名，遇到 `AttributeError` 先查官方文档或 `dir(Descriptors)` |
| `MolFromSmiles` 返回 `None`，后面全崩 | 必须先判断 `if mol is None: 跳过并计数` |
| 输出列顺序乱了 | 用 `df[列名列表]` 明确指定顺序 |
| 图例太长把图挤变形 | 用 `name_en`，别用 `name_cn + name_en` 拼起来 |
| `Draw` 模块导入成功但画图很慢 | 正常，RDKit 画图是 CPU 密集的；批量画图建议分块 |

卡超过 40 分钟就带着代码和报错来问我。
