# 基础练习包

> 配套 [12周执行计划.md](../12周执行计划.md) 第 1–2 周
> 大项目在 [projects/qsar-egfr](../projects/qsar-egfr/README.md)（第 3–6 周）

---

## 为什么有这两个练习

第 3–6 周的项目是一个完整的 QSAR 流程，但它建立在两个基础能力上：

1. **会用 pandas 处理真实脏数据** —— 真实数据永远是脏的，这不是运气问题
2. **会用 RDKit 把分子变成数字** —— 这是所有 AIDD 工作的第一步

直接上大项目，你会在这两件事上反复卡壳，然后误以为"我不适合做这个"。
所以先把基础练熟。

---

## 两个练习

| 练习 | 主题 | 难度 | 用时 | 交付物 |
|---|---|---|---|---|
| [week1-pandas](week1-pandas/README.md) | 清洗药物活性数据 | 入门 | 3–4 小时 | `clean_activities.csv` |
| [week2-rdkit](week2-rdkit/README.md) | 计算分子理化性质 + Lipinski | 入门–中等 | 4–5 小时 | `drugs_properties.csv` + 分子图 |

每个练习都有：
- `README.md` —— 任务说明、精确规则、思考题
- `template.py` —— 拆好步骤的骨架，你填空
- `check.py` —— **自动批改**。告诉你哪一步错了、可能是什么原因，但不给答案

---

## 怎么用（重要）

```powershell
conda activate chemai

# ⚠️ 用完整路径！你的终端默认位置不是这个工作区，相对路径会失败
cd C:\Users\Aurora\Documents\deepseek-harness\default-workspace\exercises\week1-pandas
# 或 week2-rdkit

python template.py            # 1. 读代码，填空，运行
python check.py               # 2. 看哪一项没过，按提示修
                              # 3. 反复 1-2 直到全绿
```

**批改脚本全绿 ≠ 学会了。** 真正的交付物是 README 里的**思考题回答**。
编程只是手段，回答问题才是在训练判断力。

---

## 三条纪律

1. **不要先看答案。** 每个练习的 `template.py` 里已经把步骤拆碎了，
   卡住时先看 `README.md` 的「常见坑」，还不行就问我。
2. **报错要读完整。** Python 的报错**从最后一行往上读**，最后一行才是真正的错误。
3. **卡超过 40 分钟就停下来问。** 死磕不是勤奋，是效率低。
   把代码和完整报错贴给我。

---

## 两个练习必须掌握的概念

### 第 1 周
- 真实数据的脏法：单位不一致、删失数据（`>`/`<`）、缺失、重复、格式混乱
- **为什么单位必须统一**（否则模型学的是单位，不是活性）
- 为什么用中位数聚合而不是平均值
- `pIC50 = 9 - log10(IC50[nM])`

### 第 2 周
- SMILES 是什么，`MolFromSmiles` 会返回 `None` 必须处理
- 分子描述符的**物理意义**（logP 管渗透性、TPSA 管极性、FractionCSP3 管三维性）
- Lipinski 五规则是**经验统计**不是硬门槛
- RDKit 绘图不支持中文（用英文图例）

---

## 环境问题

如果遇到 `DLL load failed` 之类的报错，先看
[项目 README 的踩坑记录](../projects/qsar-egfr/README.md)——
里面记录了这个环境上真实发生过的 6 个坑，包括一个查了很久的 DLL 冲突。
