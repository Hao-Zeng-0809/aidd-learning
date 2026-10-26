"""第 2 周交付物：25 个药物的 15 列性质表

════════════════════════════════════════════════════════════════
**用法（复制粘贴）：**
    conda activate chemai
    cd exercises/week2-rdkit
    python compute_all.py

运行后会生成 `drugs_properties.csv`（15 列 × 25 行）。

**任务：跑完之后，逐列读一遍下面的代码，确保你能说出每一列在算什么。**
明天我会抽问你其中任何一列。
════════════════════════════════════════════════════════════════
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")      # 关掉 RDKit 的啰嗦警告

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = Path(__file__).resolve().parent
INPUT = BASE / "data" / "drugs.csv"
OUTPUT = BASE / "drugs_properties.csv"

# ① 需要计算的 9 个理化性质：列名 → rdkit 里的函数
PROPS = {
    "MolWt":             Descriptors.MolWt,             # 分子量
    "MolLogP":           Descriptors.MolLogP,           # 脂水分配系数
    "TPSA":              Descriptors.TPSA,              # 拓扑极性表面积
    "NumHDonors":        Descriptors.NumHDonors,        # 氢键供体数
    "NumHAcceptors":     Descriptors.NumHAcceptors,     # 氢键受体数
    "NumRotatableBonds": Descriptors.NumRotatableBonds, # 可旋转键数
    "RingCount":         Descriptors.RingCount,         # 环数
    "FractionCSP3":      Descriptors.FractionCSP3,      # sp3 碳比例
    "QED":               Descriptors.qed,               # 类药性综合评分（注意小写）
}

# ② Lipinski 五规则的四个阈值（违反即 +1）
RO5 = {
    "MolWt":         (500, ">"),   # 分子量 > 500 违反
    "MolLogP":       (5,   ">"),   # logP > 5 违反
    "NumHDonors":    (5,   ">"),   # 氢键供体 > 5 违反
    "NumHAcceptors": (10,  ">"),   # 氢键受体 > 10 违反
}


def main():
    df = pd.read_csv(INPUT)
    print(f"读入 {len(df)} 个药物")

    # ── 逐分子算描述符 ────────────────────────────────────
    rows = []
    for _, r in df.iterrows():
        mol = Chem.MolFromSmiles(r["smiles"])
        if mol is None:
            print(f"  ⚠ SMILES 解析失败，跳过：{r['name_cn']}")
            continue
        rec = {"name_cn": r["name_cn"], "name_en": r["name_en"],
               "smiles": r["smiles"], "drug_class": r["drug_class"]}
        for col, fn in PROPS.items():
            rec[col] = fn(mol)
        rows.append(rec)

    out = pd.DataFrame(rows)
    print(f"  成功算出 {len(out)} 个分子的性质")

    # ── Lipinski 五规则 ───────────────────────────────────
    # violations = 违反了几条（0-4 的整数）
    violations = np.zeros(len(out), dtype=int)
    for col, (thr, op) in RO5.items():
        if op == ">":
            violations += (out[col] > thr).astype(int)
    out["violations"] = violations
    # passes_ro5：违反 ≤ 1 条通常仍可接受
    out["passes_ro5"] = (violations <= 1)

    # ── 按 README 要求的列顺序整理并保存 ──────────────────
    COLS = ["name_cn", "name_en", "smiles", "drug_class",
            "MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
            "NumRotatableBonds", "RingCount", "FractionCSP3", "QED",
            "violations", "passes_ro5"]
    out = out[COLS]
    out.to_csv(OUTPUT, index=False)

    # ── 打印摘要 ─────────────────────────────────────────
    print(f"\n✅ 已保存 {OUTPUT.name}：{len(out)} 行 × {len(out.columns)} 列")
    print(f"\n五规则通过：{out['passes_ro5'].sum()}/{len(out)} 个")
    print(f"违反 0 条：{(out['violations'] == 0).sum()} 个")
    print(f"违反 1 条：{(out['violations'] == 1).sum()} 个")
    print(f"违反 2 条：{(out['violations'] == 2).sum()} 个")

    bad = out[out["violations"] >= 1][["name_cn", "violations", "MolWt", "MolLogP"]]
    if len(bad):
        print(f"\n违反规则最多的分子：")
        print(bad.sort_values("violations", ascending=False).head(8).to_string(index=False))


if __name__ == "__main__":
    main()
