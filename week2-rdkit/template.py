"""第 2 周练习：用 RDKit 计算药物分子的理化性质

目标：
  1. drugs_properties.csv —— 25 个药物的理化性质 + Lipinski 五规则
  2. drugs_grid.png       —— 分子结构网格图

运行：
    python template.py
    python check.py        # 自检

════════════════════════════════════════════════════════════════
⚠️ 图例请用英文名。RDKit 自带字体不支持中文，中文会被静默丢弃。
════════════════════════════════════════════════════════════════
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Draw

# Windows 控制台默认 GBK，打印特殊字符会报错
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# RDKit 会为 SMILES 里的小问题刷警告，练习时先关掉，保持输出干净
RDLogger.DisableLog("rdApp.*")

BASE = Path(__file__).resolve().parent
INPUT = BASE / "data" / "drugs.csv"
OUT_CSV = BASE / "drugs_properties.csv"
OUT_PNG = BASE / "drugs_grid.png"

# 输出列的固定顺序
OUTPUT_COLS = [
    "name_cn", "name_en", "smiles", "drug_class",
    "MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
    "NumRotatableBonds", "RingCount", "FractionCSP3", "QED",
    "violations", "passes_ro5",
]


def compute_properties(mol) -> dict:
    """算一个分子的 9 个理化性质，返回 dict。"""
    # TODO: 用 Descriptors 里的函数把 9 个量算出来
    # 提示：Descriptors.MolWt(mol)、Descriptors.MolLogP(mol)、Descriptors.TPSA(mol) ...
    #       ⚠️ QED 是小写的 Descriptors.qed(mol)，不是 Descriptors.QED
    raise NotImplementedError


def count_violations(mol, props: dict) -> int:
    """计算 Lipinski 五规则违反了几条。

    规则：MW<=500, logP<=5, HBD<=5, HBA<=10
    超过阈值的算一条违反。
    """
    # TODO: 四个 if，累加违反数
    # 提示：注意是「大于」才算违反
    raise NotImplementedError


def main():
    df = pd.read_csv(INPUT)
    print(f"读入 {len(df)} 个药物")
    print(f"  类别分布：{df['drug_class'].nunique()} 类")

    records = []
    mols = []          # 画图用
    legends = []       # 画图用
    n_failed = 0

    for _, row in df.iterrows():
        # ---- 1. SMILES -> Mol ----
        mol = Chem.MolFromSmiles(row["smiles"])

        # ---- 2. 必须处理解析失败 ----
        # TODO: 如果 mol is None，n_failed += 1 并 continue
        # if mol is None:
        #     ...

        # ---- 3. 算性质 ----
        # TODO: props = compute_properties(mol)

        # ---- 4. 算五规则违反数 ----
        # TODO: viol = count_violations(mol, props)

        # ---- 5. 组装这一行 ----
        # TODO:
        # records.append({...原来 4 列..., **props,
        #                 "violations": viol,
        #                 "passes_ro5": viol <= 1})
        # mols.append(mol)
        # legends.append(row["name_en"])

    # ---- 6. 输出 CSV ----
    # out = pd.DataFrame(records)[OUTPUT_COLS]
    # out.to_csv(OUT_CSV, index=False)
    # print(f"\n✅ 已保存 {len(out)} 行 -> {OUT_CSV}")

    # ---- 7. 画分子网格图 ----
    # TODO: 4 列布局，每格 330x260
    # img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(330, 260), legends=legends)
    # img.save(OUT_PNG)
    # print(f"✅ 已保存分子图 -> {OUT_PNG}")

    # ---- 8. 打印小结（留着，方便你判断是否合理）----
    # print(f"\n解析失败 {n_failed} 个")
    # print(f"Lipinski 全通过：{int(out['passes_ro5'].sum())} / {len(out)}")
    # print(f"\n违反 >= 2 条的分子：")
    # bad = out[out['violations'] >= 2][['name_cn', 'name_en', 'violations', 'passes_ro5']]
    # print(bad.to_string(index=False) if len(bad) else "   （没有）")
    # print(f"\nlogP 最高的 3 个：")
    # print(out.nlargest(3, 'MolLogP')[['name_cn', 'MolLogP']].to_string(index=False))

    print("\n⚠️ 你还没写完 template.py 里的 TODO。")


if __name__ == "__main__":
    main()
