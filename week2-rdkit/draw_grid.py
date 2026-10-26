"""第 2 周交付物 ②：25 个药物的分子结构网格图

════════════════════════════════════════════════════════════════
用法：
    conda activate chemai
    cd exercises/week2-rdkit
    python draw_grid.py

运行后生成 `drugs_grid.png`（4 列 × 7 行，25 个分子）。

⚠️ 图例必须用**英文名**。RDKit 自带字体没有中文字形，
   写中文会被**静默丢弃**（不报错，只是图上没有）。
════════════════════════════════════════════════════════════════
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw, AllChem

RDLogger.DisableLog("rdApp.*")

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = Path(__file__).resolve().parent
INPUT = BASE / "data" / "drugs.csv"
OUTPUT = BASE / "drugs_grid.png"

MOLS_PER_ROW = 4          # 4 列 × 7 行 = 28 格，装 25 个分子
SUB_IMG_SIZE = (320, 240) # 每个格子的像素尺寸


def main():
    df = pd.read_csv(INPUT)
    print(f"读入 {len(df)} 个药物")

    mols, legends = [], []
    skipped = []
    for _, r in df.iterrows():
        mol = Chem.MolFromSmiles(r["smiles"])
        if mol is None:
            skipped.append(r["name_cn"])
            continue
        # 统一重绘成 2D 坐标，保证图上结构摆正
        AllChem.Compute2DCoords(mol)
        mols.append(mol)
        # ⚠️ 只用英文名，中文会被静默丢弃
        legends.append(str(r["name_en"]))

    print(f"  成功解析 {len(mols)} 个分子")
    if skipped:
        print(f"  ⚠ 解析失败跳过：{skipped}")

    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=MOLS_PER_ROW,
        subImgSize=SUB_IMG_SIZE,
        legends=legends[:len(mols)],
    )
    img.save(str(OUTPUT))

    print(f"\n✅ 已保存 {OUTPUT.name}")
    print(f"   尺寸：{img.size[0]} × {img.size[1]} 像素")
    print(f"   布局：{MOLS_PER_ROW} 列 × {len(mols)//MOLS_PER_ROW + (1 if len(mols)%MOLS_PER_ROW else 0)} 行")
    print(f"\n   打开看看：start {OUTPUT}")


if __name__ == "__main__":
    main()
