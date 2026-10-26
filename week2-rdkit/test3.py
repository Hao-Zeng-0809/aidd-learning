"""第 2 周 · 最小可运行版本：只算 3 个描述符

目的：先确认 RDKit 能算出数，再扩展到 15 列。
跑通这个，第 2 周就过了 80%。

用法：
    conda activate chemai
    cd exercises/week2-rdkit
    python test3.py
"""
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = Path(__file__).resolve().parent

df = pd.read_csv(BASE / "data" / "drugs.csv")
print(f"读入 {len(df)} 个药物\n")

# 三个最基础的描述符，写法完全一样，只是换函数名
df["MolWt"] = df["smiles"].map(lambda s: Descriptors.MolWt(Chem.MolFromSmiles(s)))
df["MolLogP"] = df["smiles"].map(lambda s: Descriptors.MolLogP(Chem.MolFromSmiles(s)))
df["TPSA"] = df["smiles"].map(lambda s: Descriptors.TPSA(Chem.MolFromSmiles(s)))

out = df[["name_cn", "MolWt", "MolLogP", "TPSA"]].copy()
out["MolWt"] = out["MolWt"].round(2)
out["MolLogP"] = out["MolLogP"].round(2)
out["TPSA"] = out["TPSA"].round(2)

pd.set_option("display.unicode.east_asian_width", True)
print(out.head(10).to_string(index=False))

print(f"\n✅ 3 个描述符算完了。剩下 12 个是同样的写法，只是换 Descriptors.XXX")
print(f"   下一步：打开 template.py，把它扩成 15 列")
