"""第 1 周练习：清洗药物活性数据

目标：把 data/messy_activities.csv 清洗成 clean_activities.csv

规则见 README.md。下面按 12 条规则拆成了 TODO，
你只需要把每个 TODO 补上，然后运行：

    python template.py
    python check.py      # 自检

════════════════════════════════════════════════════════════════
⚠️ 别急着看答案。卡住了先看 README 的「常见坑」，再问我。
════════════════════════════════════════════════════════════════
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Windows 控制台默认 GBK，打印特殊字符会报错，强制 UTF-8
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = Path(__file__).resolve().parent
INPUT = BASE / "data" / "messy_activities.csv"
OUTPUT = BASE / "clean_activities.csv"


def main():
    # ---- 读数据 ----
    df = pd.read_csv(INPUT)
    n_raw = len(df)
    print(f"读入 {n_raw} 行，列：{list(df.columns)}")

    # ---- 规则 1：删掉「序号」列 ----
    # TODO: 用 drop(columns=[...]) 删掉它
    df = df.drop(columns=['序号'])

    # ---- 规则 2：文本列去掉首尾空格 ----
    # 提示：用 .astype("string").str.strip()
    #       不要用 .astype(str)，那会把缺失值变成字符串 "nan"
    TEXT_COLS = ["molecule_id", "smiles", "standard_type", "standard_units", "standard_value"]
    # TODO: 对上面 5 列做 strip
    for c in TEXT_COLS:
        df[c] = df[c].astype('string').str.strip()

    # ---- 规则 3：丢弃 molecule_id 或 smiles 缺失的行 ----
    # 提示：df["col"].notna() 和 (df["col"] != "")
    # TODO:
    df = df[df['molecule_id'].notna() & df['smiles'].notna()]

    # ---- 规则 4：只保留 standard_type 是 IC50 的行（大小写不敏感）----
    # 提示：.str.upper() == "IC50"
    # TODO:
    df = df[df['standard_type'].str.upper() == 'IC50']

    # ---- 规则 5：归一化 standard_units ----
    # 允许的输入写法：nM/nm/NM/" nM "、uM/um/µM/μM ...
    # 注意 µ (U+00B5) 和 μ (U+03BC) 都是希腊字母 mu，要都替换成 u
    def normalize_unit(u):
        """把各种单位写法统一成 nM / uM / mM；认不出来返回 None。"""
        if pd.isna(u):
            return None
        # TODO:
        u = u.strip()
        u = u.replace("µ", "u").replace("μ", "u")
        u = u.lower() #之后查表
        unit_map = {"nm": "nM", "um": "uM", "mm": "mM"}
        return unit_map.get(u,None)

    df["unit"] = df["standard_units"].map(normalize_unit)
    # TODO: 丢掉 unit 为 None 的行
    df = df.dropna(subset=['unit'])

    # ---- 规则 6：丢弃删失数据（以 > 或 < 开头）----
    # 提示：.str.startswith((">", "<"), na=False)
    # TODO:
    df = df[~df['standard_value'].str.startswith(('<','>'),na=False)]

    # ---- 规则 7：从 standard_value 里提取数字部分 ----
    # "100.0 nM" -> 100.0      "n.d." -> 无法提取       "-478.6" -> 提取不到（负号在开头）
    # 提示：.str.extract(r"^([0-9]*\.?[0-9]+)", expand=False)
    #       然后 pd.to_numeric(..., errors="coerce")
    # TODO:
    df["value"] = pd.to_numeric(df['standard_value'].str.extract(r"^([0-9]*\.?[0-9]+)",expand=False),errors='coerce')
    

    # ---- 规则 8：丢弃无法转数字的行，以及 <= 0 的行 ----
    # TODO:
    df = df[df["value"].notna() & (df["value"] > 0)]

    # ---- 规则 9 + 10：换算成 nM，算 pIC50 ----
    # 换算系数：1 nM = 1, 1 uM = 1000 nM, 1 mM = 1e6 nM
    # pIC50 = 9 - log10(值_nM)
    # TODO:
    factor = {"nM": 1.0, "uM": 1e3, "mM": 1e6}
    df["value_nM"] = df['unit'].map(factor) * df['value']
    df["pIC50"] = 9 - np.log10(df['value_nM'])

    # ---- 规则 11：删除完全重复的行 ----
    # 提示：df.drop_duplicates()
    # ⚠️ 必须在聚合之前做
    # TODO:
    df = df.drop_duplicates()

    # ---- 规则 12：按 molecule_id 聚合 ----
    # 需要：smiles 取第一条、pIC50 取中位数、n_records 计数
    # 提示：
    #out = (df.groupby("molecule_id", as_index=False).agg(smiles=("smiles", "first"),pIC50=("pIC50", "median"),n_records=("pIC50", "size")))
    # TODO:
    # out = ...
    out = (df.groupby("molecule_id", as_index=False).agg(smiles=("smiles", "first"),pIC50=("pIC50", "median"),n_records=("pIC50", "size")))
      # ---- 排序并保存 ----
    out = out.sort_values("pIC50", ascending=False).reset_index(drop=True)
    out.to_csv(OUTPUT, index=False)
    print(f"\n✅ 已保存 {OUTPUT.name}")
    # ---- 打印小结（这几行留着，方便你判断是否合理）----
    print(f"\n最终分子数：{len(out)}")
    print(f"pIC50 范围：{out['pIC50'].min():.2f} ~ {out['pIC50'].max():.2f}")
    print(f"pIC50 均值：{out['pIC50'].mean():.2f}")
    print(f"n_records 总和：{out['n_records'].sum()}")
    print(out.head(5).to_string(index=False))

    # ---- 自检提示（这行自动判断你做完了没有）----
    if OUTPUT.exists():
        done = pd.read_csv(OUTPUT)
        print(f"\n✅ 已生成 {OUTPUT.name}：{len(done)} 行，列 {list(done.columns)}")
        print(done.head(5).to_string(index=False))
        print("\n👉 现在运行：  python check.py    看哪几项没过。")
    else:
        print("\n⚠️ 还没有生成 clean_activities.csv。")
        print("   继续填上面的 TODO，或者检查最后几行是不是还没取消注释。")
        print("   卡住了？运行 python check.py 也会有提示。")


if __name__ == "__main__":
    main()
