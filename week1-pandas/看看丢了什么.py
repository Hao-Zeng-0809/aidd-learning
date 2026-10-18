"""看清每一步到底丢掉了什么 —— 第 4 题的工具

════════════════════════════════════════════════════════════════
**这个脚本不产出 clean_activities.csv，只做诊断。**

它按 12 条规则一步步走，在每一步告诉你：
    · 丢掉了几行
    · 丢掉的**具体是什么值**（抽 3 个例子给你看）

跑完就能直接写出第 4 题的答案表。

用法：
    conda activate chemai
    cd exercises/week1-pandas
    python 看看丢了什么.py
════════════════════════════════════════════════════════════════
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = Path(__file__).resolve().parent
INPUT = BASE / "data" / "messy_activities.csv"

TEXT_COLS = ["molecule_id", "smiles", "standard_type", "standard_units", "standard_value"]

# 记录每一步的结果，最后打印成表
LOG: list[dict] = []


def show(step: str, before: pd.DataFrame, after: pd.DataFrame, col: str, why: str):
    """记录并打印一步的丢弃情况。

    step  ：规则编号和说明
    before/after：这一步前后
    col   ：丢掉的判断依据是看哪一列（用来取例子）
    why   ：为什么丢
    """
    dropped = len(before) - len(after)
    examples = ""
    if dropped > 0:
        # 找出被丢掉的那些行：在 before 里，但不在 after 里
        # 用 index 判断：after 的 index 就是保留下来的
        lost = before.loc[~before.index.isin(after.index)]
        vals = lost[col].dropna().astype(str).unique()[:3]
        if len(vals) == 0:
            # 全部是缺失值 —— 说不出来具体值，但要说清楚
            examples = f"（这 {dropped} 行的 {col} 列本身就是空的，没有值可显示）"
        else:
            examples = " | ".join(repr(v) for v in vals)
    LOG.append({"步骤": step, "丢了几行": dropped, "丢掉的值（例子）": examples, "原因": why})
    print(f"\n{'─' * 72}")
    print(f"  {step}")
    print(f"{'─' * 72}")
    print(f"  {len(before)} 行  →  {len(after)} 行   丢弃 {dropped} 行")
    if dropped:
        print(f"  丢掉的 {col} 例子：{examples}")
    else:
        print(f"  ⚠ 这一步没有丢任何行 —— 想一想为什么")


def main():
    df = pd.read_csv(INPUT)
    print(f"\n读入 {len(df)} 行")
    print("=" * 72)
    print("  逐步追踪：每一步丢掉了什么")
    print("=" * 72)

    # ── 规则 1：删序号列（丢的是列，不是行）──────────────────
    before = df.copy()
    df = df.drop(columns=["序号"])
    print(f"\n{'─' * 72}")
    print("  规则 1：删掉[序号]列")
    print(f"{'─' * 72}")
    print(f"  {len(before)} 行  →  {len(df)} 行   丢弃 0 行")
    print(f"  注意：这一步丢的是**列**（7 列 → 6 列），不是行")
    LOG.append({"步骤": "规则 1 删序号列", "丢了几行": 0, "丢掉的值（例子）": "丢的是列不是行", "原因": "导出时的行号，无意义"})

    # ── 规则 2：去空格 ────────────────────────────────────────
    before = df.copy()
    for c in TEXT_COLS:
        df[c] = df[c].astype("string").str.strip()
    print(f"\n{'─' * 72}")
    print("  规则 2：文本列去首尾空格")
    print(f"{'─' * 72}")
    print(f"  {len(before)} 行  →  {len(df)} 行   丢弃 0 行")
    print(f"  ⚠ 这一步也不丢行 —— 但它让**下一步**的重复行暴露出来")
    print(f"     （原始数据 df.duplicated() = {before.duplicated().sum()}，去空格后 = {df.duplicated().sum()}）")
    LOG.append({"步骤": "规则 2 去空格", "丢了几行": 0,
                "丢掉的值（例子）": "不丢行；但让 22 行重复记录现形",
                "原因": "格式统一"})

    # ── 规则 3：丢缺失 ────────────────────────────────────────
    before = df.copy()
    df = df[(df["molecule_id"].notna() & (df["molecule_id"] != "")) &
            (df["smiles"].notna() & (df["smiles"] != ""))]
    show("规则 3：丢弃 molecule_id 或 smiles 缺失的行", before, df, "smiles", "没有分子结构，无法建模")

    # ── 规则 4：只留 IC50 ─────────────────────────────────────
    before = df.copy()
    df = df[df["standard_type"].str.upper() == "IC50"]
    show("规则 4：只保留 standard_type == IC50", before, df, "standard_type",
         "活性类型不同（Ki/EC50/Kd），不能混用")

    # ── 规则 5：归一单位 ──────────────────────────────────────
    before = df.copy()

    def normalize_unit(u):
        if pd.isna(u):
            return None
        u = u.strip().replace("µ", "u").replace("μ", "u").lower()
        return {"nm": "nM", "um": "uM", "mm": "mM"}.get(u, None)

    df = df.copy()
    df["unit"] = df["standard_units"].map(normalize_unit)
    df = df.dropna(subset=["unit"])
    # 注意：这一步"可能"一行都不丢 —— 原因写死会自相矛盾，所以按实际情况写
    n_lost = len(before) - len(df)
    why_unit = ("这份数据里所有单位写法都能识别 → 一行都没丢（**这是加分点，要在答案里说明**）"
                if n_lost == 0 else "单位写法认不出，无法换算成 nM")
    show("规则 5：归一化单位（丢掉认不出的）", before, df, "standard_units", why_unit)

    # ── 规则 6：丢删失数据 ────────────────────────────────────
    before = df.copy()
    df = df[~df["standard_value"].str.startswith((">", "<"), na=False)]
    show("规则 6：丢弃删失数据（以 > 或 < 开头）", before, df, "standard_value",
         "只知道大于/小于，真值不确定")

    # ── 规则 7+8：提取数字 + 丢非法 ───────────────────────────
    before = df.copy()
    df = df.copy()
    df["value"] = pd.to_numeric(
        df["standard_value"].str.extract(r"^([0-9]*\.?[0-9]+)", expand=False),
        errors="coerce")

    # 分成两类看清楚
    print(f"\n{'─' * 72}")
    print("  规则 7：从 standard_value 提取数字")
    print(f"{'─' * 72}")
    na_after_extract = df["value"].isna().sum()
    print(f"  {len(df)} 行中，提取不出数字的：{na_after_extract} 行")
    print(f"  提取不出的值：")
    bad_vals = df.loc[df["value"].isna(), "standard_value"].astype(str).unique()[:6]
    for v in bad_vals:
        print(f"      {v!r}")
    print(f"\n  为什么提取不出？")
    print(f"      'n.d.' / 空值  → 本来就不是数字")
    print(f"      '-0.0213796'   → 正则开头的 ^ 要求从第一个字符就是数字，负号挡路")

    before2 = df.copy()
    df = df[df["value"].notna() & (df["value"] > 0)]
    show("规则 8：丢掉无法转数字的 + 数值 ≤ 0 的", before2, df, "standard_value",
         "n.d./空值/负数（提不出） + 零（提出来的）")

    # ── 规则 9+10：换算 + pIC50 ───────────────────────────────
    before = df.copy()
    df = df.copy()
    factor = {"nM": 1.0, "uM": 1e3, "mM": 1e6}
    df["value_nM"] = df["unit"].map(factor) * df["value"]
    df["pIC50"] = 9 - np.log10(df["value_nM"])
    show("规则 9+10：换算成 nM，算 pIC50", before, df, "standard_value", "—")

    # ── 规则 11：去重 ─────────────────────────────────────────
    before = df.copy()
    df = df.drop_duplicates()
    show("规则 11：删除完全重复的行", before, df, "standard_value",
         "同一分子同一条件的重复记录")

    # ── 规则 12：聚合（不是丢弃！）────────────────────────────
    print(f"\n{'─' * 72}")
    print("  规则 12：按 molecule_id 聚合")
    print(f"{'─' * 72}")
    print(f"  {len(df)} 条测定记录  →  按分子合并")
    out = (df.groupby("molecule_id", as_index=False)
             .agg(smiles=("smiles", "first"), pIC50=("pIC50", "median"),
                  n_records=("pIC50", "size")))
    n_multi = (out["n_records"] >= 2).sum()
    print(f"  最终 {len(out)} 个分子   丢弃 0 行")
    print(f"  ⚠️ 这不是'丢弃'，是'合并'——信息没丢，n_records 记录了每个分子有几条测定")
    print(f"     其中 {len(out) - n_multi} 个分子只有 1 条记录")
    print(f"          {n_multi} 个分子有 2 条以上记录（共 {out.loc[out['n_records'] >= 2, 'n_records'].sum()} 条）")
    LOG.append({"步骤": "规则 12 聚合", "丢了几行": 0,
                "丢掉的值（例子）": f"{len(df)} 条记录 → {len(out)} 个分子（合并非丢弃）",
                "原因": "同一分子的多次测定取中位数"})

    # ── 汇总表 ────────────────────────────────────────────────
    print(f"\n\n{'=' * 72}")
    print("  【第 4 题答案表】直接抄进你的答案里")
    print("=" * 72)
    log = pd.DataFrame(LOG)
    print(f"\n  {'步骤':<32}{'丢了几行':>8}   原因")
    print(f"  {'─' * 68}")
    for _, r in log.iterrows():
        print(f"  {r['步骤']:<32}{r['丢了几行']:>8}   {r['原因']}")
    print(f"\n  合计：{len(pd.read_csv(INPUT))} 行 → {len(out)} 个分子")

    print(f"\n  每一步丢掉的**具体值**（写进答案，不要只写'缺失值'）：")
    print(f"  {'─' * 68}")
    for _, r in log.iterrows():
        if r['丢掉的值（例子）']:
            print(f"  · {r['步骤']}")
            print(f"      {r['丢掉的值（例子）']}")

    print(f"\n{'=' * 72}")
    print("  💡 现在你能自己回答第 4 题了：")
    print("     ① 每步丢了几行（上面表里有）")
    print("     ② 丢的是什么（上面有真实例子）")
    print("     ③ 哪一步'没丢'以及为什么（规则 2、5、12）")
    print("=" * 72)
    print()


if __name__ == "__main__":
    main()
