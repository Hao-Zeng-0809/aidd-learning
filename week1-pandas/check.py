"""第 1 周练习自检脚本

用法：
    python check.py

它会逐条检查你的 clean_activities.csv，告诉你哪一步没做对。
**它不会给你答案**，只告诉你哪里不对、可能是什么原因。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = Path(__file__).resolve().parent
SUBMISSION = BASE / "clean_activities.csv"

# ── 标准答案的关键指标（由参考解得出）──────────────────────────
EXPECTED = {
    "n_raw": 108,            # 原始行数
    "n_after_dedup": 50,     # 去重后的有效记录数
    "n_molecules": 31,       # 最终分子数
    "mean_pIC50": 6.9071,
    "median_pIC50": 6.7950,
    "min_pIC50": 4.6300,
    "max_pIC50": 8.9250,
    "sum_n_records": 50,
}

# 抽检分子：(molecule_id, 期望 pIC50, 期望 n_records)
SPOT = [
    ("CHEMBL102473", 8.9250, 2),
    ("CHEMBL103569", 8.9000, 1),
    ("CHEMBL102884", 8.6100, 1),
    ("CHEMBL101240", 8.4000, 1),
    ("CHEMBL101651", 6.7950, 2),
    ("CHEMBL104528", 6.6300, 3),
    ("CHEMBL100829", 6.5700, 3),
    ("CHEMBL104802", 6.5200, 1),
    ("CHEMBL103706", 5.0800, 1),
    ("CHEMBL103432", 4.8900, 1),
    ("CHEMBL100144", 4.6300, 1),
]

REQUIRED_COLS = ["molecule_id", "smiles", "pIC50", "n_records"]

results: list[tuple[bool, str, str]] = []


def check(cond: bool, title: str, hint: str = "") -> bool:
    results.append((bool(cond), title, "" if cond else hint))
    return bool(cond)


def main():
    print("=" * 66)
    print("  第 1 周练习 · 自检")
    print("=" * 66)

    if not SUBMISSION.exists():
        print(f"\n❌ 没找到 {SUBMISSION.name}")
        print("   先运行 python template.py 生成它。")
        sys.exit(1)

    try:
        df = pd.read_csv(SUBMISSION)
    except Exception as e:
        print(f"\n❌ 读不了这个文件：{e}")
        sys.exit(1)

    if df.empty:
        print("\n❌ 文件是空的。")
        sys.exit(1)

    print(f"\n你的文件：{len(df)} 行，{len(df.columns)} 列")
    print(f"列名：{list(df.columns)}\n")

    # ── 1. 列名 ──────────────────────────────────────────────
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    extra = [c for c in df.columns if c not in REQUIRED_COLS]
    check(
        not missing and not extra,
        "输出列名正确",
        f"缺少 {missing}；多出 {extra}。要求恰好这 4 列：{REQUIRED_COLS}。"
        "想一想：中间过程的列（value、unit、value_nM）应该留在最终输出里吗？",
    )

    if missing:
        print("\n列名都不对，后面的检查没法继续。先修好列名再跑。")
        report()
        sys.exit(1)

    # ── 1b. 列顺序 ───────────────────────────────────────────
    check(
        list(df.columns) == list(REQUIRED_COLS),
        "输出列顺序正确",
        f"你的顺序是 {list(df.columns)}，要求是 {list(REQUIRED_COLS)}。"
        "用 df[REQUIRED_COLS] 显式指定列顺序，不要依赖中间过程留下的顺序。",
    )

    # ── 2. 缺失值 ────────────────────────────────────────────
    n_na = int(df[REQUIRED_COLS].isna().sum().sum())
    check(
        n_na == 0,
        "没有缺失值",
        f"有 {n_na} 个空单元格。检查是否漏了 dropna 或过滤条件。",
    )

    # ── 3. 类型 ──────────────────────────────────────────────
    check(
        pd.api.types.is_numeric_dtype(df["pIC50"]),
        "pIC50 是数值类型",
        "pIC50 现在不是数字（可能是字符串）。检查 pd.to_numeric 那一步。",
    )
    if pd.api.types.is_numeric_dtype(df["pIC50"]):
        pass

    # ── 4. 分子数 ────────────────────────────────────────────
    n = len(df)
    ok_n = n == EXPECTED["n_molecules"]
    if n > EXPECTED["n_molecules"]:
        hint = (f"你得到 {n} 个，期望 {EXPECTED['n_molecules']} 个 —— **多了**。"
                "最可能的原因：忘了删完全重复的行（规则 11），"
                "或者没过滤掉删失数据（`>` / `<`）、非数字值、≤0 的值。")
    elif n < EXPECTED["n_molecules"]:
        hint = (f"你得到 {n} 个，期望 {EXPECTED['n_molecules']} 个 —— **少了**。"
                "可能是过滤太狠：检查单位归一化是否漏了 `µM`/`um` 这类写法，"
                "或者只保留了 `IC50` 却把 `ic50`/` IC50` 也误删了。")
    else:
        hint = ""
    check(ok_n, f"分子数 = {EXPECTED['n_molecules']}", hint)

    # ── 5. 记录总数 ──────────────────────────────────────────
    if pd.api.types.is_numeric_dtype(df["n_records"]):
        total = int(df["n_records"].sum())
        check(
            total == EXPECTED["sum_n_records"],
            f"n_records 总和 = {EXPECTED['sum_n_records']}",
            f"你得到 {total}，期望 {EXPECTED['sum_n_records']}。"
            "这个数字等于去重后的有效记录总数 —— 对不上说明去重时机或过滤条件有问题。",
        )

    # ── 6. pIC50 统计量 ─────────────────────────────────────
    if pd.api.types.is_numeric_dtype(df["pIC50"]):
        for key, label, tol in [
            ("mean_pIC50", "均值", 0.02),
            ("median_pIC50", "中位数", 0.02),
        ]:
            got = float(df["pIC50"].mean() if key == "mean_pIC50" else df["pIC50"].median())
            check(
                abs(got - EXPECTED[key]) < tol,
                f"pIC50 {label} ≈ {EXPECTED[key]:.4f}",
                f"你得到 {got:.4f}，期望 {EXPECTED[key]:.4f}。"
                "偏差不大时通常是漏了几条记录；偏差很大时检查单位换算"
                "（uM 是不是忘了乘 1000？pIC50 公式是不是 9 - log10？）。",
            )

        lo, hi = float(df["pIC50"].min()), float(df["pIC50"].max())
        check(
            abs(lo - EXPECTED["min_pIC50"]) < 0.02 and abs(hi - EXPECTED["max_pIC50"]) < 0.02,
            f"pIC50 范围 = [{EXPECTED['min_pIC50']:.2f}, {EXPECTED['max_pIC50']:.2f}]",
            f"你得到 [{lo:.2f}, {hi:.2f}]。"
            "范围偏大 → 可能有没过滤掉的异常值（负值/零/删失数据）。",
        )

    # ── 7. 排序 ──────────────────────────────────────────────
    if pd.api.types.is_numeric_dtype(df["pIC50"]):
        check(
            df["pIC50"].is_monotonic_decreasing,
            "按 pIC50 降序排列",
            "README 要求降序。用 sort_values('pIC50', ascending=False)。",
        )

    # ── 8. smiles 有效性（加分项）───────────────────────────
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        bad = [s for s in df["smiles"].astype(str) if Chem.MolFromSmiles(s) is None]
        check(
            len(bad) == 0,
            "所有 SMILES 都能被 RDKit 解析",
            f"有 {len(bad)} 个 SMILES 解析失败，例如 {bad[:2]}。"
            "提示：strip 的时候有没有把结构里的空格也处理干净？",
        )
    except ImportError:
        pass

    # ── 9. 抽检 ──────────────────────────────────────────────
    print("抽检具体分子：")
    sub = df.set_index("molecule_id")
    spot_pass = 0
    for mol_id, exp_p, exp_n in SPOT:
        if mol_id not in sub.index:
            print(f"  ❌ {mol_id:<12} 不在你的结果里（可能是被误过滤了）")
            continue
        row = sub.loc[mol_id]
        got_p = float(row["pIC50"])
        got_n = int(row["n_records"]) if pd.api.types.is_numeric_dtype(df["n_records"]) else -1
        ok = abs(got_p - exp_p) < 0.02 and (got_n == exp_n or got_n == -1)
        if ok:
            spot_pass += 1
            print(f"  ✅ {mol_id:<12} pIC50={got_p:.4f}  n_records={got_n}")
        else:
            print(f"  ❌ {mol_id:<12} 你得到 pIC50={got_p:.4f} n={got_n}，"
                  f"期望 pIC50={exp_p:.4f} n={exp_n}")

    check(
        spot_pass == len(SPOT),
        f"抽检 {len(SPOT)} 个分子全部正确",
        f"通过 {spot_pass}/{len(SPOT)}。pIC50 对但 n_records 不对 → 看去重那一步；"
        "两者都不对 → 看单位换算和过滤条件。",
    )

    report()


def report():
    print("\n" + "=" * 66)
    passed = sum(1 for ok, _, _ in results if ok)
    total = len(results)
    print(f"  结果：{passed} / {total} 通过")
    print("=" * 66)

    failed = [(t, h) for ok, t, h in results if not ok]
    if not failed:
        print("\n🎉 全部通过！")
        print("\n下一步：")
        print("  1. 回答 README 里的 5 个思考题，写进你的仓库 README")
        print("  2. 把这份代码 push 到 GitHub（这是你的第一个作品）")
        print("  3. 去把练习数据删了重跑一遍，确认你能独立写出来")
        return

    print(f"\n还有 {len(failed)} 项没过：\n")
    for i, (title, hint) in enumerate(failed, 1):
        print(f"  {i}. ❌ {title}")
        if hint:
            print(f"     💡 {hint}")
        print()

    print("建议顺序：从上往下修，修一个跑一次 check.py。")
    print("卡超过 40 分钟就带着你的代码和报错来问我。")


if __name__ == "__main__":
    main()
