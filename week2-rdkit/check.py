"""第 2 周练习自检脚本

用法：
    python check.py

用 RDKit 重新独立计算一遍标准答案，逐列和你的结果对比。
只告诉你哪里不对，不给完整答案。
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = Path(__file__).resolve().parent
INPUT = BASE / "data" / "drugs.csv"
SUB_CSV = BASE / "drugs_properties.csv"
SUB_PNG = BASE / "drugs_grid.png"

EXPECTED_COLS = [
    "name_cn", "name_en", "smiles", "drug_class",
    "MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
    "NumRotatableBonds", "RingCount", "FractionCSP3", "QED",
    "violations", "passes_ro5",
]

# 需要数值对比的列（容差）
NUM_COLS = {
    "MolWt": 0.05, "MolLogP": 0.05, "TPSA": 0.05,
    "NumHDonors": 0, "NumHAcceptors": 0, "NumRotatableBonds": 0,
    "RingCount": 0, "FractionCSP3": 1e-6, "QED": 1e-6,
    "violations": 0,
}

results: list[tuple[bool, str, str]] = []


def check(cond, title, hint=""):
    results.append((bool(cond), title, "" if cond else hint))
    return bool(cond)


def build_reference() -> pd.DataFrame:
    """用 RDKit 独立算一遍标准答案。"""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors
    RDLogger.DisableLog("rdApp.*")

    src = pd.read_csv(INPUT)
    rows = []
    for _, r in src.iterrows():
        mol = Chem.MolFromSmiles(r["smiles"])
        if mol is None:
            continue
        props = {
            "MolWt": Descriptors.MolWt(mol),
            "MolLogP": Descriptors.MolLogP(mol),
            "TPSA": Descriptors.TPSA(mol),
            "NumHDonors": Descriptors.NumHDonors(mol),
            "NumHAcceptors": Descriptors.NumHAcceptors(mol),
            "NumRotatableBonds": Descriptors.NumRotatableBonds(mol),
            "RingCount": Descriptors.RingCount(mol),
            "FractionCSP3": Descriptors.FractionCSP3(mol),
            "QED": Descriptors.qed(mol),
        }
        v = 0
        if props["MolWt"] > 500:
            v += 1
        if props["MolLogP"] > 5:
            v += 1
        if props["NumHDonors"] > 5:
            v += 1
        if props["NumHAcceptors"] > 10:
            v += 1
        rows.append({**{k: r[k] for k in ["name_cn", "name_en", "smiles", "drug_class"]},
                     **props, "violations": v, "passes_ro5": v <= 1})
    return pd.DataFrame(rows)


def main():
    print("=" * 68)
    print("  第 2 周练习 · 自检")
    print("=" * 68)

    # ── CSV ────────────────────────────────────────────────
    if not SUB_CSV.exists():
        print(f"\n❌ 没找到 {SUB_CSV.name}，先运行 python template.py")
        sys.exit(1)
    try:
        sub = pd.read_csv(SUB_CSV)
    except Exception as e:
        print(f"\n❌ 读不了 {SUB_CSV.name}：{e}")
        sys.exit(1)

    ref = build_reference()
    print(f"\n你的 CSV：{len(sub)} 行 × {len(sub.columns)} 列")
    print(f"标准答案：{len(ref)} 行 × {len(ref.columns)} 列\n")

    # 1. 列
    missing = [c for c in EXPECTED_COLS if c not in sub.columns]
    extra = [c for c in sub.columns if c not in EXPECTED_COLS]
    check(not missing and not extra, "列名完全正确（15 列）",
          f"缺少 {missing}；多出 {extra}。")

    # 1b. 列顺序（原 4 列在前，再按计算顺序接 11 个新列）
    check(list(sub.columns) == list(EXPECTED_COLS), "列顺序正确",
          f"你的顺序是 {list(sub.columns)}。\n"
          f"        要求是 {list(EXPECTED_COLS)}。\n"
          "        用 df[EXPECTED_COLS] 显式指定顺序，不要依赖中间过程留下的顺序。")

    # 2. 行数
    check(len(sub) == len(ref), f"行数 = {len(ref)}",
          f"你得到 {len(sub)} 行。少了 → 可能有 SMILES 解析失败被跳过了（但这 25 个都能解析）；"
          "多了 → 可能重复写入了。")

    if missing or len(sub) != len(ref):
        report(ref)
        return

    # 3. 对齐后逐列比
    sub = sub.sort_values("name_en").reset_index(drop=True)
    ref = ref.sort_values("name_en").reset_index(drop=True)

    if not sub["name_en"].equals(ref["name_en"]):
        check(False, "药物名称与顺序一致",
              "name_en 对不上。检查是否漏了分子或改了名字。")
    else:
        check(True, "25 个药物齐全")

        # 数值列逐列对比
        for col, tol in NUM_COLS.items():
            got = pd.to_numeric(sub[col], errors="coerce")
            if got.isna().any():
                bad = sub.loc[got.isna(), "name_en"].tolist()[:3]
                check(False, f"{col} 全部是数值",
                      f"有 {int(got.isna().sum())} 个非数值，例如 {bad}。")
                continue
            diff = (got - ref[col]).abs()
            n_bad = int((diff > tol + 1e-9).sum())
            if n_bad == 0:
                check(True, f"{col} 全部正确")
            else:
                worst = sub.loc[diff.idxmax(), "name_en"]
                check(False, f"{col} 全部正确",
                      f"{n_bad} 个对不上。偏差最大的是 {worst}"
                      f"（你 {got[diff.idxmax()]:.3f}，期望 {ref[col][diff.idxmax()]:.3f}）。"
                      + (" 检查 Descriptors 里用的函数名对不对。"
                         if col in ("MolWt", "MolLogP", "TPSA") else
                         " 检查是不是漏了某个分子或算错了规则。"))

    # 4. passes_ro5 与 violations 一致性
    if "violations" in sub.columns and "passes_ro5" in sub.columns:
        try:
            v = pd.to_numeric(sub["violations"])
            p = sub["passes_ro5"].astype(str).str.lower().isin(["true", "1", "yes"])
            check(bool((p == (v <= 1)).all()), "passes_ro5 与 violations 逻辑一致",
                  "passes_ro5 应该是 violations <= 1。注意违反 1 条通常仍算通过。")
        except Exception:
            check(False, "passes_ro5 可解析", "passes_ro5 应该是 True/False。")

    # 5. 图片
    if not SUB_PNG.exists():
        check(False, "生成了 drugs_grid.png", "还没有这个文件，检查画图那一步。")
    else:
        size = SUB_PNG.stat().st_size
        ok_png = size > 20_000
        dims = ""
        try:
            from PIL import Image
            with Image.open(SUB_PNG) as im:
                dims = f"{im.size[0]}×{im.size[1]}"
                ok_png = ok_png and im.size[0] > 600 and im.size[1] > 400
        except Exception:
            pass
        check(ok_png, "分子网格图看起来正常",
              f"文件只有 {size:,} 字节{'（尺寸 ' + dims + '）' if dims else ''}，"
              "可能没画成功或图太小。用 MolsToGridImage（有 s），molsPerRow=4。")

    report(ref)


def report(ref: pd.DataFrame):
    print("\n" + "=" * 68)
    passed = sum(1 for ok, _, _ in results if ok)
    print(f"  结果：{passed} / {len(results)} 通过")
    print("=" * 68)

    failed = [(t, h) for ok, t, h in results if not ok]
    if failed:
        print(f"\n还有 {len(failed)} 项没过：\n")
        for i, (t, h) in enumerate(failed, 1):
            print(f"  {i}. ❌ {t}")
            if h:
                print(f"     💡 {h}")
        print("\n修一个跑一次 check.py。")
    else:
        print("\n🎉 全部通过！")

    # ── 数据观察（不是答案，是让你做思考题用的素材）──
    print("\n" + "─" * 68)
    print("  从你自己的结果里能看到的（思考题素材）")
    print("─" * 68)

    if "violations" in ref.columns:
        bad = ref[ref["violations"] >= 2].sort_values("violations", ascending=False)
        print(f"\n违反 Lipinski ≥ 2 条的分子（共 {len(bad)} 个）：")
        if len(bad):
            print(bad[["name_cn", "name_en", "MolWt", "MolLogP",
                       "NumHDonors", "NumHAcceptors", "violations"]].to_string(index=False))
        else:
            print("   （没有）")

        one = ref[ref["violations"] == 1]
        print(f"\n违反 1 条的分子（共 {len(one)} 个，通常仍算通过）：")
        print("   " + "、".join(one["name_cn"].tolist()) if len(one) else "   （没有）")

    print(f"\nlogP 最高的 3 个：")
    print(ref.nlargest(3, "MolLogP")[["name_cn", "MolLogP", "drug_class"]].to_string(index=False))
    print(f"\nlogP 最低的 3 个：")
    print(ref.nsmallest(3, "MolLogP")[["name_cn", "MolLogP", "drug_class"]].to_string(index=False))
    print(f"\nFractionCSP3 最高的 3 个（最「立体」的分子）：")
    print(ref.nlargest(3, "FractionCSP3")[["name_cn", "FractionCSP3", "drug_class"]].to_string(index=False))

    print("\n💡 现在去回答 README 里的 5 个思考题。**那才是这个练习的重点。**")


if __name__ == "__main__":
    main()
