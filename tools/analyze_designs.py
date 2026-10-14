#!/usr/bin/env python
"""analyze_designs.py — D2-S1 汇总 + D4 短名单（03_TOOL_SPECS.md · T3）

产出：
  out/D2_S1_sequence_scores.csv   S1 内在打分（逐变体）
  out/D4_variant_shortlist.csv    10 条 + WT 行（含 WT_percentile）

WT 分位定义（诚实声明）：D1 fasta 的 ligand_confidence 只对采样序列存在，
WT 没有该字段；因此 WT 与 300 条变体在**同一打分流程**下的可比量是
S3 的 score_EGT（score.py，自回归、10 解码序、seed=111，13 个口袋位点平均 log-prob）。
WT_percentile = 变体中 score_EGT ≤ WT.score_EGT 的比例。
"""
import os, re, csv, sys
from collections import defaultdict

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SITES = [14, 50, 53, 107, 139, 219, 223, 231, 322, 326, 346, 350, 353]
SITE_ROLE = {14: "笼", 53: "笼", 111: "笼", 219: "笼", 326: "笼", 350: "笼",
             139: "极性衬里", 223: "极性衬里", 322: "极性衬里", 346: "极性衬里", 353: "极性衬里",
             50: "疏水衬里", 231: "疏水衬里", 107: "唯一负电(D)"}

# ---------- 1. 解析 D1 fasta ----------
variants = []
wt_seq = None
for T in ("0.1", "0.2", "0.3"):
    fa = os.path.join(PKG, f"out/D1_EGT_T{T}/seqs/mfsT1_EGT.fa")
    name, meta = None, {}
    for ln in open(fa):
        ln = ln.strip()
        if ln.startswith(">"):
            m = re.search(r"id=(\d+)", ln)
            name = int(m.group(1)) if m else 0
            if name:
                meta = {k: float(v) for k, v in
                        re.findall(r"(overall_confidence|ligand_confidence|seq_rec)=([\d.]+)", ln)}
        elif ln:
            if name == 0:
                wt_seq = ln
            else:
                variants.append(dict(
                    variant_id=f"T{T}/id{name}", temperature=T, seq=ln,
                    overall_confidence=meta["overall_confidence"],
                    ligand_confidence=meta["ligand_confidence"],
                    seq_rec=meta["seq_rec"]))
assert wt_seq, "fasta 首条应为 WT 序列"
for v in variants:
    v["pocket_fingerprint"] = "".join(v["seq"][i - 1] for i in SITES)
    v["mutations"] = [f"{wt_seq[i - 1]}{i}{v['seq'][i - 1]}" for i in SITES
                      if v["seq"][i - 1] != wt_seq[i - 1]]

# ---------- 2. D2_S1 ----------
with open(os.path.join(PKG, "out/D2_S1_sequence_scores.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["variant_id", "temperature", "overall_confidence", "ligand_confidence",
                "seq_rec", "pocket_fingerprint", "mutations"])
    for v in variants:
        w.writerow([v["variant_id"], v["temperature"], v["overall_confidence"],
                    v["ligand_confidence"], v["seq_rec"], v["pocket_fingerprint"],
                    ";".join(v["mutations"])])
print(f"[ok] D2_S1: {len(variants)} 条")

# ---------- 3. 合并 S2 / S3 ----------
s2 = {r["variant_id"]: r for r in
      csv.DictReader(open(os.path.join(PKG, "out/D2_S2_geometry_check.csv"), encoding="utf-8"))}
s3 = {r["variant_id"]: r for r in
      csv.DictReader(open(os.path.join(PKG, "out/D2_S3_selectivity_bias.csv"), encoding="utf-8"))}
for v in variants:
    g, s = s2[v["variant_id"]], s3[v["variant_id"]]
    v.update(n_pass=int(g["n_pass"]),
             score_EGT=float(s["score_EGT"]),
             dscore=float(s["dscore_vs_hercynine"]))
wt_s3 = s3["WT"]
wt_score = float(wt_s3["score_EGT"])
wt_dscore = float(wt_s3["dscore_vs_hercynine"])

# ---------- 4. WT 分位 ----------
n = len(variants)
wt_pct_score = 100.0 * sum(v["score_EGT"] <= wt_score for v in variants) / n
wt_pct_dscore = 100.0 * sum(v["dscore"] <= wt_dscore for v in variants) / n
wt_pct_ligconf = None  # fasta ligand_confidence 无 WT 值，见文件头声明
print(f"[ok] WT score_EGT={wt_score:.4f} → 分位 {wt_pct_score:.1f}% ；"
      f"WT Δscore={wt_dscore:.4f} → 分位 {wt_pct_dscore:.1f}%")

# ---------- 5. 综合排序 + 短名单 ----------
def z(xs):
    mu = sum(xs) / len(xs)
    sd = (sum((x - mu) ** 2 for x in xs) / len(xs)) ** 0.5 or 1.0
    return [(x - mu) / sd for x in xs]

zl = z([v["ligand_confidence"] for v in variants])
zd = z([v["dscore"] for v in variants])
for v, a, b in zip(variants, zl, zd):
    v["composite"] = a + b

eligible = [v for v in variants if v["n_pass"] == 3]
best_by_fp = {}
for v in sorted(eligible, key=lambda x: -x["composite"]):
    best_by_fp.setdefault(v["pocket_fingerprint"], v)
shortlist = sorted(best_by_fp.values(), key=lambda x: -x["composite"])[:10]
for rank, v in enumerate(shortlist, 1):
    v["rank"] = rank


def phenotype(v):
    if not v["mutations"]:
        return "口袋指纹 = WT；预期表型同 WT 基线"
    parts = []
    for m in v["mutations"]:
        pos = int(re.search(r"\d+", m).group())
        parts.append(f"{m}({SITE_ROLE.get(pos, '?')})")
    bias = "Δscore>0 → 计算级预测：对 EGT 上下文的偏好较 hercynine 增强" if v["dscore"] > 0 \
        else "Δscore≤0 → 无 EGT 偏好增强信号"
    return ("三段锚定保留（R21/F111/T227 未动，约束级）；" + bias +
            "（🖥 设计偏置量，非亲和力；待 HPLC-MS 胞内外定量检验）；改动: " + " ".join(parts))


out_csv = os.path.join(PKG, "out/D4_variant_shortlist.csv")
FIELDS = ["rank", "variant_id", "pocket_fingerprint", "mutations", "ligand_confidence",
          "overall_confidence", "seq_rec", "score_EGT", "dscore_vs_hercynine",
          "criterion1", "criterion2", "criterion3", "n_pass", "is_WT", "WT_percentile",
          "expected_phenotype"]
with open(out_csv, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(FIELDS)
    wt_g = s2["WT"]
    w.writerow(["", "WT", "".join(wt_seq[i - 1] for i in SITES), "(WT)", "", "", "",
                f"{wt_score:.4f}", f"{wt_dscore:.4f}",
                wt_g["criterion1"], wt_g["criterion2"], wt_g["criterion3"], wt_g["n_pass"],
                1, f"score_EGT 分位 {wt_pct_score:.1f}% ｜ Δscore 分位 {wt_pct_dscore:.1f}%",
                "WT 基线：三段约束全过（T227–S 3.34 Å / R21–羧基 3.30 Å / 笼 4.19 Å）"])
    for v in shortlist:
        g = s2[v["variant_id"]]
        w.writerow([v["rank"], v["variant_id"], v["pocket_fingerprint"],
                    ";".join(v["mutations"]), v["ligand_confidence"], v["overall_confidence"],
                    v["seq_rec"], f"{v['score_EGT']:.4f}", f"{v['dscore']:.4f}",
                    g["criterion1"], g["criterion2"], g["criterion3"], v["n_pass"],
                    0, "", phenotype(v)])
print(f"[ok] {out_csv}: 10 条 + WT 行（指纹去重后候选池 {len(best_by_fp)} 条）")
