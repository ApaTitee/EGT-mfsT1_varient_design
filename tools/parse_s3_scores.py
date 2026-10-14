#!/usr/bin/env python
"""parse_s3_scores.py — 解析 score.py 的 .pt，产出 D2_S3_selectivity_bias.csv

score 定义（写进结果说明，防误用）：
  score(variant | context) = 在 13 个可设计口袋位点上，"该变体实际氨基酸"的
  自回归 log-prob（10 个解码顺序取平均）再对位点取平均。
  Δscore = score(EGT 上下文) − score(hercynine 上下文) > 0 表示更偏好"有硫"上下文。
⛔ Δscore 是设计偏置量，不是亲和力差；不得换算 ΔΔG。

用法: python tools/parse_s3_scores.py  （在包根目录执行）
"""
import os, csv
import numpy as np
import torch

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SITES = [f"A{r}" for r in (14, 50, 53, 107, 139, 219, 223, 231, 322, 326, 346, 350, 353)]


def score_pt(path):
    d = torch.load(path, map_location="cpu", weights_only=False)
    names = d["residue_names"]              # {idx: "A14"} 或 list
    if isinstance(names, dict):
        pos = {v: int(k) for k, v in names.items()}
    else:
        pos = {v: i for i, v in enumerate(names)}
    lp = d["log_probs"]                     # [n_repeats, L, 21]
    nat = d["native_sequence"]              # [L]
    vals = []
    for s in SITES:
        i = pos[s]
        vals.append(lp[:, i, nat[i]].mean())
    return float(np.mean(vals))


rows = []
for ctx in ("EGT", "hercynine"):
    d = os.path.join(PKG, "out/D2_S3", f"scores_{ctx}")
    for f in sorted(os.listdir(d)):
        if f.endswith(".pt"):
            rows.append((f[:-3], ctx, score_pt(os.path.join(d, f))))

table = {}
for vid, ctx, s in rows:
    table.setdefault(vid, {})[ctx] = s

out_csv = os.path.join(PKG, "out/D2_S3_selectivity_bias.csv")
with open(out_csv, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["variant_id", "score_EGT", "score_hercynine", "dscore_vs_hercynine"])
    w.writerow(["WT", f"{table['WT']['EGT']:.4f}", f"{table['WT']['hercynine']:.4f}",
                f"{table['WT']['EGT'] - table['WT']['hercynine']:.4f}"])
    for vid in sorted(table):
        if vid == "WT":
            continue
        e, h = table[vid]["EGT"], table[vid]["hercynine"]
        w.writerow([vid.replace("_id", "/id").replace("_", "/"), f"{e:.4f}", f"{h:.4f}", f"{e - h:.4f}"])
print(f"[ok] {out_csv}: {len(table)} 条（含 WT）")
