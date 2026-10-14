#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_site_table.py — D0：生成「设计位点表」（档 4 的第一步，纯结构、无外部依赖）

输入：一个含配体的复合物结构（Chai-1 cif 或已转好的 PDB）
输出：D0_site_table.csv
      字段：chain,resnum,resname,aa1,TM,min_dist_to_ligand,role,layer,
            allowed_AA,fixed,is_cross_bundle,note

判定规则（全部写在下方常量里，便于复核）：
  · 核心锚点（core）：R21（羧基盐桥）/ T227（读硫 CH 供体）/ F111（笼主接触）
    → 只允许同功能替代；默认在 LigandMPNN 里 **固定**
  · 次级衬里（shell）：可采样位点
  · 观察项（observe）：D107（口袋内唯一负电）
  · 跨束界面（cross-bundle）：该残基 Cα 在 8 Å 内**同时**碰到 TM1–6 与 TM7–12 的残基
    （= 改它更可能影响交替开放循环，见档 4 任务书 §D0）

用法： python make_site_table.py <complex.cif|.pdb> <out.csv> [--cutoff 9.0]
"""
import sys, os, csv, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_cif, read_pdb, by_res, tm_of, RING_ATOMS, write_pdb

CORE = {21: "R", 227: "T", 111: "F"}
SHELL = {14: "Y", 53: "Y", 219: "H", 326: "F", 350: "F",
         50: "I", 231: "L", 223: "Q", 346: "N", 322: "S", 139: "Q", 353: "S"}
OBSERVE = {107: "D"}
# TM 扭结处的结构 Gly：不改（它们是折叠所需，且不是识别化学的一部分）
GLY_KINK = {135: "G", 138: "G", 319: "G", 349: "G"}

AROMATIC_SHELL = [14, 53, 219, 326, 350]
POLAR_SHELL = [223, 346, 322, 139, 353]
HYDROPHOBIC_SHELL = [50, 231]

ALLOWED = {}
ALLOWED.update({k: "RK" for k in CORE if CORE[k] == "R"})
ALLOWED.update({k: "TSV" for k in CORE if CORE[k] == "T"})       # 读硫：Thr 首选，Ser/Val 为同几何备选
ALLOWED.update({k: "FYW" for k in CORE if CORE[k] == "F"})
ALLOWED.update({k: "FYWH" for k in AROMATIC_SHELL})
ALLOWED.update({k: "QNSTDEHKR" for k in POLAR_SHELL})
ALLOWED.update({k: "ILVMAF" for k in HYDROPHOBIC_SHELL})
ALLOWED.update({k: "DENQ" for k in OBSERVE})

LAYER = {}
LAYER.update({k: "core" for k in CORE})
LAYER.update({k: "shell" for k in SHELL})
LAYER.update({k: "observe" for k in OBSERVE})
LAYER.update({k: "fixed_backbone" for k in GLY_KINK})
ROLE = {}
ROLE.update({21: "carboxyl_partner", 227: "read_sulfur", 111: "cage_major"})
ROLE.update({k: "cage" for k in AROMATIC_SHELL})
ROLE.update({k: "lining_polar" for k in POLAR_SHELL})
ROLE.update({k: "lining_hydrophobic" for k in HYDROPHOBIC_SHELL})
ROLE.update({107: "only_negative"})

LIG_NAMES = ("LIG", "LIG2", "LW8")


def load_any(path):
    if path.lower().endswith((".cif", ".mmcif")):
        prot, lig = read_cif(path)
    else:
        prot, lig = read_pdb(path, het=None)
        lig = [a for a in lig if a["resname"] in LIG_NAMES]
    return prot, lig


def cross_bundle_flags(res, prot):
    """返回 {resnum: bool}：Cα 是否同时靠近 TM1–6 与 TM7–12"""
    ca = {}
    for a in prot:
        if a["atom"] == "CA":
            ca[a["resnum"]] = a["xyz"]
    lo = {n: p for n, p in ca.items() if tm_of(n) in [f"TM{i}" for i in range(1, 7)]}
    hi = {n: p for n, p in ca.items() if tm_of(n) in [f"TM{i}" for i in range(7, 13)]}
    out = {}
    for n, p in ca.items():
        near_lo = any(np.linalg.norm(p - q) <= 8.0 for m, q in lo.items() if m != n)
        near_hi = any(np.linalg.norm(p - q) <= 8.0 for m, q in hi.items() if m != n)
        out[n] = bool(near_lo and near_hi)
    return out


def main(complex_path, out_csv, cutoff=9.0):
    prot, lig = load_any(complex_path)
    if not lig:
        sys.exit("[ERROR] 结构里没有配体（HETATM）—— D0 需要复合物")
    res = by_res(prot)
    xb = cross_bundle_flags(res, prot)

    # 口袋中心 = 档 1–3 的锚点（Y14/Q57/F111/Q139/H219/A260/G361）Cα 质心
    ANCHORS = (14, 57, 111, 139, 219, 260, 361)
    cen_pts = [a["xyz"] for (num, _), ats in res.items() if num in ANCHORS
               for a in ats if a["atom"] == "CA"]
    if not cen_pts:
        sys.exit("[ERROR] 找不到任何锚点残基的 Cα —— 残基编号可能不对")
    cen = np.mean(cen_pts, axis=0)
    print(f"口袋中心取自 {len(cen_pts)} 个锚点 Cα")

    rows = []
    for (num, rn), ats in sorted(res.items()):
        d_lig = min(float(np.linalg.norm(a["xyz"] - b["xyz"])) for a in ats for b in lig
                    if a.get("el") != "H" and b.get("el") != "H")
        d_cen = min(float(np.linalg.norm(a["xyz"] - cen)) for a in ats if a["atom"] == "CA") \
            if any(a["atom"] == "CA" for a in ats) else 999
        in_pocket = d_cen <= cutoff
        layer = LAYER.get(num, "outside" if not in_pocket else "pocket_unassigned")
        if layer == "outside":
            continue
        rows.append(dict(
            chain=ats[0]["chain"], resnum=num, resname=rn,
            aa1="",
            TM=tm_of(num), min_dist_to_ligand=round(d_lig, 2),
            dist_to_pocket_center=round(d_cen, 2),
            role=ROLE.get(num, "unassigned"), layer=layer,
            allowed_AA=ALLOWED.get(num, ""),
            fixed=int(layer in ("core", "fixed_backbone")),
            is_cross_bundle=int(xb.get(num, False)),
            note=""))
    # 补 aa1
    from common import AA3to1
    for r in rows:
        r["aa1"] = AA3to1.get(r["resname"], "X")

    cols = ["chain", "resnum", "resname", "aa1", "TM", "min_dist_to_ligand",
            "dist_to_pocket_center", "role", "layer", "allowed_AA", "fixed",
            "is_cross_bundle", "note"]
    with open(out_csv, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print(f"写出 {out_csv}（{len(rows)} 个位点）")
    for r in rows:
        tag = {"core": "★", "shell": "·", "observe": "◇", "fixed_backbone": "="}.get(r["layer"], "?")
        print(f"  {tag} {r['resname']}{r['resnum']:<4} {r['TM']:<5} d_lig={r['min_dist_to_ligand']:<5} "
              f"{r['role']:<20} allowed={r['allowed_AA']:<12} cross_bundle={r['is_cross_bundle']}")
    fixed = " ".join(f"A{r['resnum']}" for r in rows if r["fixed"])
    redes = " ".join(f"A{r['resnum']}" for r in rows if not r["fixed"])
    print("\n[LigandMPNN 参数建议]")
    print(f"  --fixed_residues \"{fixed}\"")
    print(f"  --redesigned_residues \"{redes}\"")
    return out_csv


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("complex_path"); ap.add_argument("out_csv")
    ap.add_argument("--cutoff", type=float, default=9.0)
    a = ap.parse_args()
    main(a.complex_path, a.out_csv, a.cutoff)
