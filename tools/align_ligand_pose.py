#!/usr/bin/env python
"""align_ligand_pose.py — D2-S3 前置（03_TOOL_SPECS.md · T1）

把 hercynine 复合物（另一次独立 Chai-1 预测）的配体位姿，经口袋 Cα Kabsch 叠合
搬到 EGT 复合物的骨架坐标系，使 S3 反向筛选的两个配体上下文共用同一骨架。

用法:
  python tools/align_ligand_pose.py <ref_complex.cif> <donor_complex.cif> <out.pdb> \
      --sites out/D0_site_table.csv [--resid out/D2_S3_align_residuals.csv]
"""
import os, sys, csv, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from common import read_cif, kabsch, apply_rt, write_pdb

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref"); ap.add_argument("donor"); ap.add_argument("out")
    ap.add_argument("--sites", required=True)
    ap.add_argument("--resid", default=None)
    a = ap.parse_args()

    # 口袋位点 = D0 表里 layer in (core, shell)
    site_resnums = []
    with open(a.sites, encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            if row["layer"] in ("core", "shell"):
                site_resnums.append(int(row["resnum"]))

    ref_prot, _ = read_cif(a.ref)
    don_prot, don_lig = read_cif(a.donor)

    def ca_map(prot):
        return {x["resnum"]: x["xyz"] for x in prot
                if x["atom"] == "CA" and x["chain"] == "A"}
    R_ca, D_ca = ca_map(ref_prot), ca_map(don_prot)
    common = [r for r in site_resnums if r in R_ca and r in D_ca]
    P = np.array([D_ca[r] for r in common]); Q = np.array([R_ca[r] for r in common])
    R, t, rmsd = kabsch(P, Q)

    resid = [(r, float(np.linalg.norm(P[i] @ R.T + t - Q[i]))) for i, r in enumerate(common)]
    print(f"pocket_CA_RMSD={rmsd:.2f} Å (n={len(common)})")
    if a.resid:
        with open(a.resid, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh); w.writerow(["resnum", "resname", "ca_dist_after_fit"])
            rn = {x["resnum"]: x["resname"] for x in ref_prot if x["atom"] == "CA"}
            for r, d in resid:
                w.writerow([r, rn.get(r, "?"), f"{d:.3f}"])
        print(f"[ok] 逐位点残差 -> {a.resid}")
    if rmsd > 1.5:
        print(f"[warn] RMSD > 1.5 Å —— 需在结果说明中注明")

    moved = [dict(x, xyz=apply_rt(x["xyz"][None, :], R, t)[0]) for x in don_lig]
    write_pdb(a.out, ref_prot, moved, lig_chain="B", lig_resname="LIG")
    print(f"[ok] 写出 {a.out}（参考蛋白 + 变换后的供体配体，{len(moved)} 个配体原子）")

if __name__ == "__main__":
    main()
