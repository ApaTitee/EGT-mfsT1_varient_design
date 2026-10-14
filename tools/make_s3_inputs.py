#!/usr/bin/env python
"""make_s3_inputs.py — 为 D2-S3 生成两套配体上下文的 PDB 集合 + pdb_path_multi JSON

- pdbs_EGT/      : 变体 packed 蛋白 + 原 EGT 配体（packed 自带）+ WT
- pdbs_hercynine/: 同一蛋白 + T1 叠合后的 hercynine 位姿（共用 EGT 骨架坐标系）

用法: python tools/make_s3_inputs.py  （在包根目录执行）
"""
import os, sys, glob, json, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_pdb, write_pdb

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HERC_PDB = os.path.join(PKG, "out/backbone_EGT_with_hercynine_pose.pdb")
OUT = os.path.join(PKG, "out/D2_S3")

_, herc_lig = read_pdb(HERC_PDB, het="LIG")
assert len(herc_lig) == 14, f"hercynine 应为 14 重原子，得到 {len(herc_lig)}"

egt_dir = os.path.join(OUT, "pdbs_EGT"); her_dir = os.path.join(OUT, "pdbs_hercynine")
os.makedirs(egt_dir, exist_ok=True); os.makedirs(her_dir, exist_ok=True)

# 收集 (variant_id, packed_pdb)
variants = [("WT", os.path.join(PKG, "out/mfsT1_EGT.pdb"))]
for T in ("0.1", "0.2", "0.3"):
    for p in sorted(glob.glob(os.path.join(PKG, f"out/D1_EGT_T{T}/packed/*.pdb")),
                    key=lambda x: int(x.rsplit("_", 2)[1])):
        i = int(p.rsplit("_", 2)[1])
        variants.append((f"T{T}_id{i}", p))

egt_list, her_list = [], []
for vid, pdb in variants:
    dst = os.path.join(egt_dir, vid + ".pdb")
    shutil.copy(pdb, dst)
    egt_list.append(dst)
    prot, _ = read_pdb(pdb)
    hp = os.path.join(her_dir, vid + ".pdb")
    write_pdb(hp, prot, herc_lig, lig_chain="B", lig_resname="LIG")
    her_list.append(hp)

with open(os.path.join(OUT, "list_EGT.json"), "w") as fh:
    json.dump(egt_list, fh, indent=1)
with open(os.path.join(OUT, "list_hercynine.json"), "w") as fh:
    json.dump(her_list, fh, indent=1)
print(f"[ok] {len(variants)} 个变体 × 2 套上下文 -> {OUT}")
