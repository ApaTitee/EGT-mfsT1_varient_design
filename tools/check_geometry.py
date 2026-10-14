#!/usr/bin/env python
"""check_geometry.py — D2-S2 三段约束几何复核（03_TOOL_SPECS.md · T2）

输入：一个或多个 PDB（蛋白 + HETATM 配体；LigandMPNN packed/*.pdb 或 WT 复合物）
输出：CSV，每条记录给出三段约束的距离与 PASS/FAIL/NA。

判据（阈值写死，与 03_TOOL_SPECS.md 一致）：
  ① 三甲铵 N → 芳香笼：到 6 个笼残基环原子最近 ≤5.0 Å 且 ≤5 Å 的不同笼残基 ≥2
  ② 羧基 C/O → 21 位正电 N（ARG:NH1/NH2/NE；LYS:NZ）≤4.2 Å（spec 4.5，见下）
  ③ 硫酮 S → CH 供体（227:THR→CG2/SER→CB；备选 50:CD1、231:CD1/CD2）≤4.2 Å
  ③′ 配体无 S 时改用咪唑 C2，阈值同 ③
"""
import os, sys, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from common import read_pdb, by_res, lig_sites, RING_ATOMS

CAGE_SITES = [14, 53, 111, 219, 326, 350]          # D0 表中 role in (cage, cage_major)
THR_CAGE, THR_COO = 5.0, 4.5                        # 判据①②阈值
THR_S = 4.2                                         # 判据③阈值
POS_PARTNER = {"ARG": ["NH1", "NH2", "NE"], "LYS": ["NZ"]}
CH_DONOR = {227: {"THR": ["CG2"], "SER": ["CB"]},
            50:  {"ILE": ["CD1"], "LEU": ["CD1"], "VAL": ["CG1", "CG2"]},
            231: {"LEU": ["CD1", "CD2"], "ILE": ["CD1"], "VAL": ["CG1", "CG2"]}}


def check(pdb_path, variant_id):
    prot, lig = read_pdb(pdb_path)
    if not lig:
        return dict(variant_id=variant_id, note="no_ligand", criterion1="NA",
                    criterion2="NA", criterion3="NA", n_pass=0)
    res = by_res(prot)
    sites = lig_sites(lig)

    # ① 三甲铵 → 芳香笼
    n_q = sites["N_quat"]
    if n_q is not None:
        mind, hits = 1e9, set()
        for rn in CAGE_SITES:
            key = next((k for k in res if k[0] == rn), None)
            if not key:
                continue
            ring = [a for a in res[key] if a["atom"] in RING_ATOMS.get(a["resname"], [])]
            for a in ring:
                d = float(np.linalg.norm(a["xyz"] - lig[n_q]["xyz"]))
                mind = min(mind, d)
                if d <= THR_CAGE:
                    hits.add(rn)
        c1 = "PASS" if (mind <= THR_CAGE and len(hits) >= 2) else "FAIL"
        d1, n1 = round(mind, 2), len(hits)
    else:
        c1, d1, n1 = "NA", "", 0

    # ② 羧基 → 21 位正电伙伴
    key21 = next((k for k in res if k[0] == 21), None)
    c_atoms = [lig[i] for i in ([sites["C_carbox"]] + sites["O_carbox"]) if i is not None]
    if key21 and c_atoms and key21[1] in POS_PARTNER:
        don = [a for a in res[key21] if a["atom"] in POS_PARTNER[key21[1]]]
        mind2 = min(float(np.linalg.norm(c["xyz"] - d["xyz"])) for c in c_atoms for d in don) if don else 1e9
        c2 = "PASS" if mind2 <= THR_COO else "FAIL"
        d2 = round(mind2, 2)
    else:
        c2, d2 = "NA", ""

    # ③ / ③′ 硫端 → CH 供体
    probe = sites["S"] if sites["S"] is not None else sites["C2_imidazole"]
    if probe is not None:
        cand = []
        for rn, tbl in CH_DONOR.items():
            key = next((k for k in res if k[0] == rn), None)
            if key and key[1] in tbl:
                cand += [a for a in res[key] if a["atom"] in tbl[key[1]]]
        mind3 = min((float(np.linalg.norm(a["xyz"] - lig[probe]["xyz"])) for a in cand), default=1e9)
        c3 = "PASS" if mind3 <= THR_S else "FAIL"
        d3 = round(mind3, 2)
        if sites["S"] is None:
            c3 = c3  # ③′：同阈值；note 里标注
    else:
        c3, d3 = "NA", ""

    n_pass = sum(1 for c in (c1, c2, c3) if c == "PASS")
    note = [] if sites["S"] is not None else ["no_sulfur_probe=C2_imidazole"]
    return dict(variant_id=variant_id, n_quat_to_cage_min=d1, cage_contacts_le5=n1,
                carboxyl_to_R21=d2, sulfur_to_donor=d3,
                criterion1=c1, criterion2=c2, criterion3=c3,
                n_pass=n_pass, note=";".join(note))


FIELDS = ["variant_id", "n_quat_to_cage_min", "cage_contacts_le5", "carboxyl_to_R21",
          "sulfur_to_donor", "criterion1", "criterion2", "criterion3", "n_pass", "note"]

if __name__ == "__main__":
    # 用法: check_geometry.py out.csv pdb1 pdb2 ...   （variant_id 取文件名去扩展名）
    out_csv, pdbs = sys.argv[1], sys.argv[2:]
    rows = [check(p, os.path.splitext(os.path.basename(p))[0]) for p in pdbs]
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    npass = sum(r["n_pass"] == 3 for r in rows)
    print(f"[ok] {out_csv}: {len(rows)} 条，{npass} 条三段全过")
