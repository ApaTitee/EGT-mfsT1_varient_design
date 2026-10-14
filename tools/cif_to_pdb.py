#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""cif_to_pdb.py — 把 AF3 / Chai-1 的 mmCIF 转成 PDB（蛋白质链 A + 配体链 B）

用途：LigandMPNN 与 Prody 只吃 PDB。
用法：
    python cif_to_pdb.py <in.cif> <out.pdb> [--lig LIG2] [--lig-name LIG]
不带参数时按 config.yaml 里的默认列表批量转换（见 README）。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_cif, write_pdb


def convert(src, dst, lig_name="LIG"):
    prot, lig = read_cif(src)
    if not prot:
        raise SystemExit(f"[ERROR] {src} 里没有 ATOM 记录")
    # 只保留链 A 的蛋白质（AF3/Chai-1 都是单链）
    prot = [a for a in prot if a["chain"] == "A"] or prot
    for a in lig:
        a["resname"] = lig_name
    write_pdb(dst, prot, lig)
    print(f"  {os.path.basename(src)} -> {os.path.basename(dst)}  "
          f"({len(prot)} protein atoms, {len(lig)} ligand atoms)")
    return dst


if __name__ == "__main__":
    if len(sys.argv) >= 3:
        convert(sys.argv[1], sys.argv[2],
                lig_name=(sys.argv[4] if len(sys.argv) > 4 else "LIG"))
    else:
        print(__doc__)
