#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fix_ligand_bond_orders.py — 用 RDKit 把 Chai-1 配体的**键级**修好，写出 SDF

为什么必须做：Chai-1 的 mmCIF **不含 `_chem_comp_bond`**（已核：0 条），
只有坐标 + 由距离可推的连键关系 → 硫酮 C=S 会被当成单键 C–S。
本脚本：
  1) 从 cif 取配体坐标与距离图（≤1.9 Å 成键）；
  2) 用给定 SMILES 建"键级模板"（**必须用 Kekulé 式**，见 ligand_smiles.txt）；
  3) 在两张图之间做小图同构搜索（≤20 重原子，回溯，秒级）得到原子映射；
  4) 按模板键级写出 SDF（3D；含正确 C=S 双键）。

用法： python fix_ligand_bond_orders.py <in.cif> <smiles> <out.sdf>
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_cif, lig_graph

from rdkit import Chem
from rdkit.Geometry import Point3D


def clean_sym(s):
    s = (s or "C").strip()
    return s[0].upper() + s[1:].lower() if len(s) > 1 else s.upper()


def isomorphism(src_adj, tgt_adj, src_el, tgt_el):
    """把 src 映射到 tgt（同构回溯）。返回 src_idx -> tgt_idx；失败返回 None"""
    n = len(src_adj)
    if n != len(tgt_adj):
        return None
    order = sorted(range(n), key=lambda i: -len(src_adj[i]))
    mapping, used = {}, set()

    def bt(k):
        if k == n:
            return True
        i = order[k]
        for j in range(n):
            if j in used or tgt_el[j] != src_el[i] or len(tgt_adj[j]) != len(src_adj[i]):
                continue
            if all((mapping[m] in tgt_adj[j]) for m in src_adj[i] if m in mapping):
                mapping[i] = j; used.add(j)
                if bt(k + 1):
                    return True
                used.discard(j); del mapping[i]
        return False

    return dict(mapping) if bt(0) else None


def main(cif, smiles, out, cutoff=1.9):
    prot, lig = read_cif(cif)
    if not lig:
        raise SystemExit("[ERROR] cif 里没有 HETATM 配体")
    adj = lig_graph(lig, cutoff)
    src_el = [clean_sym(a.get("el") or a["atom"]).upper() for a in lig]

    tmpl = Chem.MolFromSmiles(smiles)
    if tmpl is None:
        raise SystemExit("[ERROR] SMILES 解析失败：" + smiles +
                         "\n        提示：硫酮请用 Kekulé 式（见 ligand_smiles.txt）")
    # ⚠️ 关键：先把模板凯库勒化并清掉芳香标记。否则从模板搬过来的 AROMATIC 键
    #    会被加到未标记芳香的 RWMol 原子上 → SDWriter 报 KekulizeException。
    try:
        Chem.Kekulize(tmpl, clearAromaticFlags=True)
    except Exception as e:
        raise SystemExit(f"[ERROR] 模板无法凯库勒化（{e}）—— 请改用 Kekulé 式 SMILES")
    M = tmpl.GetNumAtoms()
    t_adj = {i: [b.GetOtherAtomIdx(i) for b in tmpl.GetAtomWithIdx(i).GetBonds()]
             for i in range(M)}
    t_el = [a.GetSymbol().upper() for a in tmpl.GetAtoms()]
    t_bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondType()) for b in tmpl.GetBonds()]
    print(f"  模板：{M} 重原子，{len(t_bonds)} 根键")
    print(f"  位姿：{len(lig)} 个重原子，{sum(len(v) for v in adj.values())//2} 根距离键")

    mp = isomorphism(adj, t_adj, src_el, t_el)
    if mp is None:
        raise SystemExit("[ERROR] 配体骨架与 SMILES 不同构 —— 检查 SMILES 或 --cutoff")
    inv = {v: k for k, v in mp.items()}          # 模板 idx -> 位姿 idx
    print("  同构成立：", {f"pose{i}": f"tmpl{j}" for i, j in sorted(mp.items())})

    em = Chem.RWMol()
    for i in range(len(lig)):
        em.AddAtom(Chem.Atom(clean_sym(src_el[i])))
    for i, j, bt in t_bonds:
        em.AddBond(inv[i], inv[j], bt)
    m = em.GetMol()
    conf = Chem.Conformer(len(lig))
    for i, a in enumerate(lig):
        conf.SetAtomPosition(i, Point3D(*[float(v) for v in a["xyz"]]))
    m.AddConformer(conf, assignId=True)
    Chem.AssignStereochemistryFrom3D(m)

    w = Chem.SDWriter(out); w.write(m); w.close()
    dbl = [(b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol())
           for b in m.GetBonds() if abs(b.GetBondTypeAsDouble() - 2.0) < 1e-6]
    print(f"  写出 {out}：{m.GetNumAtoms()} 重原子，{len(dbl)} 个双键 -> {dbl}")
    return out


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(__doc__); raise SystemExit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
