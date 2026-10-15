#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
common.py — 档 4 工具集的共享底层（纯 numpy / 标准库；RDKit 仅个别脚本需要）

提供：
  · mmCIF / PDB 读取（自动识别 AF3 与 Chai-1 两种列布局）
  · PDB 写出（固定宽度、含 element 列，供 Prody / LigandMPNN 读）
  · Kabsch 叠合 / 距离工具 / 元素与残基工具
  · 配体小分子图（按距离建键、识别季铵 N / 羧基 C / 硫酮 S）

⛔ 所有坐标输出都不做任何"修正"：只搬运与变换，不推测几何。
"""
import os, math, json
import numpy as np

# ---------------------------------------------------------------- 常量
AA3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
       "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
AA3to1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
          "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
          "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
          "TYR": "Y", "VAL": "V"}
AA1to3 = {v: k for k, v in AA3to1.items()}

# mfsT1 的 12 段 TM（DeepTMHMM 2，见 docs/01_family_and_topology.md）
TM_RANGES = [(8, 28), (42, 64), (72, 89), (95, 113), (131, 152), (161, 179),
             (214, 234), (251, 270), (280, 297), (307, 326), (344, 361), (377, 394)]
RING_ATOMS = {"PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
              "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
              "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
              "TRP": ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
              "ASN": ["CG", "OD1", "ND2"]}


def tm_of(resnum):
    for i, (a, b) in enumerate(TM_RANGES, 1):
        if a <= resnum <= b:
            return f"TM{i}"
    return "loop"


# ---------------------------------------------------------------- 读取
def read_cif(path):
    """返回 (prot, lig)；自动识别 AF3 / Chai-1 两种 _atom_site 列布局"""
    prot, lig = [], []
    for ln in open(path, encoding="utf-8", errors="ignore"):
        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
            continue
        p = ln.split()
        if p[15].isdigit():                       # AF3
            atom, comp, resnum, chain, el = p[3], p[5], int(p[8]), p[6], p[2]
            xyz = p[10:13]
        else:                                     # Chai-1
            atom, comp, resnum, chain, el = p[3], p[5], int(p[7]), p[9], p[2]
            xyz = p[10:13]
        d = dict(atom=atom, resname=comp, resnum=resnum, chain=chain, el=el,
                 xyz=np.array([float(v) for v in xyz]))
        (lig if ln.startswith("HETATM") else prot).append(d)
    return prot, lig


def read_pdb(path, chain=None, het=None, keep_h=False):
    """读 PDB；chain=None 表示所有链；het=None 表示只读 ATOM"""
    prot, lig = [], []
    for ln in open(path, encoding="utf-8", errors="ignore"):
        rec = ln[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        ch = ln[21]
        el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
        if not keep_h and el == "H":
            continue
        d = dict(atom=ln[12:16].strip(), resname=ln[17:20].strip(),
                 resnum=int(ln[22:26]), icode=ln[26].strip(), chain=ch, el=el,
                 xyz=np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
        if rec == "ATOM":
            if chain and ch != chain:
                continue
            prot.append(d)
        else:
            if het is not None and d["resname"] != het:
                continue
            if chain and ch != chain:
                continue
            lig.append(d)
    return prot, lig


def by_res(atoms):
    out = {}
    for a in atoms:
        out.setdefault((a["resnum"], a["resname"]), []).append(a)
    return out


def seq_of(prot):
    """按链返回 {(chain,resnum): aa1} 与有序列表"""
    order, seq = [], {}
    for a in prot:
        if a["atom"] != "CA":
            continue
        k = (a["chain"], a["resnum"])
        if k not in seq:
            order.append(k)
            seq[k] = AA3to1.get(a["resname"], "X")
    return order, seq


# ---------------------------------------------------------------- 写出
def write_pdb(path, prot, lig=(), lig_chain="B", lig_resname=None):
    """写 PDB（固定宽度；element 列必填，Prody 依赖）"""
    lines = []
    serial = 1

    def fmt(rec, a, chain, resname, resnum, charge=""):
        nonlocal serial
        el = a["el"].rjust(2)
        nm = a["atom"]
        nm = (" " + nm if len(nm) < 4 else nm)[:4]
        lines.append(
            f"{rec:<6}{serial:>5} {nm:<4}{'':1}{resname:>3} {chain}{resnum:>4}{'':1}   "
            f"{a['xyz'][0]:>8.3f}{a['xyz'][1]:>8.3f}{a['xyz'][2]:>8.3f}"
            f"{1.00:>6.2f}{0.00:>6.2f}          {el:>2}{charge:>2}")
        serial += 1

    last_chain = None
    for a in prot:
        if last_chain is not None and a["chain"] != last_chain:
            lines.append("TER")
        fmt("ATOM", a, a["chain"], a["resname"], a["resnum"])
        last_chain = a["chain"]
    if lig:
        lines.append("TER")
        rn = lig_resname or (lig[0]["resname"][:3] if lig[0].get("resname") else "LIG")
        for a in lig:
            fmt("HETATM", a, lig_chain, rn, 1)
    lines.append("END")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return path


# ---------------------------------------------------------------- 几何
def kabsch(P, Q):
    """把 P 叠合到 Q：返回 (R, t, rmsd)，使 (P @ R.T + t) ≈ Q"""
    P = np.asarray(P, float); Q = np.asarray(Q, float)
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = qc - R @ pc
    rmsd = float(np.sqrt(((P @ R.T + t - Q) ** 2).sum() / len(P)))
    return R, t, rmsd


def apply_rt(xyz, R, t):
    return np.asarray(xyz, float) @ R.T + t


def contacts(prot, lig, cutoff=4.5, exclude_h=True):
    """返回 [(resnum, resname, dist, atom_prot, atom_lig)]，按距离升序"""
    out = []
    for a in prot:
        if exclude_h and a["el"] == "H":
            continue
        for b in lig:
            if exclude_h and b["el"] == "H":
                continue
            d = float(np.linalg.norm(a["xyz"] - b["xyz"]))
            if d <= cutoff:
                out.append((a["resnum"], a["resname"], d, a["atom"], b["atom"]))
    return sorted(out, key=lambda x: x[2])


# ---------------------------------------------------------------- 配体图
def lig_graph(lig, cutoff=1.9):
    n = len(lig)
    adj = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if float(np.linalg.norm(lig[i]["xyz"] - lig[j]["xyz"])) < cutoff:
                adj[i].append(j); adj[j].append(i)
    return adj


def lig_sites(lig, adj=None):
    """识别官能团：季铵 N / 羧基 C（+两个 O）/ 硫（硫酮或硫醇）/ 咪唑 C2"""
    adj = adj or lig_graph(lig)
    info = dict(N_quat=None, C_carbox=None, O_carbox=[], S=None, C2_imidazole=None)
    for i, a in enumerate(lig):
        nm, el = a["atom"], a.get("el", a["atom"][0])
        nbr = adj[i]
        if el == "N" or nm.startswith("N"):
            me = [j for j in nbr if lig[j]["atom"].startswith("C") and len(adj[j]) == 1]
            if len(me) == 3:
                info["N_quat"] = i
        if el == "C" or nm.startswith("C"):
            # 咪唑 C2 = 与 2 个 N 相连的环碳（重原子邻居恰好为 2 个 N）
            n_n = sum(1 for j in nbr if lig[j]["atom"].startswith("N"))
            n_c = sum(1 for j in nbr if lig[j]["atom"].startswith("C"))
            if n_n == 2 and n_c == 0:
                info["C2_imidazole"] = i
        if nm.startswith("C") and sum(1 for j in nbr if lig[j]["atom"].startswith("O")) == 2:
            info["C_carbox"] = i
            info["O_carbox"] = [j for j in nbr if lig[j]["atom"].startswith("O")]
        if (el == "S" or nm.startswith("S")) and info["S"] is None:
            info["S"] = i
    return info


def dist_min(lig, i, prot_ats):
    if i is None or not prot_ats:
        return None
    return float(min(np.linalg.norm(a["xyz"] - lig[i]["xyz"]) for a in prot_ats))


# ---------------------------------------------------------------- 杂项
def load_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def die(msg):
    raise SystemExit("[ERROR] " + msg)
