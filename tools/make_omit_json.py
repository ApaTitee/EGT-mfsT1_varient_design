#!/usr/bin/env python
"""从 D0 位点表生成 LigandMPNN 的 --omit_AA_per_residue JSON（把每个可设计位点约束在其 allowed_AA 集合内）。

用法: python tools/make_omit_json.py out/D0_site_table.csv out/D1_inputs/omit_AA_per_residue.json
"""
import csv, json, sys

ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"  # LigandMPNN 字母表（含 X，一并 omit）

def main(site_table: str, out_json: str) -> None:
    omit = {}
    designed = []
    with open(site_table, encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            allowed = (row["allowed_AA"] or "").strip()
            if not allowed:
                continue  # fixed 层不参与设计
            if row["layer"] not in ("shell", "observe"):
                continue  # core 锚点在 run 时用 --fixed_residues 固定，不在这里约束
            res = f"{row['chain']}{row['resnum']}"
            omit[res] = "".join(a for a in ALPHABET if a not in allowed)
            designed.append(f"{res}:{allowed}")
    with open(out_json, "w") as fh:
        json.dump(omit, fh, indent=1)
    print(f"[ok] 写出 {out_json}，约束 {len(omit)} 个位点：")
    for d in designed:
        print("   ", d)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
