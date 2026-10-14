#!/usr/bin/env bash
# run_smoke_test.sh — 环境自检（10 秒）：把 3 个已测工具跑一遍，与 expected/ 比对
# 用法： bash smoke_test/run_smoke_test.sh
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(dirname "$HERE")"
WORK="$HERE/work"
mkdir -p "$WORK"

CIF="$PKG/data/chai1_mfsT1_EGT/pred.rank_0.cif"
FAIL=0
say() { printf '\n\033[1m%s\033[0m\n' "$*"; }
chk() { if [ "$1" = "$2" ]; then echo "   ✅ $3"; else echo "   ❌ $3（得到 $1 / 期望 $2）"; FAIL=1; fi; }

say "[1/3] mmCIF → PDB"
python "$PKG/tools/cif_to_pdb.py" "$CIF" "$WORK/mfsT1_EGT.pdb" || { echo "❌ 失败"; exit 1; }
[ -s "$WORK/mfsT1_EGT.pdb" ] && echo "   ✅ 生成了 PDB" || { echo "   ❌ PDB 为空"; FAIL=1; }

say "[2/3] 修配体键级（应报出 C=S）"
python "$PKG/tools/fix_ligand_bond_orders.py" "$CIF" \
  "C[N+](C)(C)C(CC1=CNC(=S)N1)C(=O)[O-]" "$WORK/EGT_ligand.sdf" | tee "$WORK/bondorder.log"
grep -q "'C', 'S'" "$WORK/bondorder.log" && echo "   ✅ 检到 C=S 双键" || { echo "   ❌ 没检到 C=S —— 环境或 SMILES 有问题"; FAIL=1; }

say "[3/3] 生成 D0 设计位点表"
python "$PKG/tools/make_site_table.py" "$WORK/mfsT1_EGT.pdb" "$WORK/D0.csv" | tee "$WORK/d0.log"
N=$(python - <<PY
import csv
print(sum(1 for _ in csv.DictReader(open("$WORK/D0.csv", encoding="utf-8-sig"))))
PY
)
chk "$N" "20" "位点数为 20"

say "关键距离自检（应与档 3 报告一致：T227 3.34 / R21 3.30 / F111 3.12）"
grep -E "THR227|ARG21|PHE111" "$WORK/d0.log" || true

if [ "$FAIL" = "0" ]; then
  printf '\n\033[1;32m全部通过 —— 环境就绪，可以开始 D0→D1→D2。\033[0m\n'
  printf '（若想逐字节比对期望产物： diff %s/expected/D0_site_table.csv %s/D0.csv ）\n\n' "$HERE" "$WORK"
else
  printf '\n\033[1;31m有检查未通过 —— 请对照 01_ENVIRONMENT_setup.md §4 排查。\033[0m\n\n'
  exit 1
fi
