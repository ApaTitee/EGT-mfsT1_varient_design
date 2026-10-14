#!/usr/bin/env bash
# run_ligandmpnn.sh — D1 序列设计（03_TOOL_SPECS.md · T4）
# 用法: conda activate egt4 && bash tools/run_ligandmpnn.sh
# 参数从 out/D0_site_table.csv 的约定自动生成（见 make_omit_json.py）。
set -eu
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LMPNN="${LMPNN_DIR:-$HOME/work/LigandMPNN}"
PDB="$PKG/out/mfsT1_EGT.pdb"                 # Chai-1 rank_0 复合物（蛋白链A + EGT 链B，HETATM）
OMIT="$PKG/out/D1_inputs/omit_AA_per_residue.json"
SEED=111

mkdir -p "$PKG/out/D1_inputs"
[ -f "$OMIT" ] || python "$PKG/tools/make_omit_json.py" "$PKG/out/D0_site_table.csv" "$OMIT"

for T in 0.1 0.2 0.3; do
  OUT="$PKG/out/D1_EGT_T${T}"
  echo "==== T=$T -> $OUT ===="
  ( cd "$LMPNN" && python run.py \
    --model_type ligand_mpnn \
    --pdb_path "$PDB" \
    --out_folder "$OUT" \
    --redesigned_residues "A14 A50 A53 A107 A139 A219 A223 A231 A322 A326 A346 A350 A353" \
    --fixed_residues "A21 A111 A135 A138 A227 A319 A349" \
    --omit_AA_per_residue "$OMIT" \
    --ligand_mpnn_use_atom_context 1 \
    --ligand_mpnn_use_side_chain_context 1 \
    --ligand_mpnn_cutoff_for_score 6.0 \
    --batch_size 2 --number_of_batches 50 \
    --temperature "$T" \
    --seed $SEED \
    --save_stats 1 \
    --pack_side_chains 1 \
    --number_of_packs_per_design 1 \
    --pack_with_ligand_context 1 \
  )
done
echo "[done] 3 个温度 × 100 条 = 300 条序列（含 packed 全原子模型）"
