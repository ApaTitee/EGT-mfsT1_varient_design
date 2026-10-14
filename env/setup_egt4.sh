#!/usr/bin/env bash
# setup_egt4.sh — 一键搭建档 4 的 WSL 环境（按 2026-10-14 实测写定；不需要 sudo）
# 用法： bash env/setup_egt4.sh            （在包根目录执行）
#       SKIP_LIGANDMPNN=1 bash env/setup_egt4.sh   （只建环境，不 clone LigandMPNN）
set -u
PKG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LMPNN_DIR="${LMPNN_DIR:-$HOME/work/LigandMPNN}"
ENV_NAME="${ENV_NAME:-egt4}"

echo "=============================================="
echo " 档 4 环境搭建"
echo " 包目录   : $PKG_DIR"
echo " LigandMPNN: $LMPNN_DIR"
echo " conda env: $ENV_NAME"
echo "=============================================="

# ---------- 0. 依赖检查 ----------
for c in conda git python; do
  command -v "$c" >/dev/null || { echo "[ERROR] 找不到 $c"; exit 1; }
done
echo "[ok] conda / git 都在"
echo "[info] 代理变量： http_proxy=${http_proxy:-<未设>}"

# ---------- 1. 建环境 ----------
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "[skip] conda 环境 $ENV_NAME 已存在"
else
  echo "[1/4] 建 conda 环境（python 3.11）..."
  conda create -y -n "$ENV_NAME" python=3.11 || exit 1
fi

# 在脚本里直接调用该环境的 python/pip，避免依赖 conda activate
EPY="$(conda run -n "$ENV_NAME" which python)"
echo "[info] 环境 python = $EPY"
run_py() { conda run --no-capture-output -n "$ENV_NAME" python "$@"; }
run_pip() { conda run --no-capture-output -n "$ENV_NAME" pip "$@"; }

# ---------- 2. PyTorch ----------
if run_py -c "import torch" >/dev/null 2>&1; then
  echo "[skip] torch 已装： $(run_py -c 'import torch;print(torch.__version__)')"
else
  echo "[2/4] 安装 PyTorch（cu128，适配 RTX 5060 / sm_120）..."
  run_pip install torch --index-url https://download.pytorch.org/whl/cu128 || {
    echo "[warn] cu128 安装失败 → 回退 CPU 版（本任务量级足够）"
    run_pip install torch || exit 1
  }
fi
run_py - <<'PY'
import torch
print("[check] torch", torch.__version__, "| cuda_ok:", torch.cuda.is_available(),
      "|", (torch.cuda.get_device_name(0) if torch.cuda.is_available() else "CPU 模式"))
PY

# ---------- 3. 本包工具需要的包 ----------
echo "[3/4] 安装 numpy / scipy / matplotlib / pandas / biopython / rdkit / gemmi ..."
run_pip install numpy scipy matplotlib pandas biopython rdkit gemmi || exit 1
run_py -c "import numpy,scipy,matplotlib,rdkit,gemmi;print('[check] 依赖 ok')"

# ---------- 4. LigandMPNN ----------
if [ "${SKIP_LIGANDMPNN:-0}" = "1" ]; then
  echo "[skip] 按 SKIP_LIGANDMPNN=1 跳过 LigandMPNN"
else
  echo "[4/4] 取 LigandMPNN 源码与权重（约 1 GB，视网速 3–15 分钟）..."
  if [ ! -d "$LMPNN_DIR" ]; then
    mkdir -p "$(dirname "$LMPNN_DIR")"
    git clone --depth 1 https://github.com/dauparas/LigandMPNN.git "$LMPNN_DIR" || exit 1
  else
    echo "[skip] 已存在 $LMPNN_DIR"
  fi
  ( cd "$LMPNN_DIR" && run_pip install -r requirements.txt ) || {
    echo "[warn] requirements.txt 安装失败 —— 常见原因见 01_ENVIRONMENT_setup.md §4（代理/版本）"; }
  if [ ! -d "$LMPNN_DIR/model_params" ] || [ -z "$(ls -A "$LMPNN_DIR/model_params" 2>/dev/null)" ]; then
    ( cd "$LMPNN_DIR" && bash get_model_params.sh "./model_params" ) || {
      echo "[warn] 权重下载失败 —— 重试：cd $LMPNN_DIR && bash get_model_params.sh ./model_params"; }
  fi
  ls -la "$LMPNN_DIR/model_params" 2>/dev/null | head
  ( cd "$LMPNN_DIR" && run_py run.py --help >/dev/null 2>&1 ) && echo "[check] LigandMPNN run.py --help 正常" \
    || echo "[warn] run.py --help 未通过 —— 先查 01_ENVIRONMENT_setup.md §4"
fi

echo
echo "=============================================="
echo " 完成。用法："
echo "   conda activate $ENV_NAME"
echo "   bash smoke_test/run_smoke_test.sh          # 环境自检"
echo "=============================================="
