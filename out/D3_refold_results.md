# D3 · 结构复核结果（Boltz-2 本地推理，替代 AF3 服务器）

> 日期：2026-10-15 ｜ 执行：WSL 侧 ｜ ⛔ 全部为计算级结果（🖥/⚠️）
> **路线变更说明**：AF3 服务器不接受 EGT（CCD `LW8` 不在其配体白名单内，且不支持 SMILES），
> 改用**本地 Boltz-2**（v2.2.1，独立预测器，SMILES 输入 EGT 硫酮式）。由 LigandMPNN/Chai-1 之外的
> 第三个模型复核，独立性满足 D3 原意。

## 0. 可行性实测（8 GB 显存约束）

| 项 | 实测值 |
|---|---|
| 环境 | conda `boltz2`，torch 2.14.1+cu130（sm_120 可用），需 `--no_kernels`（cuequivariance ops 与 cu13 不兼容，回退纯 PyTorch 路径） |
| 参数 | recycling 3 / diffusion_samples 5 / sampling 200 步 / **max_msa_seqs 512**（AF3 同款 MSA，原深 15672） |
| **峰值显存** | **5.6 GB / 8 GB** ✅ |
| 单目标耗时 | ~4–5 分钟（5 个扩散样本） |

## 1. 逐模型结果（`D3_boltz2_per_model.csv` 全量）

判据与 D2-S2 完全同源（`tools/check_geometry.py`，阈值不变）。口袋 Cα RMSD = Boltz-2 模型 vs Chai-1 rank_0 设计骨架（15 个口袋位点，Kabsch）。

### WT + EGT（对照）
pTM 0.90–0.92；口袋 Cα RMSD 0.86–1.77 Å；①笼 5/5 PASS（4.5–4.8 Å）、②R21 4/5 PASS（2.1–3.5 Å）、**③S→T227 0/5（8.3–9.5 Å）**。

### 变体（top 5）

| 变体 | 突变 | pTM | 口袋Cα RMSD (Å) | ①②通过 | 三段全过 |
|---|---|---|---|---|---|
| rank1 T0.1/id74 | I50L;Q139H;Q223E;N346T;F350Y | 0.81–0.87 | 1.69–4.16 | 2/5 | **1/5**（m0：笼 4.42 / R21 2.24 / S→T227 **3.18**） |
| rank2 T0.1/id52 | I50L;Q139N;Q223E;N346T;F350Y | 0.77–0.84 | 2.07–3.43 | 1/5 | 0/5 |
| rank3 T0.1/id19 | I50L;Q139N;Q223E;F350Y | 0.74–0.80 | 1.88–2.87 | 3/5 | 0/5 |
| rank4 T0.2/id89 | I50L;Q139T;Q223E;N346T;F350Y | 0.75–0.80 | 2.07–2.75 | 2/5 | 0/5 |
| rank5 T0.2/id49 | I50L;Q139H;Q223E;F350Y | 0.76–0.84 | 2.39–3.14 | 4/5 | 0/5（m4 判据① NA：配体图识别失败，记为缺失） |

## 2. 结论

1. **折叠层面全部通过** 🖥：WT 与 5 个变体在 Boltz-2 下均折叠为高置信 MFS 结构（pTM 0.74–0.92），口袋区与 Chai-1 设计骨架偏差 0.9–4.2 Å（WT 最小）——**设计序列没有破坏折叠**。
2. **两段识别化学跨预测器复现** 🖥：① 三甲铵→芳香笼与 ② 羧基→R21 盐桥在 WT（5/5、4/5）与变体中大部分模型重现。
3. **⭐ 第三段"读硫"（S→T227 CH···S）未被 Boltz-2 复现** ⚠️：30 个模型（含 WT）中仅 rank1-m0 一例通过（3.18 Å）；其余 8–17 Å。即咪唑硫酮端的朝向在两个预测器间不一致 —— **"读硫"几何证据是预测器依赖的**，其假设强度应下调；T227A 湿实验的判别价值反而因此更高（它是唯一能裁决两个预测器分歧的手段）。
4. 变体的 pTM 与口袋 RMSD 整体略逊于 WT（WT 本身就是 Boltz-2 最"顺手"的序列）——短名单排序结论（基于 LigandMPNN 似然）与 Boltz-2 的偏好**不一致**，这是两个模型先验差异的体现，如实记录。

## 3. 通过率（按手册 §4 定义：≥4/5 模型通过全部适用项）

- 含配体完整三段标准：**0/5 变体通过**（WT 对照也不通过③）。
- 若按"折叠 + 两段"标准（pTM≥0.7 且 ①② ≥4/5）：WT 通过；变体中 T0.2/id49 通过，其余 1–3/5 部分通过。

## 4. 未做与边界

1. AF3 服务器路径未走（LW8 白名单限制；apo JSON 已备在 `out/D3_af3_jobs/`，如日后想补折叠级对照可用）。
2. Boltz-2 用 `--no_kernels`（精度无影响，仅速度）；MSA 截到 512 条（显存约束），与满深 MSA 的差异未评估。
3. 单 seed / 单组参数；③ 判据在 rank1-m0 的个例通过不构成"复现"（1/30）。
4. hercynine 未做 Boltz-2 对照跑（可做：同样的 YAML 换 SMILES）。
5. Chai-1 的 ipTM 0.31 vs Boltz-2 的 ligand_iptm 0.95：两预测器的界面置信度口径不同，**均不作为界面证据**，判据一律看几何（本报告即如此）。

## 5. 产物

| 文件 | 内容 |
|---|---|
| `out/D3_boltz2/yamls/boltz_*_EGT.yaml` | 6 个输入（WT + top5，EGT 硫酮 SMILES + AF3 MSA） |
| `out/D3_boltz2/full/boltz_results_*/predictions/` | 每目标 5 个 mmCIF + confidence json + pae/plddt npz |
| `out/D3_boltz2/D3_boltz2_per_model.csv` | ⭐ 30 行逐模型几何 + 置信度 + 口袋 RMSD |
| `out/D3_boltz2/full_run.log` `variants_run.log` `vram_monitor*.log` | 运行与显存记录 |

## 6. 复跑

```bash
conda activate boltz2   # pip install boltz cuequivariance-torch（后者可省，需 --no_kernels）
boltz predict out/D3_boltz2/yamls/boltz_WT_EGT.yaml --out_dir out/D3_boltz2/full \
  --recycling_steps 3 --diffusion_samples 5 --sampling_steps 200 \
  --max_msa_seqs 512 --accelerator gpu --devices 1 --output_format mmcif --no_kernels
# 逐模型分析：tools/cif_to_pdb.py + tools/check_geometry.py + common.kabsch（见 D4_results_notes.md §6）
```
