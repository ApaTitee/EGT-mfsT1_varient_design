# D3 · 结构复核结果（Boltz-2 本地推理，替代 AF3 服务器）

> 日期：2026-10-15 ｜ 执行：WSL 侧 ｜ ⛔ 全部为计算级结果（🖥/⚠️）
> **路线变更说明**：AF3 服务器不接受 EGT（CCD `LW8` 不在其配体白名单内，且不支持 SMILES），
> 改用**本地 Boltz-2**（v2.2.1，独立预测器，SMILES 输入 EGT 硫酮式）。由 LigandMPNN/Chai-1 之外的
> 第三个模型复核，独立性满足 D3 原意。
> **v2 补全**：hercynine 对照已跑（§1b），`tools/common.py` 的咪唑 C2 识别 bug 已修（见 §4 附注）。

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

## 1b. hercynine 对照（`D3_boltz2_hercynine_per_model.csv` 全量）

同一批序列 + hercynine（SMILES 咪唑式），判据③自动切换为 ③′（咪唑 C2 → T227 供体，同阈值 4.2 Å）。

| 变体 | ①②通过 | ③′通过 | C2→T227 (Å) | pTM | 口袋Cα RMSD (Å) |
|---|---|---|---|---|---|
| WT | 4/5 | 0/5 | 7.70–8.54 | 0.89–0.93 | 0.98–1.82 |
| rank1 T0.1/id74 | 1/5 | 0/5 | 11.54–15.61 | 0.72–0.79 | 2.22–4.17 |
| rank2 T0.1/id52 | 4/5 | 0/5 | 9.43–12.95 | 0.84–0.91 | 2.85–4.02 |
| rank3 T0.1/id19 | 4/5 | 0/5 | 6.28–16.71 | 0.68–0.79 | 2.11–3.88 |
| rank4 T0.2/id89 | 3/5 | 0/5 | 8.36–16.16 | 0.74–0.85 | 2.57–4.12 |
| rank5 T0.2/id49 | 5/5 | 0/5 | 11.47–15.66 | 0.73–0.77 | 2.64–3.80 |

### ⭐ EGT vs hercynine 逐模型对照（硫端探针 → T227 距离，Å）

| 变体 | EGT（S→T227） | hercynine（C2→T227） |
|---|---|---|
| WT | 8.26 / 9.49 / 9.08 / 9.45 / 8.51 | 8.54 / 7.70 / 7.90 / 7.78 / 8.42 |
| rank1 id74 | **3.18** / 12.69 / 5.42 / 12.42 / 8.38 | 12.22 / 11.54 / 14.68 / 14.33 / 15.61 |
| rank2 id52 | 16.33 / 13.70 / 12.47 / 17.29 / 12.38 | 12.95 / 11.80 / 11.27 / 11.87 / 9.43 |
| rank3 id19 | 11.66 / 12.02 / 10.68 / 10.72 / 13.36 | 16.71 / 14.06 / 10.78 / 9.42 / 6.28 |
| rank4 id89 | 5.89 / 5.33 / 8.65 / 5.41 / 7.28 | 9.16 / 10.16 / 8.36 / 16.16 / 10.21 |
| rank5 id49 | 12.67 / 11.01 / 10.84 / 13.00 / 10.63 | 11.95 / 15.66 / 11.47 / 15.21 / 11.79 |

**结论：Boltz-2 对“有没有那个硫”不敏感**——WT 的 EGT 与 hercynine 硫端距离分布几乎重合（8.3–9.5 vs 7.7–8.5 Å），
与 Chai-1 的鲜明对比（3.34 vs 9.21 Å）完全不同。变体层面也无一致的“EGT 更近”模式
（id74 唯一的一例 3.18 Å 出现在 EGT 侧，但 id19 最近的一例反而在 hercynine 侧 6.28 Å）。

## 2. 结论

1. **折叠层面全部通过** 🖥：WT 与 5 个变体在 Boltz-2 下均折叠为高置信 MFS 结构（pTM 0.74–0.92），口袋区与 Chai-1 设计骨架偏差 0.9–4.2 Å（WT 最小）——**设计序列没有破坏折叠**。
2. **两段识别化学跨预测器复现** 🖥：① 三甲铵→芳香笼与 ② 羧基→R21 盐桥在 WT（5/5、4/5）与变体中大部分模型重现。
3. **⭐ 第三段"读硫"（S→T227 CH···S）未被 Boltz-2 复现，且 Boltz-2 对硫不敏感** ⚠️：EGT 侧 30 个模型中仅 rank1-m0 一例通过（3.18 Å）；hercynine 对照补齐后显示两个配体上下文的硫端距离分布无差别（WT：8.3–9.5 vs 7.7–8.5 Å）——即 **Chai-1 观察到的"硫端朝向 T227"在独立预测器下既不重现、也不区分硫的有无**。“读硫”假设目前只剩 Chai-1 单侧几何证据，强度应下调；T227A 湿实验是唯一裁决手段，价值反而更高。
4. 变体的 pTM 与口袋 RMSD 整体略逊于 WT（WT 本身就是 Boltz-2 最"顺手"的序列）——短名单排序结论（基于 LigandMPNN 似然）与 Boltz-2 的偏好**不一致**，这是两个模型先验差异的体现，如实记录。

## 3. 通过率（按手册 §4 定义：≥4/5 模型通过全部适用项）

- 含配体完整三段标准：**0/5 变体通过**（EGT 与 hercynine 两侧均如此；WT 对照也不过③）。
- 若按"折叠 + 两段"标准（pTM≥0.7 且 ①② ≥4/5）：EGT 侧 WT 与 T0.2/id49 通过；hercynine 侧 WT、T0.1/id52、T0.1/id19、T0.2/id49 通过。

## 4. 未做与边界

1. AF3 服务器路径未走（LW8 白名单限制；apo JSON 已备在 `out/D3_af3_jobs/`，如日后想补折叠级对照可用）。
2. Boltz-2 用 `--no_kernels`（精度无影响，仅速度）；MSA 截到 512 条（显存约束），与满深 MSA 的差异未评估。
3. 单 seed / 单组参数；③ 判据在 rank1-m0 的个例通过不构成"复现"（1/30）。
4. ~~hercynine 未做 Boltz-2 对照跑~~ → **已补**（§1b；6 目标 × 5 模型，同参数）。
5. Chai-1 的 ipTM 0.31 vs Boltz-2 的 ligand_iptm 0.95：两预测器的界面置信度口径不同，**均不作为界面证据**，判据一律看几何（本报告即如此）。

### 附：工具修复记录（v2）

`tools/common.py` 的 `lig_sites()` 原实现把咪唑 C2 检测误嵌套在 N 原子分支内（永远触发不到）。
本次修复为"C 原子且重原子邻居恰为 2 个 N"。回归验证：① EGT 侧全部已交付数字不变（WT 4.19/3.30/3.34）；
② Chai-1 hercynine rank_0 的 C2→T227 = **9.21 Å**，与档 3 报告逐字一致。

## 5. 产物

| 文件 | 内容 |
|---|---|
| `out/D3_boltz2/yamls/boltz_*_{EGT,hercynine}.yaml` | 12 个输入（WT + top5 × 两配体，AF3 MSA） |
| `out/D3_boltz2/full/boltz_results_*/predictions/` | 每目标 5 个 mmCIF + confidence json + pae/plddt npz（共 60 模型） |
| `out/D3_boltz2/D3_boltz2_per_model.csv` | ⭐ EGT 侧 30 行逐模型几何 + 置信度 + 口袋 RMSD |
| `out/D3_boltz2/D3_boltz2_hercynine_per_model.csv` | ⭐ hercynine 侧 30 行（判据③′ = 咪唑 C2） |
| `out/D3_boltz2/full_run.log` `variants_run.log` `vram_monitor*.log` | 运行与显存记录 |

## 6. 复跑

```bash
conda activate boltz2   # pip install boltz cuequivariance-torch（后者可省，需 --no_kernels）
boltz predict out/D3_boltz2/yamls/boltz_WT_EGT.yaml --out_dir out/D3_boltz2/full \
  --recycling_steps 3 --diffusion_samples 5 --sampling_steps 200 \
  --max_msa_seqs 512 --accelerator gpu --devices 1 --output_format mmcif --no_kernels
# 逐模型分析：tools/cif_to_pdb.py + tools/check_geometry.py + common.kabsch（见 D4_results_notes.md §6）
```
