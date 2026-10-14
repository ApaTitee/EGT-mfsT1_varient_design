# 档 4 · 结果说明：mfsT1 口袋的变体设计与计算筛选

> 日期：2026-10-15（WSL 侧执行）｜ 实际完成档位：**完整版**（D0 + D1 + D2-S1/S2/S3 + D3 + D4；D3 由本地 Boltz-2 完成，见 `D3_refold_results.md`）
> 体例按 `docs/01–09`：结论 / 证据等级 / 未做清单 / 复跑方式。
> ⛔ 本档全部产出为**计算级假设**（🖥/⚠️）。Δscore 是**设计偏置量，不是亲和力差**，只用于在序列空间排序，不可换算 ΔΔG。

---

## 1. 结论（一句话）

在 Chai-1 EGT 复合物骨架上，用 LigandMPNN（配体原子上下文 + 位点级氨基酸约束）采样 300 条口袋变体，经三段约束几何复核与 hercynine 反向筛选，得到 **10 条去重短名单**；**WT 天然口袋在全部变体中处于 score_EGT 分布的 2.7% 分位、Δscore 分布的 0.0% 分位** —— 按模型口径，"天然口袋离（模型认为的）最优还有明显距离"（🖥）。

## 2. 三个可上页数字

| 数字 | 值 | 含义 |
|---|---|---|
| ① 设计臂 N | **10 条**（79 条去重指纹中综合排序前 10；采样 300 条 = 3 温度 × 100） | 湿实验设计臂的规模 |
| ② top 10 的 Δscore 区间 | **+0.279 ~ +0.415**（WT = +0.054） | 全部短名单成员对"有硫上下文"的偏置都强于 WT（🖥 偏置量） |
| ③ WT 分位 | **score_EGT 2.7% 分位 ｜ Δscore 0.0% 分位**（300 条变体中） | "天然口袋离最优还有多远"的定量（模型口径） |

## 3. 分档结果

### D0 位点定义 —— 完成
- `out/D0_site_table.csv`：20 位点 = 3 核心锚点（R21/T227/F111，fixed）+ 12 次级衬里（shell）+ 1 观察位（D107，约束 DENQ 参与采样）+ 4 结构 Gly（fixed）。
- **跨束界面位点已单列**（`is_cross_bundle` 列）：I50、G135、Q139、F326、N346、F350 = 6 个。
- 关键距离复算与档 3 逐条一致（T227–S 3.34 / R21–羧基 3.30 / F111–N 3.12）🖥。

### D1 序列设计 —— 完成
- 环境：torch 2.11.0+cu128（sm_120 实测 GPU 可用）；LigandMPNN 官方 main；捆绑 openfold 的 `np.int/np.bool` 已打补丁（numpy 1.26 兼容）。
- **配体键级已修正**：`out/ligand_EGT.sdf` 含 C=S 双键（hercynine 无，正确）；验收用 🖥。
- 采样：T=0.1/0.2/0.3 各 100 条 + 侧链全原子 packing（配体上下文开启）。**每个设计位点用 `omit_AA_per_residue` 约束在 D0 允许集合内**（`out/D1_inputs/omit_AA_per_residue.json`）。
- 多样性：去重口袋指纹 T0.1=12 / T0.2=36 / T0.3=68，合计 79 条唯一指纹。

### D2 计算筛选 —— 完成（S1/S2/S3）
- **S1**（`out/D2_S1_sequence_scores.csv`）：300 条的 overall/ligand_confidence。
- **S2**（`out/D2_S2_geometry_check.csv`）：301 条（含 WT）逐条三段约束判定，WT 复现档 3 数字。
  ⚠️ **诚实记录：S2 本轮无区分度（301/301 全过）**。原因有二（均为方法固有）：① 允许替代集合本身保守（笼位 FYWH 等），② 所有变体共用同一配体位姿与骨架。S2 在本轮的角色是**一致性检查**（packing 未破坏三段几何），不是筛选器。任务书 §D2 的"简化路线"警示与实际一致。
- **S3**（`out/D2_S3_selectivity_bias.csv`）：hercynine 位姿经口袋 Cα Kabsch 搬到 EGT 骨架（**RMSD 0.69 Å**，n=15，<1.5 Å 阈值，`out/D2_S3_align_residuals.csv` 逐位点残差），同一序列在两种配体上下文下用 score.py 各自打分（自回归、10 解码序、seed 111）。
  - score 定义：13 个口袋位点上实际氨基酸的 log-prob（对解码序与位点取平均）。
  - **Δscore = score(EGT 上下文) − score(hercynine 上下文) > 0 表示更偏好"有硫"上下文。⛔ 它是设计偏置量，不是亲和力差，不得换算 ΔΔG。**
  - 结果：300 条变体 Δscore 全部为正（+0.11 ~ +0.42），WT = +0.054 为全场最低。

### D4 汇总 —— 完成
- `out/D4_variant_shortlist.csv`：**10 条 + WT 行**（按 z(ligand_confidence)+z(Δscore) 综合排序、指纹去重）；WT 分位已写入。
- `out/D4_expected_phenotype.csv`：每条变体的预期表型与判别实验（沿用 `docs/08` §3 框架：读硫端 / 羧基锚 / 笼）。
- **复现设计信号**：I50L（300/300）、F350Y（≈全部）、Q223E、Q139H/N/T、N346T/D 高频出现 —— 模型稳定偏好这些替代（🖥，其化学含义待检验）。

## 4. 未做清单（诚实清单）

1. ~~**D3 结构复核未做**~~ → **已完成**（2026-10-15，本地 Boltz-2，WT + top5 × 5 模型）：折叠全通过；两段化学（笼 + R21）跨预测器复现；**第三段"读硫"未被复现（1/30）**，详见 `D3_refold_results.md`。AF3 服务器路径因 LW8 白名单限制未走。
2. **Δscore 未做 seed 稳健性**：全部打分用 seed=111；不同随机解码序带来的涨落未定量。
3. **S2 为约束类型级 + packing 距离级，但没有独立侧链采样方法（如 Rosetta）交叉验证**。
4. **空腔体积、对接、自由能类计算均未做**（沿用档 3 的决定：ipTM 0.31/0.35 已说明置信度边界）。
5. **AF3 骨架（apo）第二轮对照设计未做**（任务书列为可选项）。
6. Δscore 只覆盖 EGT vs hercynine 一对；终点"选择性比"里的另外两个对照分子（甘氨酸甜菜碱、组氨酸）未进入打分。

## 5. 边界（沿用任务书 §5，仍然成立）

- 骨架是 AF3/Chai-1 **预测结构**，设计可信度上限 = 模型可信度上限。
- "读硫"只有几何证据；唯一证伪途径是湿实验（T227A 对 EGT 与 hercynine 无差别 → 假设被证伪）。
- 改的是口袋化学，**不承诺转运速率**；mfsT1 底物身份仍是 putative。

## 6. 复跑方式

```bash
conda activate egt4
bash smoke_test/run_smoke_test.sh                 # 环境自检
cp smoke_test/work/D0.csv out/D0_site_table.csv   # D0（或与本次产物直接比对）
bash tools/run_ligandmpnn.sh                      # D1（3 温度 × 100，约 40 分钟/GPU）
python tools/check_geometry.py ...                # D2-S2（本轮由会话内脚本批量跑 301 条）
python tools/align_ligand_pose.py data/chai1_mfsT1_EGT/pred.rank_0.cif \
  data/chai1_mfsT1_hercynine/pred.rank_0.cif out/backbone_EGT_with_hercynine_pose.pdb \
  --sites out/D0_site_table.csv --resid out/D2_S3_align_residuals.csv   # T1
python tools/make_s3_inputs.py                    # S3 输入（301 × 2 套上下文）
cd ~/work/LigandMPNN && python score.py --model_type ligand_mpnn --seed 111 \
  --autoregressive_score 1 --use_sequence 1 --batch_size 1 --number_of_batches 10 \
  --ligand_mpnn_use_atom_context 1 --ligand_mpnn_use_side_chain_context 1 \
  --ligand_mpnn_cutoff_for_score 6.0 \
  --pdb_path_multi <pkg>/out/D2_S3/list_{EGT,hercynine}.json \
  --out_folder <pkg>/out/D2_S3/scores_{EGT,hercynine}                    # S3 打分（两次）
python tools/parse_s3_scores.py                   # Δscore 表
python tools/analyze_designs.py                   # S1 + 短名单 + WT 分位
```

## 7. 产物清单（`out/`，回传 Windows 侧）

| 文件 | 内容 |
|---|---|
| `D0_site_table.csv` | 20 位点三层分类 + 跨束界面标注 |
| `mfsT1_EGT.pdb` / `ligand_EGT.sdf` / `ligand_hercynine.sdf` | 设计输入骨架 + 修好键级的配体（EGT 含 C=S） |
| `D1_inputs/omit_AA_per_residue.json` | 位点级氨基酸约束 |
| `D1_EGT_T{0.1,0.2,0.3}/`（seqs/packed/backbones/stats） | 300 条采样 + 全原子模型 |
| `D2_S1_sequence_scores.csv` | S1 内在打分 |
| `D2_S2_geometry_check.csv` | 301 条三段约束判定 |
| `backbone_EGT_with_hercynine_pose.pdb` + `D2_S3_align_residuals.csv` | T1 叠合产物（RMSD 0.69 Å） |
| `D2_S3/`（pdbs_*, scores_*, list_*.json, *.log） | 反向筛选原始数据 |
| `D2_S3_selectivity_bias.csv` | 301 条 Δscore |
| `D4_variant_shortlist.csv` | ⭐ 10 条 + WT 行（含 WT 分位） |
| `D4_expected_phenotype.csv` | 预期表型 + 判别实验 |
| `D1_run.log` | D1 运行日志 |
| `D3_boltz2/` | D3 完整版：Boltz-2 输入 YAML、30 个模型、逐模型几何表（`D3_boltz2_per_model.csv`）、显存记录 |
| `D3_refold_results.md` | ⭐ D3 结构复核报告（含"读硫段未被 Boltz-2 复现"的关键发现） |
| `D3_af3_jobs/` + `D3_full_version_manual.md` | AF3 服务器路径的任务 JSON（apo 后备）与手册（v2） |

**本轮实际做到：完整版**（D3 由本地 Boltz-2 完成，替代 AF3 服务器——后者不接受 CCD `LW8`）。下一步建议：① Δscore 换 2–3 个 seed 复查稳健性；② Boltz-2 补 hercynine 对照跑；③ 短名单第 1–3 条可先做湿实验引物设计（注意：D3 显示"读硫"段为预测器依赖，T227A 实验的裁决价值更高）。
