# `out/` — 产物目录（要产出什么，与验收标准一一对应）

> 本目录默认只放本 README。跑完档 4 后，所有产物写在这里，**整个 `out/` 目录回传 Windows 侧**。

## 按档分层的产物清单

### 最小版（1.5–2.5 天）

| 文件 | 内容 | 验收要点 |
|---|---|---|
| `D0_site_table.csv` | 20 个口袋位点的三层分类 + 允许替代集合 + **跨束界面标注** | 跨束界面位点已单列（`is_cross_bundle=1`） |
| `D1_*/seqs/*.fa` | LigandMPNN 采样序列（**≥200 条**） | 每个温度一个目录；头部含 `ligand_confidence` |
| `D1_*/packed/*.pdb` | 这些序列的**全原子侧链**模型 | 用于 D2-S2 的真实距离判据 |
| `ligand_EGT.sdf` / `ligand_hercynine.sdf` | 修好键级的配体 | **EGT 的 SDF 里必须有 C=S 双键** |
| `D2_S2_geometry_check.csv` | 逐条变体的三段约束判定 | 每条都给出 PASS/FAIL 与依据距离 |
| `D4_variant_shortlist.csv` | **10 条 + WT 行** | 含 `is_WT` 与 **`WT_percentile`** |
| `D4_results_notes.md` | 结论 / 证据等级 / **未做清单** / 复跑命令 | 按 `docs/01–09` 的文件体例写 |

### 标准版（3–4 天，+ 反向筛选）

| 文件 | 内容 | 验收要点 |
|---|---|---|
| `backbone_EGT_with_hercynine_pose.pdb` | 把 hercynine 位姿套到 EGT 骨架上（同一骨架构象） | 报告 **口袋 Cα RMSD**（>1.5 Å 需注明） |
| `D2_S3_selectivity_bias.csv` | 每个变体的 `Δscore = score(EGT 上下文) − score(hercynine 上下文)` | ⛔ 必须写明"**这是设计偏置量，不是亲和力差**" |
| 更新后的 `D4_variant_shortlist.csv` | 加入 Δscore 列并参与排序 | — |

### 完整版（4–6 天，+ 结构复核）

| 文件 | 内容 | 验收要点 |
|---|---|---|
| `D3_refold_results.md` | top 3–5 条用 AF3 服务器（+EGT，CCD `LW8`）复折叠 | 通过率 + 每条的 pTM 与口袋是否仍成笼 |

---

## ⛔ 三条纪律（写进 `D4_results_notes.md` 的显著位置）

1. 全部产出为**计算级假设**。**禁止**出现"结合""亲和力""ΔΔG""识别机制已获得"。
2. **Δscore 是设计偏置量**，只用于排序，不换算能量。
3. 骨架本身是 AF3/Chai-1 **预测结构**（非实验结构）→ 设计结果的可信度上限 = 模型可信度上限。

## 回传时请一并带上

- `D4_results_notes.md` 里写清：**这一轮实际做到了哪一档**、停了在哪一步、下一步建议。
- 若某一步失败/跳过，**如实记录**（不要用推理填补）——Windows 侧要靠这份说明决定 PPT 上能写什么。
