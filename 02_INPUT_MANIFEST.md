# 02 · 输入资源清单（每个文件是什么、用在哪一步）

> 版本：2026-10-14。所有数字的**验证等级**与出处同时给出。等级：✅ 全文原句｜🖥 计算（可复跑）｜⚠️ 我方推断。

---

## 1. `data/` — 结构模型与判据

### 1.1 `data/af3_mfsT1_apo/`（AF3 服务器 · apo 单体）

| 文件 | 用途 | 关键数字 |
|---|---|---|
| `fold_mfst1_apo_stage2_model_0..4.cif` | **apo 骨架**（5 个模型）；D1 的备用骨架 / D3 复折叠的对照 | model_0：全局 pLDDT **84.7**，TM 区 **82.2–96.1** |
| `fold_mfst1_apo_stage2_summary_confidences_*.json` | 置信度（pTM / pLDDT / PAE） | **pTM 0.82**（model_4 = 0.73） |
| `fold_mfst1_apo_stage2_job_request.json` | 任务参数留档（seed 957414333） | — |
| ⚠️ `model0_plus_NO3_from4U4T.pdb` | **无效中间产物**：全局叠合把 NarK 的 NO3 转到本模型上，只留档备查，**勿用于任何分析** | — |

**来源**：档 2（`docs/07_stage2_results_stage3_plan.md` §1）。等级 🖥。

### 1.2 `data/chai1_mfsT1_EGT/` 与 `data/chai1_mfsT1_hercynine/`（Chai-1 网页 · 复合物）

每个目录各 5 个模型（rank_0–4）：

| 文件 | 用途 |
|---|---|
| `pred.rank_*.cif` | ⭐ **主输入**：蛋白链 A + 配体 `LIG2` 链 B（EGT 或 hercynine）。**这是设计的骨架 + 配体上下文** |
| `scores.rank_*.json` | 置信度：`ptm` / `iptm` / `per_chain_pair_iptm` / `has_inter_chain_clashes` |
| `pae.rank_*.npy` | PAE 矩阵（785 KB ×5）。用于界面置信度分析（可选） |

**关键数字**（🖥）：pTM **0.868**（EGT）/ **0.882**（hercynine）；ipTM **0.31 / 0.35**；配体平均 pLDDT 56–61。
⚠️ **ipTM 偏低是小分子单配体的常态 → 不得据此下"界面不可信"的结论**；判据一律看几何。

**来源**：档 3（`docs/08_stage3_results.md`）。

### 1.3 判据与距离表

| 文件 | 用途 | 说明 |
|---|---|---|
| `data/chai1_criteria_data.json` | 档 3 五条判据的原始计算（笼距离 / R21 盐桥 / S→T227 / hercynine 对照 / 模型一致性） | 🖥 可复跑 |
| `data/chai1_5model_distance_table.csv` | ⭐ 5 个模型 × 两个配体 × 三个关键距离 + pTM/ipTM | 由 `tools/` 同套逻辑重算，与档 3 报告**逐条一致** |
| `data/pocket_mapping_table.csv` | 档 1 的 NarK 口袋 8 残基 → mfsT1 映射（**Arg 对双失**：R89→Q57、R305→A260） | 🖥 |
| `data/pocket_contacts_selfcheck_7TXL_4U4T.csv` | 抽取流程的自检（7TXL 复现文献 10/10 识别残基） | 🖥 |

### 1.4 `data/evidence/` — 档 1–3 的**原始证据留档**（小体积，使包可自证）

| 文件 | 支撑哪个结论 | 类型 |
|---|---|---|
| `interproscan_mfsT1.tsv` | 家族判定：MFS 四库一致（SUPERFAMILY 5.4e-48 / cd06174 1.7e-23 / PF07690 2.1e-19 / Gene3D 双 bundle）；**PANTHER PTHR43124（PbuE 嘌呤外排泵）1.1e-13** | 🖥 原始输出 |
| `topology/`（5 个文件 + 拓扑图 PNG） | 拓扑：**12 TM、N-in/C-in**、中央胞质环 180–213（DeepTMHMM 2） | 🖥 原始输出 |
| `mfsT1_vs_NNP.aln` + `msa_input.fasta` | NNP 样仅为 **profile 级**（NS1/NS2 不保留：2/12、2/13） | 🖥 |
| `narK_P10903.fasta` / `narU_P37758.fasta` / `SLC22A4_Q9H015.fasta` | identity：vs NarK **16.4%** / NarU 18.4% / 人 SLC22A4 20.2% | 🖥 可复算 |
| `ding2024_SI_fulltext.txt` | **Table S3**：mfsT1 FC **3.45×**（p=0.033）；备选 gene3725/NGR_a01730 FC **4.74×** | ✅ SI 全文 |
| `europepmc_search_records.json` | 空白/占位复核的**检索式原文**（7 条） | 🖥 |
| `pubchem_egt.json` / `pubchem_hercynine.json` | 配体身份核验（CID 5351619 / 5459798） | 🖥 |

> 这些文件**不参与 D0–D4 的计算**，作用是：当某个上游结论被质疑时，**不必回到 Windows 侧**就能重算。

---

## 2. `ref/` — 参考结构与序列

| 文件 | 内容 | 为什么需要 |
|---|---|---|
| `mfsT1_WP_030133796.1.faa` | 靶点序列（428 aa） | 序列层判定 / 比对 |
| `mfsT1_CDS_NZ_JMDW01000005.1_319880-321166.fna` | CDS（1287 bp，起始 GTG） | 出处留档 |
| `7TXK.pdb` | ⭐ **EgtUC^CTT + ET 复合物，1.82 Å**（= PPT 引的那个分辨率） | 识别化学的**晶体参照**（三段式的"标准答案"） |
| `7TXL.pdb` | 野生型 SpEgtUC + ET，2.44 Å | 分辨率口径的另一半（**不要与 7TXK 混淆**） |
| `4U4T.pdb` / `4U4V.pdb` | NarK（inward-open） | 同源模板与构象参照 |
| `4U4W.pdb` | NarK（occluded） | 三态参照的第三态 |
| `4JR9.pdb` | NarK（按 2015 原文重判 = **部分闭合中间态**） | 构象态判读 |
| `ligand_smiles.txt` | ⭐ EGT（**Kekulé 硫酮式**）与 hercynine 的 SMILES + PubChem CID | `fix_ligand_bond_orders.py` 必需；**硫酮必须用 Kekulé 式**（芳香式在 RDKit 里无法凯库勒化） |

⚠️ **NarK 只有 occluded + inward-open；外向态至今无实验结构**（✅ 全文）→ 涉及构象态时不要写"两态"。

---

## 3. `docs/` — 档 1–3 的证据文件（只读参考）

| 文件 | 内容 |
|---|---|
| `01_family_and_topology.md` | Q1 家族（MFS 四库一致）+ Q4 拓扑（**12 TM，N-in/C-in**，中央环 180–213） |
| `02_NNP_motif_and_pocket_mapping.md` | ⭐ Q2/Q3：NS1/NS2 不保留（2/12、2/13）；Arg 对双失；7 位点×同号 TM 互证 |
| `02a_selfcheck_EgtUC_contacts.md` / `02b_NarK_pocket_residues.md` | 方法自检（7TXL 10/10）与 NarK 口袋 8 残基 |
| `03_homolog_structure_availability.md` | RCSB 0 命中；无 HOLD；NarK 仅 16.4% identity |
| `04_genomic_neighborhood_annotation.md` | ⭐ 不与 egt 基因簇共定位；Table S3（mfsT1 FC **3.45×**；备选 gene3725/NGR_a01730 FC **4.74×**） |
| `05_gap_and_precedent_check.md` | 外排空白成立；全球仅 Ding 2024 一篇 |
| `06_stage1_conclusion_and_gate.md` | 档 1 判定 **GO** |
| `07_stage2_results_stage3_plan.md` | ⭐ AF3 分析：口袋 24 衬里残基 / 6 芳香成笼 / R21+D107 / 半径 ~9.5 Å |
| `08_stage3_results.md` | ⭐ 五条判据全文 + **可检验突变预测 T227A / R21A / F111A** |
| `09_交接文档_...任务书.md` | 关键数字速查表（带等级与出处） |
| `EGT_key_numbers.csv` | 滴度 / 外排占比 / 亲和力 / 结构分辨率等（含来源 Zotero key） |

---

## 4. `tools/` — 脚本

| 脚本 | 状态 | 说明 |
|---|---|---|
| `common.py` | ✅ 已测 | 读 cif/pdb、写 PDB、Kabsch、配体图与官能团识别、TM 区间表 |
| `cif_to_pdb.py` | ✅ 已测 | mmCIF → PDB（自动识别 AF3 / Chai-1 两种列布局） |
| `fix_ligand_bond_orders.py` | ✅ 已测 | RDKit 修键级 → SDF（含 **C=S**） |
| `make_site_table.py` | ✅ 已测 | **D0 位点表**（三层 + 允许替代 + 跨束界面） |
| `render_structure.py` | ✅ 已测 | 纯 matplotlib 结构渲染器（出图用，不依赖 PyMOL） |

---

## 5. ⚠️ 使用这些数据时的三条硬约束

1. **不要把"mfsT1 转运 EGT"当事实**：原文自称 *putative*，无结合/转运实验，且不与 egt 基因簇共定位。严格表述 = "**唯一被命名、且被表型支持的 EGT 外排候选；其分子身份与转运机制均未被表征**"。
2. **结构是模型，不是实验结构**：mfsT1 无任何实验结构。所有几何都是 🖥+⚠️。
3. **`7TXK` 与 `7TXL` 不要混**：7TXK = EgtUC^CTT = **1.82 Å**；7TXL = 野生型 = **2.44 Å**（两者都是 ET 复合物）。
