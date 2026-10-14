# D3 · 升级完整版操作手册：AF3 服务器复折叠复核

> 目标：对短名单 top 3–5（+WT 对照）用 **AlphaFold 3 服务器**复折叠，
> 检验设计序列在独立预测器下是否仍维持口袋化学 → 把本档从「标准版」升级为「完整版」。
> ⛔ 产出仍为计算级假设（🖥/⚠️）；AF3 的 ipTM 预计仍低（小分子单配体常态），**不得据此下界面置信度结论**，判据一律看几何。
>
> **⚠️ 2026-10-15 修正（v2）**：
> 1. **JSON 方言已改为 `alphafoldserver`**（服务器只认 `proteinChain`/`ligand{ligand,count}` 等键；此前 v1 误用了本地版 AF3 的 `alphafold3` 方言，服务器报 `Unknown molecule type(s): protein`）。
> 2. **配体白名单风险**：AF3 服务器的 `ligand` 只接受 UI 白名单内的 CCD 码（官方 README 列为 ADP/ATP/AMP/GTP/GDP/FAD/NAD/NAP/NDP/HEM/HEC/PLM/OLA/MYR/CIT/CLA/CHL/BCL/BCB 等约 19 种），**`LW8`（EGT）大概率不在其列**，服务器也不接受 SMILES / 自定义 CCD。因此每个变体提供**两个 JSON**：先试含配体版，被拒则用 apo 版。

---

## 1. 任务文件（`out/D3_af3_jobs/`，每个变体两份）

| 含配体版（先试） | apo 后备版 | 变体 |
|---|---|---|
| `af3job_WT.json` | `af3job_WT_apo.json` | ⭐ WT 对照（**必跑**） |
| `af3job_T0p1_id74.json` | `af3job_T0p1_id74_apo.json` | rank 1：I50L;Q139H;Q223E;N346T;F350Y |
| `af3job_T0p1_id52.json` | `af3job_T0p1_id52_apo.json` | rank 2：I50L;Q139N;Q223E;N346T;F350Y |
| `af3job_T0p1_id19.json` | `af3job_T0p1_id19_apo.json` | rank 3：I50L;Q139N;Q223E;F350Y |
| `af3job_T0p2_id89.json` | `af3job_T0p2_id89_apo.json` | rank 4（做 top 5 时用） |
| `af3job_T0p2_id49.json` | `af3job_T0p2_id49_apo.json` | rank 5（做 top 5 时用） |

JSON 要点（已按服务器方言写好，勿手改）：顶层为**任务数组**；`proteinChain{sequence, count:1}`（428 aa，突变已逐位核验）；`dialect: "alphafoldserver"`；`version: 1`；`modelSeeds: []` = 服务器自选种子；含配体版加 `{"ligand": {"ligand": "CCD_LW8", "count": 1}}`。

## 2. 提交步骤（约 10 分钟 + 排队等待）

1. 打开 https://alphafoldserver.com/ → 登录（Google 账号）。
2. **Upload job → Upload JSON file**，**先只传 `af3job_WT.json` 试配体**：
   - ✅ 若服务器接受 `CCD_LW8` → 继续上传其余含配体版；
   - ❌ 若报"ligand not allowed / unsupported"之类 → **改用 `*_apo.json` 系列**（此时 D3 退化为"折叠 + 口袋骨架"复核，仍有效，见 §4 注）。
3. 服务器有**每日任务配额**（历史上 ~20 jobs/天，本轮 ≤6 个含配体版或 ≤6 个 apo 版均够）。
4. 每个任务返回 5 个模型。下载每个任务的 **full results zip**。

## 3. 回传与落盘

把每个 zip 解压到（目录名 = 任务名，保持 ASCII）：
```
out/D3_refold/mfsT1_WT_EGT/           ← WT 含配体（或 mfsT1_WT_apo/）
out/D3_refold/mfsT1_T0p1_id74_EGT/    ← rank 1
...
```

## 4. 分析（WSL 侧，回来后执行）

```bash
conda activate egt4
# ① cif → pdb（配体若存在会一并转出为 HETATM）
python tools/cif_to_pdb.py out/D3_refold/<task>/<job>_model_0.cif out/D3_refold/<task>/model_0.pdb
# ② 三段约束几何复核（与 D2-S2 同一判据、同一阈值；apo 版无配体 → 三判据记 NA，只做 ③′）
python tools/check_geometry.py tmp.csv out/D3_refold/<task>/model_0.pdb
# ③ 与 Chai-1 设计骨架比口袋 Cα RMSD（tools/common.py 提供 kabsch）
```

**逐条检查表**（写进 `out/D3_refold_results.md`）：

| 检查项 | 通过标准 | 含配体版 | apo 版 |
|---|---|---|---|
| pTM | ≥ 0.7（对照 apo 0.82） | ✔ | ✔ |
| 口袋成笼 | 6 个芳香残基环仍朝向腔体（判据① PASS） | ✔ | ✔（改为检查 6 个笼位残基侧链取向） |
| R21 盐桥 | 判据② PASS（≤4.5 Å） | ✔ | NA |
| T227 读硫 | 判据③ PASS（≤4.2 Å） | ✔ | NA |
| 口袋 Cα RMSD vs Chai-1 骨架 | ≤ 2 Å（跨预测器放宽） | ✔ | ✔ |
| 位姿三段式 | 与 Chai-1 方向一致 | ✔ | NA |

> **apo 版的定位**：若 LW8 被拒，apo 复折叠仍能回答"设计序列是否仍折叠成带芳香笼的 MFS"（pTM + 笼位取向 + 口袋 RMSD），但**不能复核配体三段约束** —— 这部分在结果说明里如实降级。
> **可选替代**（若坚持要配体复合物复核）：本地跑 **Boltz-2**（开源、接受 SMILES）或再用 **Chai-1**；注意 Chai-1 与设计骨架同源，独立性弱于 AF3，Boltz-2 则是独立预测器。本机 8 GB 显存跑 428 aa 需实测。

**通过率定义**：一条变体 = 5 个模型中 ≥4 个同时满足上表该版本全部适用项，记"通过复折叠"。

## 5. 结果说明模板（`out/D3_refold_results.md`）

```markdown
# D3 · 结构复核结果
- 任务提交日期 / 服务器版本 / 配体：LW8 是否被接受（若否，apo 版并注明降级）
- 每变体：pTM 区间 / ipTM（注明：低 ipTM 为小分子单配体常态，不作界面置信度依据）
  / 各检查项通过模型数 / 口袋 Cα RMSD / 通过与否
- 通过率：x / N（+WT 对照表现）
- 未做与边界：（如实记录）
```

## 6. 完成定义（完整版验收）

- `out/D3_refold_results.md` 存在，含每条变体的通过判定与通过率；
- 结果回写 `out/D4_results_notes.md`（把"未做清单"第 1 条划掉，补 D3 结论；若 LW8 被拒，注明复核范围为 apo 级）；
- 若某变体未通过：**如实记录**，在短名单 CSV 中加注（不删除、不粉饰）。

---

## 附：v1 → v2 变更记录

| 日期 | 变更 |
|---|---|
| 2026-10-15 v1 | 初版（误用本地 AF3 方言 `alphafold3`，服务器报 `Unknown molecule type(s): protein`） |
| 2026-10-15 v2 | 改用 `alphafoldserver` 方言（`proteinChain`/`ligand{ligand,count}`）；新增 LW8 白名单风险说明与 `*_apo.json` 后备系列 |
