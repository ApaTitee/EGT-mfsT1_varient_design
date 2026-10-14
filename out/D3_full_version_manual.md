# D3 · 升级完整版操作手册：AF3 服务器复折叠复核

> 目标：对短名单 top 3–5（+WT 对照）用 **AlphaFold 3 服务器**复折叠（含 EGT 配体），
> 检验设计序列在独立预测器下是否仍维持口袋化学 → 把本档从「标准版」升级为「完整版」。
> ⛔ 产出仍为计算级假设（🖥/⚠️）；AF3 的 ipTM 预计仍低（小分子单配体常态），**不得据此下界面置信度结论**，判据一律看几何。

---

## 1. 任务文件（已生成，`out/D3_af3_jobs/`）

| 文件 | 任务名 | 说明 |
|---|---|---|
| `af3_job_WT.json` | `mfsT1_WT_EGT` | ⭐ WT 对照（**必跑**，是"变体是否退化"的参照） |
| `af3_job_rank1_T0p1_id74.json` | `mfsT1_T0p1_id74_EGT` | rank 1：I50L;Q139H;Q223E;N346T;F350Y |
| `af3_job_rank2_T0p1_id52.json` | `mfsT1_T0p1_id52_EGT` | rank 2：I50L;Q139N;Q223E;N346T;F350Y |
| `af3_job_rank3_T0p1_id19.json` | `mfsT1_T0p1_id19_EGT` | rank 3：I50L;Q139N;Q223E;F350Y |
| `af3_job_rank4_T0p2_id89.json` | `mfsT1_T0p2_id89_EGT` | rank 4（做 top 5 时用） |
| `af3_job_rank5_T0p2_id49.json` | `mfsT1_T0p2_id49_EGT` | rank 5（做 top 5 时用） |

JSON 要点（已按 AF3 服务器格式写好，勿手改）：
- 顶层是**任务数组** `[ {...} ]`；
- 蛋白链 `id: ["A"]`（428 aa，突变已逐位核验）；配体链 `id: ["B"]`，**CCD 码 `LW8`**（= EGT；AF3 服务器无 SMILES 入口，只能用 19 种白名单外的 CCD 码时走 `ccdCodes`）；
- `modelSeeds: []` = 服务器自选种子。

## 2. 提交步骤（约 10 分钟 + 排队等待）

1. 打开 https://alphafoldserver.com/ → 登录（Google 账号）。
2. **Upload job**（或 Add job → "Upload JSON file"），逐个上传上表 JSON。
   - ⚠️ 服务器有**每日任务配额**（历史上 ~20 jobs/天，够本轮 4–6 个任务）。
   - 也可以不走 JSON：手动粘贴序列 + 配体处选 CCD 码 `LW8`，效果相同。
3. 每个任务返回 5 个模型。下载每个任务的 **full results zip**。

## 3. 回传与落盘

把每个 zip 解压到：
```
out/D3_refold/mfsT1_WT_EGT/           ← WT
out/D3_refold/mfsT1_T0p1_id74_EGT/    ← rank 1
...（每个任务一个目录，目录名 = 任务名）
```
（目录名保持 ASCII；zip 内的 `*_model_*.cif` / `*_summary_confidences_*.json` 保持原名。）

## 4. 分析（WSL 侧，回来后执行）

对每个任务、每个模型（至少 model_0）：

```bash
conda activate egt4
# ① cif → pdb（配体一并转出）
python tools/cif_to_pdb.py out/D3_refold/<task>/<job>_model_0.cif out/D3_refold/<task>/model_0.pdb
# ② 三段约束几何复核（与 D2-S2 同一判据、同一阈值）
python tools/check_geometry.py tmp.csv out/D3_refold/<task>/model_0.pdb
# ③ 与 Chai-1 设计骨架比口袋 Cα RMSD（Kabsch，tools/common.py 提供 kabsch）
```

**逐条检查表**（写进 `D3_refold_results.md`）：

| 检查项 | 通过标准 | 含义 |
|---|---|---|
| pTM | ≥ 0.7（对照 apo 0.82） | 折叠整体未崩 |
| 口袋成笼 | 6 个芳香残基环仍朝向腔体（判据① PASS） | 设计位点未破坏笼 |
| R21 盐桥 | 判据② PASS（≤4.5 Å） | 羧基锚仍在 |
| T227 读硫 | 判据③ PASS（≤4.2 Å） | 硫端 CH···S 仍在 |
| 口袋 Cα RMSD vs Chai-1 骨架 | ≤ 2 Å（参考值；跨预测器放宽） | 口袋骨架保守 |
| 位姿三段式 | 三甲铵→笼 / 羧基→R21 / S→T227 方向与 Chai-1 一致 | 识别化学雏形跨预测器复现 |

**通过率定义**：一条变体 = 5 个模型中 ≥4 个同时满足上表全部几何项，记该变体"通过复折叠"。

## 5. 结果说明模板（`out/D3_refold_results.md`）

```markdown
# D3 · 结构复核结果
- 任务提交日期 / 服务器版本：
- 每变体：pTM 区间 / ipTM（注明：低 ipTM 为小分子单配体常态，不作界面置信度依据）
  / 判据①②③通过模型数 / 口袋 Cα RMSD / 通过与否
- 通过率：x / N（+WT 对照表现）
- 未做与边界：（如实记录）
```

## 6. 完成定义（完整版验收）

- `out/D3_refold_results.md` 存在，含每条变体的通过判定与通过率；
- 结果回写 `out/D4_results_notes.md`（把"未做清单"第 1 条划掉，补 D3 结论）；
- 若某变体未通过：**如实记录**，在短名单 CSV 中加注（不删除、不粉饰）。
