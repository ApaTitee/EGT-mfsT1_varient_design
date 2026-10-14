# 03 · 待写工具规格（需要时再写，规格已定死）

> 这里只列**本包未提供**、但 D1–D4 会用到的脚本。规格写得足够具体，照写即可；
> 也可由 WSL 侧会话自行实现，只要**输出字段一致**（下游脚本与 PPT 都按这些字段读）。

---

## T1 · `align_ligand_pose.py`（D2-S3 的前置，⚠️ 方法论关键）

**问题**：hercynine 的复合物是**另一次独立的 Chai-1 预测**，骨架构象与 EGT 的那次略有差异。
直接比较"同一序列在两次预测里的分数"会把**骨架差异**混进 Δscore。

**做法**（把 hercynine 的配体位姿搬到 EGT 的骨架上）：
1. 读 `data/chai1_mfsT1_EGT/pred.rank_0.cif`（参考，含 EGT）与
   `data/chai1_mfsT1_hercynine/pred.rank_0.cif`（供体，含 hercynine）。
2. 取两者**共有的口袋残基 Cα**（默认 = `smoke_test/expected/D0_site_table.csv` 里 `layer in (core, shell)` 的位点），
   用 Kabsch 把**供体**叠合到**参考**上（`common.kabsch` 已提供）。
3. 把供体的**配体坐标**用同一变换搬到参考坐标系。
4. 写出 `out/backbone_EGT_with_hercynine_pose.pdb` = 参考蛋白 + 变换后的 hercynine 配体（resname `LIG`，链 B）。
5. **打印口袋 Cα RMSD**（这一步的"可信度指标"；若 > 1.5 Å 要在结果说明里注明）。

**输出**：
```csv
resnum,resname,ca_dist_after_fit   # 逐位点残差，供检查
```
+ 一行汇总：`pocket_CA_RMSD=0.xx Å (n=NN)`。

**调用方式**：`python tools/align_ligand_pose.py <ref_complex.cif> <donor_complex.cif> <out.pdb> --sites out/D0_site_table.csv`

---

## T2 · `check_geometry.py`（D2-S2 的判据计算）

**输入**：一个 PDB（蛋白 + HETATM 配体），来源可以是
① `cif_to_pdb.py` 转出的 Chai-1 复合物（当 WT 基线），或
② LigandMPNN `--pack_side_chains 1` 产出的 `packed/*.pdb`（全原子设计模型）。

**要算的三段约束**（阈值写死，便于复现）：

| # | 判据 | 计算 | 通过阈值 |
|---|---|---|---|
| ① | 三甲铵 → 芳香笼 | 季铵 N 到 6 个笼残基**环原子**（`common.RING_ATOMS`）的最近距离；并统计 ≤5 Å 的接触数 | 最近距离 **≤ 5.0 Å** 且 ≤5 Å 的**不同笼残基 ≥ 2 个** |
| ② | 羧基 → 正电伙伴 | 羧基 C（或两个 O）到 21 位残基 N 原子（ARG 的 NH1/NH2/NE；若换成 LYS 则用 NZ）的最近距离 | **≤ 4.5 Å** |
| ③ | 硫酮 S → CH 供体 | S 到候选 CH 供体碳（227 的 CG2，若换成 SER 用 CB；备选 50 CD1 / 231 CD1,CD2 / 136 / 323 CE）的最近距离 | **≤ 4.2 Å** |
| ③′ | hercynine 版 | 配体无 S 时，自动改用**咪唑 C2**（连 2 个 N 的环碳，`common.lig_sites` 会给出） | 同一阈值 |

**配体识别**：用 `common.lig_sites()` 自动判官能团（季铵 N / 羧基 C / S / 咪唑 C2）。
**残基替换**：若设计把 227 换成 T/S/V，按上表取对应原子；把 21 换成 K/R，取 NZ 或 NH*。

**输出** `out/D2_S2_geometry_check.csv`：
```csv
variant_id,n_quat_to_cage_min,cage_contacts_le5,carboxyl_to_R21,sulfur_to_donor,
criterion1,criterion2,criterion3,n_pass,note
```
`criterion*` 取 `PASS/FAIL/NA`。

---

## T3 · `analyze_designs.py`（D2-S1 + D4 短名单）

**输入**：
- LigandMPNN 的 `seqs/*.fa`（头部自带 `overall_confidence` 与 `ligand_confidence`）
- `out/D2_S2_geometry_check.csv`
- （标准版起）`out/D2_S3_selectivity_bias.csv`

**要点**：
1. 解析 fasta：头部形如
   `>{name}, id=3, T=0.1, seed=111, overall_confidence=0.83, ligand_confidence=0.91, seq_rec=0.42`
   → 取 `overall_confidence`（全部重设计位点的平均似然）与 **`ligand_confidence`**（靠近配体位点的平均似然，⭐ 本档主排序键）。
2. 从序列里**只抽出可设计位点**（用 `D0_site_table.csv` 的 `resnum` 列表）形成"口袋指纹"字符串（便于去重与人工核对）。
3. 与 S2 几何表按 `variant_id` 合并（建议 `variant_id = <T>/<seed>/<seq_id>`）。
4. **WT 分位**：用同一套流程对 **WT 序列**跑一次（`--redesigned_residues` 留空即可，序列即 WT），
   得到 WT 的 `ligand_confidence` → 计算 `WT_percentile = 在全体变体中的百分位`。
   ⭐ 这是"天然口袋离最优还有多远"的定量，**必须算**。
5. 输出 `out/D4_variant_shortlist.csv`：按 `ligand_confidence` 与 `Δscore` 综合排序取 **前 10 条 + WT 行**。

**输出字段**：
```csv
variant_id,pocket_fingerprint,ligand_confidence,overall_confidence,seq_rec,
dscore_vs_hercynine,criterion1,criterion2,criterion3,n_pass,rank,is_WT,expected_phenotype
```
`expected_phenotype` 按 `docs/08_stage3_results.md` §3 的框架填（读硫端 / 羧基锚 / 笼 → 对 EGT 与 hercynine 的预期影响）。

---

## T4 · `run_ligandmpnn.sh`（D1 的命令封装）

把 `run.py` 的调用写成一个脚本，**参数从 `D0_site_table.csv` 自动生成**（`make_site_table.py` 末尾已经打印出建议的
`--fixed_residues` / `--redesigned_residues` 两串，可直接粘贴）。

**`run.py` 的关键参数（已核对官方 README）**：

| 参数 | 本档取值 | 说明 |
|---|---|---|
| `--model_type` | `ligand_mpnn` | 原子上下文条件模型（能读配体原子） |
| `--pdb_path` | 转好的 `mfsT1_EGT.pdb` | 蛋白链 A + 配体链 B |
| `--out_folder` | `out/D1_EGT_T{0.1,0.2,0.3}` | **每个温度一个目录** |
| `--fixed_residues` | `A21 A111 A135 A138 A227 A319 A349` | 核心锚点 + 结构 Gly（见 D0 输出） |
| `--redesigned_residues` | 其余 13 个次级衬里位点 | 只在这里采样 |
| `--ligand_mpnn_use_atom_context` | `1` | ⭐ 让配体原子参与条件（default 即 1，显式写出防误改） |
| `--ligand_mpnn_use_side_chain_context` | `1` | 把固定残基的侧链当额外上下文 |
| `--ligand_mpnn_cutoff_for_score` | `6.0` | 影响 `ligand_confidence` 的统计范围 |
| `--batch_size` / `--number_of_batches` | `2` / `50`（GPU） | 目标 **100 条/温度**；OOM 时降到 `1` |
| `--temperature` | `0.1 / 0.2 / 0.3` | 三档多样性 |
| `--save_stats` | `1` | 保存 `log_probs` 等（`stats/*.pt`），供更细的打分 |
| `--pack_side_chains` | `1` | ⭐ **让模型顺便把这些序列的侧链摆出来**（`packed/*.pdb`）→ 使 D2-S2 能做**真实距离判据**而非规则判据 |
| `--number_of_packs_per_design` | `1` | 先 1 个 pack（足够做几何复核） |
| `--pack_with_ligand_context` | `1` | 摆侧链时考虑配体 |
| `--seed` | `111` | 固定，便于复跑 |

> ⚠️ **先 `python run.py --help` 核对 flag 名称**（版本会演进），再跑正式任务。

**输出目录约定**：`seqs/*.fa`、`backbones/*.pdb`、`packed/*.pdb`、`stats/*.pt`。

---

## T5 · `run_score.sh`（D2-S3 反向筛选的命令封装）

用官方 `score.py` 给"**给定序列 + 给定配体上下文**"打分：

```bash
# 以 EGT 骨架为参考，分别用「EGT 配体」与「hercynine 位姿」两套上下文给同一序列打分
python score.py --model_type ligand_mpnn --seed 111 \
  --single_aa_score 1 --use_sequence 1 \
  --pdb_path <backbone_with_EGT.pdb>       --out_folder out/D2_S3/EGT \
  --batch_size 1 --number_of_batches 10
python score.py --model_type ligand_mpnn --seed 111 \
  --single_aa_score 1 --use_sequence 1 \
  --pdb_path <backbone_with_hercynine.pdb> --out_folder out/D2_S3/hercynine \
  --batch_size 1 --number_of_batches 10
```
⚠️ 两个 PDB 必须**共用同一套骨架构象**（这就是 T1 存在的理由）。
⚠️ 输入 PDB 的**序列必须是该变体的序列**——`packed/*.pdb` 里已经是设计序列，可直接用（只需替换配体）。
⚠️ `--number_of_batches ≥ 10` 是官方建议（自回归打分依赖解码顺序）。

**Δscore 定义**（写进结果说明，防误用）：
```
Δscore = score(该序列 | 骨架 + EGT 上下文) − score(该序列 | 骨架 + hercynine 上下文)
      > 0 表示该序列在"有硫"的上下文里更被偏好
```
⛔ **它是设计偏置量，不是亲和力差；不得换算成 ΔΔG，不得声称选择性已被预测。**
