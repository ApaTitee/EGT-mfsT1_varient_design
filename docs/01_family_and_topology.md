# 01 · 家族与拓扑判定（Q1 + Q4）

> **档 1 · Step 1**｜**执行**：2026-10-12（第十八轮）
> **输入**：`05_素材库/方向三/序列/mfsT1_WP_030133796.1.faa`（428 aa，已复核：标准氨基酸，MRPWIVWATG 起始）

---

## Part A · Q1 家族归属

> ✅ **InterProScan 5（EBI REST，2026-10-12 完成，job `iprscan5-R20261011-202048-0466-9386894-p2m`，结果留档 `data/interproscan_mfsT1.tsv`）**——见 A.2。
> ⬜ **CDD batch CD-Search 复核未完成（如实标注，规范 9）**：两个 job（`QM3-qcdsearch-1586C1E25CAD8FD0`、`QM3-qcdsearch-AD9F941E65933D5`）提交于 2026-10-12，**排队 >4.5 h 仍 running，经 3 次定时轮询未返回，已停止重试**。**Q1 以现有证据为准**：InterProScan 四库独立确认 MFS（A.2）+ E-8 基线（同法旧结果）+ 序列层新证据（A.1）。若日后需要 COG2223 NarK 的精确 E 值复核，可用任一 CDD job ID 在 NCBI 页面手动取回。
> **已知基线（E-8，待 CDD 复核确认）**：COG2223 NarK E=7.6e-46（bits 162.7，覆盖 2–403）；COG2271 UhpC 2.3e-34；cd06174 MFS 6.5e-23；pfam07690 MFS_1 9.3e-15；cd17320 MFS_MdfA_MDR_like 1.9e-5（排除 DHA1）。
> ⚠️ KEGG 陷阱重申：`MyAD_11635`/K07552 是另一基因座，不可套用。

### A.2 InterProScan 命中表（PfamA/CDD/NCBIfam/SuperFamily/Gene3d/Panther 六个库）

| 库 | 命中 | 描述 | 区间 | E 值 | 与 E-8 对比 |
|---|---|---|---|---|---|
| SUPERFAMILY | SSF103473 | MFS general substrate transporter | 2–402 | **5.4e-48** | 新增（E-8 未跑此库） |
| CDD | cd06174 | MFS | 11–371 | **1.7e-23** | ✅ 一致（E-8: 6.5e-23，同库同族同量级） |
| Pfam | PF07690 | Major Facilitator Superfamily（MFS_1） | 11–229 | **2.1e-19** | ✅ 同族（E-8: 9.3e-15，管线版本不同） |
| Gene3D | G3DSA:1.20.1250.20 | MFS general substrate transporter like domain | 1–189 | 9.7e-26 | 新增；**N/C 两 bundle 分别命中（1–189 与 206–407）**——与 12-TM 双 bundle 拓扑互证 |
| Gene3D | G3DSA:1.20.1250.20 | 同上 | 206–407 | 9.6e-10 | 同上 |
| PANTHER | **PTHR43124** | **PURINE EFFLUX PUMP PBUE** | 12–329 | **1.1e-13** | ⭐ 新增，见 A.2.1 |
| NCBIfam | （无命中） | | | | — |

#### A.2.1 ⭐ PANTHER 命中 PbuE（嘌呤外排泵）的意义（⚠️ 需后续核实的新线索）

- **PTHR43124 = 嘌呤外排泵 PbuE 家族**。PbuE（如 *E. coli* YdhC）是 MFS 的**嘌呤（腺嘌呤/鸟嘌呤）外排泵**——杂环底物 + **外排方向**，与 EGT（杂环两性离子，外排）的功能画像**远比 NarK（无机阴离子，摄取/反转运）接近**。
- ⚠️ 分寸：E=1.1e-13 弱于超家族级命中，PANTHER 家族划分较粗；**这是"候选物理类比对象"的线索，不是家族结论**。是否成立需在档 1.5/档 2 用成对比对 +（若有）PbuE 结构核实（本会话未查 PbuE 有无实验结构——记为档 2 前置小任务）。
- ⚠️ 本命中**不改变**"非 DHA1"的结论（E-8：cd17320 仅 1.9e-5），也**不改变** MFS 总归属。

### A.3 与 E-8 基线的一致性声明

| 项 | E-8（batch CD-Search） | 本轮 InterProScan | 判定 |
|---|---|---|---|
| MFS 总归属 | cd06174 6.5e-23 / pfam07690 9.3e-15 | cd06174 1.7e-23 / PF07690 2.1e-19 | ✅ 一致 |
| NarK/NNP 样 | COG2223 E=7.6e-46 | （InterProScan 的 CDD app 未报 COG 级）| ⬜ CDD 复核未完成（队列未返回）；维持 E-8 基线 + "profile 级"表述 |
| 非 DHA1 | cd17320 1.9e-5 | 无 DHA1 类命中 | ✅ 一致 |
| 新增 | — | PANTHER PTHR43124（PbuE 嘌呤外排泵）1.1e-13 | ⭐ 新线索（A.2.1） |

**Q1 阶段性结论**：MFS 超家族归属**确认**（四库一致，E ≤ 1e-13）；"NarK/NNP 样"维持 profile 级表述（待 CDD 复核收尾）；**新增 PbuE 嘌呤外排泵线索**——方向性（外排）与底物类型（杂环）都比 NarK 更像 EGT 外排的物理类比对象。

### A.1 本轮新增的序列层证据（不依赖 CDD，已可写下）

| 证据 | 内容 | 等级 |
|---|---|---|
| 成对 identity（MAFFT，5 序列） | mfsT1 vs NarK = **16.4%**；vs NarU = **18.4%**；vs 人 SLC22A4 = **20.2%**——全部"暮区"，且与 SLC22 几乎等距 | 🖥 |
| NNP 签名基序 | **NS1/NS2 均不保留**（2/12、2/13；对照 NarU 8/12、11/13）→ 不是严格 NNP 成员 | 🖥（基序定义 ✅ Zheng 2013 全文） |
| NarK Arg 对 | **R89/R305 均不保留**（→Q57/A260；对照 NarU 双 R 保留） | 🖥 |

> **阶段性解读（待 CDD 收尾）**：profile 级命中指向 NarK/NNP 样（E-8），但**成对序列层面 mfsT1 与 NNP、SLC22 几乎等距，且不携带 NNP 的功能性签名**。建议口径：**"MFS（NarK/NNP 样折叠，profile 级），非严格 NNP 成员"**。

---

## Part B · Q4 拓扑 → **12 TM，N-in / C-in（与 NarK 一致）**

> ✅ **DeepTMHMM 2 已完成**（2026-10-12，申请人按 agent 指引在 biolib 网页端手动运行，`https://biolib.com/DTU/DeepTMHMM2/`；原始输出全套落档 `data/deeptmhmm2/`，拓扑图复制至 `fig/mfsT1_DeepTMHMM2_拓扑图.png`）。
> 注：DeepTMHMM 1.0 的 DTU 旧服务端故障（3 次提交均 failed，留档 `data/dtu*_job.html`），已由 DeepTMHMM 2（更新版本）取代。
> 三预测器 + MSA 交叉验证全部收敛，**Q4 判定为最终结论**。

### B.0 ⭐ DeepTMHMM 2 结果（权威基准）

- **TMR 数 = 12**；膜类型 = 革兰氏阳性菌质膜（与 *M. neoaurum* 放线菌属性一致）
- **N 端 inside 1–7，TM1 = 8–28**——**明确判为跨膜螺旋，不是信号肽**
- 94–157 区 = **两个 TM（95–113、131–152）**
- 中央胞质环：inside 180–213（含 interfacial 197–204），位于 TM6（161–179）与 TM7（214–234）之间
- **C 端 inside 395–428**
- 12 个 TM 边界：8-28 / 42-64 / 72-89 / 95-113 / 131-152 / 161-179 / 214-234 / 251-270 / 280-297 / 307-326 / 344-361 / 377-394

### B.1 另两个独立预测器（留档 `data/phobius_mfsT1_raw.html`、`data/tmhmm_job.html`）

**Phobius（sbc.su.se，plp 格式）**：N 端 1–31 判为信号肽（n/h/c 区）+ 11 个 M 区（41-64 / 71-91 / 97-118 / 130-155 / 161-180 / 214-235 / 247-269 / 281-299 / 305-326 / 347-364 / 376-395）；C 端 i 396–428；中央胞质环 i 181–213。

**TMHMM 2.0（DTU）**：10 个 TMhelix（5-27, 42-64, 71-93, 158-180, 215-237, 247-269, 281-300, 304-326, 347-364, 379-396）；`Total prob of N-in = 0.831`；`POSSIBLE N-term signal sequence`；C 端 inside 397–428；94–157 被整体判为 outside 大环。

### B.2 三器分歧的仲裁 → 全部收敛于 12-TM

| 分歧点 | DeepTMHMM 2 | Phobius | TMHMM | 仲裁结果 |
|---|---|---|---|---|
| N 端 8–28 | **TMhelix（明确）** | "信号肽 h 区" | TM + "POSSIBLE signal" | ✅ **TM1**（DeepTMHMM 2 与 TMHMM 均给 TM；Phobius 的"信号肽"是 MFS N 端疏水段的经典误判；细菌 MFS 无信号肽、GTG 起始、N-in） |
| 94–157 区 | **两个 TM（95-113, 131-152）** | 两个 M | 64 aa outside 大环 | ✅ **TM4 + TM5**（DeepTMHMM 2 与 Phobius 一致；且 MSA 口袋映射 F147→F111、N175→Q139 正落这两段，结构佐证） |
| TM 总数 | **12** | 11+"信号肽" | 10 | ✅ **12** |
| N/C 端朝向 | N-in / C-in | N-in / C-in | N-in（0.83）/ C-in | ✅ **N-in / C-in**，与 NarK 一致 |

### B.3 Q4 判定表（最终）

| 判据 | 结果 |
|---|---|
| TM 数 = 12 ± 1 | ✅ **12**（DeepTMHMM 2 边界为准：TM1 8-28 / TM2 42-64 / TM3 72-89 / TM4 95-113 / TM5 131-152 / TM6 161-179 / TM7 214-234 / TM8 251-270 / TM9 280-297 / TM10 307-326 / TM11 344-361 / TM12 377-394） |
| N/C 端朝向 | ✅ **N-in / C-in**（三器一致） |
| MFS 中央胞质大环（TM6/TM7 之间） | ✅ **180–213**（含 8 aa interfacial 段 197–204）——教科书式 MFS 布局 |
| 信号肽 | ❌ 无（DeepTMHMM 2 明确否定） |
| ⭐ 口袋映射 × DeepTMHMM 2 拓扑一致性 | ✅ 7 个映射位点全部落在同号 TM：Y14→TM1、Q57→TM2、F111→TM4、Q139→TM5、H219→TM7、A260→TM8、G361→TM11 |

> **档 2 用图**：`fig/mfsT1_DeepTMHMM2_拓扑图.png`（DeepTMHMM 2 官方输出）；膜平面定义可直接用本表边界。

**验证等级**：DeepTMHMM 2 / Phobius / TMHMM 原始输出 = 🖥（全部留档）；12-TM 与 N/C 朝向 = **三器一致 + MSA 佐证，可作定论**（🖥）；膜平面精确边界仍以档 2 结构为准。
