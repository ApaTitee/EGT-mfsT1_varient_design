# 05 · 空白与占位复核（Q7 + Q8）

> **档 1 · Step 6**｜**执行**：2026-10-12（第十八轮）｜**检索平台**：Europe PMC REST API（`europepmc.org/webservices/rest/search`）
> **留档**：`data/evidence/europepmc_search_records.json`（检索式原文 + 命中数 + 前 15–20 条标题）

---

## 1. 检索式与命中数（原样留档）

| # | 检索式（原文） | 命中数 |
|---|---|---|
| Q1 | `(ergothioneine OR EGT OR "L-ergothioneine") AND (export* OR efflux OR secret* OR excret*) AND (transporter OR permease OR pump OR channel)` | **931** |
| Q2 | `"mfsT1"` | **6** |
| Q3 | `(ergothioneine) AND (MFS OR "major facilitator") AND (efflux OR export)` | **39** |
| Q4 | `(ergothioneine) AND (Mycobacterium OR Mycolicibacterium OR Corynebacterium)` | **381** |
| Q1b（收窄补检） | `TITLE_ABS:"ergothioneine" AND TITLE_ABS:(efflux OR export* OR secret* OR excret*) AND TITLE_ABS:(transporter OR permease OR pump OR channel)` | **8** |
| Q1c（收窄补检） | `TITLE_ABS:"ergothioneine" AND TITLE_ABS:"efflux"` | **7** |
| Q4b（收窄补检） | `TITLE_ABS:"ergothioneine" AND TITLE_ABS:(efflux OR export OR secretion) AND TITLE_ABS:(Corynebacterium OR Mycolicibacterium OR Mycobacterium)` | **0** |

> Q1/Q3/Q4 宽式命中量大且绝大多数为"全文任意位置命中"的噪声（综述、医理研究、其他代谢物），故补做 TITLE_ABS 级收窄（Q1b/Q1c/Q4b）用于判定。

## 2. Q7 · "外排侧空白"是否仍然成立？→ **成立**

逐条核查（题名/摘要级 ◻，重点条目摘要级 ✅）：

| 核查点 | 结果 |
|---|---|
| Q1b 的 8 条 | 全部为人/动物 **OCTN1 摄取**、*C. elegans* 摄取（2013）、综述——**无**微生物 EGT 外排蛋白鉴定 |
| Q1c 的 7 条 | 无 EGT 外排蛋白（命中词为其他代谢物的 efflux，如 gliotoxin/heme） |
| Q4b | **0 命中**——分枝杆菌/棒杆菌属无 EGT 外排/分泌蛋白的题录级报道 |
| 2025 *C. glutamicum* 渗透压增产（PMID 41359152，摘要 ✅） | 渗透压**生理效应**增产 459 mg/L，**未鉴定任何外排蛋白**；摘要自述"目前研究的都是 EGT **import** 系统" |
| 2026 酵母 transport×antioxidant 综述（PMID 42792589，摘要 ✅） | 碳源摄取×抗氧化合成的泛综述，**不涉及 EGT 外排蛋白** |
| 2026 EGT 机制综述（PMID 42707330，摘要 ✅） | 全部围绕 **OCTN1（摄取）**，无外排侧新蛋白 |

> ✅ **结论（Q7）**：截至 2026-10-12，**没有任何团队鉴定出 EGT 外排/分泌蛋白**（*M. smegmatis* 2013 仍只有"分泌现象"；2024 之后的新工作集中在合成通路与宿主工程，未触外排蛋白身份）。**立论不变。**

## 3. Q8 · 占位检查（逐格点名，规范 11）

**我核实的是：mfsT1 本身是否已有团队在做工程化/结构/设计——不是"别的 MFS 工程化"、也不是"别的宿主 EGT 外排"。**

Q2（`"mfsT1"`）6 条命中逐条判定：

| 条目 | 判定 |
|---|---|
| Ding 2024（PMID 38253966） | **mfsT1 唯一来源报道**（已知，即我方基准文献） |
| 综述 2024（PMID 39717282, "Ergothioneine biosynthesis: present state and future prospect"） | 综述引用 Ding 2024，非独立研究 |
| 2011/2015/2019 的 3 条（米糠/小麦胚芽提取物） | ⛔ **时间性排除**：均早于 mfsT1 命名（2024），为全文检索误命中（非题名/摘要） |
| 2015 小麦胚芽（PMID 25866747） | 同上，误命中 |

> ✅ **结论（Q8）**：**全球范围内 mfsT1 全长报道仍只有 Ding 2024 一篇**；无第二团队在做 mfsT1 的工程化、结构或设计。**占位检查通过。**

## 4. 边界与诚实声明

- 检索限于 Europe PMC（覆盖 PubMed + PMC 全文 + 预印本）；未单独检 Google Scholar——按 PLAN 约定 Europe PMC 为主平台，结论在该范围内成立。
- Q1b/Q1c/Q4b 为收窄补检（非 PLAN 原 4 式），判定级为 ◻（题录/摘要级）；三篇重点文献已读摘要 ✅。
- 无检索服务不可用项（本轮 Europe PMC 全部成功）。

**验证等级**：检索式与命中数 = 🖥 原始返回留档；逐条判定 = ◻ 题录/摘要级（三篇重点 = 摘要级 ✅）。
