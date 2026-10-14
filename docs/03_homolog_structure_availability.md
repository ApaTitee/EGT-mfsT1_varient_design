# 03 · 同源结构可得性（Q5）

> **档 1 · Step 4**｜**执行**：2026-10-12（第十八轮）
> **问题**：有没有已存在的同源实验结构覆盖口袋？（若有 → HOLD 分支：跳过档 2 AF2）

---

## 1. 三条检索路径的结果

| 路径 | 状态 | 结果 |
|---|---|---|
| **RCSB sequence search**（MMseqs2，identity_cutoff 0.1、evalue 10，experimental only） | ✅ 完成 | **0 命中**（对照实验：同参数下 NarK P10903 命中 7 个结构含 4JR9/4JRE/4U4T → API 正常，0 命中为真阴性；查询 JSON 与对照留档 `data/rcsb_q.json`、`data/rcsb_nark_test.json`） |
| **NCBI BLAST blastp vs pdb**（RID `CRNN9ZV1014`） | ⬜ **未完成**（队列 >4.5 h 未返回，3 次轮询后停止重试，规范 9） | 未获结果；**不影响判定**——RCSB 0 命中（含阳性对照）+ 实测最高 identity 16.4%（NarK）已封顶"高 identity 模板"的可能性；如需收尾可用 RID 在 NCBI 页面手动取回 |
| **FoldSeek**（search.foldseek.com，序列模式 3diaa = ProstT5+Foldseek） | ⚠️ **未完成** | 服务端对同一序列两次返回 `ERROR`（ticket 留档 `data/foldseek_ticket.json`）；按规范 9 如实标注，不以推理填补 |

## 2. 判定

> **不存在可直接用作骨架的高 identity 实验结构。**
> - RCSB（MMseqs2）在 10% identity 截断下 0 命中——这是最灵敏的序列-结构库初筛之一。
> - 我方 MSA 实测：mfsT1 vs 最近实验结构 NarK = **16.4% identity**（384 同列非 gap）——处于"暮区"，**即便 BLAST 返回 NarK 系命中，16% identity 也不足以做可靠的同源建模骨架**（口袋几何不能直接移植）。
> - 唯一可能的"改判"情形：BLAST 返回某个 **identity ≥30% 且覆盖 TM 口袋段**的意外命中 → 届时升级为 HOLD 分支。**概率极低，但若出现必须执行 HOLD**。

## 3. 模板推荐（给档 2）

| 用途 | 推荐 | 说明 |
|---|---|---|
| 建模方法 | **AF3 / ColabFold（共进化信号）**，不做经典同源建模 | identity 太低，共进化是主要信号源；NarK 只作**后验比对模板** |
| 后验比对模板 | **4U4T**（nitrate-bound inward-open，2.40 Å）+ **4U4V**（apo inward-open） | 本库已有；口袋残基清单见 `02b`；映射见 `02` |
| （建议补下） | **4U4W**（nitrate-bound occluded，2.40 Å） | 本库尚无，档 2 时从 RCSB 下载——inward/occluded 两态比对更完整 |

**验证等级**：RCSB 0 命中 = 🖥（含阳性对照）；identity 数字 = 🖥（MAFFT 实算）；BLAST/FoldSeek 状态 = 如实标注。
