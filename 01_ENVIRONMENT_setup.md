# 01 · 环境与工具链（按 2026-10-14 WSL 实测写定）

> **实测来源**：申请人执行的 `WSL环境核查输出.md`（Ubuntu 24.04.5 / WSL2）。以下命令**已按实测结果写死**，不留"二选一"。
> ⚠️ **不需要 sudo**：下面全部是用户态操作（已确认 `sudo` 需要密码，所以刻意全部避开）。

---

## 0. 实测到的环境（结论摘要）

| 项 | 实测值 | 对档 4 的意义 |
|---|---|---|
| 发行版 / 内核 | **Ubuntu 24.04.5 LTS**，WSL2，`6.18.54.1-microsoft-standard-WSL2` | ✅ 一切常规 Linux 工具链可用 |
| CPU / 内存 | **16 核** / **15 GiB**（可用 13 GiB） | ✅ 即使全程 CPU 也够 |
| 磁盘 | `/` = **1007 GB，可用 904 GB** | ✅ 远超所需的 ~20 GB |
| **GPU** | **RTX 5060 Laptop, 8151 MiB, compute_cap 12.0**（sm_120），驱动 617.42 | ✅ GPU 可见 → **优先走 cu128 torch**；⚠️ 8 GB 显存适中，OOM 时降 `batch_size` |
| Python | `python3` = **3.13.12**（`/home/apatite/miniconda3`），conda **26.1.1** | ⚠️ base 的 3.13 偏新 → **另建 3.11 环境** |
| 已有包 | torch / pymol / rdkit **均未安装**；numpy/scipy/matplotlib **不全** | 需在本环境内装齐（见 §2） |
| 网络 | pypi / conda-forge / github **全部可达**；走本机代理 **`http://127.0.0.1:15732`**（`http_proxy` 等已设） | ✅ pip/conda 会自动继承；若某次下载卡住见 §4 |
| 工具 | `unzip 7z zip wget curl rsync git gcc g++ make` **都在**；`cmake` / `aria2c` **缺**（不需要） | ✅ 压缩包可用 `.zip` |
| 图形库 | **`libGL.so.1` 存在** | ✅ PyMOL 无头渲染不需要 apt |
| `/mnt/c/...` 中文路径 | **可以访问** | ⚠️ 但仍建议解压到 **WSL 原生盘**（见 §1） |

---

## 1. 工作目录：放 WSL 原生盘，不要放 `/mnt/c`

```bash
mkdir -p ~/work && cd ~/work
# 假设压缩包已拷到 ~/  （例如： cp /mnt/c/.../EGT-D4_*.zip ~/ ）
unzip -O UTF-8 ~/EGT-D4_*.zip -d ~/work/
cd ~/work/EGT-D4
```

**为什么不用 `/mnt/c`**：① 跨文件系统 I/O 明显更慢（模型权重读写、多轮采样会放大）；② 中文+空格路径在部分脚本里会出问题；③ `/mnt/c` 上的 Linux 权限位是伪造的。

---

## 2. 建环境（一次性，约 10–20 分钟，主要时间花在 torch 与权重下载）

```bash
# ① 建一个干净环境（3.11 —— 与 LigandMPNN 官方 README 一致，wheel 生态最稳）
conda create -y -n egt4 python=3.11
conda activate egt4

# ② PyTorch：先试 cu128（sm_120 需要 torch ≥ 2.7 的 cu128 轮子）
pip install torch --index-url https://download.pytorch.org/whl/cu128
python -c "import torch;print('torch',torch.__version__,'cuda_ok',torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else '')"
#   ↳ 期望：cuda_ok True 且打印出 RTX 5060
#   ↳ 若为 False：pip install torch   （改 CPU 版，本任务量级完全够，见 §3 末）

# ③ 本包自带工具需要的包
pip install numpy scipy matplotlib pandas biopython rdkit gemmi
```

### 2.1 LigandMPNN 本体 + 权重（约 1 GB）

```bash
cd ~/work
git clone --depth 1 https://github.com/dauparas/LigandMPNN.git
cd LigandMPNN
pip install -r requirements.txt          # 含 prody 等；PDB 读写靠它
bash get_model_params.sh "./model_params"
ls -la model_params/ | head             # 应能看到 ligandmpnn_v_32_*.pt 等
python run.py --help | head -40         # ✅ 冒烟：能打印参数即环境通了
```

> ⚠️ **先跑 `--help` 再跑真任务**：确认 flag 名称与你环境里的版本一致（本包文档按官方 `main` 分支的 README 写，但版本会演进）。

### 2.2 PyMOL（可选，只在想画"漂亮结构图"时装）

```bash
conda install -y -c conda-forge pymol-open-source
python -c "import pymol;print('pymol ok')"
```
> ⚠️ Windows 侧的 PyMOL pip wheel 曾因缺 DLL 而不可用；**WSL/conda-forge 是它的正常安装路径**，且 `libGL.so.1` 已存在。
> 不装也不影响主线：本包 `tools/render_structure.py` 是纯 matplotlib 的结构渲染器，已产出过正式图件。

---

## 3. 环境自检（跑完 §2 立刻做，10 秒）

```bash
cd ~/work/EGT-D4
bash smoke_test/run_smoke_test.sh
```
它会：① 把 Chai-1 的 EGT 复合物 cif 转 PDB；② 修配体键级（应报出 **C=S 双键**）；③ 生成 `D0_site_table.csv`。
把输出与 `smoke_test/expected/` 下的期望产物对比 —— **三者都对上，环境即就绪**。

---

## 4. 已知坑与对策（都与本次实测有关）

| 坑 | 表现 | 对策 |
|---|---|---|
| **代理** | 下载卡住 / 超时；报 TLS 或连接被重置 | 你的 `http_proxy/https_proxy` 已指向 `127.0.0.1:15732`，pip 会自动用。若仍失败：`pip install --proxy http://127.0.0.1:15732 <pkg>`；conda 则 `conda config --set proxy_servers.http http://127.0.0.1:15732` |
| **sm_120 不被支持** | `torch.cuda.is_available()` 为 False，或运行时报 `no kernel image is available` | 说明 torch 版本过旧 → 重装 cu128（torch ≥ 2.7）；仍不行就用 CPU 版 |
| **显存不够（8 GB）** | `CUDA out of memory` | `--batch_size 1`；仍 OOM 则 `--pack_side_chains 0` 分两轮跑 |
| **conda 环境未激活就跑脚本** | `ModuleNotFoundError: torch/rdkit` | 每个新终端先 `conda activate egt4` |
| **中文路径** | Prody / 脚本读文件报编码错 | 工作目录尽量用纯 ASCII（若解压目录名含中文且遇问题 → 包名已是 `EGT-D4`，无需改） |
| **CPU 路线的可行性** | —— | 428 aa 单链、200–500 条序列，**16 核 CPU 跑数十分钟量级**；GPU 不是必需品，只是加速 |

---

## 5. 与 Windows 侧的差异（避免踩同样的坑）

- **在 Windows 侧踩到过的**：conda 装 PyMOL 时中断事务 → 把 base 环境的 `certifi` 删了（已修）。**WSL 侧不要把两次安装并行跑**，装完一个再装下一个。
- **控制台编码**：Linux 下默认 UTF-8，**不需要** Windows 那套 `PYTHONUTF8=1`（设了也无害）。
- **本包工具全部是纯 numpy + 标准库**（`fix_ligand_bond_orders.py` 额外用 RDKit），**不依赖 PyMOL**、不依赖 GPU，随时可跑。
