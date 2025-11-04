# LGA 跨机器 Segmentation Fault 检修报告（含加固与堆分配改造方案）

## 摘要
- 问题本质：不同机器默认栈大小差异导致的运行时“栈溢出（stack overflow）”，继而触发 Segmentation fault；与架构、GLIBC 或编译器无关。
- 根因机制：LGA 源码在多个函数中分配了“超大”的局部数组/结构体（位于函数栈上），典型如 `MAXRES`（9501）和 `MAXATOMS`（95001）规模的二维或一维数组；在默认 8MB 栈限制的机器上必然溢出。
- 已实施修复：将最关键的超大局部对象迁移到静态存储（BSS），绕开栈限制；该方法最小侵入、行为等价。
- 建议加固：继续“迁移中等体量局部数组至静态/堆”，提供了两种可选实现（A. 静态迁移；B. 堆分配重构），并给出具体补丁思路与接口设计。

---

## 1. 现象与复现
- 现象：`setup.sh` 按现有方式编译得到的 `lga` 二进制，仅在一台机器上可以运行；其他“看似相同架构和环境”的机器运行即 Segmentation fault。
- 复现实验：
  - 在“问题机器”重新编译得到的新 `lga`：在“好机器”上可运行，在“问题机器”上仍崩溃。
  - 这说明二进制本身“没问题”，问题出在运行环境差异——典型就是“默认栈限制”。
- 快速验证：在问题机器运行 `ulimit -s`，常见输出 8192（8MB）；在好机器往往是 `unlimited`。若先执行 `ulimit -s unlimited` 再运行 `./lga ...`，通常可暂时规避崩溃。

---

## 2. 根因与触发条件
- 代码中多个“超大局部对象”导致栈溢出：
  - `main` 里大结构体（已修复）：
    - `LGA_package_src/SRC/lga.c:335` `static pdata_struct pdata;`
    - `LGA_package_src/SRC/lga.c:336` `static part_struct  part;`
  - `part_struct` 内部自身极大（`LGA_package_src/SRC/tools.h`）：
    - `printlga[MAXRES+1]` 且 `stext.line[MAXRES+1]`，约 9502×9502 ≈ 90MB。
    - `opt_sup/opt_sup1/opt_sup2/opt_set_nb`，总计 ~8MB。
    - `pairs[2][MAXGAP2+1]`（MAXGAP2=65536），~1MB。
  - 关键算法路径中的局部大数组（示例）：
    - `tools_2.c: rmsd_isp` 中的 `calphas[MAXITR+1]`（含 2×(MAXRES+1) 的 long 组合），约 18.5MB。
    - `tools_2.c: rmsd_dist` 的 `opt1/opt2`（21×MAXRES），约 3.2MB。
    - `tools_2.c: lcs_run/gdt_run/lcs_gdt_analysis` 中的多组 `calphas/opt_lcs125/gdt_chars/l_atoms`，合计 ~2MB。
    - `tools_3.c` 与 `tools_4.c` 中还有多处 `MAXRES` 级局部数组（见 §4.1/§4.2）。
- 触发条件：一旦函数栈帧尝试容纳上述数组，在默认 8MB 栈上即会溢出（进程创建或函数进入时即触发）。

---

## 3. 已实施的最小修复（已提交）
- 目标：以“最小侵入”方式消除最易触发的栈溢出点。
- 具体修改：
  1) `LGA_package_src/SRC/lga.c`
     - 将 `pdata_struct pdata;` 与 `part_struct part;` 改为 `static`，将其存储迁移到 BSS 段。
  2) `LGA_package_src/SRC/tools_2.c`
     - `rmsd_isp`：`iter_struct calphas[MAXITR+1]` 改为 `static`。
     - `rmsd_dist`：`long opt1[21][MAXRES], opt2[21][MAXRES]` 改为 `static`。
     - `lcs_run` / `gdt_run`：`long calphas[2][MAXRES+1]` 改为 `static`。
     - `lcs_gdt_analysis`：`opt_lcs125/gdt_chars/l_atoms` 改为 `static`。
- 取舍：
  - `static` 简单直接、修改少，立即消除栈开销；代价是增加 BSS/数据段驻留内存（单进程 CLI 可接受），函数非线程安全（对 LGA 单线程程序非问题）。

---

## 4. 建议的进一步加固

### 4.1 方案 A：继续将“中等体量局部数组”改为 static（快速稳妥）
- 目标：覆盖剩余仍在栈上的中等数组，避免在极端“超大输入 + 小栈”的组合下触发边缘问题。
- 推荐改动清单（示例）：
  - `LGA_package_src/SRC/tools_4.c`
    - `rms_window(...)`：
      - `long  calphas[2][MAXRES+1], rwa1[MAXRES+1], rwa2[MAXRES+1];` → 前者改 `static`，后两者可改 `static`（总计约几十 KB）。
    - `print_lga(...)`：
      - `char span1[MAXRES+1], span2[MAXRES+1];`、`char span1_fix[MAXRES+1], span2_fix[MAXRES+1];`、`char sub[MAXRES+1];` → 改 `static`。
    - `check_all_atoms(...)`：
      - `eval_list eset[MAXRES+1], esup[MAXRES+1];`、`long list0[MAXRES+1], list1[MAXRES+1];` → 改 `static`。
    - `best_fit(...)`：
      - `long sup[2][MAXRES+1], sup1[2][MAXRES+1], sup2[2][MAXRES+1];` → 改 `static`。
  - `LGA_package_src/SRC/tools_3.c`
    - `align_search(...)`：
      - `long tmp[MAXRES+1];`、`float xa[MAXRES+1][3], xb[MAXRES+1][3];` → 改 `static`。
    - `alignment_run(...)`：
      - `long calphas[2][MAXRES+1];` → 改 `static`。
    - `alignment_best_match(...)`：
      - `long l_atoms[2][MAXRES+1], ind[3][MAXRES+1], out[2][MAXRES+1];` → 改 `static`。
- 示例补丁（片段，仅示意）：
  ```diff
  --- a/LGA_package_src/SRC/tools_4.c
  +++ b/LGA_package_src/SRC/tools_4.c
   void rms_window(...)
   {
  -  long  calphas[2][MAXRES+1], rwa1[MAXRES+1], rwa2[MAXRES+1];
  +  static long  calphas[2][MAXRES+1];
  +  static long  rwa1[MAXRES+1];
  +  static long  rwa2[MAXRES+1];
   }
  
   long align_search(...)
   {
  -  long  tmp[MAXRES+1];
  -  float xa[MAXRES+1][3], xb[MAXRES+1][3];
  +  static long  tmp[MAXRES+1];
  +  static float xa[MAXRES+1][3];
  +  static float xb[MAXRES+1][3];
   }
  ```
- 优点：
  - 实现量小、易代码审阅；立即消除栈风险。
- 局限：
  - 增加 BSS/数据段驻留；函数非线程安全；不同调用间共享缓冲（但对当前 CLI 单线程行为影响可忽略）。

### 4.2 方案 B：堆分配重构（推荐长期方案）
- 目标：将所有“超大/中等”缓冲迁移到堆，集中管理，按需分配与释放；减少全局驻留内存，提高可维护性与可移植性。
- 设计要点：
  1) 新建工作区结构 `lga_workspace`，收纳所有大缓冲；并在 `part_struct` 中增加一个指针 `lga_workspace *ws;`（或独立为模块全局单例）。
  2) 提供初始化/释放接口 `lga_ws_init(...)` / `lga_ws_free(...)`，在 `main` 中成对调用。
  3) 以“一维数组 + 索引宏”替代原二维数组，避免指针数组碎片，提升局部性。
  4) 渐进式迁移：先迁 `tools_4.c` 的窗口缓冲与列表，再迁 `tools_3.c` 的 `xa/xb/tmp`、`ind/out/l_atoms`，最后考虑把 `tools_2.c` 的静态缓冲也迁到堆。

- 接口与数据结构（示例）：
  ```c
  // LGA_package_src/SRC/tools.h
  typedef struct {
    // 标准索引宏：行主序
    // 2×(MAXRES+1) 的 long 矩阵
    long *calphas2;   // size = 2*(MAXRES+1)

    // 窗口缓冲
    long *rwa1;       // size = MAXRES+1
    long *rwa2;       // size = MAXRES+1

    // 对齐/匹配相关
    long *l_atoms2;   // size = 2*(MAXRES+1)
    long *ind3;       // size = 3*(MAXRES+1)
    long *out2;       // size = 2*(MAXRES+1)

    // GDT/LCS 相关
    long *opt1;       // size = 21*MAXRES
    long *opt2;       // size = 21*MAXRES
    long *opt_lcs125; // size = 3*(MAXRES+1)
    long *gdt_chars;  // size = 21*(MAXRES+1)

    // 几何缓冲：xa/xb 为 (MAXRES+1)×3 的 float
    float *xa3;       // size = (MAXRES+1)*3
    float *xb3;       // size = (MAXRES+1)*3

    // 迭代缓冲（rmsd_isp）：(MAXITR+1) × [2×(MAXRES+1)]
    long *iter_calphas; // size = (MAXITR+1)*2*(MAXRES+1)
  } lga_workspace;

  int  lga_ws_init(lga_workspace **ws);
  void lga_ws_free(lga_workspace **ws);

  // 在 part_struct 中新增：
  // lga_workspace *ws;
  ```

- 索引宏（示例）：
  ```c
  #define IDX2(LD, i, j)    ((i)*(LD) + (j))             // 行主序
  #define L2(a, i, j)       ((a)[ IDX2((MAXRES+1), (i), (j)) ])
  #define L2_21(a, i, j)    ((a)[ IDX2((MAXRES),   (i), (j)) ])  // 21×MAXRES 等场景
  #define F2(a, i, j)       ((a)[ IDX2(3, (i), (j)) ])           // (MAXRES+1)×3
  ```

- 初始化/释放实现（示意）：
  ```c
  // LGA_package_src/SRC/tools_*.c 或新增 workspace.c
  int lga_ws_init(lga_workspace **pws) {
    lga_workspace *ws = calloc(1, sizeof(*ws));
    if (!ws) return -1;
    size_t n = MAXRES + 1;
    size_t n2 = 2*n, n3 = 3*n;
    size_t m21 = 21*MAXRES, m3 = 3*n, mit = (MAXITR+1)*n2;

    ws->calphas2     = calloc(n2, sizeof(long));
    ws->rwa1         = calloc(n,  sizeof(long));
    ws->rwa2         = calloc(n,  sizeof(long));
    ws->l_atoms2     = calloc(n2, sizeof(long));
    ws->ind3         = calloc(3*n,sizeof(long));
    ws->out2         = calloc(n2, sizeof(long));
    ws->opt1         = calloc(m21,sizeof(long));
    ws->opt2         = calloc(m21,sizeof(long));
    ws->opt_lcs125   = calloc(m3, sizeof(long));
    ws->gdt_chars    = calloc(21*n,sizeof(long));
    ws->xa3          = calloc(n3, sizeof(float));
    ws->xb3          = calloc(n3, sizeof(float));
    ws->iter_calphas = calloc(mit,sizeof(long));

    if (!ws->calphas2 || !ws->rwa1 || !ws->rwa2 || !ws->l_atoms2 || !ws->ind3 || !ws->out2 ||
        !ws->opt1 || !ws->opt2 || !ws->opt_lcs125 || !ws->gdt_chars || !ws->xa3 || !ws->xb3 ||
        !ws->iter_calphas) { lga_ws_free(&ws); return -1; }

    *pws = ws;
    return 0;
  }

  void lga_ws_free(lga_workspace **pws) {
    if (!pws || !*pws) return;
    lga_workspace *ws = *pws;
    free(ws->calphas2);   free(ws->rwa1);      free(ws->rwa2);
    free(ws->l_atoms2);   free(ws->ind3);      free(ws->out2);
    free(ws->opt1);       free(ws->opt2);      free(ws->opt_lcs125);
    free(ws->gdt_chars);
    free(ws->xa3);        free(ws->xb3);
    free(ws->iter_calphas);
    free(ws);
    *pws = NULL;
  }
  ```

- 入口处接入（`main` 示例）：
  ```c
  // lga.c 内 main 的开始处（在 clean_part(&part) 之后）
  if (part.ws == NULL) {
    if (lga_ws_init(&part.ws) != 0) {
      fprintf(stderr, "ERROR: failed to allocate LGA workspace\n");
      return 1;
    }
  }
  // 在 main 结束前或所有使用点后
  lga_ws_free(&part.ws);
  ```

- 函数内替换示例（`tools_4.c:rms_window`）：
  ```c
  void rms_window(..., part_struct *part, ...) {
    lga_workspace *ws = part->ws;
    long *calphas2 = ws->calphas2;   // 2×(MAXRES+1)
    long *rwa1     = ws->rwa1;       // MAXRES+1
    long *rwa2     = ws->rwa2;       // MAXRES+1

    // 写 calphas[0][rwj] = atoms[0][ ... ]  改为：
    L2(calphas2, 0, rwj) = atoms[0][rwa1[rwi]];
    L2(calphas2, 1, rwj) = atoms[1][rwa2[rwi]];
    // ... 其余类似位置按宏替换
  }
  ```

- 函数内替换示例（`tools_3.c:align_search`）：
  ```c
  long align_search(..., part_struct *part, ...) {
    lga_workspace *ws = part->ws;
    float *xa3 = ws->xa3;  // (MAXRES+1)×3
    float *xb3 = ws->xb3;  // (MAXRES+1)×3

    F2(xa3, i, 0) = pdb[0].atom[atoms[0][j]].R.x;
    F2(xa3, i, 1) = pdb[0].atom[atoms[0][j]].R.y;
    F2(xa3, i, 2) = pdb[0].atom[atoms[0][j]].R.z;
    // ... 其余访问按 F2 宏替换
  }
  ```

- 迁移建议顺序与影响评估：
  1) 先迁移 `tools_4.c` 的 `rms_window/print_lga/check_all_atoms/best_fit` 中的数组（调用频度较低，改动局部）。
  2) 再迁移 `tools_3.c` 的 `align_search/alignment_run/alignment_best_match` 中的数组（与比对核心相关）。
  3) 最后评估是否将 `tools_2.c` 中已改为 `static` 的缓冲也迁到堆，以获得完全的可重入性与更低的常驻内存。

- 优点：
  - 彻底消除对栈大小的依赖；可按需分配与释放；利于多实例并发；更利于将来参数（如 MAXRES）扩展。
- 注意：
  - 需要谨慎替换二维数组访问为索引宏；建议分函数分步骤提交并配套小规模验证用例。

---

## 5. 验证与回归测试建议
- 快速检查：
  - `./lga -h` 不再 Segfault。
  - 对原始能复现崩溃的用例再次运行，确认通过。
- 运行前后对比：
  - 对比两机 `ulimit -s`；在问题机临时执行 `ulimit -s unlimited` 应也可通过（用于侧证栈问题）。
- 性能与内存：
  - `static` 方案：二进制 BSS 增大，常驻内存 ~120–140MB（视编译器对齐略有差异）；单进程 CLI 可接受。
  - `heap` 方案：在用时分配，释放即回收；峰值内存与 `static` 相当，但更灵活。

---

## 6. 风险与注意事项
- `static` 的局限：
  - 不可重入/不可线程安全；但对当前单线程 CLI 影响不大。
  - 多次调用共享同一缓冲，注意不要跨调用持久化指针到这些缓冲。
- `heap` 改造：
  - 替换二维数组为一维 + 宏索引，需小心越界；建议添加断言。
  - `init/free` 必须成对；在异常路径也要释放；建议集中在 `main` 或外层调度处管理生命周期。

---

## 7. 编译与使用建议
- 继续使用现有 Makefile 构建即可：
  ```bash
  cd LGA_package_src/SRC
  make -f Makefile.linux_x86_64_msse2 clean && make -f Makefile.linux_x86_64_msse2
  ```
- 如遇数值稳定性问题（与本次崩溃无关）：可对比 `-ffast-math` 的影响，或添加 `-fno-strict-aliasing`。
- 集群作业脚本建议在运行前加：`ulimit -s unlimited`，作为防御式设置。

---

## 8. 附：关键代码位置（便于审阅）
- `LGA_package_src/SRC/lga.c:335` `static pdata_struct pdata;`
- `LGA_package_src/SRC/lga.c:336` `static part_struct  part;`
- `LGA_package_src/SRC/tools_2.c:21`  `static iter_struct  calphas[MAXITR+1];`
- `LGA_package_src/SRC/tools_2.c:166` `static long  opt1[21][MAXRES], opt2[21][MAXRES];`
- `LGA_package_src/SRC/tools_2.c:278` `static long  calphas[2][MAXRES+1];`
- `LGA_package_src/SRC/tools_2.c:408` `static long  calphas[2][MAXRES+1];`
- `LGA_package_src/SRC/tools_2.c:472` `static long  opt_lcs125[3][MAXRES+1], gdt_chars[21][MAXRES+1], l_atoms[2][MAXRES+1];`
- `LGA_package_src/SRC/tools_3.c:12`  `long  tmp[MAXRES+1];  float xa[MAXRES+1][3], xb[MAXRES+1][3];`
- `LGA_package_src/SRC/tools_3.c:208` `long  calphas[2][MAXRES+1];`
- `LGA_package_src/SRC/tools_4.c:502` `char span1[MAXRES+1], span2[MAXRES+1]; ...`
- `LGA_package_src/SRC/tools_4.c:983` `eval_list eset[MAXRES+1], esup[MAXRES+1];`
- `LGA_package_src/SRC/tools_4.c:1986` `long  calphas[2][MAXRES+1], rwa1[MAXRES+1], rwa2[MAXRES+1];`
- `LGA_package_src/SRC/tools_4.c:2131` `long sup[2][MAXRES+1], sup1[2][MAXRES+1], sup2[2][MAXRES+1];`

---

## 9. 结论
- 问题根因是“栈溢出”，而非平台/库差异；采用“静态迁移 +（可选）堆分配改造”的组合可彻底消除此类问题。

# 检修思路
总体思路

先界定问题类型：同一架构与环境、同一二进制在A机正常、B机崩溃，强烈指向“运行时资源/限制差异”，首选怀疑栈溢出（stack overflow），而非编译器/库/指令集。
用最小代价验证假设：在问题机用 ulimit -s 检查栈；临时 ulimit -s unlimited 再运行。若问题消失，基本锁定为栈问题。
代码审计定位栈压力源：搜索所有以 MAXRES、MAXATOMS 等上限定义的“巨大局部数组/结构体”（函数内分配），按调用路径逐层排查，优先修复最常走的路径。
修复策略分两阶段：
阶段一（最快的止血）：把巨大局部对象从“栈”迁移到“静态存储区（BSS）”，或静态全局，保证功能不变。
阶段二（长期可维护）：设计“堆分配工作区（workspace）”，把所有剩余的大缓冲统一迁移到“堆”，集中管理生命周期，保持算法不变，但消除对栈大小的依赖。
渐进式硬化：优先修复崩溃所处的代码路径（-4 对齐路径 + -stral 报告路径），再补齐其余调用点的防御性初始化，最后消除编译/链接细节问题（如头文件缺失）。
锁定根因

大量超大局部（栈上）对象：
main 中大结构体 pdata_struct 与 part_struct（包含约100MB级数组，如 printlga[MAXRES+1] 每行是 char line[MAXRES+1]）。
tools_2.c（GDT/LCS/ISP 算法）里的临时大缓冲（如 iter_struct calphas[MAXITR+1]、opt1/opt2、calphas[2][MAXRES+1] 等）。
tools_3.c（-4 路径，结构对齐核心）里的 tmp[MAXRES+1]、xa/xb[MAXRES+1][3]、calphas[2][MAXRES+1]、l_atoms/ind/out 等。
tools_4.c（打印与 STRAL）里的 calphas/rwa1/rwa2、sup/sup1/sup2、以及 print_lga 内的 check_mol2 mol[2] 与 stext straline[200]（stext 的一行就是 char line[MAXRES+1]，很大）。
为什么“只在 -4 + -stral 崩”？因为这条路径叠加了多层 MAXRES 级局部数组；虽然你的数据只有 35 residues，但数组按上限（~9501）分配，和输入大小无关。
修复决策与取舍

第一刀（最小侵入、最快止血）：把最显眼的大对象从栈搬出：
lga.c 中 pdata_struct/part_struct 改为 static，BSS 存储，避免进程启动即栈溢出。
tools_2.c 中 ISP/GDT/LCS 的大临时数组改为 static，消除在这条路径上的栈风险。
第二刀（长治久安、跨路径统一）：实现“堆分配工作区”：
设计 lga_workspace，集中管理所有大数组（一次 calloc，在程序结束前 free）。
在 -4 路径与 -stral 输出涉及的函数内，将局部大数组替换为指向 ws 的堆缓冲。
为兼容不同构建/目录（SRC vs SRC2）或不同初始化路径，给核心函数加“防御性初始化”（如果 part->ws == NULL，则就地 lga_ws_init），避免 NULL 解引用。
链接与头文件修补：
tools.h 自身包含 <stdio.h> 以提供 FILE* 声明，解决 workspace.c 单独编译时的 unknown type name 'FILE'。
Makefile 把 workspace.o 加入 OBJS，确保链接。
具体实施（关键文件位置）

栈 → 静态（BSS）
LGA_package_src/SRC/lga.c:335-336：pdata_struct、part_struct 变 static。
LGA_package_src/SRC/tools_2.c：rmsd_isp 的 calphas[MAXITR+1]、rmsd_dist 的 opt1/opt2、lcs_run/gdt_run/lcs_gdt_analysis 的 calphas/opt_lcs125/gdt_chars/l_atoms 改为 static。
引入堆工作区
定义与 API：LGA_package_src/SRC/tools.h:168、245 新增 lga_workspace、lga_ws_init/lga_ws_free、以及 part_struct.ws 指针。
实现：LGA_package_src/SRC/workspace.c 中集中 calloc/free 所有缓冲：
calphas2, l_atoms2, out2, ind3, tmp, rwa1, rwa2, xa3, xb3, sup, sup1, sup2, mol(2个), straline(200个).
链接：LGA_package_src/SRC/Makefile.linux_x86_64* 把 workspace.o 加入 OBJS。
替换 -4 路径局部数组为堆缓冲
LGA_package_src/SRC/tools_3.c
align_search：tmp/xa/xb → ws->tmp/ws->xa3/ws->xb3；并加防御性 lga_ws_init。
alignment_run：calphas[2] → ws->calphas2；并加防御性 lga_ws_init。
alignment_best_match：l_atoms/ind/out → ws->l_atoms2/ws->ind3/ws->out2；并加防御性 lga_ws_init。
LGA_package_src/SRC/tools_4.c
rms_window：calphas/rwa1/rwa2 → ws->calphas2/ws->rwa1/ws->rwa2；并加防御性 lga_ws_init。
best_fit：sup/sup1/sup2 → ws->sup/ws->sup1/ws->sup2；并加防御性 lga_ws_init。
print_lga：check_mol2 mol[2] 与 stext straline[200] → ws->mol 与 ws->straline；并加防御性 lga_ws_init。
主流程初始化与释放（兜底）
LGA_package_src/SRC/lga.c：在进入循环之前 lga_ws_init(&part.ws)，退出时 lga_ws_free(&part.ws)。
验证与定位技巧

栈问题快速验证：
ulimit -s 对比两机；问题机 8MB 左右常见。
ulimit -s unlimited 再运行，若不崩，几乎可判定为栈。
逐路径隔离：
./lga -h OK → 程序初始化没问题。
仅 -4，不带 -stral OK → 输出/窗口计算路径问题（print_lga/rms_window）。
-4 -stral 崩 → 锁定 print_lga 分支，继续查内部使用的临时数组。
栈轨迹/核心转储：
ulimit -c unlimited，运行复现。
gdb ./lga core，bt/frame/list 定位行号；必要时用 -g -fno-omit-frame-pointer 重编译生成更清晰的栈信息。
编译/链接核查：
nm -C lga | rg lga_ws_init 检查是否链接进工作区；strings lga | rg .stral 等确认路径。
对计算结果与性能的影响

结果：不改变算法流程、不改变计算顺序或参数，计算结果应完全一致。若出现差异，通常是之前“未初始化的栈内存”被静态/堆零清掉后暴露了原本隐患。
性能：数组仍为连续内存，访问模式不变；仅改变分配位置（栈→BSS/堆）。对运行时间影响可以忽略。
内存：BSS（静态）与堆的总占用与原实际使用规模相当，但变得明确、稳定；再无“栈峰值”风险。堆模式还便于将来按需分配/释放。
为什么问题只在某些机器出现

默认栈大小是系统/用户/容器级的运行时配置，不属于“软件环境”的一部分。两台机器即使“架构/库/编译器”相同，栈限制不同也会造成截然不同的运行结果。
这也是为什么“问题机上编译出的 lga 在好机器能跑”。二进制没问题，崩溃是运行时资源限制使然。
通用工程经验与最佳实践

避免在函数栈上分配巨型数组/结构。对巨型缓冲统一使用堆或静态。优选堆分配工作区（workspace）这种集中管理手法。
老代码中常见的 MAXRES/MAXATOMS 上限数组应评估其真实需求：若可能，按输入规模动态分配；若必须保留上限，尽量放入堆。
在输出/报表路径（如 -stral）同样注意大对象：stext line[MAXRES+1] 这种一行就近万个字符，200 行就是约 19 MB，需要堆化或复用缓冲。
给关键路径加防御性初始化，使模块在不同构建/入口场景下都能安全运行（例如在函数内检测 ws 是否为 NULL）。
将“资源限制检测/提示”加入启动或日志（例如打印 ulimit -s），方便部署时快速定位。
调试期可用 -g、-fno-omit-frame-pointer、ASan (-fsanitize=address) 等，但上线时去除以确保性能与兼容性。
你可以复盘的步骤（从零到解）

观察现象 → 怀疑栈 → ulimit 验证。
全局搜“栈大户”：
rg -n "\[MAXRES\+1\]"、rg -n "MAXATOMS"、rg -n "typedef.*struct.*\[MAXRES\+1\]" 快速定位大数组。
优先修主干路径（-4）：
tools_3.c 栈数组 → 堆；tools_4.c 报表/窗口缓冲 → 堆。
给函数加防御性初始化，避免不同构建目录/入口遗漏全局初始化。
修头文件/Makefile 细节，保证编译顺畅。
回归运行：-h、-4、-4 -stral，通过后再做压力测试。
后续扩展（可选）

若要线程并发，将 lga_workspace 改为“每线程一份”，在外层以线程上下文传入；接口上不再依赖全局 part->ws。
进一步细化工作区：将 MAXRES 改为基于输入动态尺寸的“最小必要大小”，减少常驻/峰值内存。
对 -ffast-math 等优化做一次对比基线验证，确保数值稳定性满足你的业务（非这次问题根因，但长期建议）。
以上就是完整的思路与落地步骤。你已经完成了从“定位资源类崩溃”到“工程化消除风险”的关键链路，这套方法也适用于大量传统 C 科学计算代码的移植与稳健化。
