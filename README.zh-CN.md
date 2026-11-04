# MSA-Kit

[English](README.md) | **简体中文**

用同一个项目理解多序列比对、完成第一份结果，再逐步开发和评测自己的算法。

MSA-Kit 是 DNA/RNA 多序列比对程序和可复用的 C++ 模板，当前版本 **0.3.0**。先用本页的合成数据认识 MSA，再编译并运行程序；想开发算法时，继续阅读开发指南。运行下面的练习不要求先会 C++。

**阅读顺序：**[什么是 MSA](#what-is-msa) → [FASTA 输入](#fasta-input) → [准备环境](#setup) → [构建](#build-program) → [运行与读结果](#first-result) → [自己的输入](#your-own-input) → [结果检查](#check-result) → [开发算法](#develop-an-algorithm) → [网站评测](#easymsa-benchmark)。

<a id="what-is-msa"></a>

## 1. 我们到底要比对什么？

生物序列可以写成一串字母。DNA 常用 A、C、G、T，RNA 用 U 代替 T。**两两比对**比较两条序列；**多序列比对（Multiple Sequence Alignment，MSA）**通常比较三条或更多序列，把字母排成列，以便观察可能对应的位置。MSA-Kit 也接受一条或两条记录，用于小型测试。

第一次练习请使用仓库自带的合成数据。处理生物数据时，先选择同一基因或可比区域、方向一致的相关序列。把无关序列放在同一个文件里，不会让对齐列自动具有生物学意义。MSA-Kit 不负责寻找同源序列、反向互补或测序组装。

先看三条很短的输入：

```text
名称      序列          长度
center    AC            2
q1        ATC           3
q2        AGGC          4
```

本项目的参考算法会产生：

```text
列号      1 2 3 4
center    A - - C
q1        A T - C
q2        A G G C
```

要**竖着看列**：A 都在第 1 列，C 都在第 4 列，gap 给中间的字母留出空间。现在每行都有 4 个位置；逐行去掉 `-`，还能得到各自的输入序列。

gap 表示这条序列在当前对齐位置没有碱基。gap 的位置表达的是序列之间的一种对应关系；仅凭这个玩具结果，不能断定发生过哪一次进化插入或删除。不同算法可能选择不同的 gap 位置。

| 符号或术语 | 在这里的含义 |
| --- | --- |
| `A`、`C`、`G`、`T`、`U` | 碱基字母。 |
| `N` | 具体种类未确定的碱基，仍占一个序列位置。 |
| `-`（gap，空位） | 该序列在当前对齐位置没有碱基。 |
| 匹配 | 在所选评分规则下，比较的字母相符。 |
| 错配 | 在该规则下，比较的字母不同。 |
| 对齐列 | 对齐后所有行的同一位置，可包含碱基和 gap。 |

N 和 gap 含义不同，不能为了让行看起来对齐而互相替换。

<a id="fasta-input"></a>

## 2. 读懂输入文件：FASTA

打开 [examples/tiny.fasta](examples/tiny.fasta)：

```fasta
>center
AC
>q1
ATC
>q2
AGGC
```

每个 `>` 开始一条新记录，同行后面的文字是**标题**，可以是名称及描述。后续序列行属于这条记录，直到下一个标题。序列可以分多行写，读取时会拼接。本例把三条记录放在一个 FASTA 文件中，作为一次 MSA 的输入。

初次练习请使用互不相同的简单名称、UTF-8 纯文本，并只在序列行写碱基字母。输入序列可以不等长。请提供**未比对输入**：程序拒绝解析后序列中已有的 `-`、`.`。把 Excel 或 Word 文件改名为 `.fasta`，不会把内容转换成 FASTA。

输出仍是 FASTA，但序列中可以有新增 gap，所有对齐行等长。输入 FASTA、输出的 aligned FASTA、Benchmark 指标表是三类不同文件。

<a id="setup"></a>

## 3. 准备 Linux 终端

以下命令在 Linux 中执行。Windows 用户可先进入 WSL 的 Ubuntu 等 Linux 环境及其终端；项目依赖 Linux 文件操作，不能作为原生 Windows 程序构建。

需要 Git 获取项目、C++17 编译器和 CMake 编译程序、zlib 开发文件支持 gzip、Python >=3.8 运行测试。普通 C++ 可执行程序不依赖 Python，独立的网站适配器需要 Python。

在 Ubuntu/Debian 上，缺少这些包时，可由管理员安装：

```bash
sudo apt-get update
sudo apt-get install -y git build-essential cmake zlib1g-dev python3
```

安装系统包和下载仓库需要网络。准备完成后，MSA-Kit 的普通构建使用附带的依赖，不再下载库。在共享服务器上，请让管理员提供缺少的包。

获取项目，然后进入仓库根目录：

```bash
git clone https://github.com/malabz/MSA-Kit.git
cd MSA-Kit
```

已有副本时，直接进入该目录即可。你应能看到 `CMakeLists.txt`、`src/`、`examples/`。下文命令都从这里执行。代码块中的内容就是命令，不需要额外输入终端提示符 `$`。

<a id="build-program"></a>

<a id="build-and-test"></a>
<a id="quick-start"></a>

## 4. 构建程序

```bash
cmake -S . -B build-learning -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build-learning --parallel 2
./build-learning/msakit --version
```

`-S .` 指当前源码目录；`-B build-learning` 把生成文件放进单独构建目录。第二条命令进行编译，`./build-learning/msakit` 就是生成的程序，版本输出应为 `msakit 0.3.0`。

编译进度属于正常输出。配置或编译失败时，先解决错误再往下做。编译成功表示得到了程序，还没有处理你的序列。

<a id="first-result"></a>

## 5. 运行自带示例，读懂结果

```bash
./build-learning/msakit -i examples/tiny.fasta -o build-learning/tiny.aligned.fasta --force
cat build-learning/tiny.aligned.fasta
python3 -c 'from pathlib import Path; assert Path("examples/tiny.expected.fasta").read_text() == Path("build-learning/tiny.aligned.fasta").read_text(); print("Example: PASS")'
```

- `-i` 后面是输入文件。
- `-o` 后面是输出文件。
- `--force` 允许重复练习时覆盖这个示例输出；希望保护已有文件时不加它。
- `cat` 显示保存的结果；比对程序的进度写到终端错误流，与 FASTA 分开。
- 最后一条 Python 命令与参考核心示例答案比较，兼容 LF 和 CRLF 换行；相同则打印 `Example: PASS`。

预期文件内容：

```fasta
>center
A--C
>q1
AT-C
>q2
AGGC
```

它就是第 1 节的四列对齐。名称和顺序保留，center 仍只有原来的 A、C；算法增加的是 gap，没有增加碱基。这份具体布局属于附带的参考算法，更换算法后可以得到另一种合法布局。

继续前试着回答：一共有几列？q1 去掉 gap 后是什么？产生这个文件是否说明生物学比对准确？

<details>
<summary>查看答案</summary>

一共四列，q1 去掉 gap 后是 ATC。这组合成练习展示参考程序的输出；判断生物学质量，还需要用合适证据或可用的参考比对评估位置对应关系。

</details>

<a id="your-own-input"></a>

## 6. 自己准备一个小文件

再创建三条短的合成 DNA 序列：

```bash
cat > build-learning/my-input.fasta <<'FASTA'
>seq_a
ACGT
>seq_b
AGT
>seq_c
ACGTT
FASTA
./build-learning/msakit -i build-learning/my-input.fasta -o build-learning/my-alignment.fasta --force
cat build-learning/my-alignment.fasta
```

逐行去掉结果中的 gap，应按顺序得到 ACGT、AGT、ACGTT，而且所有对齐行等长。两行 `FASTA` 是 Shell 写文件块的边界，不会写进文件。处理真实数据时使用新的输出路径，保留原输入。

使用模糊碱基前先了解清洗：R（A 或 G）等字符是合法的碱基模糊符号，但普通 MSA-Kit CLI 只保留 A/C/G/T/U/N，将对应小写大写化，其余返回的序列字节转为 N。例如 `acRUn` 会变成 `ACNUN`。普通 CLI 按这个**清洗后的序列**检查保真。先使用 A/C/G/T/U/N 练习；使用自己的模糊字符数据前，阅读 [输入准备](docs/guide.zh-CN.md#3-进入核心前序列会怎样变化)。网站适配器还会单独恢复原始模糊字符，接入教程会解释这一点。

<a id="check-result"></a>

## 7. 怎样判断自己读对了结果？

先问四个问题：

1. 输入的记录是否都在，名称和顺序是否符合预期？
2. 所有对齐行是否非空、等长？
3. 每行去掉 `-` 后，能否还原对应的清洗后输入？
4. 字符是否合法，是否没有整列都为 gap 的位置？

MSA-Kit 会在成功写出前检查这些性质，说明程序保留了要求的数据、结果结构合法。生物学质量还要进一步判断：这些列是否对应相关序列中有意义的位置？结构合法的结果也可能不是好的生物学比对。

运行附带的软件检查：

```bash
(cd build-learning && ctest --output-on-failure)
```

原模板注册了四个 CTest 入口，完整成功运行应报告四个都通过。这些测试使用小型合成数据，检查软件行为和参考算法，并不据此证明你的数据集具有很高生物学准确率。

<a id="develop-an-algorithm"></a>

## 8. 从这个程序开始开发自己的算法

完成上面的练习后，打开 [src/align.cpp](src/align.cpp)，找到 **Reference star alignment — replaceable core** 标记。附带算法选择第一条序列为中心，用 Needleman–Wunsch 做全局两两对齐，再合并中心的插入槽。

本例 A、C 之间需要容纳最长插入 GG，所以槽宽为 2；q1 的 T 左对齐得到 AT-C。核心使用完整动态规划表，适合学习和小型检查，长输入可能有较大时间和内存代价。

开始编写算法需要 C++ 知识，但无需重新实现文件读取或命令行。框架把准备好的记录交给 `align()`，你按输入位置返回对齐字符串，外围负责读写及验证。接着阅读 [开发指南](docs/guide.zh-CN.md)，跟踪同一组示例如何经过核心接口，再替换算法。

修改代码后重新编译，并运行通用约定检查：

```bash
(cd build-learning && ctest -L contract --output-on-failure)
```

contract 检查可复用规则；reference 检查原星比对的评分和具体布局。更换算法时，为新策略添加测试，并调整算法专属答案。

<a id="easymsa-benchmark"></a>

## 9. 在网站评测一个程序

想先体验流程，可打开 [EasyMSA 提交页](https://easymsa.cn/#/benchmark/submit)，选择 MSA-Kit 示例，直接查看已完成的真实成绩，不启动新的评测任务。

有了自己的程序，再按 [Benchmark 接入教程](benchmark/README.zh-CN.md) 准备公开仓库及 `benchmark/build.sh`、`benchmark/run.sh`。网站固定源码提交、构建、预检小样本，再运行支持的数据集。附带适配器在检查核心输出后恢复完整标题及原始 IUPAC 字符。

| 结果项 | 怎样理解 |
| --- | --- |
| 时间 | 此次测量的运行耗时，低表示更快。 |
| 峰值内存 | 记录的最大内存使用量，低表示更省内存。 |
| Q／TC | 基于参考比对的质量指标，在报告定义下越高越好，有参考时查看。 |
| 归一化 SP | 报告中的对齐评分指标，在相同评分规则下越高越好；不是准确率百分比或概率。 |
| 超时／内存超限／输出无效 | 实际运行状态，不是质量得分为零。 |

在相同数据集上比较方法，并查看执行与评分条件。速度和质量回答不同问题。参考比对是质量评估使用的已有目标，不代表每个算法都必须输出一模一样的结果。Q/TC 需要参考，缺少指标与零分不同。简单星比对即使通过接口测试，在大数据上仍可能超时或内存不足。

“下载数据”提供未比对的 FASTA；“导出成绩”提供测量得到的 CSV/JSON。你自己生成的 aligned FASTA 是程序的另一项输出。

<a id="troubleshooting"></a>

## 10. 遇到问题先看这里

| 现象 | 先检查什么 |
| --- | --- |
| 找不到 cmake、g++ 或 python3 | 回到第 3 节，检查 Linux 环境和依赖。 |
| CMake 提示缺少 zlib | 需要开发包，只有运行时库不够。 |
| 找不到 `./build-learning/msakit` | 当前目录是否正确、编译是否完成。 |
| 打不开输入文件 | 检查 `-i` 和相对于仓库根目录的路径。 |
| 输出文件已存在 | 换新名称，或在确定要覆盖时加 `--force`。 |
| 提示已有 gap | 输入应是原始未比对序列，而不是上次的对齐结果。 |
| 字母变成了 N | 看第 6 节的 CLI 清洗规则。 |
| 大文件失败或很慢 | 先用小子集复现，查看条数和长度；完整 DP 教学核心的时间/内存代价较高。 |
| 换算法后只有 reference 失败 | 检查算法专属答案，并运行可复用的 contract 组。 |

可继续阅读 [EMBL-EBI 的 MSA 入门](https://www.ebi.ac.uk/training/online/courses/guide-to-sequence-analysis-tools/sequence-alignment/multiple-sequence-alignment/) 与 [NCBI 的核酸 FASTA 说明](https://www.ncbi.nlm.nih.gov/genbank/fastaformat)。NCBI 页面还包含 GenBank 提交专属规则；MSA-Kit 的精确解析行为以 [设计规范](docs/design.zh-CN.md#4-input-contract) 为准。


## 文档导航

| 你要继续做什么 | 文档 |
| --- | --- |
| 理解框架，替换算法核心 | [开发指南](docs/guide.zh-CN.md) |
| 编写脚本并提交网站评测 | [Benchmark 接入](benchmark/README.zh-CN.md) |
| 查询精确输入输出、CLI 和错误码 | [设计规范](docs/design.zh-CN.md) |
| 编写与扩展软件测试 | [测试指南](docs/testing.zh-CN.md) |
| 查看已执行的环境和测试结果 | [验证记录](docs/validation.zh-CN.md) |
| 理解小型合成测试数据 | [夹具说明](tests/data/msa/README.zh-CN.md) |
| 查看来源、固定依赖及许可 | [第三方声明](THIRD_PARTY_NOTICES.zh-CN.md) |

## 能力、依赖与来源

程序支持普通/gzip 文件及 stdin 的 FASTA、FASTQ 和固定 kseq 识别的混合记录；FASTQ 质量长度检查后丢弃。输出为 aligned FASTA，命名文件以 `.gz` 结尾时压缩，stdout 始终普通 FASTA。默认每行 80 字符，`--line-width 0` 不折行。标题由 kseq name/comment 重建，保留输入顺序及重复记录；日志在 stderr。命名输出具备覆盖保护、同文件与链接身份检查、临时写入完成后原子发布。

使用 Linux、C++17、CMake >=3.18 和必需系统 zlib。附带 CLI11 2.6.1、spdlog 1.16.0、doctest 2.4.12、固定 kseq；正常构建不下载依赖。测试需要 Python >=3.8；普通 C++ 可执行程序运行不依赖 Python，独立网站适配器需要 Python。`BUILD_TESTING=OFF` 只构建程序，可用 `cmake --install build-learning --prefix <prefix>` 安装。

附带单线程参考核心不支持蛋白或 profile，不保证生物学最优 MSA。设计参考 [HAlign-5 b0d5608811d44e30bc323da54f0583ceaaf684d2](https://github.com/pinglu-zhang/HAlign-5/tree/b0d5608811d44e30bc323da54f0583ceaaf684d2)，采用 [MIT 许可](LICENSE)，复制时保留上游和依赖声明。[依赖清单](extern/provenance.json) 记录版本及校验值。
