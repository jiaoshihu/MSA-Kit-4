# MSA-Kit

**English** | [简体中文](README.zh-CN.md)

Use one project to understand multiple sequence alignment, produce your first result, and then develop and evaluate your own algorithm.

MSA-Kit is a DNA/RNA alignment program and a reusable C++ template, currently **version 0.3.0**. Work through the synthetic example on this page, build and run the program, then continue to the developer guide when you want to change its algorithm. You do not need C++ knowledge to run these exercises.

**Reading path:** [What MSA does](#what-is-msa) → [FASTA input](#fasta-input) → [Setup](#setup) → [Build](#build-program) → [Run and read the result](#first-result) → [Your own input](#your-own-input) → [Result checks](#check-result) → [Develop an algorithm](#develop-an-algorithm) → [Website evaluation](#easymsa-benchmark).

<a id="what-is-msa"></a>

## 1. What are we trying to align?

A biological sequence is represented as a string of letters. DNA commonly uses A, C, G and T; RNA uses U in place of T. A **pairwise alignment** compares two sequences. **Multiple sequence alignment (MSA)** usually compares three or more, arranging their letters into columns so you can examine corresponding positions. MSA-Kit also accepts one or two records for small tests.

For a first exercise, use the supplied synthetic data. For biological work, begin with related sequences from the same gene or comparable region, in the same orientation. Merely putting unrelated sequences in a file does not make the resulting columns biologically meaningful. MSA-Kit does not select homologous sequences, reverse-complement them, or perform read assembly for you.

Consider three short inputs:

```text
Name      Sequence      Length
center    AC            2
q1        ATC           3
q2        AGGC          4
```

The reference algorithm in this project produces:

```text
Column    1 2 3 4
center    A - - C
q1        A T - C
q2        A G G C
```

Read this **by column**: the A letters occupy column 1, the C letters column 4, and gaps make room for letters between them. Every aligned row now has four positions. Removing `-` from each row recovers its input sequence.

A gap is a position where a row has no nucleotide in this alignment. Its placement expresses a proposed correspondence between sequences; this toy result does not by itself establish a particular evolutionary insertion or deletion. Different algorithms can produce different gap placements.

| Symbol or term | Meaning here |
| --- | --- |
| `A`, `C`, `G`, `T`, `U` | Nucleotide letters. |
| `N` | An unspecified nucleotide; it still occupies a sequence position. |
| `-` (gap) | No nucleotide for this row at this alignment position. |
| Match | Two compared letters agree under the chosen scoring rule. |
| Mismatch | Two compared letters differ under that rule. |
| Alignment column | One position across all aligned rows, including gaps. |

`N` and `-` have different meanings. Do not replace one with the other to make rows look aligned.

<a id="fasta-input"></a>

## 2. Read the input file: FASTA

Open [examples/tiny.fasta](examples/tiny.fasta). It contains:

```fasta
>center
AC
>q1
ATC
>q2
AGGC
```

Each `>` starts a new record. The rest of that line is its **header** (a name, optionally followed by a description). The following sequence lines belong to that record until the next header. A record may span several sequence lines; the lines are joined when read. Three records in one FASTA file form the input for this example.

For your first files, choose distinct simple names, use plain-text UTF-8, and put only nucleotide letters on sequence lines. Input rows may have different lengths. Supply **unaligned** input: this program rejects existing `-` and `.` in parsed sequences. A file extension such as `.fasta` names the format; it does not convert a spreadsheet or a Word document into FASTA.

The output will also be FASTA, but its sequence rows contain inserted gaps and have equal lengths. Input FASTA, output aligned FASTA and a table of benchmark scores are three different things.

<a id="setup"></a>

## 3. Prepare a Linux terminal

Run the commands below in Linux. On Windows, use a Linux environment such as Ubuntu in WSL and enter its terminal first. The project uses Linux file operations and does not build as a native Windows application.

You need Git to obtain the project, a C++17 compiler and CMake to build it, zlib development files for gzip support, and Python 3.8 or newer for tests. The standalone C++ program does not need Python; the separate website adapter does.

On Ubuntu/Debian, if these packages are missing, an administrator can install them with:

```bash
sudo apt-get update
sudo apt-get install -y git build-essential cmake zlib1g-dev python3
```

Package installation and downloading the repository use the network. Once those are ready, MSA-Kit's normal build uses its bundled dependencies without downloading more libraries. On a shared server, ask its administrator to provide missing packages rather than changing the server environment yourself.

Get the project, then work from its root directory:

```bash
git clone https://github.com/malabz/MSA-Kit.git
cd MSA-Kit
```

If you already have a copy, enter that directory instead. You should see `CMakeLists.txt`, `src/` and `examples/`. Commands below assume you are in this directory. Text inside a code block is the command; do not type a shell prompt such as `$` before it.

<a id="build-program"></a>

<a id="build-and-test"></a>
<a id="quick-start"></a>

## 4. Build the program

```bash
cmake -S . -B build-learning -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build-learning --parallel 2
./build-learning/msakit --version
```

`-S .` selects the current source directory. `-B build-learning` keeps generated files in a separate build directory. The second command compiles the program. `./build-learning/msakit` is the resulting executable; `--version` should print `msakit 0.3.0`.

Compiler progress messages are normal. If configuration or compilation fails, resolve that error before proceeding. A successful build gives you a program; it has not yet aligned your data.

<a id="first-result"></a>

## 5. Run and understand the supplied example

```bash
./build-learning/msakit -i examples/tiny.fasta -o build-learning/tiny.aligned.fasta --force
cat build-learning/tiny.aligned.fasta
python3 -c 'from pathlib import Path; assert Path("examples/tiny.expected.fasta").read_text() == Path("build-learning/tiny.aligned.fasta").read_text(); print("Example: PASS")'
```

- `-i` identifies the input file.
- `-o` identifies the output file.
- `--force` permits replacing this example's existing output when you repeat the exercise. Leave it out when you want existing files protected.
- `cat` displays the saved result; progress messages from the aligner go to the terminal's error stream, separately from FASTA.
- The final Python command compares the saved result with the reference-core example, treating LF and CRLF line endings equivalently. It prints `Example: PASS` when they agree.

Expected file contents:

```fasta
>center
A--C
>q1
AT-C
>q2
AGGC
```

This is the four-column alignment from section 1. The original names and their order are retained. `center` still has its original A and C; the algorithm inserted gaps, not new nucleotides. This exact expected layout is for the supplied reference algorithm. A replacement algorithm can return a different valid layout.

Before continuing, try answering: how many alignment columns are there? What is q1 after removing its gap? Does producing this file establish biological accuracy?

<details>
<summary>Check your answers</summary>

There are four columns. q1 becomes ATC. This synthetic exercise demonstrates the reference program's output; biological quality requires judging the proposed correspondence, with suitable evidence or a reference where available.

</details>

<a id="your-own-input"></a>

## 6. Try a small file of your own

Create three more short synthetic DNA sequences:

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

Read each aligned row and remove its gaps by eye. You should recover ACGT, AGT and ACGTT in that order, and all aligned rows should have the same length. The `FASTA` lines delimit the shell's file-writing block; they are not part of the file. Use a new output path for real data, and keep your original inputs.

Before using ambiguity symbols, learn the cleaning rule. Codes such as R (A or G) are valid nucleotide ambiguity symbols, but the ordinary MSA-Kit CLI keeps only A/C/G/T/U/N, uppercases their lowercase forms, and maps other returned sequence bytes to N. For example, `acRUn` becomes `ACNUN`. The CLI's preservation check uses this **cleaned sequence**. Start with A/C/G/T/U/N, then read [input preparation](docs/guide.md#3-input-preparation) before using your own ambiguous data. The website adapter separately restores original ambiguity symbols, as explained in its tutorial.

<a id="check-result"></a>

## 7. Check that you understand the result

Ask these questions about an alignment:

1. Are all input records present, with the expected names and order?
2. Are all aligned rows equally long and nonempty?
3. Does removing `-` recover each corresponding cleaned input sequence?
4. Are the characters valid, with no column containing gaps in every row?

MSA-Kit checks these properties before writing a successful result. They establish that the program preserved the required data and produced a well-formed alignment. Biological quality asks a further question: do the proposed columns correspond to meaningful positions in related sequences? A structurally valid result can still be a poor biological alignment.

Run the supplied software checks:

```bash
(cd build-learning && ctest --output-on-failure)
```

The supplied project registers four CTest entries. A successful complete run reports all four passing. The tests use small synthetic data; they check software behavior and the reference algorithm, rather than certifying biological accuracy on your dataset.

<a id="develop-an-algorithm"></a>

## 8. Develop an algorithm inside this program

After completing the exercise, open [src/align.cpp](src/align.cpp) and find **Reference star alignment — replaceable core**. The supplied algorithm chooses the first sequence as the center, uses Needleman–Wunsch global pairwise alignment, and merges insertions around that center.

In the example, the slot between A and C needs room for GG, so its width is two. q1's T is left-aligned there to give AT-C. The core stores full dynamic-programming tables: it is useful for learning and small checks, but long inputs can have substantial time and memory cost.

Writing an algorithm requires C++ knowledge, while readers and the CLI can be reused. The framework supplies prepared records to `align()`, you return aligned strings in input order, and the outer layers validate and write them. Continue with the [developer guide](docs/guide.md), which traces the same example through this interface before you replace the algorithm.

After source changes, rebuild and run the reusable contract checks:

```bash
(cd build-learning && ctest -L contract --output-on-failure)
```

`contract` tests reusable rules; `reference` tests the supplied star algorithm's scores and exact layouts. Add tests for your replacement strategy and adapt algorithm-specific expectations.

<a id="easymsa-benchmark"></a>

## 9. Evaluate a program on the website

To see the process before running a full experiment, open [EasyMSA's submission page](https://easymsa.cn/#/benchmark/submit) and choose the MSA-Kit example. It displays already completed real results without starting a new benchmark job.

When you have a program to evaluate, follow the [Benchmark integration tutorial](benchmark/README.md). You provide a public repository with `benchmark/build.sh` and `benchmark/run.sh`; the website fixes a source commit, builds it, checks small inputs, then runs supported datasets. The supplied adapter preserves complete titles and restores original IUPAC characters after validating the core output.

| Result | How to read it |
| --- | --- |
| Time | Runtime for that measurement; lower means faster. |
| Peak memory | Maximum recorded memory use; lower means less memory. |
| Q / TC | Reference-based alignment quality measures; higher is better under the report's definition. Read them where a reference is available. |
| Normalized SP | A reported alignment-scoring measure; higher is better under the same scoring rules. It is not an accuracy percentage or a probability. |
| Timeout / memory limit / invalid output | Actual run outcomes, not zero-quality scores. |

Compare methods on the same dataset and inspect their execution/scoring conditions. Speed and quality answer different questions. A reference alignment is an existing target used for quality assessment, not the expected output of every algorithm. Q/TC need a reference; their absence is different from a zero score. The template's simple star core can time out or exhaust memory on large inputs even after passing interface checks.

Use **Download data** for unaligned dataset FASTA; use **Export results** for measured CSV/JSON. Your own aligned FASTA is a separate output of the alignment program.

<a id="troubleshooting"></a>

## 10. When something goes wrong

| What you see | First action |
| --- | --- |
| `cmake`, `g++` or `python3`: command not found | Return to the Linux environment and prerequisites in section 3. |
| CMake reports a missing zlib | Install/provide the zlib development package, not only the runtime library. |
| `./build-learning/msakit`: not found | Check your current directory and whether the build completed. |
| Input file cannot be opened | Check `-i` and the path relative to the repository root. |
| Output already exists | Choose a new output name, or use `--force` if replacement is intended. |
| An existing gap is rejected | Supply the original unaligned sequences, not a previous aligned output. |
| A letter became N | Review the CLI cleaning rule in section 6. |
| A large input fails or takes too long | Reproduce with a small subset and inspect lengths/counts; the full-DP teaching core has substantial memory/time cost. |
| A changed algorithm fails only `reference` tests | Review those algorithm-specific expectations and run the reusable `contract` group. |

For further reading, use [EMBL-EBI's introduction to MSA](https://www.ebi.ac.uk/training/online/courses/guide-to-sequence-analysis-tools/sequence-alignment/multiple-sequence-alignment/) and [NCBI's nucleotide FASTA guide](https://www.ncbi.nlm.nih.gov/genbank/fastaformat). The latter also includes GenBank-specific submission rules; MSA-Kit's exact parser behavior is documented in its [design specification](docs/design.md#4-input-contract).


## Documentation map

| What you want to do next | Document |
| --- | --- |
| Understand the framework and replace the core | [Developer guide](docs/guide.md) |
| Prepare scripts and submit to the website | [Benchmark integration](benchmark/README.md) |
| Look up precise input/output, CLI and errors | [Design specification](docs/design.md) |
| Write and extend software tests | [Testing guide](docs/testing.md) |
| Inspect actual environments and recorded checks | [Validation record](docs/validation.md) |
| Understand small synthetic test files | [Fixture notes](tests/data/msa/README.md) |
| Inspect sources, dependency versions and licensing | [Third-party notices](THIRD_PARTY_NOTICES.md) |

## Capabilities, dependencies and sources

The program reads FASTA, FASTQ and mixed records recognized by pinned kseq, through plain/gzip files or stdin. FASTQ quality lengths are checked, then qualities are discarded. Output is aligned FASTA; named files ending in `.gz` are compressed and stdout remains plain FASTA. Default sequence width is 80; `--line-width 0` disables wrapping. Headers are rebuilt from kseq name/comment fields, preserving record order and duplicates. Diagnostics go to stderr. Named output uses overwrite protection, same-file/link identity checks, checked temporary writes and atomic publication.

Linux, C++17, CMake >=3.18 and system zlib are required. CLI11 2.6.1, spdlog 1.16.0, doctest 2.4.12 and pinned kseq are bundled; normal builds do not download dependencies. Tests require Python >=3.8. The C++ executable needs no Python at runtime; the separate website adapter does. `BUILD_TESTING=OFF` builds just the program; optional installation uses `cmake --install build-learning --prefix <prefix>`.

The supplied single-threaded reference core supports neither protein nor profile alignment and does not guarantee optimal biological MSA. The design draws on [HAlign-5 b0d5608811d44e30bc323da54f0583ceaaf684d2](https://github.com/pinglu-zhang/HAlign-5/tree/b0d5608811d44e30bc323da54f0583ceaaf684d2). MSA-Kit uses the [MIT license](LICENSE); preserve upstream and dependency notices when copying it. [The dependency manifest](extern/provenance.json) records versions and checksums.
