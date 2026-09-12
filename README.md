# CPM-LDPC distance verification

Exponent arrays, verification programs, and archived search records for
*Classical CPM-LDPC Codes Attaining the Minimum-Distance Bound* by Kenta Kasai.
These are binary CPM-LDPC codes, also called QC-LDPC codes.

## Code definitions and results

For a J-by-L exponent array E and a positive integer P, let I(e) be the P-by-P
binary matrix with I(e)[s,t] = 1 exactly when s = t + e (mod P), for
0 <= s,t < P. Form H from the blocks I(E[j][ell]). The code consists of all
binary vectors c satisfying Hc = 0. Its length is n = LP, column weight is J,
row weight is L, and dimension is k = n - rank(H) over GF(2).
All array entries, including the zero first row and column, are given in
[codes.json](codes.json). The full computation metadata are in
[joint_examples.json](joint_examples.json).

| J | L | P | n | k | Minimum distance d |
|---:|---:|---:|---:|---:|:---|
| 3 | 4 | 24 | 96 | 26 | 24 |
| 3 | 5 | 45 | 225 | 92 | 24 |
| 3 | 6 | 71 | 426 | 215 | 24 |
| 3 | 7 | 111 | 777 | 446 | 24 |
| 3 | 8 | 159 | 1272 | 797 | 24 |
| 4 | 5 | 23 | 115 | 26 | 30 |
| 4 | 6 | 29 | 174 | 61 | 28 |
| 4 | 7 | 43 | 301 | 132 | 32 |
| 4 | 8 | 73 | 584 | 295 | 30 <= d <= 46 |

The distance of the J=4, L=8 example remains between 30 and 46; its exact
value is not known. The values of P in this table are not claimed to be minimal.
The first entry, J=3, L=4, P=24, is the [96,26,24] example used to illustrate
the parity-check matrix, a weight-24 codeword, and the graph of paired positions.

## Get the published version

```sh
git clone --branch v1.1.3 https://github.com/kasaikenta/cpm-ldpc-distance.git
cd cpm-ldpc-distance
```

Python 3.10+ and a C++17 compiler named `c++` are required. Python uses only the
standard library. The repository includes the three compressed record archives,
so there is no separate data download or Git LFS setup.

## Check the data and implementation

```sh
python3 verification/check_records.py
python3 verification/example34.py
python3 verification/independent24.py
python3 verification/enumerate_j4_l5.py
python3 verification/review_checks.py
```

- `check_records.py` checks the nine arrays, binary ranks, explicit codewords,
  archive hashes, and completion and coverage of the recorded search partitions.
  It reads archives directly without unpacking thousands of files.
- `example34.py` checks the first table entry and generates the matrix bitmap
  and TikZ graph in `figures/`. The 72 highlighted matrix entries correspond
  exactly to the 36 graph edges, with matching colors for each block row.
- `independent24.py` independently enumerates all 67,108,863 nonzero words of
  the [96,26,24] example; it finds 668 words of weight 24.
- `enumerate_j4_l5.py` independently enumerates all 67,108,863 nonzero words of
  the [115,26,30] example; it finds 23 words of weight 30.
- `review_checks.py` compares the search program with exhaustive enumeration
  for 54 small codes and checks explicit cofactor codewords.

Checking stored search records does not repeat the large searches. The records
are solver logs, not independently checkable nonexistence proofs. The algorithm,
pruning arguments, and precise scope of each check are explained in the
[verification documentation](verification/README.md).

## Repeat a distance search

```sh
python3 verification/reproduce.py --J 3 --L 4
python3 verification/reproduce.py --J 3 --L 5 --threads 4 --seconds 60
python3 verification/reproduce.py --J 4 --L 6 --threads 4 --seconds 60
```

Repeat the same command to resume saved search frontiers. Completed partitions
are skipped. A time-limited run reports its incomplete scope; a lower bound is
reported only after all assigned partitions finish. Large cases can require
substantial computation. J=3, L=4 and J=4, L=5 use independent full enumeration instead.

The search program uses position index `t*L+ell`. This is a fixed permutation
of the block-column order `ell*P+t` in the definition of H and preserves distance.
The legacy input-format identifier `CPM_CSS_V1` does not change the code being
checked: every command above selects classical mode with no stabilizer rows.

## Files and version

- `codes.json`: compact full arrays and current parameters.
- `joint_examples.json`: table data and computation metadata.
- `verification/`: Python checks and the C++ search and enumeration sources.
- `supporting_material/`: baseline J=3 and J=4 records and the J=4 refinements.
- `figures/`: matrix and codeword PNG bitmaps, a complete TikZ graph, and the
  support and graph data for the J=3, L=4 example.
  Generate these files with `python3 verification/example34.py`; only typesetting
  the TikZ source requires LaTeX and TikZ.
- `SHA256SUMS`: hashes of the published source, data, and documentation files.

The tag `v1.1.3` identifies the verification package for the table above.
The earlier tag `v1.0.0` retains the original eight-entry package.
Historical archive filenames retain their original distance targets; the table
and refinement records give the current results.

## 日本語の案内

論文に掲載した9符号の完全な指数配列、距離の検証プログラム、探索記録を公開しています。
符号の定義とパラメータは上の表と `codes.json` にあります。
上記の検証コマンドを実行してください。
`example34.py` は J=3, L=4 の例の行列ビットマップとTikZグラフを生成します。
行列の強調した72個の1とグラフの36辺の対応も確認します。
探索そのものを再実行する場合は `reproduce.py` を使い、同じコマンドで再開できます。
J=4, L=8 の最小距離は未確定で、現在の範囲は 30 <= d <= 46 です。

The codeword strip and graph now share 24 distinct symbols (four colors and six shapes),
so each nonzero component can be matched with its graph vertex. Vertex coordinates,
edges, codeword support, and matrix highlights are unchanged.
The generated TikZ strip is `figures/example34_codeword.tex`; both panels use
`figures/example34_symbols.tex` for colors. The binary PNG remains a reference.
