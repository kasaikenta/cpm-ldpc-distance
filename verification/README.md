# Reproducing the distance checks

These files support the distance computations in the bilingual paper.
Public repository: https://github.com/kasaikenta/cpm-ldpc-distance
All codes are classical binary codes.
The input magic `CPM_CSS_V1` is a legacy file-format name. Every invocation here
uses `--mode classical --side X` and no stabilizer rows. No quotient by a
stabilizer space or quantum-code equivalence is used.

## Commands (from the repository root or extracted distribution directory)

```
python3 verification/independent24.py
python3 verification/check_records.py
python3 verification/example34.py
python3 verification/enumerate_j4_l5.py
python3 verification/review_checks.py
python3 verification/reproduce.py --J 4 --L 5 --threads 4
```

Python 3.10+ and a C++17 compiler are required. The independent24 command enumerates all 2^26-1 nonzero words of the first table
[96,26,24] example and checks that exactly 668 words have weight 24.
The check_records command reads the
J3/J4 archives in `supporting_material/` without extracting thousands of files.
It checks all nine table entries, all referenced run-record hashes, each root
and partition, the input arrays, ranks, and actual upper-bound codewords. The J=3 table entries have exact distance 24.
The J=4 entries use the fixed (L,P)=(5,23),(6,29),(7,43),(8,73).
Their exact distances are established for L=5 (30), L=6 (28), and L=7 (32).
For L=8, the current bounds are 30 <= d <= 46. These values are recorded in
`joint_examples.json` and `J4_distance_refinement_verification.zip`.
The latter archive contains the complete search records supporting each current
lower bound, explicit upper-bound words, the L=5 enumeration record, and the
search source. All four J=4 arrays have no 4-cycles.
It does **not** rerun the lower-bound search. The record names and archived
input tables provide experiment provenance; the JSON run records do not embed
the full exponent arrays and are not independently checkable nonexistence proofs.
Their historical filenames contain “certificate” but denote solver run records.
For the J=3,L=4 table entry, the record checker verifies the matrix, rank, and
a weight-24 word; use independent24.py to re-establish the lower bound by enumeration.

The example34 command generates the matrix and codeword bitmaps and the TikZ graph
in figures/. Its 24 support positions give 36 graph edges. Each edge corresponds
to two highlighted matrix entries in the same check row, for 72 highlighted 1s.
The three block-row colors agree between the matrix and graph. All other 1s are
light gray. This command also checks the syndrome, all pairings, and connectivity.
The figure generator uses only the Python standard library; typesetting its graph
requires LaTeX with TikZ.

The review_checks command independently enumerates every nonzero kernel word of 54
small codes, including J=2,3,4, prime/composite lift sizes, repeated exponents,
and short cycles. It compares 216 unpartitioned solver calls (with and without
QC symmetry) and 54 partitioned calls with these exact answers. It also checks
12 cofactor words of weights 6,24,120 for the explicit construction, using sparse
supports at both threshold+1 and threshold+2. These finite checks test the
implementation and examples; the general lower bound uses the paper's proof.

The reproduce command reconstructs H from the displayed E and searches through
`d_lower-2` for the selected example, using the even-weight property.
`--J` may be 3 or 4. For J=3, `--L` may be 4,5,6,7,8; for J=4, it may be 5,6,7,8.
For J=3,L=4, it runs the independent enumeration of all 67,108,863 nonzero words,
obtaining minimum weight 24 and exactly 668 words of that weight.
For J=4,L=5, it instead runs
the short independent enumeration of all 67,108,863 nonzero words, obtaining
minimum weight 30 and exactly 23 words of that weight.

For other cases, use, for example:

```
python3 verification/reproduce.py --J 4 --L 6 --threads 4 --seconds 60
```

The time budget is per invocation. Repeating the same command restores the
actual DFS frontiers. Completed partitions are skipped. The output directory
contains per-partition checkpoints and a `summary.json` with the precise
completed and incomplete scope. A global lower bound is reported only when
all roots and partitions are complete. Large examples may need substantial
computation. Do not change `--shards` or the input arrays while resuming.

The serial rerun source has SHA256
`0dd279be2e41cd81ce7d70b767854012bd106116d57fd455dc9262bbb1aa3620`.
The refined archive also includes a parallel variant that uses the same
branching and pruning rules and stores each worker's remaining DFS frontier.
The two variants share a checkpoint representation; the input, cutoff, and
partition settings are fingerprinted. Short timeout-and-resume tests and
small-code exhaustive comparisons are retained with the experiment records.

The parallel checkpoint code was corrected to retain a discovered word in every
periodic snapshot. A forced-exit test confirms that the word survives restart.
The reused historical checkpoints were audited: all had zero kernel events, so
the lost-witness case did not occur in these runs. The record checker also rejects
any complete no-witness record with a nonzero kernel-event count. The refined
archive includes a compact audit summary; the full checkpoint manifests remain
with the experiment records.

For the last 58 original L=7, weight-30 partitions, the remaining assigned DFS
frontiers were subdivided into 928 fragments for load balancing across nodes.
The subdivision preserves every saved frontier node exactly once. The parent's
accumulated counters are carried by one fragment, so they are not double counted.
An original partition is marked complete only after all its fragments finish;
the native search executable then emits its usual completion record from the
verified empty frontier. The additional archive, when present, contains the
parent states, initial fragments, completed fragment states, and merged records.
The checker validates the complete correspondence and ties merged records to the
main certificate list. A small-code test gives the same complete search counters
before and after subdivision.

## Why the classical search is complete

A variable has integer index t*L+l. Under a common cyclic translation, choose a
support point in the first occupied block column l and translate its t to zero.
This places a representative at root l with every earlier block column empty.
The L roots cover all nonzero supports up to this symmetry for every P. Further
representatives can remain within a root; uniqueness is unnecessary.

At a partial support, an odd check needs at least one unused incident variable.
Choose its candidates in a deterministic order. The branch choosing candidate v
forbids the earlier candidates. Every completion belongs to the branch for its
first candidate. At syndrome zero a nonzero word has already been found, so its
supersets need not be searched to decide whether a word exists below the cutoff.
A check with no available variable and odd parity is impossible. With exactly
one available variable, its value is forced by the check parity.

All other classical pruning rules are necessary lower bounds on added weight:

1. If s checks are odd and a column has at most J ones, at least ceil(s/J)
   variables must be added.
2. A full CPM block row contains one neighbor of each variable, so its number
   of odd checks is also a lower bound.
3. Odd checks with pairwise disjoint available-variable sets require distinct
   added variables. A greedy packing gives a valid lower bound.
4. The implementation also builds an integer-scaled feasible dual of set cover:
   assign nonnegative weights to odd checks, ensuring that each available
   variable meets checks of total weight at most J. The total assigned weight
   divided by J, rounded up, lower-bounds the required variables. This follows
   by summing these inequalities over any completing set of variables.
5. For two CPM block rows with a and b odd checks, let s be the maximum matching
   size in the bipartite graph connecting these odd checks by available variables.
   Any completion needs at least a+b-s variables. To see this, take a minimum
   subset of completion variables covering all a+b marked vertices. Edges
   covering only one marked vertex count separately. The edges between marked
   vertices form disjoint stars in a minimal edge cover; choosing one edge from
   each nontrivial star yields a matching of size at least a+b-q if q edges were
   used. Thus s>=a+b-q. The implementation checks b-s<=q-a after a>=b.
6. All codewords have even weight, because summing the P checks in any block
   row gives the all-one parity equation. Round a total-weight lower bound up
   to an even integer.

The additional stabilizer-based functions in the shared source are inactive
when the stabilizer row list is empty. The test and record check assert this.

The partitioned search numbers frontier states in deterministic order and assigns
state h to partition h modulo the partition count. Branches closed before that
frontier are handled during frontier generation. Checking every partition for
all roots covers the full search. Completed record sets are validated without
assuming that every root uses the same partition count.

## Bounds of this verification

Hashes detect changed files; they do not establish mathematical correctness.
Finite oracle tests do not prove a program correct on all inputs. Confidence in
the reported distances rests on the algorithm's completeness argument, the
archived complete runs, the independently checked upper-bound words, and the
available rerun source. No proof-assistant or SAT nonexistence certificate is
claimed. Source and data are available in the public repository linked above.
The v1.1.4 tag identifies the published verification package for this table.

## Current and historical archives

The baseline records are `J3_L5to8_distance24_certificates.zip` and
`J4_L5to8_distance26_verification.zip` in `supporting_material/`.
The additional J=4 bounds and exact distances are verified from
`J4_distance_refinement_verification.zip` in the same directory.
Older J=4 distance-24 archives are not included in this repository and are
not the source of the current table. Use the commands above for this version.

The figure generator uses 24 distinct color/shape pairs for support positions.
The TikZ codeword strip and graph use identical symbols; vertex coordinates,
edges, and the highlighted matrix entries are unchanged.
`figures/example34_symbols.tex` defines their four colors.
The paper uses `figures/example34_codeword.tex`; the original binary PNG is retained as a reference.

In v1.1.4 the TikZ codeword strip spans 16 cm, with 1.5 mm symbols.
The graph coordinates use a scale of 0.82 cm per unit, retaining 2 mm vertex symbols and 7 pt labels.
