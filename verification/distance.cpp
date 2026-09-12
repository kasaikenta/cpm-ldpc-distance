#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <csignal>
#include <cstdio>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mutex>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

constexpr int MAX_N = 8192;
constexpr int MAX_M = 4096;
constexpr int MAX_N_WORDS = (MAX_N + 63) / 64;
constexpr int MAX_M_WORDS = (MAX_M + 63) / 64;

template <std::size_t W>
struct Bits {
    std::array<std::uint64_t, W> a{};

    bool test(int i) const { return (a[i >> 6] >> (i & 63)) & 1ULL; }
    void set(int i) { a[i >> 6] |= 1ULL << (i & 63); }
    void clear(int i) { a[i >> 6] &= ~(1ULL << (i & 63)); }
    void xor_with(const Bits &b, int words) {
        for (int i = 0; i < words; ++i) a[i] ^= b.a[i];
    }
    void or_with(const Bits &b, int words) {
        for (int i = 0; i < words; ++i) a[i] |= b.a[i];
    }
    bool intersects(const Bits &b, int words) const {
        for (int i = 0; i < words; ++i) if (a[i] & b.a[i]) return true;
        return false;
    }
    bool zero(int words) const {
        for (int i = 0; i < words; ++i) if (a[i]) return false;
        return true;
    }
    int popcount(int words) const {
        int s = 0;
        for (int i = 0; i < words; ++i) s += __builtin_popcountll(a[i]);
        return s;
    }
    int first(int words) const {
        for (int i = 0; i < words; ++i) {
            if (a[i]) return 64 * i + __builtin_ctzll(a[i]);
        }
        return -1;
    }
    int highest(int words) const {
        for (int i = words - 1; i >= 0; --i) {
            if (a[i]) return 64 * i + 63 - __builtin_clzll(a[i]);
        }
        return -1;
    }
};

using Support = Bits<MAX_N_WORDS>;
using Syndrome = Bits<MAX_M_WORDS>;

struct CPMInstance {
    std::string name;
    bool sparse_binary = false;
    int n = 0;
    int p = 0;
    int l = 0;
    int jx = 0;
    int jz = 0;
    std::vector<std::vector<int>> ex;
    std::vector<std::vector<int>> ez;
    std::vector<std::vector<int>> hx_rows;
    std::vector<std::vector<int>> hz_rows;
};

struct Options {
    std::string input;
    std::string output;
    std::string mode = "css";
    std::string side = "X";
    int max_weight = -1;
    int root_start = 0;
    int root_end = -1;
    int threads = 1;
    std::string checkpoint;
    std::string progress;
    int checkpoint_interval_sec = 60;
    int soft_deadline_sec = 0;
    int h_rt_sec = 0;
    int split_weight = 5;
    int shard_count = 1;
    int shard_index = 0;
    int target_tasks = 500000;
    int check_group_size = 0;
    bool partition_only = false;
    bool qc_reduce = true;
    bool stop_on_witness = true;
    bool quiet = false;
};

struct SearchStats {
    std::uint64_t states = 0;
    std::uint64_t branches = 0;
    std::uint64_t forced = 0;
    std::uint64_t weight_prunes = 0;
    std::uint64_t lower_bound_prunes = 0;
    std::uint64_t infeasible_prunes = 0;
    std::uint64_t kernel_circuits = 0;
    std::uint64_t stabilizer_circuits = 0;

    SearchStats &operator+=(const SearchStats &b) {
        states += b.states; branches += b.branches; forced += b.forced;
        weight_prunes += b.weight_prunes;
        lower_bound_prunes += b.lower_bound_prunes;
        infeasible_prunes += b.infeasible_prunes;
        kernel_circuits += b.kernel_circuits;
        stabilizer_circuits += b.stabilizer_circuits;
        return *this;
    }
};

struct RootResult {
    int root = -1;
    bool complete = false;
    bool found = false;
    Support witness;
    int witness_weight = 0;
    SearchStats stats;
    double elapsed = 0.0;
};

struct DFSNode {
    Support support;
    Support forbidden;
    Syndrome syndrome;
    int weight = 0;
    bool partitioned = false;
};

struct CheckpointState {
    std::uint64_t fingerprint = 0;
    int root = -1;
    int max_weight = 0;
    int split_weight = 0;
    int shard_count = 1;
    int shard_index = 0;
    std::uint64_t frontier_counter = 0;
    std::uint64_t resume_count = 0;
    std::int64_t started_at_epoch = 0;
    RootResult result;
    std::vector<DFSNode> stack;
};

volatile std::sig_atomic_t signal_stop_requested = 0;

void request_signal_stop(int) {
    signal_stop_requested = 1;
}

CPMInstance read_instance(const std::string &path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open input: " + path);
    std::string magic;
    f >> magic;
    CPMInstance x;
    if (magic == "SPARSE_CSS_V1") {
        int mx = 0, mz = 0;
        f >> x.name >> x.n >> mx >> mz;
        if (!f || x.n <= 0 || mx < 0 || mz < 0)
            throw std::runtime_error("invalid sparse CSS header");
        if (x.n > MAX_N || std::max(mx, mz) > MAX_M)
            throw std::runtime_error("sparse instance exceeds compiled MAX_N or MAX_M");
        x.sparse_binary = true;
        x.p = 1;
        x.l = x.n;
        x.jx = mx;
        x.jz = mz;
        auto read_rows = [&](int count, std::vector<std::vector<int>> &rows) {
            rows.assign(count, {});
            for (auto &row : rows) {
                int weight = -1;
                f >> weight;
                if (!f || weight < 0) throw std::runtime_error("invalid sparse row weight");
                row.resize(weight);
                for (int &column : row) f >> column;
                if (!f) throw std::runtime_error("truncated sparse row");
                if (!std::is_sorted(row.begin(), row.end()) ||
                    std::adjacent_find(row.begin(), row.end()) != row.end())
                    throw std::runtime_error("sparse row must be sorted and duplicate-free");
                for (int column : row)
                    if (column < 0 || column >= x.n)
                        throw std::runtime_error("sparse column out of range");
            }
        };
        read_rows(mx, x.hx_rows);
        read_rows(mz, x.hz_rows);
        return x;
    }
    if (magic != "CPM_CSS_V1")
        throw std::runtime_error("expected CPM_CSS_V1 or SPARSE_CSS_V1");
    f >> x.name >> x.p >> x.l >> x.jx >> x.jz;
    if (!f || x.p <= 0 || x.l <= 0 || x.jx < 0 || x.jz < 0)
        throw std::runtime_error("invalid header");
    x.ex.assign(x.jx, std::vector<int>(x.l));
    x.ez.assign(x.jz, std::vector<int>(x.l));
    for (auto &row : x.ex) for (int &v : row) f >> v;
    for (auto &row : x.ez) for (int &v : row) f >> v;
    if (!f) throw std::runtime_error("truncated exponent arrays");
    for (const auto &row : x.ex) for (int v : row)
        if (v < 0 || v >= x.p) throw std::runtime_error("X exponent out of range");
    for (const auto &row : x.ez) for (int v : row)
        if (v < 0 || v >= x.p) throw std::runtime_error("Z exponent out of range");
    if (x.p * x.l > MAX_N || std::max(x.jx, x.jz) * x.p > MAX_M)
        throw std::runtime_error("instance exceeds compiled MAX_N or MAX_M");
    x.n = x.p * x.l;
    return x;
}

std::vector<std::vector<int>> lifted_rows(
    const std::vector<std::vector<int>> &e, int p, int l
) {
    std::vector<std::vector<int>> rows;
    rows.reserve(e.size() * p);
    for (const auto &block_row : e) {
        for (int a = 0; a < p; ++a) {
            std::vector<int> row;
            row.reserve(l);
            for (int t = 0; t < l; ++t) {
                int g = (a - block_row[t]) % p;
                if (g < 0) g += p;
                row.push_back(l * g + t);
            }
            rows.push_back(std::move(row));
        }
    }
    return rows;
}

bool kernel_has_even_weight(const std::vector<std::vector<int>> &checks, int n) {
    // Every vector in ker(H) has even weight iff the all-one vector belongs
    // to row(H).  Testing odd column degree is only a sufficient condition;
    // it misses the CPM case with an even number of block rows, where the P
    // rows of any one block-row family already sum to the all-one vector.
    int words = (n + 63) / 64;
    std::vector<Support> basis(n);
    std::vector<unsigned char> has_pivot(n, 0);
    for (const auto &row : checks) {
        Support x;
        for (int column : row) x.set(column);
        while (!x.zero(words)) {
            int pivot = x.highest(words);
            if (has_pivot[pivot]) {
                x.xor_with(basis[pivot], words);
            } else {
                basis[pivot] = x;
                has_pivot[pivot] = 1;
                break;
            }
        }
    }

    Support all_one;
    for (int column = 0; column < n; ++column) all_one.set(column);
    while (!all_one.zero(words)) {
        int pivot = all_one.highest(words);
        if (!has_pivot[pivot]) return false;
        all_one.xor_with(basis[pivot], words);
    }
    return true;
}

class SearchEngine {
public:
    SearchEngine(
        int n_, const std::vector<std::vector<int>> &checks_,
        const std::vector<std::vector<int>> &stabilizers_, int max_weight_,
        std::atomic<bool> &global_stop_, int qc_base_columns_ = 0,
        bool kernel_even_weight_ = false, int check_group_size_ = 0
    ) : n(n_), max_weight(max_weight_), checks(checks_),
        stabilizers(stabilizers_), global_stop(global_stop_),
        qc_base_columns(qc_base_columns_),
        kernel_even_weight(kernel_even_weight_),
        check_group_size(check_group_size_) {
        m = static_cast<int>(checks.size());
        n_words = (n + 63) / 64;
        m_words = (m + 63) / 64;
        col_masks.assign(n, Syndrome{});
        basis.assign(n, Support{});
        has_pivot.assign(n, 0);
        stabilizer_rows_by_variable.assign(n, {});
        stabilizer_min_variable.assign(stabilizers.size(), -1);
        dual_load.assign(n, 0);
        dual_stamp.assign(n, 0);
        max_col_weight = 0;
        if (check_group_size > 0 && m % check_group_size != 0)
            throw std::runtime_error(
                "check row count is not divisible by --check-group-size");
        if (check_group_size > 0)
            variable_check_in_group.assign(
                n, std::vector<int>(m / check_group_size, -1));
        std::vector<int> degree(n, 0);
        for (int r = 0; r < m; ++r) {
            for (int v : checks[r]) {
                if (v < 0 || v >= n) throw std::runtime_error("column out of range");
                col_masks[v].set(r);
                ++degree[v];
                if (!variable_check_in_group.empty()) {
                    const int group = r / check_group_size;
                    if (variable_check_in_group[v][group] >= 0)
                        throw std::runtime_error(
                            "CPM variable has two checks in one row group");
                    variable_check_in_group[v][group] = r;
                }
            }
        }
        if (!variable_check_in_group.empty()) {
            for (const auto &by_group : variable_check_in_group)
                if (std::find(by_group.begin(), by_group.end(), -1) !=
                    by_group.end())
                    throw std::runtime_error(
                        "CPM variable misses a check row group");
        }
        for (int d : degree) max_col_weight = std::max(max_col_weight, d);
        for (int row = 0; row < static_cast<int>(stabilizers.size()); ++row) {
            if (!stabilizers[row].empty())
                stabilizer_min_variable[row] = *std::min_element(
                    stabilizers[row].begin(), stabilizers[row].end());
            for (int variable : stabilizers[row])
                stabilizer_rows_by_variable[variable].push_back(row);
        }
        build_stabilizer_basis();
    }

    RootResult run_root(int root) {
        RootResult result;
        result.root = root;
        current = &result;
        auto t0 = std::chrono::steady_clock::now();
        Support support, forbidden = initial_forbidden(root);
        support.set(root);
        Syndrome syndrome = col_masks[root];
        dfs(root, support, forbidden, syndrome, 1);
        result.complete = !global_stop.load(std::memory_order_relaxed) || result.found;
        result.elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
        current = nullptr;
        return result;
    }

    RootResult run_subtree(int root, const DFSNode &node) {
        RootResult result;
        result.root = root;
        current = &result;
        auto t0 = std::chrono::steady_clock::now();
        dfs(root, node.support, node.forbidden, node.syndrome, node.weight);
        result.complete =
            !global_stop.load(std::memory_order_relaxed) || result.found;
        result.elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
        current = nullptr;
        return result;
    }

    void initialize_checkpoint(CheckpointState &state, int root) const {
        state.root = root;
        state.result.root = root;
        DFSNode node;
        node.support.set(root);
        node.forbidden = initial_forbidden(root);
        node.syndrome = col_masks[root];
        node.weight = 1;
        node.partitioned = state.shard_count == 1;
        state.stack.push_back(node);
    }

    // Process exactly one explicit DFS node.  Keeping the frontier explicit
    // makes the search state small, deterministic, and exactly restartable.
    void checkpoint_step(
        CheckpointState &state, bool partition_only = false,
        std::vector<DFSNode> *assigned = nullptr
    ) {
        if (state.stack.empty() || state.result.found) return;
        DFSNode node = state.stack.back();
        state.stack.pop_back();

        if (partition_only && node.partitioned) {
            assigned->push_back(std::move(node));
            return;
        }

        if (!node.partitioned && state.shard_count > 1) {
            if (node.weight < state.split_weight) {
                // Still above the deterministic partition frontier.
            } else {
                std::uint64_t frontier_id = state.frontier_counter++;
                if (static_cast<int>(frontier_id % state.shard_count) !=
                    state.shard_index) return;
                node.partitioned = true;
                if (partition_only) {
                    assigned->push_back(std::move(node));
                    return;
                }
            }
        }

        ++state.result.stats.states;
        if (!unit_propagate(
                state.root, node.support, node.forbidden, node.syndrome,
                node.weight, state.result.stats)) {
            ++state.result.stats.infeasible_prunes;
            return;
        }
        if (!node.partitioned && state.shard_count > 1 &&
            node.weight >= state.split_weight) {
            std::uint64_t frontier_id = state.frontier_counter++;
            if (static_cast<int>(frontier_id % state.shard_count) !=
                state.shard_index) return;
            node.partitioned = true;
            if (partition_only) {
                assigned->push_back(std::move(node));
                return;
            }
        }
        if (node.syndrome.zero(m_words)) {
            ++state.result.stats.kernel_circuits;
            if (stabilizers.empty() || !in_stabilizer_span(node.support)) {
                state.result.found = true;
                state.result.witness = node.support;
                state.result.witness_weight = node.weight;
                state.stack.clear();
            } else {
                ++state.result.stats.stabilizer_circuits;
            }
            return;
        }
        if (node.weight >= max_weight) {
            ++state.result.stats.weight_prunes;
            return;
        }

        BoundInfo info = analyze(
            state.root, node.support, node.forbidden, node.syndrome);
        if (!info.feasible) {
            ++state.result.stats.infeasible_prunes;
            return;
        }
        int completion_lower = info.lower;
        if (kernel_even_weight && ((node.weight + completion_lower) & 1))
            ++completion_lower;
        if (node.weight + completion_lower > max_weight) {
            ++state.result.stats.lower_bound_prunes;
            return;
        }
        if (!tight_pair_matching_possible(
                state.root, node.support, node.forbidden, node.syndrome,
                max_weight - node.weight)) {
            ++state.result.stats.lower_bound_prunes;
            return;
        }

        std::vector<DFSNode> children;
        children.reserve(info.branch_vars.size());
        Support local_forbidden = node.forbidden;
        for (int v : info.branch_vars) {
            ++state.result.stats.branches;
            DFSNode child;
            child.support = node.support;
            child.support.set(v);
            child.forbidden = local_forbidden;
            child.syndrome = node.syndrome;
            child.syndrome.xor_with(col_masks[v], m_words);
            child.weight = node.weight + 1;
            child.partitioned = node.partitioned;
            if (reducible_by_incident_stabilizer(child.support, v)) {
                ++state.result.stats.infeasible_prunes;
            } else {
                children.push_back(child);
            }
            local_forbidden.set(v);
        }
        // LIFO reverse push preserves the original recursive DFS order.
        for (auto it = children.rbegin(); it != children.rend(); ++it)
            state.stack.push_back(*it);
    }

private:
    int n = 0, m = 0, n_words = 0, m_words = 0;
    int max_weight = 0, max_col_weight = 0;
    const std::vector<std::vector<int>> &checks;
    const std::vector<std::vector<int>> &stabilizers;
    std::vector<Syndrome> col_masks;
    // A full support vector is comparatively large.  Keep the elimination
    // basis on the heap: macOS worker threads have a small default stack.
    std::vector<Support> basis;
    std::vector<unsigned char> has_pivot;
    std::vector<std::vector<int>> stabilizer_rows_by_variable;
    std::vector<int> stabilizer_min_variable;
    std::vector<unsigned char> dual_load;
    std::vector<std::uint64_t> dual_stamp;
    std::vector<std::vector<int>> variable_check_in_group;
    std::uint64_t dual_epoch = 0;
    std::atomic<bool> &global_stop;
    int qc_base_columns = 0;
    bool kernel_even_weight = false;
    // For a full CPM lift, check rows are consecutive P-row groups and every
    // variable touches exactly one row in every group.  A completion must
    // therefore contain at least as many variables as the largest group
    // syndrome weight.  Zero disables this CPM-specific exact lower bound.
    int check_group_size = 0;
    RootResult *current = nullptr;

    Support initial_forbidden(int root) const {
        Support forbidden;
        if (qc_base_columns <= 0) return forbidden;
        // Root t represents supports whose earliest occupied base-column
        // type is t.  A common cyclic shift moves one coordinate of that
        // type to residue zero (residue-major index t).  Forbid all lift
        // coordinates of earlier types so the L roots are disjoint.
        for (int variable = 0; variable < n; ++variable) {
            if (variable % qc_base_columns < root) forbidden.set(variable);
        }
        return forbidden;
    }

    void build_stabilizer_basis() {
        for (const auto &row : stabilizers) {
            Support x;
            for (int v : row) x.set(v);
            while (!x.zero(n_words)) {
                int p = x.highest(n_words);
                if (has_pivot[p]) x.xor_with(basis[p], n_words);
                else { basis[p] = x; has_pivot[p] = 1; break; }
            }
        }
    }

    bool in_stabilizer_span(Support x) const {
        while (!x.zero(n_words)) {
            int p = x.highest(n_words);
            if (!has_pivot[p]) return false;
            x.xor_with(basis[p], n_words);
        }
        return true;
    }

    bool reducible_by_incident_stabilizer(
        const Support &support, int last_variable
    ) const {
        for (int row : stabilizer_rows_by_variable[last_variable]) {
            int overlap = 0;
            for (int variable : stabilizers[row])
                overlap += support.test(variable) ? 1 : 0;
            const int row_weight = static_cast<int>(stabilizers[row].size());
            if (2 * overlap > row_weight ||
                (2 * overlap == row_weight &&
                 !support.test(stabilizer_min_variable[row])))
                return true;
        }
        return false;
    }

public:
    void filter_reducible_frontier(CheckpointState &state) const {
        std::vector<DFSNode> kept;
        kept.reserve(state.stack.size());
        std::vector<unsigned char> overlap(stabilizers.size(), 0);
        std::vector<int> touched;
        for (const DFSNode &node : state.stack) {
            bool reducible = false;
            touched.clear();
            touched.reserve(node.weight * 4);
            for (int word = 0; word < n_words; ++word) {
                std::uint64_t bits = node.support.a[word];
                while (bits) {
                    int offset = __builtin_ctzll(bits);
                    int variable = 64 * word + offset;
                    for (int row : stabilizer_rows_by_variable[variable]) {
                        if (overlap[row]++ == 0) touched.push_back(row);
                    }
                    bits &= bits - 1;
                }
            }
            for (int row : touched) {
                const int row_weight =
                    static_cast<int>(stabilizers[row].size());
                if (2 * overlap[row] > row_weight ||
                    (2 * overlap[row] == row_weight &&
                     !node.support.test(stabilizer_min_variable[row]))) {
                    reducible = true;
                    break;
                }
            }
            for (int row : touched) overlap[row] = 0;
            if (reducible) ++state.result.stats.infeasible_prunes;
            else kept.push_back(node);
        }
        state.stack.swap(kept);
    }

private:

    std::vector<int> available_on_check(
        int check, int root, const Support &support, const Support &forbidden
    ) const {
        std::vector<int> out;
        for (int v : checks[check]) {
            if (v > root && !support.test(v) && !forbidden.test(v)) out.push_back(v);
        }
        return out;
    }

    // For two row groups with a>=b odd checks, q future variables give q-a
    // degree slots beyond the mandatory one at each of the a larger-side odd
    // checks.  If a maximum matching between the two odd sets has size s,
    // at least b-s smaller-side odd checks must use such extra slots.  Thus
    // b-s <= q-a is a necessary completion condition.
    bool tight_pair_matching_possible(
        int root, const Support &support, const Support &forbidden,
        const Syndrome &syndrome, int remaining_budget
    ) const {
        if (variable_check_in_group.empty()) return true;
        const int groups = m / check_group_size;
        std::vector<std::vector<int>> odd(groups);
        for (int group = 0; group < groups; ++group) {
            const int first = group * check_group_size;
            for (int check = first; check < first + check_group_size; ++check)
                if (syndrome.test(check)) odd[group].push_back(check);
        }

        for (int ga = 0; ga < groups; ++ga) {
            for (int gb = ga + 1; gb < groups; ++gb) {
                int small = ga, large = gb;
                if (odd[small].size() > odd[large].size())
                    std::swap(small, large);
                if (odd[small].empty()) continue;
                const int allowed_deficit = remaining_budget -
                    static_cast<int>(odd[large].size());
                if (allowed_deficit < 0) return false;
                if (static_cast<int>(odd[small].size()) <= allowed_deficit)
                    continue;

                std::vector<int> large_index(check_group_size, -1);
                const int large_first = large * check_group_size;
                for (int index = 0;
                     index < static_cast<int>(odd[large].size()); ++index)
                    large_index[odd[large][index] - large_first] = index;
                std::vector<int> matched_by(odd[large].size(), -1);

                auto augment = [&](auto &&self, int small_index,
                                   std::vector<unsigned char> &seen) -> bool {
                    const int check = odd[small][small_index];
                    for (int variable : checks[check]) {
                        if (variable <= root || support.test(variable) ||
                            forbidden.test(variable))
                            continue;
                        const int target =
                            variable_check_in_group[variable][large];
                        const int target_index =
                            large_index[target - large_first];
                        if (target_index < 0 || seen[target_index]) continue;
                        seen[target_index] = 1;
                        if (matched_by[target_index] < 0 ||
                            self(self, matched_by[target_index], seen)) {
                            matched_by[target_index] = small_index;
                            return true;
                        }
                    }
                    return false;
                };

                int matched = 0;
                for (int index = 0;
                     index < static_cast<int>(odd[small].size()); ++index) {
                    std::vector<unsigned char> seen(odd[large].size(), 0);
                    if (augment(augment, index, seen)) ++matched;
                }
                if (static_cast<int>(odd[small].size()) - matched >
                    allowed_deficit)
                    return false;
            }
        }
        return true;
    }

    bool unit_propagate(
        int root, Support &support, Support &forbidden,
        Syndrome &syndrome, int &weight, SearchStats &stats
    ) const {
        while (true) {
            bool changed = false;
            for (int check = 0; check < m; ++check) {
                int only = -1;
                int count = 0;
                for (int variable : checks[check]) {
                    if (variable > root && !support.test(variable) &&
                        !forbidden.test(variable)) {
                        only = variable;
                        if (++count == 2) break;
                    }
                }
                const bool odd = syndrome.test(check);
                if (count == 0) {
                    if (odd) return false;
                    continue;
                }
                if (count != 1) continue;
                ++stats.forced;
                if (odd) {
                    support.set(only);
                    syndrome.xor_with(col_masks[only], m_words);
                    ++weight;
                    if (weight > max_weight ||
                        reducible_by_incident_stabilizer(support, only))
                        return false;
                } else {
                    forbidden.set(only);
                }
                changed = true;
                break;
            }
            if (!changed) return true;
        }
    }

    struct BoundInfo {
        bool feasible = true;
        int lower = 0;
        int branch_check = -1;
        std::vector<int> branch_vars;
    };

    BoundInfo analyze(
        int root, const Support &support, const Support &forbidden,
        const Syndrome &syndrome
    ) {
        BoundInfo info;
        int unsat = syndrome.popcount(m_words);
        if (unsat == 0) return info;
        info.lower = (unsat + std::max(1, max_col_weight) - 1) /
                     std::max(1, max_col_weight);
        if (check_group_size > 0 && m % check_group_size == 0) {
            for (int first = 0; first < m; first += check_group_size) {
                int group_unsat = 0;
                const int last = first + check_group_size;
                for (int check = first; check < last; ++check)
                    group_unsat += syndrome.test(check) ? 1 : 0;
                info.lower = std::max(info.lower, group_unsat);
            }
        }

        struct CandidateSet { int check; std::vector<int> vars; };
        std::vector<CandidateSet> sets;
        sets.reserve(unsat);
        for (int wi = 0; wi < m_words; ++wi) {
            std::uint64_t z = syndrome.a[wi];
            while (z) {
                int b = __builtin_ctzll(z);
                int c = 64 * wi + b;
                z &= z - 1;
                if (c >= m) continue;
                auto vars = available_on_check(c, root, support, forbidden);
                if (vars.empty()) { info.feasible = false; return info; }
                if (info.branch_check < 0 || vars.size() < info.branch_vars.size()) {
                    info.branch_check = c;
                    info.branch_vars = vars;
                }
                sets.push_back({c, std::move(vars)});
            }
        }

        std::sort(sets.begin(), sets.end(), [](const auto &a, const auto &b) {
            return a.vars.size() < b.vars.size();
        });

        // Feasible integer-scaled dual of the set-cover relaxation.  Every
        // unsatisfied check must be hit by a future variable.  A variable has
        // capacity max_col_weight; the accumulated dual value divided by that
        // capacity lower-bounds the number of future variables.
        if (++dual_epoch == 0) {
            std::fill(dual_stamp.begin(), dual_stamp.end(), 0);
            ++dual_epoch;
        }
        int dual_units = 0;
        const int dual_scale = std::max(1, max_col_weight);
        for (const auto &s : sets) {
            int increment = dual_scale;
            for (int v : s.vars) {
                int load = dual_stamp[v] == dual_epoch ? dual_load[v] : 0;
                increment = std::min(increment, dual_scale - load);
            }
            if (increment <= 0) continue;
            dual_units += increment;
            for (int v : s.vars) {
                if (dual_stamp[v] != dual_epoch) {
                    dual_stamp[v] = dual_epoch;
                    dual_load[v] = 0;
                }
                dual_load[v] += static_cast<unsigned char>(increment);
            }
        }
        info.lower = std::max(
            info.lower, (dual_units + dual_scale - 1) / dual_scale);

        Support used;
        int packing = 0;
        for (const auto &s : sets) {
            bool disjoint = true;
            for (int v : s.vars) if (used.test(v)) { disjoint = false; break; }
            if (disjoint) {
                ++packing;
                for (int v : s.vars) used.set(v);
            }
        }
        info.lower = std::max(info.lower, packing);

        std::sort(info.branch_vars.begin(), info.branch_vars.end(), [&](int a, int b) {
            int sa = 0, sb = 0;
            for (int wi = 0; wi < m_words; ++wi) {
                sa += __builtin_popcountll(col_masks[a].a[wi] & syndrome.a[wi]);
                sb += __builtin_popcountll(col_masks[b].a[wi] & syndrome.a[wi]);
            }
            return sa != sb ? sa > sb : a < b;
        });
        return info;
    }

    void dfs(
        int root, Support support, Support forbidden,
        Syndrome syndrome, int weight
    ) {
        if (global_stop.load(std::memory_order_relaxed)) return;
        ++current->stats.states;

        if (!unit_propagate(
                root, support, forbidden, syndrome, weight,
                current->stats)) {
            ++current->stats.infeasible_prunes;
            return;
        }

        if (syndrome.zero(m_words)) {
            ++current->stats.kernel_circuits;
            if (stabilizers.empty() || !in_stabilizer_span(support)) {
                current->found = true;
                current->witness = support;
                current->witness_weight = weight;
                global_stop.store(true, std::memory_order_relaxed);
            } else {
                ++current->stats.stabilizer_circuits;
            }
            // A minimum-weight non-stabilizer kernel vector contains no proper
            // nonzero kernel subset, so closed circuits need not be extended.
            return;
        }
        if (weight >= max_weight) { ++current->stats.weight_prunes; return; }

        BoundInfo info = analyze(root, support, forbidden, syndrome);
        if (!info.feasible) { ++current->stats.infeasible_prunes; return; }
        int completion_lower = info.lower;
        if (kernel_even_weight && ((weight + completion_lower) & 1))
            ++completion_lower;
        if (weight + completion_lower > max_weight) {
            ++current->stats.lower_bound_prunes;
            return;
        }
        if (!tight_pair_matching_possible(
                root, support, forbidden, syndrome, max_weight - weight)) {
            ++current->stats.lower_bound_prunes;
            return;
        }

        Support local_forbidden = forbidden;
        for (int v : info.branch_vars) {
            if (global_stop.load(std::memory_order_relaxed)) return;
            ++current->stats.branches;
            Support child_support = support;
            child_support.set(v);
            Syndrome child_syndrome = syndrome;
            child_syndrome.xor_with(col_masks[v], m_words);
            if (reducible_by_incident_stabilizer(child_support, v)) {
                ++current->stats.infeasible_prunes;
            } else {
                dfs(root, child_support, local_forbidden, child_syndrome, weight + 1);
            }
            local_forbidden.set(v);
        }
    }
};

std::string json_escape(const std::string &s) {
    std::ostringstream o;
    for (char c : s) {
        if (c == '\\' || c == '"') o << '\\' << c;
        else if (c == '\n') o << "\\n";
        else o << c;
    }
    return o.str();
}

std::vector<int> support_list(const Support &s, int n) {
    std::vector<int> out;
    for (int i = 0; i < n; ++i) if (s.test(i)) out.push_back(i);
    return out;
}

void print_int_array(std::ostream &o, const std::vector<int> &v) {
    o << '[';
    for (std::size_t i = 0; i < v.size(); ++i) { if (i) o << ','; o << v[i]; }
    o << ']';
}

std::string utc_now() {
    std::time_t now = std::time(nullptr);
    std::tm tm{};
    gmtime_r(&now, &tm);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return buf;
}

std::string utc_from_epoch(std::int64_t epoch) {
    std::time_t when = static_cast<std::time_t>(epoch);
    std::tm tm{};
    gmtime_r(&when, &tm);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return buf;
}

std::string temporary_path(const std::string &path) {
    auto tick = std::chrono::high_resolution_clock::now()
                    .time_since_epoch().count();
    return path + ".tmp." + std::to_string(tick);
}

void atomic_write_text(const std::string &path, const std::string &text) {
    if (path.empty()) return;
    std::string tmp = temporary_path(path);
    {
        std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
        if (!out) throw std::runtime_error("cannot write temporary file: " + tmp);
        out << text;
        out.flush();
        if (!out) throw std::runtime_error("failed writing temporary file: " + tmp);
    }
    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        std::remove(tmp.c_str());
        throw std::runtime_error("cannot atomically rename output: " + path);
    }
}

template <typename T>
void write_pod(std::ostream &out, const T &value) {
    out.write(reinterpret_cast<const char *>(&value), sizeof(value));
}

template <typename T>
void read_pod(std::istream &in, T &value) {
    in.read(reinterpret_cast<char *>(&value), sizeof(value));
    if (!in) throw std::runtime_error("truncated checkpoint");
}

void write_words(std::ostream &out, const std::uint64_t *words, int count) {
    out.write(reinterpret_cast<const char *>(words), count * sizeof(std::uint64_t));
}

void read_words(std::istream &in, std::uint64_t *words, int count) {
    in.read(reinterpret_cast<char *>(words), count * sizeof(std::uint64_t));
    if (!in) throw std::runtime_error("truncated checkpoint bitset");
}

std::uint64_t run_fingerprint(
    const std::string &input, const Options &opt, int root, int effective_max_weight
) {
    std::ifstream in(input, std::ios::binary);
    if (!in) throw std::runtime_error("cannot fingerprint input: " + input);
    std::uint64_t h = 1469598103934665603ULL;
    char buf[8192];
    while (in) {
        in.read(buf, sizeof(buf));
        for (std::streamsize i = 0; i < in.gcount(); ++i) {
            h ^= static_cast<unsigned char>(buf[i]);
            h *= 1099511628211ULL;
        }
    }
    std::ostringstream key;
    key << "qc-base-type-root-v2" << '|'
        << opt.mode << '|' << opt.side << '|' << root << '|'
        << effective_max_weight << '|' << opt.split_weight << '|'
        << opt.shard_index << '|' << opt.shard_count << '|'
        << (opt.qc_reduce ? 1 : 0);
    for (unsigned char c : key.str()) {
        h ^= c;
        h *= 1099511628211ULL;
    }
    // Preserve fingerprints of pre-option checkpoints.  An explicitly
    // enabled sparse grouping must nevertheless be part of the resume key.
    if (opt.check_group_size > 0) {
        const std::string grouping =
            "|check_group_size=" + std::to_string(opt.check_group_size);
        for (unsigned char c : grouping) {
            h ^= c;
            h *= 1099511628211ULL;
        }
    }
    return h;
}

void write_checkpoint_atomic(
    const std::string &path, const CheckpointState &state,
    int n_words, int m_words
) {
    if (path.empty()) return;
    std::string tmp = temporary_path(path);
    {
        std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
        if (!out) throw std::runtime_error("cannot write checkpoint: " + tmp);
        const char magic[8] = {'C','P','M','C','H','K','2','\0'};
        out.write(magic, sizeof(magic));
        std::uint32_t version = 3;
        write_pod(out, version);
        write_pod(out, state.fingerprint);
        write_pod(out, state.root);
        write_pod(out, state.max_weight);
        write_pod(out, state.split_weight);
        write_pod(out, state.shard_count);
        write_pod(out, state.shard_index);
        write_pod(out, n_words);
        write_pod(out, m_words);
        write_pod(out, state.frontier_counter);
        write_pod(out, state.resume_count);
        write_pod(out, state.started_at_epoch);
        write_pod(out, state.result.elapsed);
        const SearchStats &s = state.result.stats;
        write_pod(out, s.states); write_pod(out, s.branches);
        write_pod(out, s.forced); write_pod(out, s.weight_prunes);
        write_pod(out, s.lower_bound_prunes); write_pod(out, s.infeasible_prunes);
        write_pod(out, s.kernel_circuits); write_pod(out, s.stabilizer_circuits);
        std::uint8_t found = state.result.found ? 1 : 0;
        write_pod(out, found);
        write_pod(out, state.result.witness_weight);
        write_words(out, state.result.witness.a.data(), n_words);
        std::uint64_t stack_size = state.stack.size();
        write_pod(out, stack_size);
        for (const DFSNode &node : state.stack) {
            write_pod(out, node.weight);
            std::uint8_t partitioned = node.partitioned ? 1 : 0;
            write_pod(out, partitioned);
            write_words(out, node.support.a.data(), n_words);
            write_words(out, node.forbidden.a.data(), n_words);
            write_words(out, node.syndrome.a.data(), m_words);
        }
        out.flush();
        if (!out) throw std::runtime_error("failed writing checkpoint: " + tmp);
    }
    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        std::remove(tmp.c_str());
        throw std::runtime_error("cannot atomically rename checkpoint: " + path);
    }
}

bool load_checkpoint(
    const std::string &path, CheckpointState &state,
    int expected_n_words, int expected_m_words
) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    char magic[8];
    in.read(magic, sizeof(magic));
    if (!in || std::string(magic) != "CPMCHK2")
        throw std::runtime_error("bad checkpoint magic: " + path);
    std::uint32_t version = 0;
    read_pod(in, version);
    if (version != 2 && version != 3)
        throw std::runtime_error("unsupported checkpoint version");
    std::uint64_t fingerprint = 0;
    int root = 0, max_weight = 0, split_weight = 0;
    int shard_count = 0, shard_index = 0, n_words = 0, m_words = 0;
    read_pod(in, fingerprint); read_pod(in, root); read_pod(in, max_weight);
    read_pod(in, split_weight); read_pod(in, shard_count); read_pod(in, shard_index);
    read_pod(in, n_words); read_pod(in, m_words);
    if (fingerprint != state.fingerprint || root != state.root ||
        max_weight != state.max_weight || split_weight != state.split_weight ||
        shard_count != state.shard_count || shard_index != state.shard_index ||
        n_words != expected_n_words || m_words != expected_m_words)
        throw std::runtime_error("checkpoint does not match requested search");
    read_pod(in, state.frontier_counter);
    read_pod(in, state.resume_count);
    if (version >= 3) read_pod(in, state.started_at_epoch);
    read_pod(in, state.result.elapsed);
    SearchStats &s = state.result.stats;
    read_pod(in, s.states); read_pod(in, s.branches);
    read_pod(in, s.forced); read_pod(in, s.weight_prunes);
    read_pod(in, s.lower_bound_prunes); read_pod(in, s.infeasible_prunes);
    read_pod(in, s.kernel_circuits); read_pod(in, s.stabilizer_circuits);
    std::uint8_t found = 0;
    read_pod(in, found);
    state.result.found = found != 0;
    read_pod(in, state.result.witness_weight);
    read_words(in, state.result.witness.a.data(), n_words);
    std::uint64_t stack_size = 0;
    read_pod(in, stack_size);
    if (stack_size > 10000000ULL) throw std::runtime_error("unreasonable checkpoint stack");
    state.stack.assign(static_cast<std::size_t>(stack_size), DFSNode{});
    for (DFSNode &node : state.stack) {
        read_pod(in, node.weight);
        std::uint8_t partitioned = 0;
        read_pod(in, partitioned);
        node.partitioned = partitioned != 0;
        read_words(in, node.support.a.data(), n_words);
        read_words(in, node.forbidden.a.data(), n_words);
        read_words(in, node.syndrome.a.data(), m_words);
    }
    state.result.root = state.root;
    return true;
}

std::string progress_json(
    const CPMInstance &inst, const Options &opt, const CheckpointState &state,
    const std::string &status, const std::string &reason
) {
    const int effective_check_group_size = opt.check_group_size > 0
        ? opt.check_group_size : (inst.sparse_binary ? 0 : inst.p);
    bool complete = !state.result.found && state.stack.empty();
    std::ostringstream out;
    out << "{\n"
        << "  \"schema\": \"cpm-distance-progress-v4\",\n"
        << "  \"experiment_id\": \"" << json_escape(inst.name) << "_w"
        << opt.max_weight << "\",\n"
        << "  \"task_id\": \"" << opt.side << "_root" << state.root
        << "_shard" << state.shard_index << "of" << state.shard_count << "\",\n"
        << "  \"status\": \"" << status << "\",\n"
        << "  \"requested_scope\": {\"side\":\"" << opt.side
        << "\",\"root\":" << state.root << ",\"max_weight\":"
        << opt.max_weight << ",\"shard_index\":" << state.shard_index
        << ",\"shard_count\":" << state.shard_count
        << ",\"check_group_size\":" << effective_check_group_size << "},\n"
        << "  \"completed_scope\": "
        << (complete ? "[{\"side\":\"" + opt.side + "\",\"root\":" +
                         std::to_string(state.root) + ",\"shard_index\":" +
                         std::to_string(state.shard_index) + "}]" : "[]") << ",\n"
        << "  \"current_cursor_or_frontier\": {\"frontier_nodes_seen\":"
        << state.frontier_counter << ",\"pending_nodes\":" << state.stack.size()
        << "},\n"
        << "  \"units_processed\": " << state.result.stats.states << ",\n"
        << "  \"witness_found\": " << (state.result.found ? "true" : "false") << ",\n"
        << "  \"best_weight\": ";
    if (state.result.found) out << state.result.witness_weight;
    else out << "null";
    out << ",\n  \"started_at\": \""
        << utc_from_epoch(state.started_at_epoch) << "\",\n"
        << "  \"updated_at\": \"" << utc_now() << "\",\n"
        << "  \"elapsed_sec\": " << state.result.elapsed << ",\n"
        << "  \"h_rt_sec\": " << opt.h_rt_sec << ",\n"
        << "  \"soft_deadline_sec\": " << opt.soft_deadline_sec << ",\n"
        << "  \"last_checkpoint_at\": \"" << utc_now() << "\",\n"
        << "  \"resume_count\": " << state.resume_count << ",\n"
        << "  \"termination_reason\": \"" << json_escape(reason) << "\"\n"
        << "}\n";
    return out.str();
}

std::string shard_certificate_json(
    const CPMInstance &inst, const Options &opt, const CheckpointState &state,
    int n, int check_rows, int stabilizer_rows, bool kernel_even_weight,
    int effective_max_weight, const std::string &reason
) {
    const int effective_check_group_size = opt.check_group_size > 0
        ? opt.check_group_size : (inst.sparse_binary ? 0 : inst.p);
    bool complete = !state.result.found && state.stack.empty();
    std::ostringstream out;
    out << "{\n"
        << "  \"schema\": \"cpm-distance-shard-certificate-v4\",\n"
        << "  \"instance\": \"" << json_escape(inst.name) << "\",\n"
        << "  \"mode\": \"" << opt.mode << "\",\n"
        << "  \"side\": \"" << opt.side << "\",\n"
        << "  \"n\": " << n << ", \"P\": " << inst.p
        << ", \"base_columns\": " << inst.l << ",\n"
        << "  \"check_rows\": " << check_rows
        << ", \"stabilizer_rows\": " << stabilizer_rows << ",\n"
        << "  \"check_group_size\": " << effective_check_group_size << ",\n"
        << "  \"max_weight\": " << opt.max_weight
        << ", \"effective_max_weight\": " << effective_max_weight << ",\n"
        << "  \"kernel_even_weight\": " << (kernel_even_weight ? "true" : "false")
        << ", \"qc_root_reduction\": " << (opt.qc_reduce ? "true" : "false") << ",\n"
        << "  \"qc_root_partition\": "
        << (opt.qc_reduce ? "\"earliest_base_column_type_v2\"" : "null") << ",\n"
        << "  \"coset_representative_filter\": "
        << "\"single_stabilizer_weight_lexmax_canonical_v2\",\n"
        << "  \"root\": " << state.root << ", \"split_weight\": "
        << state.split_weight << ", \"shard_index\": " << state.shard_index
        << ", \"shard_count\": " << state.shard_count << ",\n"
        << "  \"frontier_nodes_seen\": " << state.frontier_counter << ",\n"
        << "  \"complete\": " << (complete ? "true" : "false")
        << ", \"found\": " << (state.result.found ? "true" : "false") << ",\n"
        << "  \"witness\": ";
    if (state.result.found) {
        out << "{\"weight\":" << state.result.witness_weight << ",\"support\":";
        print_int_array(out, support_list(state.result.witness, n));
        out << "}";
    } else out << "null";
    const SearchStats &s = state.result.stats;
    out << ",\n  \"stats\": {\"states\":" << s.states
        << ",\"branches\":" << s.branches
        << ",\"forced\":" << s.forced
        << ",\"weight_prunes\":" << s.weight_prunes
        << ",\"lower_bound_prunes\":" << s.lower_bound_prunes
        << ",\"infeasible_prunes\":" << s.infeasible_prunes
        << ",\"kernel_circuits\":" << s.kernel_circuits
        << ",\"stabilizer_circuits\":" << s.stabilizer_circuits << "},\n"
        << "  \"elapsed_sec\": " << state.result.elapsed
        << ", \"resume_count\": " << state.resume_count << ",\n"
        << "  \"termination_reason\": \"" << json_escape(reason) << "\"\n"
        << "}\n";
    return out.str();
}

Options parse_options(int argc, char **argv) {
    Options o;
    auto value = [&](int &i) -> std::string {
        if (++i >= argc) throw std::runtime_error("missing option value");
        return argv[i];
    };
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--input") o.input = value(i);
        else if (a == "--output") o.output = value(i);
        else if (a == "--mode") o.mode = value(i);
        else if (a == "--side") o.side = value(i);
        else if (a == "--max-weight") o.max_weight = std::stoi(value(i));
        else if (a == "--root-start") o.root_start = std::stoi(value(i));
        else if (a == "--root-end") o.root_end = std::stoi(value(i));
        else if (a == "--threads") o.threads = std::stoi(value(i));
        else if (a == "--checkpoint") o.checkpoint = value(i);
        else if (a == "--progress") o.progress = value(i);
        else if (a == "--checkpoint-interval-sec")
            o.checkpoint_interval_sec = std::stoi(value(i));
        else if (a == "--soft-deadline-sec") o.soft_deadline_sec = std::stoi(value(i));
        else if (a == "--h-rt-sec") o.h_rt_sec = std::stoi(value(i));
        else if (a == "--split-weight") o.split_weight = std::stoi(value(i));
        else if (a == "--shard-count") o.shard_count = std::stoi(value(i));
        else if (a == "--shard-index") o.shard_index = std::stoi(value(i));
        else if (a == "--target-tasks") o.target_tasks = std::stoi(value(i));
        else if (a == "--check-group-size")
            o.check_group_size = std::stoi(value(i));
        else if (a == "--partition-only") o.partition_only = true;
        else if (a == "--no-qc-reduce") o.qc_reduce = false;
        else if (a == "--no-stop") o.stop_on_witness = false;
        else if (a == "--quiet") o.quiet = true;
        else if (a == "--help") {
            std::cout << "Usage: cpm_distance --input FILE --mode css|classical "
                      << "--side X|Z --max-weight W [--root-start A --root-end B] "
                      << "[--threads T] [--no-qc-reduce] [--output FILE] "
                      << "[--check-group-size G] "
                      << "[--checkpoint FILE --progress FILE --soft-deadline-sec S "
                      << "--h-rt-sec S "
                      << "--shard-count N --shard-index I --split-weight W]\n";
            std::exit(0);
        } else throw std::runtime_error("unknown option: " + a);
    }
    if (o.input.empty() || o.max_weight < 1) throw std::runtime_error("--input and --max-weight are required");
    if (o.mode != "css" && o.mode != "classical") throw std::runtime_error("--mode must be css or classical");
    if (o.side != "X" && o.side != "Z") throw std::runtime_error("--side must be X or Z");
    o.threads = std::max(1, o.threads);
    if (o.checkpoint_interval_sec < 1)
        throw std::runtime_error("--checkpoint-interval-sec must be positive");
    if (o.soft_deadline_sec < 0)
        throw std::runtime_error("--soft-deadline-sec must be nonnegative");
    if (o.h_rt_sec < 0)
        throw std::runtime_error("--h-rt-sec must be nonnegative");
    if (o.shard_count < 1 || o.shard_index < 0 || o.shard_index >= o.shard_count)
        throw std::runtime_error("invalid shard index/count");
    if (o.split_weight < 1)
        throw std::runtime_error("--split-weight must be positive");
    if (o.target_tasks < 1)
        throw std::runtime_error("--target-tasks must be positive");
    if (o.check_group_size < 0)
        throw std::runtime_error("--check-group-size must be nonnegative");
    if (o.partition_only && (o.threads != 1 || o.shard_count <= 1))
        throw std::runtime_error(
            "--partition-only requires --threads 1 and --shard-count > 1");
    if (!o.checkpoint.empty() && (o.output.empty() || o.progress.empty()))
        throw std::runtime_error("checkpoint mode requires --output and --progress");
    return o;
}

} // namespace

int main(int argc, char **argv) {
    try {
        Options opt = parse_options(argc, argv);
        CPMInstance inst = read_instance(opt.input);
        int n = inst.n;
        auto hx = inst.sparse_binary
            ? inst.hx_rows : lifted_rows(inst.ex, inst.p, inst.l);
        auto hz = inst.sparse_binary
            ? inst.hz_rows : lifted_rows(inst.ez, inst.p, inst.l);

        const std::vector<std::vector<int>> *checks = nullptr;
        const std::vector<std::vector<int>> *stabilizers = nullptr;
        static const std::vector<std::vector<int>> empty;
        if (opt.mode == "classical") {
            checks = (opt.side == "X") ? &hx : &hz;
            stabilizers = &empty;
        } else if (opt.side == "X") {
            checks = &hz; stabilizers = &hx;
        } else {
            checks = &hx; stabilizers = &hz;
        }
        const int effective_check_group_size = opt.check_group_size > 0
            ? opt.check_group_size : (inst.sparse_binary ? 0 : inst.p);

        int root_count = n;
        if (opt.qc_reduce) {
            if (inst.sparse_binary)
                throw std::runtime_error(
                    "sparse binary input requires --no-qc-reduce");
            // A simultaneous cyclic shift within every P-block is a code
            // automorphism for every CPM modulus, prime or composite.  Shift
            // a support coordinate in its earliest occupied base column to
            // residue zero; that coordinate is then the minimum coordinate.
            // Hence one root per base column is complete without any
            // primality or weight restriction.
            root_count = inst.l;
        }
        int begin = std::max(0, opt.root_start);
        int end = opt.root_end < 0 ? root_count : std::min(root_count, opt.root_end);
        if (begin >= end) throw std::runtime_error("empty root range");

        // An odd requested cutoff can be reduced by one when the all-one
        // vector lies in the selected check rowspace.
        bool kernel_even_weight = kernel_has_even_weight(*checks, n);
        int effective_max_weight = opt.max_weight;
        if (kernel_even_weight && (effective_max_weight & 1))
            --effective_max_weight;

        if (!opt.checkpoint.empty()) {
            if (end - begin != 1)
                throw std::runtime_error(
                    "checkpoint mode requires exactly one root");
            if (opt.shard_count > 1 && opt.split_weight > effective_max_weight)
                throw std::runtime_error(
                    "--split-weight must not exceed the effective maximum weight");

            const int root = begin;
            const int n_words = (n + 63) / 64;
            const int m_words = (static_cast<int>(checks->size()) + 63) / 64;
            CheckpointState state;
            state.root = root;
            state.result.root = root;
            state.max_weight = effective_max_weight;
            state.split_weight = opt.split_weight;
            state.shard_count = opt.shard_count;
            state.shard_index = opt.shard_index;
            state.fingerprint = run_fingerprint(
                opt.input, opt, root, effective_max_weight);

            std::atomic<bool> checkpoint_stop{false};
            SearchEngine engine(
                n, *checks, *stabilizers, effective_max_weight, checkpoint_stop,
                opt.qc_reduce ? inst.l : 0, kernel_even_weight,
                effective_check_group_size);
            bool resumed = load_checkpoint(
                opt.checkpoint, state, n_words, m_words);
            if (resumed) {
                ++state.resume_count;
                // Older checkpoints contain an unfiltered explicit frontier.
                // Removing reducible nodes here upgrades them safely: xor with
                // an over-half-overlap stabilizer strictly lowers weight while
                // preserving the logical coset.
                engine.filter_reducible_frontier(state);
            } else {
                state.started_at_epoch = static_cast<std::int64_t>(std::time(nullptr));
                engine.initialize_checkpoint(state, root);
            }
            if (state.started_at_epoch == 0) {
                state.started_at_epoch = static_cast<std::int64_t>(
                    std::time(nullptr) - static_cast<std::time_t>(state.result.elapsed));
            }

            signal_stop_requested = 0;
            std::signal(SIGTERM, request_signal_stop);
#ifdef SIGUSR1
            std::signal(SIGUSR1, request_signal_stop);
#endif
            auto run_started = std::chrono::steady_clock::now();
            auto last_checkpoint = run_started;
            const double elapsed_before = state.result.elapsed;
            auto update_elapsed = [&] {
                state.result.elapsed = elapsed_before +
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - run_started).count();
            };
            auto publish = [&](const std::string &status, const std::string &reason) {
                update_elapsed();
                write_checkpoint_atomic(opt.checkpoint, state, n_words, m_words);
                atomic_write_text(
                    opt.progress,
                    progress_json(inst, opt, state, status, reason));
                last_checkpoint = std::chrono::steady_clock::now();
            };

            if (opt.partition_only) {
                std::vector<DFSNode> assigned;
                auto partition_snapshot = [&] (
                    const std::string &status, const std::string &reason
                ) {
                    update_elapsed();
                    CheckpointState saved = state;
                    saved.stack.insert(
                        saved.stack.end(), assigned.begin(), assigned.end());
                    write_checkpoint_atomic(
                        opt.checkpoint, saved, n_words, m_words);
                    atomic_write_text(
                        opt.progress,
                        progress_json(inst, opt, saved, status, reason));
                    last_checkpoint = std::chrono::steady_clock::now();
                };
                std::string partition_reason;
                std::uint64_t partition_steps = 0;
                while (!state.stack.empty() && !state.result.found) {
                    engine.checkpoint_step(state, true, &assigned);
                    ++partition_steps;
                    if ((partition_steps & 0x3fffULL) != 0) continue;
                    const auto now = std::chrono::steady_clock::now();
                    if (signal_stop_requested) {
                        partition_reason = "scheduler_signal";
                        break;
                    }
                    if (opt.soft_deadline_sec > 0 &&
                        std::chrono::duration<double>(now - run_started).count()
                            >= opt.soft_deadline_sec) {
                        partition_reason = "soft_deadline";
                        break;
                    }
                    if (std::chrono::duration<double>(
                            now - last_checkpoint).count() >=
                        opt.checkpoint_interval_sec) {
                        partition_snapshot("partitioning", "partition_frontier");
                    }
                }
                if (state.result.found) {
                    state.stack.clear();
                    update_elapsed();
                    write_checkpoint_atomic(
                        opt.checkpoint, state, n_words, m_words);
                    atomic_write_text(
                        opt.progress,
                        progress_json(
                            inst, opt, state, "found", "witness_found"));
                    atomic_write_text(
                        opt.output,
                        shard_certificate_json(
                            inst, opt, state, n, checks->size(),
                            stabilizers->size(), kernel_even_weight,
                            effective_max_weight, "witness_found"));
                    return 0;
                }
                state.stack.insert(
                    state.stack.end(), assigned.begin(), assigned.end());
                update_elapsed();
                if (partition_reason.empty()) {
                    write_checkpoint_atomic(
                        opt.checkpoint, state, n_words, m_words);
                    atomic_write_text(
                        opt.progress,
                        progress_json(
                            inst, opt, state, "partitioned",
                            "shard_frontier_ready"));
                    return 0;
                }
                write_checkpoint_atomic(
                    opt.checkpoint, state, n_words, m_words);
                atomic_write_text(
                    opt.progress,
                    progress_json(
                        inst, opt, state, "timed_out", partition_reason));
                return 0;
            }

            if (opt.threads > 1) {
                if (!resumed)
                    throw std::runtime_error(
                        "parallel checkpoint mode requires an existing checkpoint");
                // Refine the saved frontier into independent exact subtrees.
                // Periodic snapshots retain every unfinished subtree in its
                // original form (including a subtree currently being worked
                // on), so partial recursive work may be repeated after an
                // interruption but exact coverage can never be lost.
                std::vector<DFSNode> frontier = state.stack;
                // Refine the leading (DFS-hard) portion of the saved frontier
                // before assigning work.  A raw stack can contain a few
                // enormous subtrees beside many tiny ones; replacing a node
                // by all of its exact children makes interrupted work finite
                // without changing coverage.  Stop close to the target rather
                // than expanding a whole round, which could overshoot by 12x.
                const std::size_t target_tasks = std::max<std::size_t>(
                    static_cast<std::size_t>(opt.target_tasks),
                    static_cast<std::size_t>(opt.threads) * 4096);
                for (int round = 0;
                     round < 10 && !frontier.empty()
                         && frontier.size() < target_tasks;
                     ++round) {
                    std::vector<DFSNode> refined;
                    refined.reserve(target_tasks + 16);
                    std::size_t index = 0;
                    for (; index < frontier.size(); ++index) {
                        const std::size_t projected =
                            refined.size() + (frontier.size() - index);
                        if (projected >= target_tasks) break;
                        const DFSNode &node = frontier[index];
                        CheckpointState fragment;
                        fragment.root = root;
                        fragment.result.root = root;
                        fragment.max_weight = effective_max_weight;
                        fragment.split_weight = opt.split_weight;
                        fragment.shard_count = 1;
                        fragment.shard_index = 0;
                        fragment.stack.push_back(node);
                        engine.checkpoint_step(fragment);
                        state.result.stats += fragment.result.stats;
                        if (fragment.result.found) {
                            state.result.found = true;
                            state.result.witness = fragment.result.witness;
                            state.result.witness_weight =
                                fragment.result.witness_weight;
                            break;
                        }
                        refined.insert(
                            refined.end(),
                            std::make_move_iterator(fragment.stack.begin()),
                            std::make_move_iterator(fragment.stack.end()));
                    }
                    if (state.result.found) break;
                    refined.insert(
                        refined.end(),
                        std::make_move_iterator(frontier.begin() + index),
                        std::make_move_iterator(frontier.end()));
                    frontier.swap(refined);
                }
                if (state.result.found) {
                    state.stack.clear();
                    update_elapsed();
                    write_checkpoint_atomic(
                        opt.checkpoint, state, n_words, m_words);
                    atomic_write_text(
                        opt.progress,
                        progress_json(
                            inst, opt, state, "found", "witness_found"));
                    atomic_write_text(
                        opt.output,
                        shard_certificate_json(
                            inst, opt, state, n, checks->size(),
                            stabilizers->size(), kernel_even_weight,
                            effective_max_weight, "witness_found"));
                    std::cout << progress_json(
                        inst, opt, state, "found", "witness_found");
                    return 0;
                }
                const int worker_count = std::min<int>(
                    opt.threads, std::max<std::size_t>(1, frontier.size()));
                std::atomic<std::size_t> next_node{0};
                std::atomic<int> finished_workers{0};
                std::atomic<bool> parallel_stop{false};
                std::mutex parallel_mutex;
                std::vector<unsigned char> task_complete(frontier.size(), 0);
                std::size_t completed_tasks = 0;
                SearchStats completed_stats;
                bool parallel_found = false;
                RootResult parallel_witness;
                auto parallel_snapshot = [&] (
                    const std::string &status, const std::string &reason
                ) {
                    CheckpointState saved = state;
                    {
                        std::lock_guard<std::mutex> lock(parallel_mutex);
                        saved.result.stats += completed_stats;
                        saved.stack.clear();
                        saved.stack.reserve(frontier.size() - completed_tasks);
                        for (std::size_t index = 0;
                             index < frontier.size(); ++index) {
                            if (!task_complete[index])
                                saved.stack.push_back(frontier[index]);
                        }
                    }
                    saved.result.elapsed = elapsed_before +
                        std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - run_started).count();
                    write_checkpoint_atomic(
                        opt.checkpoint, saved, n_words, m_words);
                    atomic_write_text(
                        opt.progress,
                        progress_json(inst, opt, saved, status, reason));
                    last_checkpoint = std::chrono::steady_clock::now();
                    return saved;
                };
                parallel_snapshot("running", "parallel_resume");
                std::vector<std::thread> workers;
                for (int tid = 0; tid < worker_count; ++tid) {
                    workers.emplace_back([&] {
                        SearchEngine worker_engine(
                            n, *checks, *stabilizers, effective_max_weight,
                            parallel_stop, opt.qc_reduce ? inst.l : 0,
                            kernel_even_weight, effective_check_group_size);
                        while (!parallel_stop.load(std::memory_order_relaxed)) {
                            std::size_t index = next_node.fetch_add(1);
                            if (index >= frontier.size()) break;
                            RootResult result = worker_engine.run_subtree(
                                root, frontier[index]);
                            std::lock_guard<std::mutex> lock(parallel_mutex);
                            if (result.complete || result.found) {
                                task_complete[index] = 1;
                                ++completed_tasks;
                                completed_stats += result.stats;
                            }
                            if (result.found && !parallel_found) {
                                parallel_found = true;
                                parallel_witness = result;
                            }
                        }
                        ++finished_workers;
                    });
                }
                std::string parallel_reason;
                while (finished_workers.load() < worker_count) {
                    if (signal_stop_requested) {
                        parallel_reason = "scheduler_signal";
                        parallel_stop.store(true);
                        break;
                    }
                    if (opt.soft_deadline_sec > 0 &&
                        std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - run_started).count()
                            >= opt.soft_deadline_sec) {
                        parallel_reason = "soft_deadline";
                        parallel_stop.store(true);
                        break;
                    }
                    if (std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - last_checkpoint)
                            .count() >= opt.checkpoint_interval_sec) {
                        parallel_snapshot("running", "periodic_checkpoint");
                    }
                    std::this_thread::sleep_for(std::chrono::milliseconds(100));
                }
                for (auto &worker : workers) worker.join();

                {
                    std::lock_guard<std::mutex> lock(parallel_mutex);
                    state.result.stats += completed_stats;
                    if (parallel_found) {
                        state.result.found = true;
                        state.result.witness = parallel_witness.witness;
                        state.result.witness_weight =
                            parallel_witness.witness_weight;
                    }
                }
                const bool found = state.result.found;
                const bool all_complete = completed_tasks == frontier.size();
                update_elapsed();
                if (found || all_complete) {
                    state.stack.clear();
                    state.result.complete = all_complete && !found;
                    parallel_reason = found ? "witness_found" : "search_exhausted";
                    const std::string status = found ? "found" : "complete";
                    write_checkpoint_atomic(
                        opt.checkpoint, state, n_words, m_words);
                    atomic_write_text(
                        opt.progress,
                        progress_json(inst, opt, state, status, parallel_reason));
                    atomic_write_text(
                        opt.output,
                        shard_certificate_json(
                            inst, opt, state, n, checks->size(),
                            stabilizers->size(), kernel_even_weight,
                            effective_max_weight, parallel_reason));
                    std::cout << progress_json(
                        inst, opt, state, status, parallel_reason);
                    return 0;
                }
                if (parallel_reason.empty()) parallel_reason = "parallel_incomplete";
                // Rebuild a durable checkpoint containing only unfinished
                // whole subtrees.  Work inside an interrupted subtree is
                // intentionally discarded, never mistaken for a proof.
                CheckpointState saved = state;
                saved.result.stats = state.result.stats;
                saved.stack.clear();
                for (std::size_t index = 0; index < frontier.size(); ++index)
                    if (!task_complete[index]) saved.stack.push_back(frontier[index]);
                write_checkpoint_atomic(
                    opt.checkpoint, saved, n_words, m_words);
                atomic_write_text(
                    opt.progress,
                    progress_json(
                        inst, opt, saved, "timed_out", parallel_reason));
                std::cout << progress_json(
                    inst, opt, saved, "timed_out", parallel_reason);
                return 0;
            }

            publish("running", resumed ? "resumed" : "started");
            std::string reason;
            std::uint64_t loop_steps = 0;
            while (!state.stack.empty() && !state.result.found) {
                engine.checkpoint_step(state);
                ++loop_steps;
                if ((loop_steps & 0x3fffULL) != 0) continue;
                auto now = std::chrono::steady_clock::now();
                if (signal_stop_requested) {
                    reason = "scheduler_signal";
                    break;
                }
                if (opt.soft_deadline_sec > 0 &&
                    std::chrono::duration<double>(now - run_started).count() >=
                        opt.soft_deadline_sec) {
                    reason = "soft_deadline";
                    break;
                }
                if (std::chrono::duration<double>(now - last_checkpoint).count() >=
                    opt.checkpoint_interval_sec) {
                    publish("running", "periodic_checkpoint");
                }
            }

            std::string status;
            if (state.result.found) {
                status = "found";
                reason = "witness_found";
            } else if (state.stack.empty()) {
                status = "complete";
                reason = "search_exhausted";
                state.result.complete = true;
            } else {
                status = "timed_out";
                if (reason.empty()) reason = "soft_deadline";
            }
            publish(status, reason);
            if (status == "complete" || status == "found") {
                atomic_write_text(
                    opt.output,
                    shard_certificate_json(
                        inst, opt, state, n, checks->size(), stabilizers->size(),
                        kernel_even_weight, effective_max_weight, reason));
            }
            std::cout << progress_json(inst, opt, state, status, reason);
            return 0;
        }

        std::atomic<bool> global_stop{false};
        std::atomic<int> next_root{begin};
        std::mutex result_mutex;
        std::vector<RootResult> results;
        results.reserve(end - begin);
        auto all_started = std::chrono::steady_clock::now();

        int worker_count = std::min(opt.threads, end - begin);
        std::vector<std::thread> workers;
        for (int tid = 0; tid < worker_count; ++tid) {
            workers.emplace_back([&] {
                SearchEngine engine(
                    n, *checks, *stabilizers, effective_max_weight, global_stop,
                    opt.qc_reduce ? inst.l : 0, kernel_even_weight,
                    effective_check_group_size);
                while (true) {
                    if (opt.stop_on_witness && global_stop.load(std::memory_order_relaxed)) break;
                    int root = next_root.fetch_add(1);
                    if (root >= end) break;
                    RootResult rr = engine.run_root(root);
                    {
                        std::lock_guard<std::mutex> lock(result_mutex);
                        results.push_back(rr);
                        if (!opt.quiet) {
                            std::cerr << "root " << root << " states=" << rr.stats.states
                                      << " found=" << rr.found << " sec=" << rr.elapsed << "\n";
                        }
                    }
                    if (rr.found && opt.stop_on_witness) break;
                    if (!opt.stop_on_witness) global_stop.store(false);
                }
            });
        }
        for (auto &t : workers) t.join();
        std::sort(results.begin(), results.end(), [](const auto &a, const auto &b) { return a.root < b.root; });

        SearchStats total;
        bool found = false;
        RootResult witness;
        std::vector<int> completed;
        for (const auto &r : results) {
            total += r.stats;
            completed.push_back(r.root);
            if (r.found && !found) { found = true; witness = r; }
        }
        bool complete = !found && static_cast<int>(results.size()) == end - begin;
        double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - all_started).count();

        std::ostringstream json;
        json << "{\n"
             << "  \"schema\": \"cpm-distance-certificate-v1\",\n"
             << "  \"instance\": \"" << json_escape(inst.name) << "\",\n"
             << "  \"mode\": \"" << opt.mode << "\",\n"
             << "  \"side\": \"" << opt.side << "\",\n"
             << "  \"n\": " << n << ",\n"
             << "  \"P\": " << inst.p << ",\n"
             << "  \"base_columns\": " << inst.l << ",\n"
             << "  \"check_rows\": " << checks->size() << ",\n"
             << "  \"stabilizer_rows\": " << stabilizers->size() << ",\n"
             << "  \"check_group_size\": "
             << effective_check_group_size << ",\n"
             << "  \"max_weight\": " << opt.max_weight << ",\n"
             << "  \"effective_max_weight\": " << effective_max_weight << ",\n"
             << "  \"kernel_even_weight\": " << (kernel_even_weight ? "true" : "false") << ",\n"
             << "  \"qc_root_reduction\": " << (opt.qc_reduce ? "true" : "false") << ",\n"
             << "  \"qc_root_partition\": "
             << (opt.qc_reduce ? "\"earliest_base_column_type_v2\"" : "null") << ",\n"
             << "  \"root_start\": " << begin << ",\n"
             << "  \"root_end\": " << end << ",\n"
             << "  \"completed_roots\": ";
        print_int_array(json, completed);
        json << ",\n  \"complete\": " << (complete ? "true" : "false")
             << ",\n  \"found\": " << (found ? "true" : "false") << ",\n"
             << "  \"witness\": ";
        if (found) {
            json << "{\"weight\":" << witness.witness_weight << ",\"support\":";
            print_int_array(json, support_list(witness.witness, n));
            json << "}";
        } else json << "null";
        json << ",\n  \"stats\": {"
             << "\"states\":" << total.states
             << ",\"branches\":" << total.branches
             << ",\"forced\":" << total.forced
             << ",\"weight_prunes\":" << total.weight_prunes
             << ",\"lower_bound_prunes\":" << total.lower_bound_prunes
             << ",\"infeasible_prunes\":" << total.infeasible_prunes
             << ",\"kernel_circuits\":" << total.kernel_circuits
             << ",\"stabilizer_circuits\":" << total.stabilizer_circuits
             << "},\n  \"elapsed_sec\": " << elapsed << "\n}\n";

        if (!opt.output.empty()) {
            std::ofstream out(opt.output);
            if (!out) throw std::runtime_error("cannot write output: " + opt.output);
            out << json.str();
        }
        std::cout << json.str();
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "cpm_distance: " << e.what() << "\n";
        return 2;
    }
}
