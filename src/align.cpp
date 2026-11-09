#include "internal.hpp"
#include "msakit/core.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <queue>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace msakit {
namespace {
using Score = double;
using Counts = std::array<double, 4>;

bool is_residue(char value) {
    return value == 'A' || value == 'C' || value == 'G' || value == 'T' ||
           value == 'U' || value == 'N';
}

std::size_t base_index(char value) {
    switch (value) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T':
        case 'U': return 3;
        default: return 4;  // N and gaps do not contribute to substitution scores.
    }
}

struct Profile {
    std::vector<std::size_t> input_rows;
    std::vector<std::string> rows;
    std::vector<Counts> columns;
};

Profile make_leaf(std::size_t row, std::string_view sequence) {
    Profile profile;
    profile.input_rows.push_back(row);
    profile.rows.emplace_back(sequence);
    profile.columns.resize(sequence.size());
    for (std::size_t column = 0; column < sequence.size(); ++column) {
        profile.columns[column].fill(0.0);
        const auto index = base_index(sequence[column]);
        if (index < 4) profile.columns[column][index] = 1.0;
    }
    return profile;
}

void rebuild_columns(Profile& profile) {
    const auto width = profile.rows.front().size();
    profile.columns.assign(width, Counts{});
    const auto denominator = static_cast<Score>(profile.rows.size());
    for (const auto& row : profile.rows) {
        for (std::size_t column = 0; column < width; ++column) {
            const auto index = base_index(row[column]);
            if (index < 4) profile.columns[column][index] += 1.0 / denominator;
        }
    }
}

Score column_score(const Counts& left, const Counts& right) {
    Score matches = 0.0;
    Score left_occupancy = 0.0;
    Score right_occupancy = 0.0;
    for (std::size_t base = 0; base < 4; ++base) {
        left_occupancy += left[base];
        right_occupancy += right[base];
        matches += left[base] * right[base];
    }
    if (left_occupancy == 0.0 || right_occupancy == 0.0) return 0.0;
    return 2.0 * matches - left_occupancy * right_occupancy;
}

// Global profile/profile alignment with a linear gap cost.  The score uses
// nucleotide frequencies in each profile column, so every sequence in both
// profiles contributes to the decision rather than only one chosen center.
Profile merge_profiles(const Profile& left, const Profile& right) {
    const auto height = detail::checked_add(left.columns.size(), 1);
    const auto width = detail::checked_add(right.columns.size(), 1);
    const auto cells = detail::checked_multiply(height, width);
    (void)detail::checked_multiply(cells, sizeof(Score));
    if (cells > std::vector<Score>().max_size() ||
        cells > std::vector<unsigned char>().max_size()) {
        throw detail::Error(detail::kCoreError, "profile alignment matrix exceeds container size");
    }

    std::vector<Score> scores(cells, 0.0);
    std::vector<unsigned char> trace(cells, 0);
    const auto index = [width](std::size_t i, std::size_t j) { return i * width + j; };
    constexpr Score gap_penalty = -1.0;
    for (std::size_t i = 1; i < height; ++i) {
        scores[index(i, 0)] = static_cast<Score>(i) * gap_penalty;
        trace[index(i, 0)] = 1;
    }
    for (std::size_t j = 1; j < width; ++j) {
        scores[index(0, j)] = static_cast<Score>(j) * gap_penalty;
        trace[index(0, j)] = 2;
    }
    for (std::size_t i = 1; i < height; ++i) {
        for (std::size_t j = 1; j < width; ++j) {
            const auto diagonal = scores[index(i - 1, j - 1)] +
                                  column_score(left.columns[i - 1], right.columns[j - 1]);
            const auto up = scores[index(i - 1, j)] + gap_penalty;
            const auto across = scores[index(i, j - 1)] + gap_penalty;
            auto best = diagonal;
            unsigned char direction = 0;
            if (up > best) {
                best = up;
                direction = 1;
            }
            if (across > best) {
                best = across;
                direction = 2;
            }
            scores[index(i, j)] = best;
            trace[index(i, j)] = direction;
        }
    }

    std::vector<unsigned char> path;
    path.reserve(detail::checked_add(left.columns.size(), right.columns.size()));
    std::size_t i = left.columns.size();
    std::size_t j = right.columns.size();
    while (i != 0 || j != 0) {
        const auto direction = trace[index(i, j)];
        if (i != 0 && j != 0 && direction == 0) {
            path.push_back(0);  // consume one column from each profile
            --i;
            --j;
        } else if (i != 0 && (j == 0 || direction == 1)) {
            path.push_back(1);  // consume left column; insert gaps in right
            --i;
        } else if (j != 0) {
            path.push_back(2);  // consume right column; insert gaps in left
            --j;
        } else {
            throw detail::Error(detail::kCoreError, "invalid profile traceback state");
        }
    }
    std::reverse(path.begin(), path.end());

    Profile merged;
    merged.input_rows.reserve(left.input_rows.size() + right.input_rows.size());
    merged.input_rows.insert(merged.input_rows.end(), left.input_rows.begin(), left.input_rows.end());
    merged.input_rows.insert(merged.input_rows.end(), right.input_rows.begin(), right.input_rows.end());
    const auto output_width = path.size();
    merged.rows.reserve(merged.input_rows.size());
    for (std::size_t row = 0; row < left.rows.size(); ++row) {
        merged.rows.emplace_back();
        merged.rows.back().reserve(output_width);
    }
    const auto left_count = merged.rows.size();
    for (std::size_t row = 0; row < right.rows.size(); ++row) {
        merged.rows.emplace_back();
        merged.rows.back().reserve(output_width);
    }

    std::size_t left_column = 0;
    std::size_t right_column = 0;
    for (const auto direction : path) {
        if (direction == 0 || direction == 1) {
            for (std::size_t row = 0; row < left_count; ++row) {
                merged.rows[row].push_back(left.rows[row][left_column]);
            }
            ++left_column;
        } else {
            for (std::size_t row = 0; row < left_count; ++row) merged.rows[row].push_back('-');
        }
        if (direction == 0 || direction == 2) {
            for (std::size_t row = 0; row < right.rows.size(); ++row) {
                merged.rows[left_count + row].push_back(right.rows[row][right_column]);
            }
            ++right_column;
        } else {
            for (std::size_t row = 0; row < right.rows.size(); ++row) {
                merged.rows[left_count + row].push_back('-');
            }
        }
    }
    rebuild_columns(merged);
    return merged;
}

std::vector<Score> kmer_vector(std::string_view sequence) {
    // A 3-mer composition gives a cheap distance for building the guide tree.
    // T and U are treated as the same nucleotide, matching the core contract.
    constexpr std::size_t alphabet = 4;
    constexpr std::size_t kmer_count = alphabet * alphabet * alphabet;
    std::vector<Score> result(kmer_count, 0.0);
    if (sequence.size() < 3) return result;
    for (std::size_t i = 0; i + 2 < sequence.size(); ++i) {
        const auto a = base_index(sequence[i]);
        const auto b = base_index(sequence[i + 1]);
        const auto c = base_index(sequence[i + 2]);
        if (a < 4 && b < 4 && c < 4) result[a * 16 + b * 4 + c] += 1.0;
    }
    Score norm = 0.0;
    for (const auto value : result) norm += value * value;
    if (norm > 0.0) {
        norm = std::sqrt(norm);
        for (auto& value : result) value /= norm;
    }
    return result;
}

Score kmer_distance(const std::vector<Score>& left, const std::vector<Score>& right) {
    Score dot = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) dot += left[i] * right[i];
    // Short sequences without a complete 3-mer are deliberately left maximally
    // distant; deterministic input order resolves ties in the guide tree.
    return 1.0 - std::clamp(dot, 0.0, 1.0);
}

struct PairDistance {
    Score value;
    std::size_t left;
    std::size_t right;
};

struct FartherFirst {
    bool operator()(const PairDistance& a, const PairDistance& b) const {
        if (a.value != b.value) return a.value > b.value;
        if (a.left != b.left) return a.left > b.left;
        return a.right > b.right;
    }
};

Profile progressive_alignment(const std::vector<SequenceRecord>& input) {
    const auto count = input.size();
    if (count == 1) return make_leaf(0, input.front().sequence);
    const auto node_count = detail::checked_multiply(count, 2);
    if (node_count > std::vector<Profile>().max_size()) {
        throw detail::Error(detail::kCoreError, "guide tree exceeds container size");
    }
    (void)detail::checked_multiply(node_count, node_count);

    std::vector<std::vector<Score>> distances(node_count, std::vector<Score>(node_count, 0.0));
    std::vector<std::size_t> cluster_sizes(node_count, 0);
    std::vector<bool> active(node_count, false);
    std::vector<Profile> profiles(node_count);
    std::vector<std::vector<Score>> kmers;
    kmers.reserve(count);
    std::priority_queue<PairDistance, std::vector<PairDistance>, FartherFirst> candidates;

    for (std::size_t i = 0; i < count; ++i) {
        profiles[i] = make_leaf(i, input[i].sequence);
        cluster_sizes[i] = 1;
        active[i] = true;
        kmers.push_back(kmer_vector(input[i].sequence));
    }
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = i + 1; j < count; ++j) {
            distances[i][j] = distances[j][i] = kmer_distance(kmers[i], kmers[j]);
            candidates.push({distances[i][j], i, j});
        }
    }

    std::size_t next_node = count;
    while (next_node < 2 * count - 1) {
        PairDistance nearest{};
        bool found = false;
        while (!candidates.empty()) {
            nearest = candidates.top();
            candidates.pop();
            if (active[nearest.left] && active[nearest.right] &&
                distances[nearest.left][nearest.right] == nearest.value) {
                found = true;
                break;
            }
        }
        if (!found) throw detail::Error(detail::kCoreError, "guide tree construction failed");

        const auto left = nearest.left;
        const auto right = nearest.right;
        profiles[next_node] = merge_profiles(profiles[left], profiles[right]);
        cluster_sizes[next_node] = cluster_sizes[left] + cluster_sizes[right];
        active[left] = false;
        active[right] = false;
        active[next_node] = true;
        for (std::size_t other = 0; other < next_node; ++other) {
            if (!active[other]) continue;
            const auto weighted =
                (static_cast<Score>(cluster_sizes[left]) * distances[left][other] +
                 static_cast<Score>(cluster_sizes[right]) * distances[right][other]) /
                static_cast<Score>(cluster_sizes[next_node]);
            distances[next_node][other] = distances[other][next_node] = weighted;
            candidates.push({weighted, std::min(other, next_node), std::max(other, next_node)});
        }
        ++next_node;
    }
    return std::move(profiles[next_node - 1]);
}
}  // namespace

// MAFFT-inspired progressive profile alignment.  This keeps the MSA-Kit core
// interface self-contained: k-mer distances build a UPGMA-style guide tree,
// then profiles are aligned progressively with a nucleotide-frequency score.
Alignment align(const std::vector<SequenceRecord>& input) {
    if (input.empty()) {
        throw detail::Error(detail::kCoreError, "the alignment core requires at least one record");
    }
    for (const auto& record : input) {
        if (record.sequence.empty() ||
            !std::all_of(record.sequence.begin(), record.sequence.end(), is_residue)) {
            throw detail::Error(detail::kCoreError,
                                "the alignment core requires nonempty cleaned ACGTUN sequences");
        }
    }

    auto profile = progressive_alignment(input);
    Alignment result;
    result.rows.resize(input.size());
    for (std::size_t row = 0; row < profile.rows.size(); ++row) {
        result.rows[profile.input_rows[row]] = std::move(profile.rows[row]);
    }
    return result;
}
}  // namespace msakit
