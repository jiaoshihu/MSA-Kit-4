#include "internal.hpp"
#include "msakit/core.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// Experimental progressive profile-profile DNA/RNA MSA.
// Guide tree: average-linkage clustering on 4-mer set distances.
// Alignment: affine-gap dynamic programming on profile column frequencies.
// All code is self-contained in src/align.cpp; no other project files need edits.
namespace msakit {
namespace {

bool is_residue(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'U' || c == 'N';
}

// 0=A 1=C 2=G 3=T/U 4=N 5=gap
int code(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': case 'U': return 3;
        case 'N': return 4;
        default: return 5;
    }
}

struct Profile {
    std::vector<std::size_t> ids;
    std::vector<std::string> rows;
};

struct Column {
    std::array<float, 6> f{};
    float occupancy = 0;
};

std::vector<Column> columns(const Profile& p) {
    const auto width = p.rows.front().size();
    std::vector<Column> result(width);
    const float inv = 1.0f / static_cast<float>(p.rows.size());
    for (const auto& s : p.rows) {
        for (std::size_t i = 0; i < width; ++i) {
            result[i].f[code(s[i])] += inv;
        }
    }
    for (auto& c : result) c.occupancy = 1.0f - c.f[5];
    return result;
}

// Expected substitution score over the cross-product of two profile columns.
// Residue-residue: +3 match, -2 mismatch, N: 0.
// Existing gap/residue: -2; gap-gap: 0.
float cross_score(const Column& a, const Column& b) {
    float value = 0;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            const float s = (i == 4 || j == 4) ? 0.0f : (i == j ? 3.0f : -2.0f);
            value += a.f[i] * b.f[j] * s;
        }
    }
    value -= 2.0f * (a.occupancy * b.f[5] + a.f[5] * b.occupancy);
    return value;
}

// Memory upper bound for three score matrices and three traceback matrices.
// Avoid exceeding the benchmark memory limit, even with large alignments.
constexpr std::size_t kMaxCells = 3000000;
constexpr float NEG = -1.0e30f;

bool merge_profiles(const Profile& a, const Profile& b, Profile& out) {
    const std::size_t n = a.rows.front().size();
    const std::size_t m = b.rows.front().size();
    if (n + 1 > std::numeric_limits<std::size_t>::max() / (m + 1)) return false;
    const std::size_t cells = (n + 1) * (m + 1);
    if (cells > kMaxCells) return false;
    const auto ca = columns(a), cb = columns(b);
    const auto idx = [m](std::size_t i, std::size_t j) { return i * (m + 1) + j; };
    std::vector<float> M(cells, NEG), X(cells, NEG), Y(cells, NEG);
    std::vector<std::uint8_t> tM(cells, 0), tX(cells, 0), tY(cells, 0);
    M[0] = 0;
    // A column aligned to new gaps in B (X), or vice versa (Y).
    // Gap open 5, extension 0.7, weighted by profile occupancy.
    for (std::size_t i = 1; i <= n; ++i) {
        const float w = std::max(0.1f, ca[i - 1].occupancy);
        X[idx(i, 0)] = (i == 1 ? M[0] - 5.0f * w : X[idx(i-1,0)] - 0.7f * w);
        tX[idx(i, 0)] = i == 1 ? 0 : 1;
    }
    for (std::size_t j = 1; j <= m; ++j) {
        const float w = std::max(0.1f, cb[j - 1].occupancy);
        Y[idx(0, j)] = (j == 1 ? M[0] - 5.0f * w : Y[idx(0,j-1)] - 0.7f * w);
        tY[idx(0, j)] = j == 1 ? 0 : 2;
    }
    for (std::size_t i = 1; i <= n; ++i) {
        for (std::size_t j = 1; j <= m; ++j) {
            const std::size_t cur = idx(i, j), diag = idx(i-1,j-1);
            float best = M[diag]; std::uint8_t from = 0;
            if (X[diag] > best) {best = X[diag]; from = 1;}
            if (Y[diag] > best) {best = Y[diag]; from = 2;}
            M[cur] = best + cross_score(ca[i-1], cb[j-1]); tM[cur] = from;

            const float wx = std::max(0.1f, ca[i-1].occupancy);
            const float xo = M[idx(i-1,j)] - 5.0f * wx;
            const float xe = X[idx(i-1,j)] - 0.7f * wx;
            X[cur] = xo >= xe ? xo : xe;
            tX[cur] = xo >= xe ? 0 : 1;

            const float wy = std::max(0.1f, cb[j-1].occupancy);
            const float yo = M[idx(i,j-1)] - 5.0f * wy;
            const float ye = Y[idx(i,j-1)] - 0.7f * wy;
            Y[cur] = yo >= ye ? yo : ye;
            tY[cur] = yo >= ye ? 0 : 2;
        }
    }
    std::size_t i = n, j = m;
    std::uint8_t state = 0;
    if (X[idx(i,j)] > M[idx(i,j)]) state = 1;
    if (Y[idx(i,j)] > (state == 1 ? X[idx(i,j)] : M[idx(i,j)])) state = 2;
    std::vector<std::uint8_t> ops;
    ops.reserve(n+m);
    while (i || j) {
        const auto cur = idx(i,j);
        if (state == 0) {
            if (!i || !j) return false;
            ops.push_back(0); state = tM[cur]; --i; --j;
        } else if (state == 1) {
            if (!i) return false;
            ops.push_back(1); state = tX[cur]; --i;
        } else {
            if (!j) return false;
            ops.push_back(2); state = tY[cur]; --j;
        }
    }
    std::reverse(ops.begin(), ops.end());
    out = Profile{};
    out.ids = a.ids;
    out.ids.insert(out.ids.end(), b.ids.begin(), b.ids.end());
    out.rows.resize(out.ids.size());
    std::size_t ia = 0, ib = 0;
    for (const auto op : ops) {
        for (std::size_t r = 0; r < a.rows.size(); ++r)
            out.rows[r].push_back(op == 2 ? '-' : a.rows[r][ia]);
        for (std::size_t r = 0; r < b.rows.size(); ++r)
            out.rows[a.rows.size()+r].push_back(op == 1 ? '-' : b.rows[r][ib]);
        if (op != 2) ++ia;
        if (op != 1) ++ib;
    }
    return true;
}

// Length-independent 4-mer fingerprints: 256 possible DNA 4-mers.
using Fingerprint = std::array<std::uint64_t,4>;
Fingerprint sketch(const std::string& s) {
    Fingerprint f{};
    if (s.size() < 4) return f;
    for (std::size_t i = 0; i + 4 <= s.size(); ++i) {
        int k = 0; bool valid = true;
        for (std::size_t j = 0; j < 4; ++j) {
            const int c = code(s[i+j]);
            if (c >= 4) {valid = false; break;}
            k = k*4 + c;
        }
        if (valid) f[static_cast<std::size_t>(k)/64] |= std::uint64_t(1) << (k%64);
    }
    return f;
}

float sequence_distance(const Fingerprint& a, const Fingerprint& b,
                        std::size_t la, std::size_t lb) {
    int inter = 0, uni = 0;
    for (int k = 0; k < 4; ++k) {
        inter += __builtin_popcountll(a[k] & b[k]);
        uni += __builtin_popcountll(a[k] | b[k]);
    }
    const float jaccard = uni == 0 ? 0.5f : 1.0f - float(inter)/float(uni);
    const float lengthdiff = float(std::max(la,lb) - std::min(la,lb))/float(std::max(la,lb));
    return 0.85f*jaccard + 0.15f*lengthdiff;
}

// Baseline star algorithm, intentionally retained for bounded-resource fallback.
// This preserves the original supported interface and output contract.
struct Projection {
    std::string at_center;
    std::vector<std::string> insertions;
};
Projection fallback_pair(std::string_view center, std::string_view query) {
    detail::check_score_bounds(center.size(), query.size());
    const auto h = detail::checked_add(center.size(),1);
    const auto w = detail::checked_add(query.size(),1);
    const auto cells = detail::checked_multiply(h,w);
    (void)detail::checked_multiply(cells,sizeof(std::int64_t));
    if (cells > std::vector<std::int64_t>().max_size())
        throw detail::Error(detail::kCoreError,"pairwise matrix exceeds container size");
    std::vector<std::int64_t> d(cells);
    const auto at = [&](std::size_t i,std::size_t j)->std::int64_t&{return d[i*w+j];};
    for (std::size_t i=1;i<h;++i) at(i,0)=-static_cast<std::int64_t>(i);
    for (std::size_t j=1;j<w;++j) at(0,j)=-static_cast<std::int64_t>(j);
    const auto sub=[](char x,char y)->std::int64_t{
        if (x=='N'||y=='N') return 0;
        return x==y||(x=='T'&&y=='U')||(x=='U'&&y=='T')?1:-1;
    };
    for (std::size_t i=1;i<h;++i) for(std::size_t j=1;j<w;++j)
        at(i,j)=std::max({at(i-1,j-1)+sub(center[i-1],query[j-1]),at(i-1,j)-1,at(i,j-1)-1});
    std::string ac,aq;
    std::size_t i=center.size(),j=query.size();
    while(i||j){
        if(i&&j&&at(i,j)==at(i-1,j-1)+sub(center[i-1],query[j-1])){
            ac.push_back(center[--i]);aq.push_back(query[--j]);
        } else if(i&&at(i,j)==at(i-1,j)-1){ac.push_back(center[--i]);aq.push_back('-');}
        else if(j&&at(i,j)==at(i,j-1)-1){ac.push_back('-');aq.push_back(query[--j]);}
        else throw detail::Error(detail::kCoreError,"invalid traceback");
    }
    std::reverse(ac.begin(),ac.end());std::reverse(aq.begin(),aq.end());
    Projection p; p.at_center.assign(center.size(),'-');p.insertions.resize(h);
    std::size_t pos=0;
    for(std::size_t k=0;k<ac.size();++k){
        if(ac[k]=='-')p.insertions[pos].push_back(aq[k]);
        else p.at_center[pos++]=aq[k];
    }
    return p;
}
Alignment star_fallback(const std::vector<SequenceRecord>& input) {
    const auto& center=input.front().sequence;
    const std::size_t count=detail::checked_add(center.size(),1);
    std::vector<std::size_t> widths(count,0);
    std::vector<Projection> pp;
    for(std::size_t r=1;r<input.size();++r){
        auto p=fallback_pair(center,input[r].sequence);
        for(std::size_t s=0;s<count;++s) widths[s]=std::max(widths[s],p.insertions[s].size());
        pp.push_back(std::move(p));
    }
    std::size_t len=center.size();
    for(auto w:widths)len=detail::checked_add(len,w);
    (void)detail::checked_multiply(input.size(),len);
    Alignment out;out.rows.resize(input.size());
    for(std::size_t r=0;r<input.size();++r){
        auto& s=out.rows[r];s.reserve(len);
        const Projection* p=r?&pp[r-1]:nullptr;
        for(std::size_t k=0;k<count;++k){
            const std::size_t used=p?p->insertions[k].size():0;
            if(p)s+=p->insertions[k];
            s.append(widths[k]-used,'-');
            if(k<center.size())s.push_back(p?p->at_center[k]:center[k]);
        }
    }
    return out;
}

} // namespace

Alignment align(const std::vector<SequenceRecord>& input) {
    if(input.empty())throw detail::Error(detail::kCoreError,"the alignment core requires at least one record");
    for(const auto& r:input) if(r.sequence.empty()||
        !std::all_of(r.sequence.begin(),r.sequence.end(),is_residue))
        throw detail::Error(detail::kCoreError,"the alignment core requires nonempty cleaned ACGTUN sequences");
    if(input.size()==1)return Alignment{{input.front().sequence}};
    // For unusually large alignments use the original reference algorithm rather
    // than allocating unbounded profile DP matrices. No silent truncation.
    if(input.size()>256) {
        std::fprintf(stderr, "[msa-progressive] fallback=star reason=sequence_count n=%zu\n", input.size());
        return star_fallback(input);
    }
    std::size_t maxlen=0;
    for(const auto& r:input)maxlen=std::max(maxlen,r.sequence.size());
    if(maxlen>2500) {
        std::fprintf(stderr, "[msa-progressive] fallback=star reason=sequence_length maxlen=%zu\n", maxlen);
        return star_fallback(input);
    }

    const std::size_t n=input.size();
    std::vector<Fingerprint> fp(n);
    for(std::size_t i=0;i<n;++i)fp[i]=sketch(input[i].sequence);
    std::vector<std::vector<float>> dist(n,std::vector<float>(n,0));
    for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j)
        dist[i][j]=dist[j][i]=sequence_distance(fp[i],fp[j],input[i].sequence.size(),input[j].sequence.size());
    // Average linkage: use mean of original pairwise distances between clusters.
    std::vector<std::vector<std::size_t>> clusters(n);
    std::vector<Profile> profiles(n);
    for(std::size_t i=0;i<n;++i){clusters[i]={i};profiles[i].ids={i};profiles[i].rows={input[i].sequence};}
    while(clusters.size()>1){
        float best=std::numeric_limits<float>::infinity();
        std::size_t aa=0,bb=1;
        for(std::size_t i=0;i<clusters.size();++i) for(std::size_t j=i+1;j<clusters.size();++j){
            double sum=0;
            for(auto x:clusters[i])for(auto y:clusters[j])sum+=dist[x][y];
            const float d=static_cast<float>(sum/(clusters[i].size()*clusters[j].size()));
            if(d<best){best=d;aa=i;bb=j;}
        }
        Profile merged;
        if(!merge_profiles(profiles[aa],profiles[bb],merged)) {
            std::fprintf(stderr, "[msa-progressive] fallback=star reason=dp_budget n=%zu\n", n);
            return star_fallback(input);
        }
        clusters[aa].insert(clusters[aa].end(),clusters[bb].begin(),clusters[bb].end());
        profiles[aa]=std::move(merged);
        clusters.erase(clusters.begin()+bb);
        profiles.erase(profiles.begin()+bb);
    }
    std::fprintf(stderr, "[msa-progressive] method=guide_tree_profile_affine n=%zu columns=%zu\n",
                 n, profiles[0].rows.front().size());
    Alignment output;
    output.rows.resize(n);
    for(std::size_t k=0;k<n;++k)output.rows[profiles[0].ids[k]]=std::move(profiles[0].rows[k]);
    return output;
}
} // namespace msakit
