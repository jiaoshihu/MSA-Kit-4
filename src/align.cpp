#include "internal.hpp"
#include "msakit/core.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <string_view>
#include <vector>

// Experimental probabilistic star variant. Pair matches receive a soft
// confidence score, then a maximum-expected-accuracy path is decoded.
namespace msakit {
namespace {
struct Projection { std::string at_center; std::vector<std::string> insertions; };
bool valid(char c) { return c=='A'||c=='C'||c=='G'||c=='T'||c=='U'||c=='N'; }
bool same(char a,char b) { return a==b || (a=='T'&&b=='U') || (a=='U'&&b=='T'); }
std::array<double,64> kmers(std::string_view s){std::array<double,64> v{};if(s.size()<3)return v;
    for(std::size_t i=0;i+2<s.size();++i){auto cv=[](char c)->int{if(c=='A')return 0;if(c=='C')return 1;if(c=='G')return 2;if(c=='T'||c=='U')return 3;return -1;};
        int a=cv(s[i]),b=cv(s[i+1]),c=cv(s[i+2]);if(a>=0&&b>=0&&c>=0)v[a*16+b*4+c]++;}return v;}
std::size_t choose_center(const std::vector<SequenceRecord>& in){std::vector<std::array<double,64>> f;for(const auto& r:in)f.push_back(kmers(r.sequence));
    std::size_t best=0;double best_sum=-1;
    for(std::size_t i=0;i<in.size();++i){double sum=0;for(std::size_t j=0;j<in.size();++j)if(i!=j){double dot=0,a=0,b=0;for(std::size_t k=0;k<64;++k){dot+=f[i][k]*f[j][k];a+=f[i][k]*f[i][k];b+=f[j][k]*f[j][k];}
            sum+=(a==0||b==0)?0:dot/std::sqrt(a*b);}if(sum>best_sum){best=i;best_sum=sum;}}return best;}

Projection align_pair(std::string_view a, std::string_view b) {
    const auto h=detail::checked_add(a.size(),1), w=detail::checked_add(b.size(),1);
    const auto n=detail::checked_multiply(h,w);
    std::vector<double> dp(n,0.0);
    std::vector<unsigned char> tr(n,0);
    auto ix=[w](std::size_t i,std::size_t j){return i*w+j;};
    // Expected-accuracy style match confidence with mild local context.
    std::vector<double> p(detail::checked_multiply(a.size(),b.size()),0.0);
    for(std::size_t i=0;i<a.size();++i) for(std::size_t j=0;j<b.size();++j) {
        double q=(a[i]=='N'||b[j]=='N')?0.35:(same(a[i],b[j])?0.90:0.04);
        if(i>0&&j>0&&same(a[i-1],b[j-1])) q=std::min(0.98,q+0.06);
        if(i+1<a.size()&&j+1<b.size()&&same(a[i+1],b[j+1])) q=std::min(0.99,q+0.04);
        p[i*b.size()+j]=q;
    }
    for(std::size_t i=1;i<h;++i) for(std::size_t j=1;j<w;++j) {
        const double m=dp[ix(i-1,j-1)]+p[(i-1)*b.size()+j-1]-0.40;
        const double u=dp[ix(i-1,j)]-0.02, l=dp[ix(i,j-1)]-0.02;
        if(m>=u&&m>=l){dp[ix(i,j)]=m;tr[ix(i,j)]=0;}
        else if(u>=l){dp[ix(i,j)]=u;tr[ix(i,j)]=1;}
        else {dp[ix(i,j)]=l;tr[ix(i,j)]=2;}
    }
    std::string x,y; x.reserve(a.size()+b.size()); y.reserve(a.size()+b.size());
    std::size_t i=a.size(),j=b.size();
    while(i||j){auto d=tr[ix(i,j)];
        if(i&&j&&d==0){x.push_back(a[--i]);y.push_back(b[--j]);}
        else if(i&&(j==0||d==1)){x.push_back(a[--i]);y.push_back('-');}
        else {x.push_back('-');y.push_back(b[--j]);}
    }
    std::reverse(x.begin(),x.end());std::reverse(y.begin(),y.end());
    Projection out;out.at_center.assign(a.size(),'-');out.insertions.resize(a.size()+1);
    std::size_t used=0;
    for(std::size_t k=0;k<x.size();++k){if(x[k]=='-')out.insertions[used].push_back(y[k]);else out.at_center[used++]=y[k];}
    return out;
}

Alignment star(const std::vector<SequenceRecord>& in) {
    if(in.empty()) throw detail::Error(detail::kCoreError,"empty input");
    for(const auto& r:in) if(r.sequence.empty()||!std::all_of(r.sequence.begin(),r.sequence.end(),valid))
        throw detail::Error(detail::kCoreError,"expected nonempty ACGTUN sequences");
    if(in.size()==1)return Alignment{{in[0].sequence}};
    const auto center=choose_center(in);const auto& c=in[center].sequence; const auto slots=c.size()+1;
    std::vector<std::size_t> widths(slots,0);std::vector<Projection> ps(in.size());
    for(std::size_t r=0;r<in.size();++r)if(r!=center){ps[r]=align_pair(c,in[r].sequence);for(std::size_t s=0;s<slots;++s)widths[s]=std::max(widths[s],ps[r].insertions[s].size());}
    Alignment out;out.rows.resize(in.size());
    for(std::size_t r=0;r<in.size();++r){auto& row=out.rows[r];const Projection* p=r==center?nullptr:&ps[r];
        for(std::size_t s=0;s<slots;++s){std::size_t used=0;if(p){row+=p->insertions[s];used=p->insertions[s].size();}
            row.append(widths[s]-used,'-');if(s<c.size())row.push_back(p?p->at_center[s]:c[s]);}}
    return out;
}
} // namespace
Alignment align(const std::vector<SequenceRecord>& input){return star(input);}
} // namespace msakit
