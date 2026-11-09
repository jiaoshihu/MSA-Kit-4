#include "internal.hpp"
#include "msakit/core.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <string_view>
#include <vector>

// Experimental unbalanced-transport-inspired star variant. A Sinkhorn plan
// supplies soft residue correspondences; a monotone DP turns them into gaps.
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
    const std::size_t n=a.size(),m=b.size();
    const auto cells=detail::checked_multiply(n,m);
    std::vector<double> kernel(cells),plan(cells),u(n,1.0),v(m,1.0);
    constexpr double eps=0.22, tau=0.45;
    const double rho=tau/(tau+eps);
    auto ix=[m](std::size_t i,std::size_t j){return i*m+j;};
    for(std::size_t i=0;i<n;++i) for(std::size_t j=0;j<m;++j){
        const double base=(a[i]=='N'||b[j]=='N')?0.45:(same(a[i],b[j])?0.0:1.0);
        const double pos=std::abs((static_cast<double>(i)+0.5)/n-(static_cast<double>(j)+0.5)/m);
        kernel[ix(i,j)]=std::exp(-(base+0.12*pos)/eps);
    }
    // Unbalanced Sinkhorn scaling: unlike balanced transport, some residue
    // mass may remain unmatched, which gives gaps a natural interpretation.
    for(int it=0;it<80;++it){
        for(std::size_t i=0;i<n;++i){double z=0;for(std::size_t j=0;j<m;++j)z+=kernel[ix(i,j)]*v[j];
            u[i]=std::pow((1.0/static_cast<double>(n))/std::max(z,1e-300),rho);}
        for(std::size_t j=0;j<m;++j){double z=0;for(std::size_t i=0;i<n;++i)z+=kernel[ix(i,j)]*u[i];
            v[j]=std::pow((1.0/static_cast<double>(m))/std::max(z,1e-300),rho);}
    }
    std::vector<double> match(cells,0.0);
    const double scale=std::sqrt((1.0/static_cast<double>(n))*(1.0/static_cast<double>(m)));
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<m;++j)
        match[ix(i,j)]=std::min(1.0,u[i]*kernel[ix(i,j)]*v[j]/scale);

    const auto h=detail::checked_add(n,1),w=detail::checked_add(m,1),dp_cells=detail::checked_multiply(h,w);
    std::vector<double> dp(dp_cells,0.0);std::vector<unsigned char> tr(dp_cells,0);
    auto di=[w](std::size_t i,std::size_t j){return i*w+j;};
    for(std::size_t i=1;i<h;++i)for(std::size_t j=1;j<w;++j){
        const double d=dp[di(i-1,j-1)]+match[ix(i-1,j-1)]-0.30;
        const double up=dp[di(i-1,j)],left=dp[di(i,j-1)];
        if(d>=up&&d>=left){dp[di(i,j)]=d;tr[di(i,j)]=0;}
        else if(up>=left){dp[di(i,j)]=up;tr[di(i,j)]=1;}else{dp[di(i,j)]=left;tr[di(i,j)]=2;}
    }
    std::string x,y;x.reserve(n+m);y.reserve(n+m);std::size_t i=n,j=m;
    while(i||j){const auto d=tr[di(i,j)];if(i&&j&&d==0){x.push_back(a[--i]);y.push_back(b[--j]);}
        else if(i&&(j==0||d==1)){x.push_back(a[--i]);y.push_back('-');}
        else{x.push_back('-');y.push_back(b[--j]);}}
    std::reverse(x.begin(),x.end());std::reverse(y.begin(),y.end());
    Projection out;out.at_center.assign(n,'-');out.insertions.resize(n+1);std::size_t used=0;
    for(std::size_t k=0;k<x.size();++k){if(x[k]=='-')out.insertions[used].push_back(y[k]);else out.at_center[used++]=y[k];}
    return out;
}

Alignment star(const std::vector<SequenceRecord>& in) {
    if(in.empty())throw detail::Error(detail::kCoreError,"empty input");
    for(const auto& r:in)if(r.sequence.empty()||!std::all_of(r.sequence.begin(),r.sequence.end(),valid))
        throw detail::Error(detail::kCoreError,"expected nonempty ACGTUN sequences");
    if(in.size()==1)return Alignment{{in[0].sequence}};
    const auto center=choose_center(in);const auto& c=in[center].sequence;const auto slots=c.size()+1;
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
