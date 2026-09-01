#include "mini_mpfr.h"
#include <gmp.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using u64 = uint64_t;
using u32 = uint32_t;
using u128 = unsigned __int128;

static constexpr int PREC = 72;
static constexpr int M = 68;
static constexpr u64 D_DEN = 1000000000000000000ULL;
static constexpr u64 GRID_DEN = 10922ULL;
static constexpr int RADIUS = 4;
static constexpr size_t N_INTERVALS = size_t(2 * RADIUS) * GRID_DEN; // 87376
static constexpr size_t N_NODES = N_INTERVALS + 1;
static constexpr int SCALE_BITS = 50;
static constexpr u64 SCALE = 1ULL << SCALE_BITS;
static constexpr int TARGET_NUM = 173639; // 1.73639
static constexpr int TARGET_DEN = 100000;

static const std::array<u64,M> dnum = {
93850787784284279ULL,93728425742790880ULL,76591171568741088ULL,76496413306980480ULL,
62712670930698392ULL,62634959976440448ULL,51405681933411200ULL,51341829408124328ULL,
42122911070025168ULL,42071884270997296ULL,34475249931136960ULL,34432415510271604ULL,
28129304910506608ULL,28095227669678320ULL,22894677772824260ULL,22867203104461428ULL,
18569496352328016ULL,18545717569651972ULL,15002733447710884ULL,14983110748912456ULL,
12067902179656702ULL,12051625579596242ULL,9657052573539880ULL,9643501383861628ULL,
7683082056180244ULL,7672247194880827ULL,6073025308358380ULL,6064358184491513ULL,
4764612193306462ULL,4757541831448853ULL,3706511845142744ULL,3700728378656680ULL,
2855481644095023ULL,2850839564382365ULL,2175555325750221ULL,2171878933904282ULL,
1636500792107330ULL,1633578463640486ULL,1213113874428726ULL,1210814871033224ULL,
884023439196695ULL,882236381121420ULL,631604803219683ULL,630218464129708ULL,
440884725651321ULL,439831511443754ULL,299372616118154ULL,298598842859409ULL,
196734122429080ULL,196180316028014ULL,124266197565460ULL,123881768868247ULL,
74739952623040ULL,74486263063135ULL,42278488554368ULL,42121786263331ULL,
22095830448095ULL,22006608368526ULL,10396698876656ULL,10350495041429ULL,
4226378023802ULL,4206165256487ULL,1378225296143ULL,1371348052322ULL,
310453649457ULL,308832577388ULL,32106640269ULL,31990196728ULL
};
static const std::array<u64,M> wnum = {
1000000000000ULL,1024544462166ULL,1052014951376ULL,1080168747263ULL,
1110137301047ULL,1143707518155ULL,1178078213459ULL,1213634093480ULL,
1252688915573ULL,1295160438244ULL,1339036607481ULL,1385620127714ULL,
1436052687332ULL,1492964075615ULL,1549920318624ULL,1613578126773ULL,
1677334426745ULL,1750903901054ULL,1828383297494ULL,1912695712185ULL,
2002651952490ULL,2098771620827ULL,2203446212683ULL,2317161226160ULL,
2441675533173ULL,2577227031199ULL,2724804092866ULL,2885234291304ULL,
3060226957132ULL,3251497724579ULL,3461769825965ULL,3693366097932ULL,
3949627244132ULL,4233165984122ULL,4547853449570ULL,4897384978629ULL,
5287191257637ULL,5722851719308ULL,6211974599272ULL,6762500981319ULL,
7385082597781ULL,8091612570413ULL,8898029323947ULL,9822356762870ULL,
10888377661604ULL,12123501810862ULL,13563810926102ULL,15252314819854ULL,
17248650068938ULL,19618807004163ULL,22478284042542ULL,25907011332280ULL,
30185385365731ULL,35307073038180ULL,42097111627793ULL,50004882353711ULL,
61675823973258ULL,74301644252826ULL,96605037642710ULL,117532151414543ULL,
166318527498026ULL,202653215957664ULL,330183940714618ULL,398817594724578ULL,
841972856358833ULL,1003088909831793ULL,3651740987208286ULL,4205535728448838ULL
};

static std::string to_string_u128(u128 x) {
  if (x == 0) return "0";
  std::string s;
  while (x) { s.push_back(char('0' + x % 10)); x /= 10; }
  std::reverse(s.begin(), s.end()); return s;
}

struct Real {
  mpfr_t v;
  Real() { mpfr_init2(v, PREC); mpfr_set_zero(v, 1); }
  Real(const Real& o) { mpfr_init2(v, PREC); mpfr_set(v,o.v,MPFR_RNDN); }
  Real& operator=(const Real& o) { if(this!=&o) mpfr_set(v,o.v,MPFR_RNDN); return *this; }
  Real(Real&& o) noexcept { mpfr_init2(v, PREC); mpfr_set(v,o.v,MPFR_RNDN); }
  Real& operator=(Real&& o) noexcept { if(this!=&o) mpfr_set(v,o.v,MPFR_RNDN); return *this; }
  ~Real(){ mpfr_clear(v); }
};
struct Interval { Real lo, hi; };

static Interval ratio_interval(u64 num, u64 den) {
  Interval x;
  mpfr_set_ui(x.lo.v,num,MPFR_RNDD); mpfr_div_ui(x.lo.v,x.lo.v,den,MPFR_RNDD);
  mpfr_set_ui(x.hi.v,num,MPFR_RNDU); mpfr_div_ui(x.hi.v,x.hi.v,den,MPFR_RNDU);
  return x;
}
static void clamp_nonnegative(Real& x) { if (mpfr_cmp_ui(x.v,0)<0) mpfr_set_zero(x.v,1); }
static u64 ceil_scaled(const Real& upper) {
  if (mpfr_cmp_ui(upper.v,0)<=0) return 0;
  if (mpfr_cmp_ui(upper.v,1)>=0) return SCALE;
  Real t;
  mpfr_mul_ui(t.v,upper.v,SCALE,MPFR_RNDU);
  u64 z=mpfr_get_ui(t.v,MPFR_RNDU);
  return std::min<u64>(z,SCALE);
}

struct GaussianData {
  std::vector<u64> Apos, Bpos; // indices 1..n
  std::vector<u64> tailUp, cdfUp; // k=0..n
};

static Interval sqrt2pi_interval() {
  Interval pi, twoPi, out;
  mpfr_const_pi(pi.lo.v,MPFR_RNDD); mpfr_const_pi(pi.hi.v,MPFR_RNDU);
  mpfr_mul_ui(twoPi.lo.v,pi.lo.v,2,MPFR_RNDD);
  mpfr_mul_ui(twoPi.hi.v,pi.hi.v,2,MPFR_RNDU);
  mpfr_sqrt(out.lo.v,twoPi.lo.v,MPFR_RNDD);
  mpfr_sqrt(out.hi.v,twoPi.hi.v,MPFR_RNDU);
  return out;
}

static void tail_density_at_k(size_t k, const Interval& sqrtD, const Interval& sqrt2pi,
                              Interval& tail, Interval& dens) {
  if (k==0) {
    mpfr_set_ui(tail.lo.v,1,MPFR_RNDN); mpfr_div_ui(tail.lo.v,tail.lo.v,2,MPFR_RNDN);
    mpfr_set(tail.hi.v,tail.lo.v,MPFR_RNDN);
    Real denLo, denHi;
    mpfr_mul(denLo.v,sqrt2pi.lo.v,sqrtD.lo.v,MPFR_RNDD);
    mpfr_mul(denHi.v,sqrt2pi.hi.v,sqrtD.hi.v,MPFR_RNDU);
    mpfr_ui_div(dens.lo.v,1,denHi.v,MPFR_RNDD);
    mpfr_ui_div(dens.hi.v,1,denLo.v,MPFR_RNDU);
    return;
  }
  // x=k/GRID_DEN; z=x/sqrt(d)
  Real xLo,xHi,zLo,zHi;
  mpfr_set_ui(xLo.v,k,MPFR_RNDD); mpfr_div_ui(xLo.v,xLo.v,GRID_DEN,MPFR_RNDD);
  mpfr_set_ui(xHi.v,k,MPFR_RNDU); mpfr_div_ui(xHi.v,xHi.v,GRID_DEN,MPFR_RNDU);
  mpfr_div(zLo.v,xLo.v,sqrtD.hi.v,MPFR_RNDD);
  mpfr_div(zHi.v,xHi.v,sqrtD.lo.v,MPFR_RNDU);
  // y=z/sqrt(2)
  Real two, sqrt2Lo,sqrt2Hi,yLo,yHi;
  mpfr_set_ui(two.v,2,MPFR_RNDN);
  mpfr_sqrt(sqrt2Lo.v,two.v,MPFR_RNDD); mpfr_sqrt(sqrt2Hi.v,two.v,MPFR_RNDU);
  mpfr_div(yLo.v,zLo.v,sqrt2Hi.v,MPFR_RNDD);
  mpfr_div(yHi.v,zHi.v,sqrt2Lo.v,MPFR_RNDU);
  // tail = erfc(y)/2; erfc decreasing
  mpfr_erfc(tail.lo.v,yHi.v,MPFR_RNDD); mpfr_div_ui(tail.lo.v,tail.lo.v,2,MPFR_RNDD);
  mpfr_erfc(tail.hi.v,yLo.v,MPFR_RNDU); mpfr_div_ui(tail.hi.v,tail.hi.v,2,MPFR_RNDU);
  // exp(-z^2/2)/(sqrt(2pi)*sqrt(d))
  Real qLo,qHi,eLo,eHi,denLo,denHi;
  mpfr_mul(qHi.v,zHi.v,zHi.v,MPFR_RNDU); mpfr_div_ui(qHi.v,qHi.v,2,MPFR_RNDU); mpfr_neg(qHi.v,qHi.v,MPFR_RNDD);
  mpfr_exp(eLo.v,qHi.v,MPFR_RNDD);
  mpfr_mul(qLo.v,zLo.v,zLo.v,MPFR_RNDD); mpfr_div_ui(qLo.v,qLo.v,2,MPFR_RNDD); mpfr_neg(qLo.v,qLo.v,MPFR_RNDU);
  mpfr_exp(eHi.v,qLo.v,MPFR_RNDU);
  mpfr_mul(denLo.v,sqrt2pi.lo.v,sqrtD.lo.v,MPFR_RNDD);
  mpfr_mul(denHi.v,sqrt2pi.hi.v,sqrtD.hi.v,MPFR_RNDU);
  mpfr_div(dens.lo.v,eLo.v,denHi.v,MPFR_RNDD);
  mpfr_div(dens.hi.v,eHi.v,denLo.v,MPFR_RNDU);
}

static GaussianData build_gaussian_data(int stage) {
  GaussianData gd;
  gd.Apos.assign(N_INTERVALS+1,0); gd.Bpos.assign(N_INTERVALS+1,0);
  gd.tailUp.assign(N_INTERVALS+1,0); gd.cdfUp.assign(N_INTERVALS+1,0);
  Interval d=ratio_interval(dnum[stage],D_DEN), sqrtD;
  mpfr_sqrt(sqrtD.lo.v,d.lo.v,MPFR_RNDD); mpfr_sqrt(sqrtD.hi.v,d.hi.v,MPFR_RNDU);
  Interval s2p=sqrt2pi_interval();
  Interval tPrev,dPrev,tCur,dCur;
  tail_density_at_k(0,sqrtD,s2p,tPrev,dPrev);
  gd.tailUp[0]=ceil_scaled(tPrev.hi);
  Real cdf0; mpfr_ui_div(cdf0.v,1,tPrev.lo.v,MPFR_RNDU); // overwritten next line; avoid special op
  mpfr_set_ui(cdf0.v,1,MPFR_RNDU); mpfr_sub(cdf0.v,cdf0.v,tPrev.lo.v,MPFR_RNDU);
  gd.cdfUp[0]=ceil_scaled(cdf0);
  // factor d/h = d*GRID_DEN
  Real facLo,facHi;
  mpfr_mul_ui(facLo.v,d.lo.v,GRID_DEN,MPFR_RNDD);
  mpfr_mul_ui(facHi.v,d.hi.v,GRID_DEN,MPFR_RNDU);
  for(size_t k=1;k<=N_INTERVALS;k++) {
    tail_density_at_k(k,sqrtD,s2p,tCur,dCur);
    gd.tailUp[k]=ceil_scaled(tCur.hi);
    Real cdfU; mpfr_set_ui(cdfU.v,1,MPFR_RNDU); mpfr_sub(cdfU.v,cdfU.v,tCur.lo.v,MPFR_RNDU);
    gd.cdfUp[k]=ceil_scaled(cdfU);
    Real qLo,qHi,ddLo,ddHi,mLo,mHi,tmp,Ahi,Bhi;
    mpfr_sub(qLo.v,tPrev.lo.v,tCur.hi.v,MPFR_RNDD); clamp_nonnegative(qLo);
    mpfr_sub(qHi.v,tPrev.hi.v,tCur.lo.v,MPFR_RNDU); clamp_nonnegative(qHi);
    mpfr_sub(ddLo.v,dPrev.lo.v,dCur.hi.v,MPFR_RNDD); clamp_nonnegative(ddLo);
    mpfr_sub(ddHi.v,dPrev.hi.v,dCur.lo.v,MPFR_RNDU); clamp_nonnegative(ddHi);
    mpfr_mul(mLo.v,facLo.v,ddLo.v,MPFR_RNDD);
    mpfr_mul(mHi.v,facHi.v,ddHi.v,MPFR_RNDU);
    mpfr_mul_ui(tmp.v,qHi.v,k,MPFR_RNDU); mpfr_sub(Ahi.v,tmp.v,mLo.v,MPFR_RNDU);
    if(mpfr_cmp_ui(Ahi.v,0)<0) throw std::runtime_error("negative A upper");
    if(k==1) mpfr_set(Bhi.v,mHi.v,MPFR_RNDU);
    else { mpfr_mul_ui(tmp.v,qLo.v,k-1,MPFR_RNDD); mpfr_sub(Bhi.v,mHi.v,tmp.v,MPFR_RNDU); }
    if(mpfr_cmp_ui(Bhi.v,0)<0) throw std::runtime_error("negative B upper");
    gd.Apos[k]=ceil_scaled(Ahi); gd.Bpos[k]=ceil_scaled(Bhi);
    tPrev=tCur; dPrev=dCur;
  }
  return gd;
}

static u64 modpow(u64 a,u64 e,u64 mod){u64 r=1;while(e){if(e&1)r=(u128)r*a%mod;a=(u128)a*a%mod;e>>=1;}return r;}
struct Prime { u32 mod, root; };
static constexpr std::array<Prime,4> PRIMES = {{{998244353u,3u},{1004535809u,3u},{469762049u,3u},{1224736769u,3u}}};

static void ntt(std::vector<u32>& a, bool invert, u32 mod, u32 root) {
  size_t n=a.size();
  for(size_t i=1,j=0;i<n;i++) { size_t bit=n>>1; for(;j&bit;bit>>=1)j^=bit; j^=bit; if(i<j)std::swap(a[i],a[j]); }
  for(size_t len=2;len<=n;len<<=1) {
    u64 wlen=modpow(root,(mod-1)/len,mod); if(invert)wlen=modpow(wlen,mod-2,mod);
    for(size_t i=0;i<n;i+=len) { u64 w=1; size_t half=len>>1;
      for(size_t j=0;j<half;j++) { u32 u=a[i+j]; u32 v=(u128)a[i+j+half]*w%mod; u64 x=u+u64(v); if(x>=mod)x-=mod; a[i+j]=u32(x); a[i+j+half]=(u>=v?u-v:u+mod-v); w=(u128)w*wlen%mod; }
    }
  }
  if(invert){u64 invn=modpow(n,mod-2,mod);for(auto&x:a)x=(u128)x*invn%mod;}
}

static u128 crt4(const std::array<u32,4>& r) {
  u128 x=r[0], mult=PRIMES[0].mod;
  for(int i=1;i<4;i++) { u64 p=PRIMES[i].mod; u64 xm=x%p, mm=mult%p; u64 diff=(r[i]+p-xm)%p; u64 t=(u128)diff*modpow(mm,p-2,p)%p; x += mult*t; mult *= p; }
  return x;
}

static std::vector<u128> exact_pair_convolution(const std::vector<u64>& L,const std::vector<u64>& KA,
                                                 const std::vector<u64>& R,const std::vector<u64>& KB) {
  size_t need=L.size()+KA.size()-1, z=1; while(z<need)z<<=1;
  std::array<std::vector<u32>,4> residues; for(auto&v:residues)v.resize(N_NODES);
  #pragma omp parallel for schedule(static)
  for(int pi=0;pi<4;pi++) {
    u32 mod=PRIMES[pi].mod, root=PRIMES[pi].root;
    std::vector<u32> fL(z),fA(z),fR(z),fB(z);
    for(size_t i=0;i<L.size();i++){fL[i]=L[i]%mod;fR[i]=R[i]%mod;}
    for(size_t i=0;i<KA.size();i++){fA[i]=KA[i]%mod;fB[i]=KB[i]%mod;}
    ntt(fL,false,mod,root); ntt(fA,false,mod,root); ntt(fR,false,mod,root); ntt(fB,false,mod,root);
    for(size_t i=0;i<z;i++) fL[i]=((u128)fL[i]*fA[i]+(u128)fR[i]*fB[i])%mod;
    ntt(fL,true,mod,root);
    for(size_t j=0;j<N_NODES;j++) residues[pi][j]=fL[j+N_INTERVALS-1];
  }
  std::vector<u128> out(N_NODES); for(size_t j=0;j<N_NODES;j++){std::array<u32,4> r;for(int p=0;p<4;p++)r[p]=residues[p][j];out[j]=crt4(r);} return out;
}

struct PowEvaluator {
  Real rlo, x, logu, expo, y;
  explicit PowEvaluator(int stage) {
    mpfr_set_ui(rlo.v,wnum[stage+1],MPFR_RNDD);
    mpfr_div_ui(rlo.v,rlo.v,wnum[stage],MPFR_RNDD);
  }
  u64 eval(u64 xint) {
    if(xint==0)return 0;
    if(xint>=SCALE)return SCALE;
    mpfr_set_ui(x.v,xint,MPFR_RNDN); mpfr_div_ui(x.v,x.v,SCALE,MPFR_RNDN);
    // x in (0,1), rlo <= r. Thus x^r <= exp(rlo*log x).
    mpfr_log(logu.v,x.v,MPFR_RNDU);
    mpfr_mul(expo.v,rlo.v,logu.v,MPFR_RNDU);
    mpfr_exp(y.v,expo.v,MPFR_RNDU);
    return ceil_scaled(y);
  }
};

static u64 curvature_eps_upper(int stage) {
  // eps = h^2/(8 sqrt(2*pi*e) d)
  Interval pi; mpfr_const_pi(pi.lo.v,MPFR_RNDD); mpfr_const_pi(pi.hi.v,MPFR_RNDU);
  Real one,eLo,eHi,cLo,cHi,tmp,den,eps;
  mpfr_set_ui(one.v,1,MPFR_RNDN); mpfr_exp(eLo.v,one.v,MPFR_RNDD); mpfr_exp(eHi.v,one.v,MPFR_RNDU);
  mpfr_mul(tmp.v,pi.lo.v,eLo.v,MPFR_RNDD); mpfr_mul_ui(tmp.v,tmp.v,2,MPFR_RNDD); mpfr_sqrt(cLo.v,tmp.v,MPFR_RNDD);
  Interval d=ratio_interval(dnum[stage],D_DEN);
  mpfr_mul(den.v,cLo.v,d.lo.v,MPFR_RNDD); mpfr_mul_ui(den.v,den.v,8,MPFR_RNDD);
  mpfr_mul_ui(den.v,den.v,GRID_DEN,MPFR_RNDD); mpfr_mul_ui(den.v,den.v,GRID_DEN,MPFR_RNDD);
  mpfr_ui_div(eps.v,1,den.v,MPFR_RNDU);
  return ceil_scaled(eps);
}

static std::vector<u32> test_convolution_mod(const std::vector<u64>& a,
                                                  const std::vector<u64>& b,
                                                  Prime prime) {
  const size_t need = a.size() + b.size() - 1;
  size_t z = 1;
  while (z < need) z <<= 1;
  if ((prime.mod - 1) % z != 0) throw std::runtime_error("self-test transform length unsupported");
  std::vector<u32> fa(z), fb(z);
  for (size_t i=0;i<a.size();i++) fa[i]=a[i]%prime.mod;
  for (size_t i=0;i<b.size();i++) fb[i]=b[i]%prime.mod;
  ntt(fa,false,prime.mod,prime.root); ntt(fb,false,prime.mod,prime.root);
  for (size_t i=0;i<z;i++) fa[i]=(u128)fa[i]*fb[i]%prime.mod;
  ntt(fa,true,prime.mod,prime.root); fa.resize(need); return fa;
}

static void self_test_ntt(){
  const std::vector<u64> a={1,2,3}, b={5,7,11,13};
  std::vector<u128> expected(a.size()+b.size()-1,0);
  for(size_t i=0;i<a.size();i++)for(size_t j=0;j<b.size();j++)expected[i+j]+=u128(a[i])*b[j];
  size_t production_need=N_INTERVALS+2*N_INTERVALS-1, production_z=1;
  while(production_z<production_need)production_z<<=1;
  for(auto prime:PRIMES){
    if((prime.mod-1)%production_z!=0)throw std::runtime_error("production NTT length unsupported");
    u64 omega=modpow(prime.root,(prime.mod-1)/production_z,prime.mod);
    if(modpow(omega,production_z,prime.mod)!=1 || (production_z>1 && modpow(omega,production_z/2,prime.mod)==1))
      throw std::runtime_error("production root has wrong order");
    auto got=test_convolution_mod(a,b,prime);
    for(size_t i=0;i<expected.size();i++)if(got[i]!=u32(expected[i]%prime.mod))
      throw std::runtime_error("NTT convolution self-test failed");
  }
  const u128 x=(u128(1)<<100)+123456789u; std::array<u32,4> residues{};
  for(int i=0;i<4;i++)residues[i]=u32(x%PRIMES[i].mod);
  if(crt4(residues)!=x)throw std::runtime_error("CRT reconstruction self-test failed");
  std::cout<<"NTT/CRT self-test: PASS\n";
}

int main(int argc,char** argv){
  std::cout.setf(std::ios::unitbuf);
  try {
    bool self_test_only=false;
    if(argc==2 && std::string(argv[1])=="--self-test-only")self_test_only=true;
    else if(argc!=1)throw std::runtime_error("usage: verifier [--self-test-only]");
    std::cout<<"MPFR "<<mpfr_get_version()<<"; grid n="<<N_INTERVALS<<"; scale=2^"<<SCALE_BITS<<"\n";
    u128 sumd=0; for(auto x:dnum)sumd+=x; if(sumd!=D_DEN)throw std::runtime_error("d sum mismatch");
    for(int i=0;i<M-1;i++)if(!(wnum[i+1]>wnum[i]))throw std::runtime_error("weights not increasing");
    u128 crtmod=1;for(auto p:PRIMES)crtmod*=p.mod;
    u128 worst=u128(2)*N_INTERVALS*SCALE*SCALE;
    std::cout<<"CRT modulus="<<to_string_u128(crtmod)<<"\n";
    if(!(worst<crtmod))throw std::runtime_error("CRT modulus too small");
    self_test_ntt();
    if(self_test_only){std::cout<<"SELF-TEST ONLY: PASS\n";return 0;}
    auto t0=std::chrono::steady_clock::now();
    GaussianData gLast=build_gaussian_data(M-1);
    std::vector<u64> U(N_NODES);
    size_t center=N_INTERVALS/2;
    for(size_t j=0;j<N_NODES;j++) { long k=long(j)-long(center); U[j]=(k>=0?gLast.cdfUp[size_t(k)]:gLast.tailUp[size_t(-k)]); }
    std::cout<<"terminal initialized; U(0)/S="<<std::setprecision(17)<<double(U[center])/double(SCALE)<<"\n";
    for(int stage=M-2;stage>=0;stage--) {
      auto st=std::chrono::steady_clock::now();
      u64 eps=curvature_eps_upper(stage+1);
      std::vector<u64> Aval(N_NODES), Pnode(N_NODES), L(N_INTERVALS),R(N_INTERVALS);
      #pragma omp parallel
      {
        PowEvaluator pe_local(stage);
        #pragma omp for schedule(static)
        for(long long jj=0;jj<(long long)N_NODES;jj++) {
          size_t j=size_t(jj);
          u128 aa=u128(U[j])+eps;
          Aval[j]=(aa>=SCALE?SCALE:u64(aa));
          Pnode[j]=(aa>=SCALE?SCALE:pe_local.eval(u64(aa)));
        }
      }
      #pragma omp parallel for schedule(static)
      for(long long kk=0;kk<(long long)N_INTERVALS;kk++) {
        size_t k=size_t(kk);
        if(Aval[k]<SCALE && Aval[k+1]<SCALE) { L[k]=Pnode[k]; R[k]=Pnode[k+1]; }
        else L[k]=R[k]=SCALE;
      }
      PowEvaluator pe_left(stage);
      u64 left=pe_left.eval(U[0]);
      GaussianData gd=build_gaussian_data(stage);
      std::vector<u64>KA(2*N_INTERVALS),KB(2*N_INTERVALS);
      for(size_t idx=0;idx<2*N_INTERVALS;idx++) {
        long dd=long(N_INTERVALS)-long(idx);
        if(dd>=1){KA[idx]=gd.Apos[size_t(dd)];KB[idx]=gd.Bpos[size_t(dd)];}
        else {size_t kk=size_t(1-dd);KA[idx]=gd.Bpos[kk];KB[idx]=gd.Apos[kk];}
      }
      auto C=exact_pair_convolution(L,KA,R,KB);
      std::vector<u64>V(N_NODES);
      for(size_t j=0;j<N_NODES;j++) {
        u128 num=C[j]+u128(left)*gd.tailUp[j]+u128(SCALE)*gd.tailUp[N_INTERVALS-j];
        u128 q=(num+SCALE-1)/SCALE; V[j]=(q>=SCALE?SCALE:u64(q));
      }
      U.swap(V);
      auto en=std::chrono::steady_clock::now();
      std::cout<<"stage "<<stage<<" eps="<<eps<<" center="<<std::setprecision(17)<<double(U[center])/double(SCALE)
               <<" sec="<<std::chrono::duration<double>(en-st).count()<<"\n";
    }
    // exact normalized a0 = w0*D_DEN / sum_i dnum_i*wnum_i
    mpz_t sumW,numA; mpz_init(sumW);mpz_init(numA);mpz_set_ui(sumW,0);
    mpz_t tmpz;mpz_init(tmpz);
    for(int i=0;i<M;i++){mpz_set_ui(tmpz,dnum[i]);mpz_mul_ui(tmpz,tmpz,wnum[i]);mpz_add(sumW,sumW,tmpz);}
    mpz_set_ui(numA,wnum[0]);mpz_mul_ui(numA,numA,D_DEN);
    Real a0lo,x,logup,nloglo,Llo,targetHi;
    mpfr_set_z(a0lo.v,numA,MPFR_RNDD); Real denz; mpfr_set_z(denz.v,sumW,MPFR_RNDU); mpfr_div(a0lo.v,a0lo.v,denz.v,MPFR_RNDD);
    mpfr_set_ui(x.v,U[center],MPFR_RNDN);mpfr_div_ui(x.v,x.v,SCALE,MPFR_RNDN);
    mpfr_log(logup.v,x.v,MPFR_RNDU);mpfr_neg(nloglo.v,logup.v,MPFR_RNDD);
    mpfr_mul(Llo.v,a0lo.v,nloglo.v,MPFR_RNDD);mpfr_mul_ui(Llo.v,Llo.v,2,MPFR_RNDD);
    mpfr_set_ui(targetHi.v,TARGET_NUM,MPFR_RNDU);mpfr_div_ui(targetHi.v,targetHi.v,TARGET_DEN,MPFR_RNDU);
    double Ld=mpfr_get_d(Llo.v,MPFR_RNDD);
    std::cout<<std::setprecision(17)<<"H0 upper="<<double(U[center])/double(SCALE)<<"\n";
    std::cout<<"certified lower endpoint L_down="<<Ld<<"\n";
    bool ok=mpfr_cmp(Llo.v,targetHi.v)>0;
    std::cout<<(ok?"PASS":"FAIL")<<": E_star > "<<TARGET_NUM<<"/"<<TARGET_DEN<<"\n";
    mpz_clear(sumW);mpz_clear(numA);mpz_clear(tmpz);
    auto t1=std::chrono::steady_clock::now();std::cout<<"total sec="<<std::chrono::duration<double>(t1-t0).count()<<"\n";
    return ok?0:2;
  } catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 1;}
}
