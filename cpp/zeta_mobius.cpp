#include <bits/stdc++.h>
using namespace std;
#define LOCAL // 提出時はコメントアウト
#define DEBUG_

typedef long long ll;
const double EPS = 1e-9;
const ll INF = ((1LL<<62)-(1LL<<31));
typedef vector<ll> vecl;
typedef pair<ll, ll> pairl;
template<typename T> using uset = unordered_set<T>;
template<typename T, typename U> using mapv = map<T,vector<U>>;
template<typename T, typename U> using umap = unordered_map<T,U>;

#define ALL(v) v.begin(), v.end()
#define REP(i, x, n) for(int i = x; i < n; i++)
#define rep(i, n) REP(i, 0, n)
#define sz(x) (ll)x.size()
ll llceil(ll a,ll b) { return (a+b-1)/b; }
template<class T> inline bool chmax(T& a, T b) { if (a < b) { a = b; return true; } return false; }
template<class T> inline bool chmin(T& a, T b) { if (a > b) { a = b; return true; } return false; }
template<class T> vector<vector<T>> genarr(ll n, ll m, T init) { return vector<vector<T>>(n,vector<T>(m,init)); }

///// DEBUG
#define DUMPOUT cerr
#define repi(itr, ds) for (auto itr = ds.begin(); itr != ds.end(); itr++)
template<typename T>istream&operator>>(istream&is,vector<T>&vec){for(T&x:vec)is>>x;return is;}
template<typename T,typename U>ostream&operator<<(ostream&os,pair<T,U>&pair_var){os<<"("<<pair_var.first<<", "<<pair_var.second<<")";return os;}
template<typename T>ostream&operator<<(ostream&os,const vector<T>&vec){os<<"{";for(int i=0;i<vec.size();i++){os<<vec[i]<<(i+1==vec.size()?"":", ");}
os<<"}";return os;}
template<typename T,typename U>ostream&operator<<(ostream&os,map<T,U>&map_var){os<<"{";repi(itr,map_var){os<<*itr;itr++;if(itr!=map_var.end())os<<", ";itr--;}
os<<"}";return os;}
template<typename T>ostream&operator<<(ostream&os,set<T>&set_var){os<<"{";repi(itr,set_var){os<<*itr;itr++;if(itr!=set_var.end())os<<", ";itr--;}
os<<"}";return os;}
void dump_func(){DUMPOUT<<endl;}
template<class Head,class...Tail>void dump_func(Head&&head,Tail&&...tail){DUMPOUT<<head;if(sizeof...(Tail)>0){DUMPOUT<<", ";}
dump_func(std::move(tail)...);}
#ifndef LOCAL
#undef DEBUG_
#endif
#ifdef DEBUG_
#define DEB
#define dump(...)                                                          \
DUMPOUT << "  " << string(#__VA_ARGS__) << ": "                            \
        << "[" << to_string(__LINE__) << ":" << __FUNCTION__ << "]"        \
        << endl                                                            \
        << "    ",                                                         \
    dump_func(__VA_ARGS__)
#else
#define DEB if (false)
#define dump(...)
#endif

//////////

#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")

// #PORT#
// name: "zeta"
// prefix: "zeta"
// description: "倍数・約数ゼータ変換 / メビウス変換"

template <class T = ll>
struct multiple_transform { // フーリエと同じ雰囲気で変換. 畳み込みは各点積へと変わる.
    static void zeta(vector<T> &f) { // 倍数ゼータ: f -> F (F[x] = Σ[x|d] f(d))
        vecl sieve(sz(f),false);
        REP(p,2,sz(f)) {
            if (sieve[p]) continue;
            for (int i = (sz(f) - 1) / p; i >= 0; i--) sieve[i*p] = true, f[i] += f[i*p];
        }
    }

    static void prod(vector<T> &f) { // 各点積
        rep(i,sz(f)) f[i] = f[i] * f[i];
    }

    static void mobius(vector<T> &f) { // 倍数メビウス F -> f
        vecl sieve(sz(f),false);
        REP(p,2,sz(f)) {
            if (sieve[p]) continue;
            for (int i = 0; i*p<sz(f); i++) sieve[i*p] = true, f[i] -= f[i*p];
        }
    }

    static void gcd_conv(vector<T> &f) { // gcd畳み込み: F[g] = Σ[gcd(x,y)=g]f(x)f(y)
        // f -> F, g -> G
        // H = F * G (各点積) = Σ[x|d1] f(d1) Σ[x|d2] f(d2) = Σ[x|d] Σ[gcd(d1,d2)=d] f(d)
        // Hを逆変換すればいいのでメビウス変換すれば h[g] = Σ[gcd(d1,d2)=g] f(g)
        zeta(f);
        prod(f);
        mobius(f);
    }
};

template <class T = ll>
struct divisor_transform { // フーリエと同じ雰囲気で変換. 畳み込みは各点積へと変わる.
    static void zeta(vector<T> &f) { // 約数ゼータ: f -> F (F[x] = Σ[d|x] f(d))
        vecl sieve(sz(f),false);
        REP(p,2,sz(f)) {
            if (sieve[p]) continue;
            for (int i = 0; i*p < sz(f); i++) sieve[i*p] = true, f[i*p] += f[i];
        }
    }

    static void prod(vector<T> &f) { // 各点積
        rep(i,sz(f)) f[i] = f[i] * f[i];
    }

    static void mobius(vector<T> &f) { // 約数メビウス F -> f
        vecl sieve(sz(f),false);
        REP(p,2,sz(f)) {
            if (sieve[p]) continue;
            for (int i = (sz(f)-1) / p; i >= 0; i--) sieve[i*p] = true, f[i*p] -= f[i];
        }
    }

    static void lcm_conv(vector<T> &f) { // lcm畳み込み: F[g] = Σ[lcm(x,y)=g]f(x)f(y)
        // f -> F, g -> G
        // H = F * G (各点積) = Σ[d1|x] f(d1) Σ[d2|x] f(d2) = Σ[d|x] Σ[lcm(d1,d2)=d] f(d)
        // Hを逆変換すればいいのでメビウス変換すれば h[g] = Σ[lcm(d1,d2)=g] f(g)
        zeta(f);
        prod(f);
        mobius(f);
    }
};



// #PORT_END#