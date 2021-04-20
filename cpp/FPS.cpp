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


// #PORT#
// name: "fps"
// prefix: "fps"
// description: "形式的べき級数"
template<int MOD> struct ModInt {
    static const int Mod = MOD; unsigned x; ModInt() : x(0) { }
    ModInt(signed sig) { x = sig < 0 ? sig % MOD + MOD : sig % MOD; }
    ModInt(signed long long sig) { x = sig < 0 ? sig % MOD + MOD : sig % MOD; }
    int get() const { return (int)x; }
    ModInt& operator+=(ModInt that) { if ((x += that.x) >= MOD) x -= MOD; return *this; }
    ModInt& operator-=(ModInt that) { if ((x += MOD - that.x) >= MOD) x -= MOD; return *this; }
    ModInt& operator*=(ModInt that) { x = (unsigned long long)x * that.x % MOD; return *this; }
    ModInt& operator/=(ModInt that) { return *this *= that.inverse(); }
    ModInt operator+(ModInt that) const { return ModInt(*this) += that; }
    ModInt operator-(ModInt that) const { return ModInt(*this) -= that; }
    ModInt operator*(ModInt that) const { return ModInt(*this) *= that; }
    ModInt operator/(ModInt that) const { return ModInt(*this) /= that; }
    ModInt inverse() const {
        long long a = x, b = MOD, u = 1, v = 0;
        while (b) { long long t = a / b; a -= t * b; std::swap(a, b); u -= t * v; std::swap(u, v); }
        return ModInt(u);
    }
    bool operator==(ModInt that) const { return x == that.x; }
    bool operator!=(ModInt that) const { return x != that.x; }
    ModInt operator-() const { ModInt t; t.x = x == 0 ? 0 : Mod - x; return t; }

    ModInt pow(ll n) { // 多分動く
        ModInt ret(1), mul((signed long long)x);
        while (n > 0) {
            if (n & 1) ret *= mul;
            mul *= mul;
            n >>= 1;
        }
        return ret;
    }
};
template<int MOD> ostream& operator<<(ostream& st, const ModInt<MOD> a) { st << a.get(); return st; };
template<int MOD> ModInt<MOD> operator^(ModInt<MOD> a, unsigned long long k) {
    ModInt<MOD> r = 1; while (k) { if (k & 1) r *= a; a *= a; k >>= 1; } return r;
}

template<typename T>
struct FormalPowerSeries {
    using Poly = vector<T>;
    using PolyMap = umap<ll,ll>;
    using Conv = function<Poly(Poly, Poly)>;
    Conv conv;
    FormalPowerSeries(Conv conv) :conv(conv) {}

    Poly pre(const Poly& as, int deg) {
        return Poly(as.begin(), as.begin() + min((int)as.size(), deg));
    }

    Poly add(Poly as, Poly bs) {
        int size_ = max(as.size(), bs.size());
        Poly cs(size_, T(0));
        for (int i = 0; i < (int)as.size(); i++) cs[i] += as[i];
        for (int i = 0; i < (int)bs.size(); i++) cs[i] += bs[i];
        return cs;
    }

    Poly sub(Poly as, Poly bs) {
        int size_ = max(as.size(), bs.size());
        Poly cs(size_, T(0));
        for (int i = 0; i < (int)as.size(); i++) cs[i] += as[i];
        for (int i = 0; i < (int)bs.size(); i++) cs[i] -= bs[i];
        return cs;
    }

    Poly mul_dence(Poly as, Poly bs) {
        return conv(as, bs);
    }

    Poly mul_sparse(Poly &as, PolyMap bs, int deg) {  // bs = [2:3] → 3x^2
        Poly ns(deg,(T)0);
        rep(i,(ll)as.size()) {
            for (auto [k,coeff] : bs) {
                if (i+k >= deg) continue;
                ns[i+k] += as[i] * coeff;
            }
        }
        return ns;
    }

    Poly mul_scalar(Poly as, T k) {
        for (auto& a : as) a *= k;
        return as;
    }

    // F(0) must not be 0
    Poly inv(Poly as, int deg) { // deg: どの次元まで係数を取得するか(恣意的に決めてOK)
        assert(as[0] != T(0));
        Poly rs({ T(1) / as[0] });
        for (int i = 1; i < deg; i <<= 1)
            rs = pre(sub(add(rs, rs), mul_dence(mul_dence(rs, rs), pre(as, i << 1))), i << 1);
        return rs;
    }

    // not zero
    Poly div(Poly as, Poly bs) {
        while (as.back() == T(0)) as.pop_back();
        while (bs.back() == T(0)) bs.pop_back();
        if (bs.size() > as.size()) return Poly();
        reverse(as.begin(), as.end());
        reverse(bs.begin(), bs.end());
        int need = as.size() - bs.size() + 1;
        Poly ds = pre(mul_dence(as, inv(bs, need)), need);
        reverse(ds.begin(), ds.end());
        return ds;
    }

    // F(0) must be 1
    Poly sqrt(Poly as, int deg) {
        assert(as[0] == T(1));
        T inv2 = T(1) / T(2);
        Poly ss({ T(1) });
        for (int i = 1; i < deg; i <<= 1) {
            ss = pre(add(ss, mul_dence(pre(as, i << 1), inv(ss, i << 1))), i << 1);
            for (T& x : ss) x *= inv2;
        }
        return ss;
    }

    Poly diff(Poly as) {
        int n = as.size();
        Poly res(n - 1);
        for (int i = 1; i < n; i++) res[i - 1] = as[i] * T(i);
        return res;
    }

    Poly integral(Poly as) {
        int n = as.size();
        Poly res(n + 1);
        res[0] = T(0);
        for (int i = 0; i < n; i++) res[i + 1] = as[i] / T(i + 1);
        return res;
    }

    // F(0) must be 1
    Poly log(Poly as, int deg) {
        return pre(integral(mul_dence(diff(as), inv(as, deg))), deg);
    }

    // F(0) must be 0
    Poly exp(Poly as, int deg) {
        Poly f({ T(1) });
        as[0] += T(1);
        for (int i = 1; i < deg; i <<= 1)
            f = pre(mul_dence(f, sub(pre(as, i << 1), log(f, i << 1))), i << 1);
        return f;
    }

    Poly partition(int n) {
        Poly rs(n + 1);
        rs[0] = T(1);
        for (int k = 1; k <= n; k++) {
            if (1LL * k * (3 * k + 1) / 2 <= n) rs[k * (3 * k + 1) / 2] += T(k % 2 ? -1LL : 1LL);
            if (1LL * k * (3 * k - 1) / 2 <= n) rs[k * (3 * k - 1) / 2] += T(k % 2 ? -1LL : 1LL);
        }
        return inv(rs, n + 1);
    }
};




#define FOR(i,n) for(int i = 0; i < (n); i++)
#define fore(i,a) for(auto &i:a)
#define ten(x) ((int)1e##x)
template<class T> T extgcd(T a, T b, T& x, T& y) { for (T u = y = 1, v = x = 0; a;) { T q = b / a; swap(x -= q * u, u); swap(y -= q * v, v); swap(b -= q * a, a); } return b; }
template<class T> T mod_inv(T a, T m) { T x, y; extgcd(a, m, x, y); return (m + x % m) % m; }
ll mod_pow(ll a, ll n, ll mod) { ll ret = 1; ll p = a % mod; while (n) { if (n & 1) ret = ret * p % mod; p = p * p % mod; n >>= 1; } return ret; }

template<ll _MOD_>
struct MathsNTTModAny {
    using mint = ModInt<_MOD_>;
    template<int mod, int primitive_root>
    class NTT {
    public:
        int get_mod() const { return mod; }
        void _ntt(vector<ll>& a, int sign) {
            const int n = sz(a);
            assert((n ^ (n & -n)) == 0); //n = 2^k

            const int g = 3; //g is primitive root of mod
            int h = (int)mod_pow(g, (mod - 1) / n, mod); // h^n = 1
            if (sign == -1) h = (int)mod_inv(h, mod); //h = h^-1 % mod

            //bit reverse
            int i = 0;
            for (int j = 1; j < n - 1; ++j) {
                for (int k = n >> 1; k > (i ^= k); k >>= 1);
                if (j < i) swap(a[i], a[j]);
            }

            for (int m = 1; m < n; m *= 2) {
                const int m2 = 2 * m;
                const ll base = mod_pow(h, n / m2, mod);
                ll w = 1;
                FOR(x, m) {
                    for (int s = x; s < n; s += m2) {
                        ll u = a[s];
                        ll d = a[s + m] * w % mod;
                        a[s] = u + d;
                        if (a[s] >= mod) a[s] -= mod;
                        a[s + m] = u - d;
                        if (a[s + m] < 0) a[s + m] += mod;
                    }
                    w = w * base % mod;
                }
            }

            for (auto& x : a) if (x < 0) x += mod;
        }
        void ntt(vector<ll>& input) {
            _ntt(input, 1);
        }
        void intt(vector<ll>& input) {
            _ntt(input, -1);
            const int n_inv = mod_inv((int)sz(input), mod);
            for (auto& x : input) x = x * n_inv % mod;
        }

        vector<ll> convolution(const vector<ll>& a, const vector<ll>& b) {
            int ntt_size = 1;
            while (ntt_size < sz(a) + sz(b)) ntt_size *= 2;

            vector<ll> _a = a, _b = b;
            _a.resize(ntt_size); _b.resize(ntt_size);

            ntt(_a);
            ntt(_b);

            FOR(i, ntt_size) {
                (_a[i] *= _b[i]) %= mod;
            }

            intt(_a);
            return _a;
        }
    };

    ll garner(vector<pair<int, int>> mr, int mod) {
        mr.emplace_back(mod, 0);

        vector<ll> coffs(sz(mr), 1);
        vector<ll> constants(sz(mr), 0);
        FOR(i, sz(mr) - 1) {
            // coffs[i] * v + constants[i] == mr[i].second (mod mr[i].first)
            ll v = (mr[i].second - constants[i]) * mod_inv<ll>(coffs[i], mr[i].first) % mr[i].first;
            if (v < 0) v += mr[i].first;

            for (int j = i + 1; j < sz(mr); j++) {
                (constants[j] += coffs[j] * v) %= mr[j].first;
                (coffs[j] *= mr[i].first) %= mr[j].first;
            }
        }

        return constants[sz(mr) - 1];
    }

    typedef NTT<167772161, 3> NTT_1;
    typedef NTT<469762049, 3> NTT_2;
    typedef NTT<1224736769, 3> NTT_3;

    vector<ll> solve(vector<ll> a, vector<ll> b, int mod = 1000000007) {
        for (auto& x : a) x %= mod;
        for (auto& x : b) x %= mod;

        NTT_1 ntt1; NTT_2 ntt2; NTT_3 ntt3;
        assert(ntt1.get_mod() < ntt2.get_mod() && ntt2.get_mod() < ntt3.get_mod());
        auto x = ntt1.convolution(a, b);
        auto y = ntt2.convolution(a, b);
        auto z = ntt3.convolution(a, b);

        const ll m1 = ntt1.get_mod(), m2 = ntt2.get_mod(), m3 = ntt3.get_mod();
        const ll m1_inv_m2 = mod_inv<ll>(m1, m2);
        const ll m12_inv_m3 = mod_inv<ll>(m1 * m2, m3);
        const ll m12_mod = m1 * m2 % mod;
        vector<ll> ret(sz(x));
        FOR(i, sz(x)) {
            ll v1 = (y[i] - x[i]) * m1_inv_m2 % m2;
            if (v1 < 0) v1 += m2;
            ll v2 = (z[i] - (x[i] + m1 * v1) % m3) * m12_inv_m3 % m3;
            if (v2 < 0) v2 += m3;
            ll constants3 = (x[i] + m1 * v1 + m12_mod * v2) % mod;
            if (constants3 < 0) constants3 += mod;
            ret[i] = constants3;
        }

        return ret;
    }

    vector<int> solve(vector<int> a, vector<int> b, int mod = 1000000007) {
        vector<ll> x(ALL(a));
        vector<ll> y(ALL(b));

        auto z = solve(x, y, mod);
        vector<int> res;
        fore(aa, z) res.push_back(aa % mod);

        return res;
    }

    vector<mint> solve(vector<mint> a, vector<mint> b, int mod = 1000000007) {
        int n = a.size();
        vector<ll> x(n);
        rep(i,n) x[i] = a[i].get();
        n = b.size();
        vector<ll> y(n);
        rep(i,n) y[i] = b[i].get();

        auto z = solve(x, y, mod);
        vector<int> res;
        fore(aa, z) res.push_back(aa % mod);

        vector<mint> res2;
        fore(x, res) res2.push_back(x);

        return res2;
    }
};


const ll MOD = 998244353;
using modint = ModInt<MOD>;

FormalPowerSeries<modint> get_fps() { // 係数情報のvectorはmodintを使うことに留意
    FormalPowerSeries<modint> fps([&](auto a, auto b) {
        MathsNTTModAny<MOD> ntt;
        return ntt.solve(a, b,MOD);
    });
    return fps;
}

// [x^N]P/Q : O(M * log N)
// http://q.c.titech.ac.jp/docs/progs/polynomial_division.html
modint bostan_mori(vector<modint> P, vector<modint> Q, ll N, MathsNTTModAny<MOD> ntt) {
    assert(Q[0] == 1);
    assert(Q.size() == P.size() + 1);

    while (N) {
        vector<modint> Q1;
        rep(i,Q.size()/2+1) {
            if (2 * i < sz(Q)) Q1.push_back(Q[2*i]); 
            if (2 * i + 1 < sz(Q)) Q1.push_back(- Q[2*i+1]);
        }
        
        ll offset = N & 1;
        auto nP = ntt.solve(P,Q1,MOD);
        auto nQ = ntt.solve(Q,Q1,MOD);

        vector<modint> nnP,nnQ;
        rep(i,nP.size()/2+1) if (2*i+offset < sz(nP)) nnP.push_back(nP[2*i+offset]);
        rep(i,nQ.size()/2+1) if (2*i < sz(nQ)) nnQ.push_back(nQ[2*i]);

        P = nnP, Q = nnQ;
        N >>= 1;
    }

    return P[0];
}

// a_n = a_n-1 * c_1 + a_n-2 * c_2 + ... + a_n-k * c_k のN項目を返す (N: 0-index)
modint get_Nth(const vector<modint> &A, const vector<modint> &C, ll N, MathsNTTModAny<MOD> ntt) {
    if (N < sz(A)) return A[N];
    assert(sz(A) == sz(C));

    vector<modint> Q(1,1);
    rep(i,sz(C)) Q.push_back(modint(-1) * C[i]);
    vector<modint> P = ntt.solve(A,Q,MOD);
    P.resize(sz(Q)-1);
    return bostan_mori(P,Q,N,ntt);
}


// #PORT_END#

int main() {
    #ifdef LOCAL
    ifstream in("../../Atcoder/input.txt");
    cin.rdbuf(in.rdbuf());
    #endif

    auto fps = get_fps();

    ll N,K;
    cin>>N>>K;

    vector<ll> L(K),R(K);
    ll max_r = 0;
    rep(i,K) {
        cin>>L[i]>>R[i];
        chmax(max_r,R[i]);
    }

    vector<modint> X(max_r+1,0);
    rep(i,K) {
        REP(j,L[i],R[i]+1) {
            X[j] = -1;
        }
    }
    X[0] = 1;

    auto h = fps.inv(X,N+1);
    cout << h[N-1] << endl;

    return 0;
}


