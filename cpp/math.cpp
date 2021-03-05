#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
#define REP(i, x, n) for (int i = x; i < n; i++)
#define rep(i, n) REP(i, 0, n)

// #PORT#
// name: "rational"
// prefix: "rational"
// description: "有理数を表すライブラリ"
struct Rational { 
    ll q,p; // = q / p

    Rational() : p(1), q(0) {}
    Rational(ll q, ll p) : p(p), q(q) {
        assert(p != 0);
        shape();
    }

    void shape() {
        if (q == 0) return;

        ll g = gcd(p,q);
        p /= g;
        q /= g;

        // (q,p): -はqの方に寄せる
        if (q < 0 && p < 0) {
            p *= -1;
            q *= -1;
        }

        if (q > 0 && p < 0) {
            p *= -1;
            q *= -1;
        }
    }

    pairl get() {
        return {q,p};
    }

    double decimal() {
        return (double) q / p;
    }

    ll floor() {
        return q / p;
    }

    ll ceil() {
        return llceil(q,p);
    }

    Rational operator+(Rational Q) const { 
        auto [s,r] = Q.get();
        return Rational(q*r+p*s,p*r);
    }

    Rational operator-(Rational Q) const { 
        auto [s,r] = Q.get();
        return Rational(q*r-p*s,p*r);
    }

    Rational operator*(Rational Q) const { 
        auto [s,r] = Q.get();
        return Rational(q*s,p*r);
    }

    Rational operator/(Rational Q) const { 
        auto [s,r] = Q.get();
        return Rational(q*r,p*s);
    }

    bool operator==(const Rational &Q) const { return p == Q.p && q == Q.q; }
    bool operator!=(const Rational &Q) const { return !(*this == Q); }
    bool operator>=(const Rational &Q) const { 
        ll r = Q.p;
        ll s = Q.q;
        return q * r >= p * s;
    }
    bool operator>(const Rational &Q) const { 
        ll r = Q.p;
        ll s = Q.q;
        return q * r > p * s;
    }
    bool operator<=(const Rational &Q) const { 
        ll r = Q.p;
        ll s = Q.q;
        return q * r <= p * s;
    }
    bool operator<(const Rational &Q) const { 
        ll r = Q.p;
        ll s = Q.q;
        return q * r < p * s;
    }
    friend ostream &operator<<(ostream &os, const Rational &x) {
        return os << x.q << "/" << x.p;
    }
};
// #PORT_END#


void dummy() {
// #PORT#
// name: "permutation"
// prefix: "perm"
// description: "順列"

    vecl perm(N,0LL);
    perm[0] = 0;
    rep(i,N-1) perm[i+1] += perm[i] + 1;

    // 上の処理はiota(ALL(perm),1); とかでもできる
    
    // do-while回るごとにpermには順列が格納されている
    do {
        f(perm);
    } while(next_permutation(ALL(perm)));

// #PORT_END#
}





int nPk_example()
{

    vector<int> nums(5);
    iota(nums.begin(), nums.end(), 0);

    do
    {
        // 処理
        // numsに所定の数字が入ってる
    } while (next_permutation(nums.begin(), nums.end()));

    return 0;
}

// #PORT#
// name: "Osak"
// prefix: "osak"
// description: "高速素因数分解"
template <ll maxim>
class Osak{
    public:
    ll spf[maxim];
    Osak() {
        ll spf[maxim];
        rep(i, maxim) {
            spf[i] = i;
        }

        ll i = 2;
        while (i * i <= maxim) {
            if (spf[i] != i) {
                i += 1;
                continue;
            }

            ll j = i * i;
            while (j <= maxim) {
                if (spf[j] != j) {
                        j += i;
                        continue;
                    }

                spf[j] = i;
                j += i;
            }
            i += 1;
        }
    }

    map<ll, vector<ll>> get(ll n) { // O(log n), Σp^eの形のmapを返す
        map<ll, vector<ll>> m;
        while (n > 1) {
            if (m.count(spf[n]) == 0) m[spf[n]] = 0;

            m[spf[n]]++;
            n /= spf[n];
        }
        return m;
    }

};

// #PORT_END#

// #PORT#
// name: "divisors"
// prefix: "div"
// description: "約数 O(√N*logN)"

// 約数 O(√N*logN)
vector<ll> divisors(ll n) {
    vector<ll> res;
    REP(i,1,int(sqrt(n))+1) {
        if (n%i==0) {
            res.push_back(i);
            
            if (i!=n/i)
                res.push_back(n/i);
        }
    }
    return res;
}

// #PORT_END#


// #PORT#
// name: "modint"
// prefix: "modint"
// description: "modint"

template <ll mod>
struct ModInt {
    ll x;

    ModInt() : x(0) {}

    ModInt(ll y) : x(y >= 0 ? y % mod : (mod - (-y) % mod) % mod) {}

    ModInt &operator+=(const ModInt &p) {
        if ((x += p.x) >= mod) x -= mod;
        return *this;
    }

    ModInt &operator-=(const ModInt &p) {
        if ((x += mod - p.x) >= mod) x -= mod;
        return *this;
    }

    ModInt &operator*=(const ModInt &p) {
        x = (int)(1LL * x * p.x % mod);
        return *this;
    }

    ModInt &operator/=(const ModInt &p) {
        *this *= p.inverse();
        return *this;
    }

    ModInt operator-() const { return ModInt(-x); }

    ModInt operator+(const ModInt &p) const { return ModInt(*this) += p; }

    ModInt operator-(const ModInt &p) const { return ModInt(*this) -= p; }

    ModInt operator*(const ModInt &p) const { return ModInt(*this) *= p; }

    ModInt operator/(const ModInt &p) const { return ModInt(*this) /= p; }

    bool operator==(const ModInt &p) const { return x == p.x; }

    bool operator!=(const ModInt &p) const { return x != p.x; }

    ModInt inverse() const {
        int a = x, b = mod, u = 1, v = 0, t;
        while (b > 0) {
            t = a / b;
            swap(a -= t * b, b);
            swap(u -= t * v, v);
        }
        return ModInt(u);
    }

    ModInt pow(ll n) const {
        ModInt ret(1), mul(x);
        while (n > 0) {
            if (n & 1) ret *= mul;
            mul *= mul;
            n >>= 1;
        }
        return ret;
    }

    friend ostream &operator<<(ostream &os, const ModInt &p) {
        return os << p.x;
    }

    friend istream &operator>>(istream &is, ModInt &a) {
        ll t;
        is >> t;
        a = ModInt<mod>(t);
        return (is);
    }

    static int get_mod() { return mod; }
};

const ll MOD = 998244353;
using modint = ModInt<MOD>;

// #PORT_END#



// 組み合わせを高速に計算

// #PORT#
// name: "comb"
// prefix: "comb"
// description: "組み合わせを高速に計算"

template <typename MODINT>
class Combinatorics {
    private:
        vector<MODINT> fact;
        vector<MODINT> invfact;

    public:
        Combinatorics(ll maxim) {
            fact.assign(maxim+1,0LL);
            invfact.assign(maxim+1,0LL);

            fact[0] = 1;
            REP(i,1,maxim+1) {
                fact[i] = fact[i-1] * i;
            }
            ll mod = MODINT().get_mod();
            invfact[maxim] = MODINT(fact[maxim]).pow(mod-2);
            REP(i,1,maxim+1) {
                ll j = maxim-i;
                invfact[j] = invfact[j+1] * (j+1);
            }
        }

        MODINT nCk(ll n,ll r) {
            if (n < 0 || n < r) return 0;
            return fact[n] * invfact[r] * invfact[n-r];
        }

        MODINT nPk(ll n,ll r) {
            if (n < 0 || n < r) return 0;
            return fact[n] * invfact[n-r];
        }

};

// #PORT_END#


// 素因数分解(vector)
vector<ll> prime_vec(ll n) {
    vector<ll> a;
    while (n % 2 == 0) {
        a.push_back(2);
        n /= 2;
    }
    ll f = 3;
    while (f*f<=n) {
        if (n % f == 0) {
            a.push_back(f);
            n /= f;
        }
        else {
            f += 2;
        }
    }

    if (n!=1) {
        a.push_back(n);
    }
    return a;
}

// 素因数分解(set)
set<ll> prime_set(ll n) {
    set<ll> a;
    while (n % 2 == 0) {
        a.insert(2);
        n /= 2;
    }
    ll f = 3;
    while (f*f<=n) {
        if (n % f == 0) {
            a.insert(f);
            n /= f;
        }
        else {
            f += 2;
        }
    }

    if (n!=1) {
        a.insert(n);
    }
    return a;
}

// 素因数分解してΣp^eの形のmapを返す / (p, e)
map<ll,ll> prime_factorize(ll n) {
    map<ll,ll> counter;
    auto add = [&](ll x) {
        if (counter.count(x) == 0) counter[x] = 0;
        counter[x] += 1;
    };
    

    while (n % 2 == 0) {
        add(2);
        n /= 2;
    }
    ll f = 3;
    while (f*f<=n) {
        if (n % f == 0) {
            add(f);
            n /= f;
        }
        else {
            f += 2;
        }
    }

    if (n!=1) {
        add(n);
    }
    return counter;
}



// N以下の素数
vector<ll> get_prime(ll N){
    bool prime[N];
    vector<ll> P;

    rep(i,N) prime[i] = 1;
    prime[0] = prime[1] = 0;
    rep(i,N){
        if(prime[i]){
            P.push_back(i);
            for(int j = i + i; j < N; j+=i){
                prime[j] = 0;
            }
        }
    }

    return P;
}




// FOR DEBUG
void printVector(const vector<ll>& vec) {
    for (int value : vec) {
        cout << value << " ";
    }
    cout << endl;
}


int main() {
    // auto d = divisors(100);
    // printVector(d);

    // auto comb = Combinatorics<13,100000000+7>();
    // cout << comb.nPk(12,12) << endl;

    // auto d = prime_vec(3000);
    // printVector(d);

    // auto d = prime_factorize(12345);

    auto d = get_prime(10000000);
    // printVector(d);
    cout << d.size() << endl;

    // for (const auto& [key, value] : d){
    //     std::cout << key << " => " << value << "\n";
    // }

    return 1;
}


// #PORT#
// name: "comb_simple"
// prefix: "comb_simple"
// description: "MOD取らずに二項係数を計算"

void nCk() {
    auto dp = genarr(x+20,x+20,0LL);
    dp[0][0] = 1;
    for (int i = 1; i < x+20; ++i) {
        dp[i][0] = 1;
        for (int j = 1; j < x+20; j++) {
            dp[i][j] = dp[i - 1][j - 1] + dp[i - 1][j];
        }
    }
}

// #PORT_END#


// #PORT#
// name: "modinv"
// prefix: "modinv"
// description: "任意のMODで逆元を取る"

// 任意のMODで逆元を取る
ll modinv(ll a, ll mod) {
    ll b = mod, u = 1, v = 0;
    while (b) {
        ll t = a/b;
        a -= t*b; swap(a, b);
        u -= t*v; swap(u, v);
    }
    u %= mod;
    if (u < 0) u += mod;
    return u;
}
 
// #PORT_END#

// #PORT#
// name: "modint_noconst"
// prefix: "modint_noconst"
// description: "可変modint"

ll MOD;

struct ModInt {
    ll x;
    ll mod = MOD;

    ModInt() : x(0) {}

    ModInt(ll y) {
        x = y >= 0 ? y % mod : (mod - (-y) % mod) % mod;
    }

    ModInt &operator+=(const ModInt &p) {
        if ((x += p.x) >= mod) x -= mod;
        return *this;
    }

    ModInt &operator-=(const ModInt &p) {
        if ((x += mod - p.x) >= mod) x -= mod;
        return *this;
    }

    ModInt &operator*=(const ModInt &p) {
        x = (int)(1LL * x * p.x % mod);
        return *this;
    }

    ModInt &operator/=(const ModInt &p) {
        *this *= p.inverse();
        return *this;
    }

    ModInt operator-() const { return ModInt(-x); }

    ModInt operator+(const ModInt &p) const { return ModInt(*this) += p; }

    ModInt operator-(const ModInt &p) const { return ModInt(*this) -= p; }

    ModInt operator*(const ModInt &p) const { return ModInt(*this) *= p; }

    ModInt operator/(const ModInt &p) const { return ModInt(*this) /= p; }

    bool operator==(const ModInt &p) const { return x == p.x; }

    bool operator!=(const ModInt &p) const { return x != p.x; }

    ModInt inverse() const {
        int a = x, b = mod, u = 1, v = 0, t;
        while (b > 0) {
            t = a / b;
            swap(a -= t * b, b);
            swap(u -= t * v, v);
        }
        return ModInt(u);
    }

    ModInt pow(ll n) const {
        ModInt ret(1), mul(x);
        while (n > 0) {
            if (n & 1) ret *= mul;
            mul *= mul;
            n >>= 1;
        }
        return ret;
    }

    friend ostream &operator<<(ostream &os, const ModInt &p) {
        return os << p.x;
    }

    friend istream &operator>>(istream &is, ModInt &a) {
        ll t;
        is >> t;
        a = ModInt(t);
        return (is);
    }
};

using modint = ModInt;

// #PORT_END#

// #PORT#
// name: "modcount"
// prefix: "modcount"
// description: "n ∈ [l,r] ≡ x (mod m)"

ll modcount(ll l, ll r, ll m, ll x) { // 閉区間[l,r]内に mod m で x となるものはいくつあるか 
    ll R = (r - x + m) / m, L = (l - x + m - 1) / m;
    return R - L;
}

// #PORT_END#