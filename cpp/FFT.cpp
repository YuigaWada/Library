#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
#define REP(i, x, n) for (ll i = x; i < n; i++)
#define rep(i, n) REP(i, 0, n)

// FFT
namespace FFT {
  using real = double;

  struct C {
    real x, y;

    C() : x(0), y(0) {}

    C(real x, real y) : x(x), y(y) {}

    inline C operator+(const C &c) const { return C(x + c.x, y + c.y); }

    inline C operator-(const C &c) const { return C(x - c.x, y - c.y); }

    inline C operator*(const C &c) const { return C(x * c.x - y * c.y, x * c.y + y * c.x); }

    inline C conj() const { return C(x, -y); }
  };

  const real PI = acosl(-1);
  ll base = 1;
  vector< C > rts = { {0, 0},
                     {1, 0} };
  vector< ll > rev = {0, 1};


  void ensure_base(ll nbase) {
    if(nbase <= base) return;
    rev.resize(1 << nbase);
    rts.resize(1 << nbase);
    for(ll i = 0; i < (1 << nbase); i++) {
      rev[i] = (rev[i >> 1] >> 1) + ((i & 1) << (nbase - 1));
    }
    while(base < nbase) {
      real angle = PI * 2.0 / (1 << (base + 1));
      for(ll i = 1 << (base - 1); i < (1 << base); i++) {
        rts[i << 1] = rts[i];
        real angle_i = angle * (2 * i + 1 - (1 << base));
        rts[(i << 1) + 1] = C(cos(angle_i), sin(angle_i));
      }
      ++base;
    }
  }

  void fft(vector< C > &a, ll n) {
    assert((n & (n - 1)) == 0);
    ll zeros = __builtin_ctz(n);
    ensure_base(zeros);
    ll shift = base - zeros;
    for(ll i = 0; i < n; i++) {
      if(i < (rev[i] >> shift)) {
        swap(a[i], a[rev[i] >> shift]);
      }
    }
    for(ll k = 1; k < n; k <<= 1) {
      for(ll i = 0; i < n; i += 2 * k) {
        for(ll j = 0; j < k; j++) {
          C z = a[i + j + k] * rts[j + k];
          a[i + j + k] = a[i + j] - z;
          a[i + j] = a[i + j] + z;
        }
      }
    }
  }

  vector< ll > multiply(const vector< ll > &a, const vector< ll > &b) {
    ll need = (ll) a.size() + (ll) b.size() - 1;
    ll nbase = 1;
    while((1 << nbase) < need) nbase++;
    ensure_base(nbase);
    ll sz = 1 << nbase;
    vector< C > fa(sz);
    for(ll i = 0; i < sz; i++) {
      ll x = (i < (ll) a.size() ? a[i] : 0);
      ll y = (i < (ll) b.size() ? b[i] : 0);
      fa[i] = C(x, y);
    }
    fft(fa, sz);
    C r(0, -0.25 / (sz >> 1)), s(0, 1), t(0.5, 0);
    for(ll i = 0; i <= (sz >> 1); i++) {
      ll j = (sz - i) & (sz - 1);
      C z = (fa[j] * fa[j] - (fa[i] * fa[i]).conj()) * r;
      fa[j] = (fa[i] * fa[i] - (fa[j] * fa[j]).conj()) * r;
      fa[i] = z;
    }
    for(ll i = 0; i < (sz >> 1); i++) {
      C A0 = (fa[i] + fa[i + (sz >> 1)]) * t;
      C A1 = (fa[i] - fa[i + (sz >> 1)]) * t * rts[(sz >> 1) + i];
      fa[i] = A0 + A1 * s;
    }
    fft(fa, sz >> 1);
    vector< ll > ret(need);
    for(ll i = 0; i < need; i++) {
      ret[i] = llround(i & 1 ? fa[i >> 1].y : fa[i >> 1].x);
    }
    return ret;
  }
};











// NTT
template< ll mod >
struct NTT {

  vector< ll > rev, rts;
  ll base, max_base, root;

  NTT() : base(1), rev{0, 1}, rts{0, 1} {
    assert(mod >= 3 && mod % 2 == 1);
    auto tmp = mod - 1;
    max_base = 0;
    while(tmp % 2 == 0) tmp >>= 1, max_base++;
    root = 2;
    while(mod_pow(root, (mod - 1) >> 1) == 1) ++root;
    assert(mod_pow(root, mod - 1) == 1);
    root = mod_pow(root, (mod - 1) >> max_base);
  }

  inline ll mod_pow(ll x, ll n) {
    ll ret = 1;
    while(n > 0) {
      if(n & 1) ret = mul(ret, x);
      x = mul(x, x);
      n >>= 1;
    }
    return ret;
  }

  inline ll inverse(ll x) {
    return mod_pow(x, mod - 2);
  }

  inline unsigned add(unsigned x, unsigned y) {
    x += y;
    if(x >= mod) x -= mod;
    return x;
  }

  inline unsigned mul(unsigned a, unsigned b) {
    return 1ull * a * b % (unsigned long long) mod;
  }

  void ensure_base(ll nbase) {
    if(nbase <= base) return;
    rev.resize(1 << nbase);
    rts.resize(1 << nbase);
    for(ll i = 0; i < (1 << nbase); i++) {
      rev[i] = (rev[i >> 1] >> 1) + ((i & 1) << (nbase - 1));
    }
    assert(nbase <= max_base);
    while(base < nbase) {
      ll z = mod_pow(root, 1 << (max_base - 1 - base));
      for(ll i = 1 << (base - 1); i < (1 << base); i++) {
        rts[i << 1] = rts[i];
        rts[(i << 1) + 1] = mul(rts[i], z);
      }
      ++base;
    }
  }


  void ntt(vector< ll > &a) {
    const ll n = (ll) a.size();
    assert((n & (n - 1)) == 0);
    ll zeros = __builtin_ctz(n);
    ensure_base(zeros);
    ll shift = base - zeros;
    for(ll i = 0; i < n; i++) {
      if(i < (rev[i] >> shift)) {
        swap(a[i], a[rev[i] >> shift]);
      }
    }
    for(ll k = 1; k < n; k <<= 1) {
      for(ll i = 0; i < n; i += 2 * k) {
        for(ll j = 0; j < k; j++) {
          ll z = mul(a[i + j + k], rts[j + k]);
          a[i + j + k] = add(a[i + j], mod - z);
          a[i + j] = add(a[i + j], z);
        }
      }
    }
  }


  vector< ll > multiply(vector< ll > a, vector< ll > b) {
    ll need = a.size() + b.size() - 1;
    ll nbase = 1;
    while((1 << nbase) < need) nbase++;
    ensure_base(nbase);
    ll sz = 1 << nbase;
    a.resize(sz, 0);
    b.resize(sz, 0);
    ntt(a);
    ntt(b);
    ll inv_sz = inverse(sz);
    for(ll i = 0; i < sz; i++) {
      a[i] = mul(a[i], mul(b[i], inv_sz));
    }
    reverse(a.begin() + 1, a.end());
    ntt(a);
    a.resize(need);
    return a;
  }
};
