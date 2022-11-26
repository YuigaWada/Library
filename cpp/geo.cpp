#include <bits/stdc++.h>
using namespace std;
#define LOCAL // 提出時はコメントアウト
#define DEBUG_

typedef long long ll;
const double EPS = 1e-9;
const ll INF = ((1LL<<62)-(1LL<<31));
typedef vector<ll> vecl;
typedef vector<ll> vl;
typedef pair<ll, ll> pl;
typedef vector<pl> vp;
template<typename T> using uset = unordered_set<T>;
template<typename T, typename U> using mapv = map<T,vector<U>>;
template<typename T, typename U> using umap = unordered_map<T,U>;

#define ALL(v) v.begin(), v.end()
#define REP(i, x, n) for(int i = x; i < n; i++)
#define rep(i, n) REP(i, 0, n)
#define sz(x) (ll)x.size()
#define pb(x) push_back(x)
#define eb(x) emplace_back(x)
#define pr(x) cout << x << endl
ll llceil(ll a,ll b) { return (a+b-1)/b; }
template<class T> inline bool chmax(T& a, T b) { if (a < b) { a = b; return true; } return false; }
template<class T> inline bool chmin(T& a, T b) { if (a > b) { a = b; return true; } return false; }
template<class T> vector<vector<T>> genarr(ll n, ll m, T init) { return vector<vector<T>>(n,vector<T>(m,init)); }

///// DEBUG
#define DUMPOUT cerr
#define repi(itr, ds) for (auto itr = ds.begin(); itr != ds.end(); itr++)
template<typename T>istream&operator>>(istream&is,vector<T>&vec){for(T&x:vec)is>>x;return is;}
template<typename T,typename U>ostream&operator<<(ostream&os,pair<T,U>&pair_var){os<<"("<<pair_var.first<<", "<<pair_var.second<<")";return os;}
template<typename T>ostream&operator<<(ostream&os,vector<T>&vec){os<<"{";for(int i=0;i<vec.size();i++){os<<vec[i]<<(i+1==vec.size()?"":", ");}
os<<"}";return os;}
template<typename T,typename U>ostream&operator<<(ostream&os,map<T,U>&map_var){os<<"{";repi(itr,map_var){os<<*itr;itr++;if(itr!=map_var.end())os<<", ";itr--;}
os<<"}";return os;}
template<typename T>ostream&operator<<(ostream&os,set<T>&set_var){os<<"{";repi(itr,set_var){os<<*itr;itr++;if(itr!=set_var.end())os<<", ";itr--;}
os<<"}";return os;}
template<typename T>ostream&operator<<(ostream&os,uset<T>&set_var){os<<"{";repi(itr,set_var){if(itr!=set_var.begin())os<<", "; os<<*itr;}
os<<"}";return os;}
template<typename T,typename U>ostream&operator<<(ostream&os,umap<T,U>&map_var){os<<"{";repi(itr,map_var){if(itr!=map_var.begin())os<<", "; os<<*itr;}
os<<"}";return os;}
void dump_func(){DUMPOUT<<endl;}
template<class Head,class...Tail>void dump_func(Head&&head,Tail&&...tail){DUMPOUT<<head;if(sizeof...(Tail)>0){DUMPOUT<<", ";}
dump_func(std::move(tail)...);}
void inp_func(){}
template<class Head,class...Tail>void inp_func(Head&&head,Tail&&...tail){cin>>head; inp_func(std::move(tail)...);}
#define ip(...) inp_func(__VA_ARGS__)
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
// name: "geo"
// prefix: "geo"
// description: "図形"

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
};

struct Argsorted {
    ll id, x, y; // x,yは変更されるのでidで識別する, (x,y)は符号情報を保持
    Rational Q;

    Argsorted() {}
    Argsorted(int idx, ll X, ll Y) : id(idx) {
        if (X == 0 && Y == 0) X = 1, Y = 0;
        else if (X == 0) X = 1, Y = Y > 0 ? INF : -INF; // "X → 0" ⇒ "Y/X → ±inf"
        Q = Rational(Y,X); // Y/X
        x = X, y = Y;
    };

    int quadrant() {
        ll _x = x, _y = y;
        if (_x == 1 && abs(_y) == INF) _x = 0;

        if (_x == 0) return _y >= 0 ? 2 : 4;
        else if (_y == 0) return _x >= 0 ? 1 : 3;
        
        if (_x * _y >= 0) return _x >= 0 && _y >= 0 ? 1 : 3;
        else return _x <= 0 && _y >= 0 ? 2 : 4;
    }
};

vector<Argsorted> argsort(const vector<pairl> &X) { // 偏角ソート / (x,y)は整数
    vector<Argsorted> res((int)X.size());
    rep(i,(int)X.size()) res[i] = Argsorted(i,X[i].first,X[i].second);
    sort(ALL(res),[](auto &x, auto &y) { // return atan2(x.y,x.x) < atan2(y.y,y.x);
        ll Xpos = x.quadrant(), Ypos = y.quadrant();
        if (Xpos == Ypos) {
            if (abs(y.y) == INF) return false;
            else if (abs(x.y) == INF) return true;

            if (y.y == 0) return false;
            else if (x.y == 0) return true;
            
            return x.Q < y.Q;
        }
        else return Xpos < Ypos;
    });
    return res;
}


ll det(pl x, pl y) {
    auto [a,b] = x;
    auto [c,d] = y;
    return a*d - b*c;
}

pl sub(pl x, pl y) { // XY
    auto [a,b] = x;
    auto [c,d] = y;
    return {c - a, d - b};
}

vl convex_hull(vp P) { // monotone chain O(NlogN) : return index
    ll size = 0;
    sort(ALL(P));

    vl ch;
    rep(i,sz(P)) {
        while (size > 1) {
            pl current = sub(P[ch[size-2]],P[ch[size-1]]);
            pl newer = sub(P[ch[size-2]],P[i]);
            if (det(current,newer) > 0) break;
            size--;
            ch.pop_back();
        }
        ch.pb(i);
        size++;
    }

    ll _size = size;
    rep(_i,sz(P)-1) {
        ll i = sz(P) - 2 - _i;
        while (size > _size) {
            pl current = sub(P[ch[size-2]],P[ch[size-1]]);
            pl newer = sub(P[ch[size-2]],P[i]);
            if (det(current,newer) > 0) break;
            size--;
            ch.pop_back();
        }
        ch.pb(i);
        size++;
    }
    ch.pop_back();
    return ch;
}

// #PORT_END#