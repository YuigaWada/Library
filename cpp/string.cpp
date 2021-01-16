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
#define contains(S,x) find(ALL(S),x) != S.end()
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
// name: "split"
// prefix: "split"
// description: "文字分割"

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        if (!item.empty()) elems.push_back(item);
    }
    return elems;
}

// #PORT_END#

// #PORT#
// name: "RollingHash"
// prefix: "rollinghash"
// description: "ロリハ"

struct RollingHash {
    const ll BASE = 9973;
    const ll NONMOD = -1;
    vecl mods;

    int N;
    void initmods() {
        // mods.push_back(1e9+7);
        // mods.push_back(999999937);
        mods.push_back(NONMOD);
    }

    vector<vecl> hashes, powers;
    RollingHash() {}
    RollingHash(const string &s) : N(sz(s)) {
        initmods();
        hashes.resize(sz(mods)), powers.resize(sz(mods));
        rep(i,sz(mods)) {
            hashes[i].assign(N+1,0LL), powers[i].assign(N+1,0LL);
            hashes[i][0] = 0, powers[i][0] = 1;
            rep(j,N) {
                if (mods[i] != NONMOD) powers[i][j + 1] = powers[i][j] * BASE % mods[i];
                else powers[i][j + 1] = powers[i][j] * BASE;

                if (mods[i] != NONMOD) hashes[i][j + 1] = (hashes[i][j] * BASE + s[j]) % mods[i];
                else hashes[i][j + 1] = (hashes[i][j] * BASE + s[j]);
            }
        }
    }

    inline long long hash(int l, int r, int i) {
        if (mods[i] != NONMOD) return ((hashes[i][r] - hashes[i][l] * powers[i][r - l]) % mods[i] + mods[i]) % mods[i];
        return hashes[i][r] - hashes[i][l] * powers[i][r - l];
    }

    inline bool match(int l1, int r1, int l2, int r2) {
        bool ret = 1;
        for (int i = 0; i < sz(mods); i++) ret &= hash(l1, r1, i) == hash(l2, r2, i);
        return ret;
    }

    inline bool match(int l1, int l2, int k) {
        return match(l1, l1 + k, l2, l2 + k);
    }

    int lcp(int i, int j) { // LCP(S[i:], S[j:]).length();
        int l = 0, r = min(N- i, N - j) + 1;
        while (l + 1 < r) {
            int mid = (l + r) / 2;
            if (match(i, i + mid, j, j + mid)) l = mid;
            else r = mid;
        }
        return l;
    }
};

// #PORT_END#