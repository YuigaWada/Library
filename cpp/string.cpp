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
// name: "AhoCorasick"
// prefix: "ahocorasick"
// description: "AhoCorasick"


using MatchedPair = vector<pair<string,pairl>>;
struct AhoCorasick {
    char e;
    vector<string> Ps;
    vector<vecl> V; // V[i]: 状態iと対応するパターン群 (iに到達=パターンiが検出)
    vecl failure; // failure[i]: 検索失敗時の状態iの遷移先
    vector<vecl> g; // g[i][c]: 状態iから遷移c+'a'(c+e)した時の遷移先

    AhoCorasick(vector<string> patterns, char e = 'a') : Ps(patterns), e(e) { // edgeに文字, nodeに状態
        build(); // Trie木を作成
        make_failure(); // BFSでfailureを埋めていく
    }

    void build() {
        g.emplace_back(vecl(26,0)); // g[0][i] = 0
        V.emplace_back(vecl());
        rep(i,sz(Ps)) {
            ll v = 0;
            for (auto _c : Ps[i]) {
                int c = _c - e;
                if (v > 0 && g[v][c] != -1) {
                    v = g[v][c];
                    continue;
                }
                else if (v == 0 && g[v][c] != 0) {
                    v = g[v][c];
                    continue;
                }

                int vs = sz(V);
                V.emplace_back(vecl());
                g.emplace_back(vecl(26,-1));
                g[v][c] = vs;
                v = vs;
            }
            assert(v < sz(V));
            V[v].emplace_back(i);
        }
    }

    void make_failure() {
        failure.assign(sz(V)+1,0LL);
        queue<ll> S;
        S.push(0);
        while (!S.empty()) {
            ll current = S.front(); S.pop();
            rep(i,26) {
                ll next = g[current][i];
                if (next <= 0) continue;

                S.push(next);
                if (current != 0) { // 最長の接尾辞をfailureに
                    ll f = failure[current];
                    while (g[f][i] == -1) f = failure[f];
                    failure[next] = g[f][i]; // fからi+'a'の遷移ができる状態fを探して,その遷移先g[f][i]をnextのfailureとする
                    for (auto p : V[failure[next]]) V[next].emplace_back(p);
                }
            }
        }
    }

    void search(const string &S, MatchedPair &X) { // マッチしたものはXに格納
        ll v = 0;
        rep(i,sz(S)) {
            int c = S[i] - e;
            while (g[v][c] == -1) v = failure[v];
            v = g[v][c];
            rep(j,sz(V[v])) {
                string x = Ps[V[v][j]];
                X.emplace_back(x,pairl{i-sz(x)+1,i});
            }
        }
    }

    void print() {
        rep(i,sz(V)) {
            rep(j,26) {
                int p = i == 0 ? 0 : -1;
                if (g[i][j] == p) continue;
                printf("%d -> %d (%c)\n",i,g[i][j],j+e);
            }
        }
        dump(failure);
    }
};

// #PORT_END#


// #PORT#
// name: "KMP"
// prefix: "kmp"
// description: "KMP"

struct KMP {
    string P;
    vecl table; // table[i]: iでfailした時, どこから検索し直すか
    KMP(string pattern) : P(pattern) { // O(|P|)
        table.assign(sz(P)+1,0LL); 
        table[0] = -1;
        for (int i = 0, p = -1; i < sz(P); i++) {
            while (p >= 0 && P[i] != P[p]) p = table[p]; // P[i]と重複している先頭部分P[0:p]を探す(ないならp=-1)
            table[i+1] = ++p; // S[i-p:i]はP[0:p]と一致しているので, i+1でfailした場合, 検索し直す先頭はP[p+1]から
        }
    }

    ll count(const string &S) { // O(|S|)
        ll res = 0;
        for (int i = 0, p = 0; i < sz(S); i++) {
            while (p >= 0 && S[i] != P[p]) p = table[p];
            if (p == sz(P) - 1) res++;
            p++;
        }
        return res;
    }

    ll find(const string &S) { // マッチしたSの先頭indexを返す, O(|S|)
        for (int i = 0, p = 0; i < sz(S); i++) {
            while (p >= 0 && S[i] != P[p]) p = table[p];
            if (p == sz(P) - 1) return i - sz(P) + 1;
            p++;
        }
        return -1;
    }
};


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