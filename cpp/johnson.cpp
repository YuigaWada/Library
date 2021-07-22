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
// name: "Johnson's algorithm"
// prefix: "johnson"
// description: "Johnson's algorithm"

struct Johnson { // Johnson's algorithm
    using Edge = tuple<int,ll,ll>; // to, cost, new_cost

    int V;
    vector<vector<Edge>> g;
    vector<ll> d; // d[x] := cost "s-x", s: 超頂点(index=N)
    bool has_negative_cycle = false;

    Johnson(int N, vector<vector<pair<int,ll>>> G) : V(N+1) {
        g.assign(V,vector<Edge>());
        rep(i,sz(G)) for (auto [b,t] : G[i]) g[i].emplace_back(b,t,INF);
        rep(i,N) g[N].emplace_back(i,0,INF); // N: 超頂点
    }

    int rebuild() { // 重みを再構成 (重みが非負整数であるように修正)
        bellman_ford();
        if (has_negative_cycle) return 0;

        rep(i,V-1) {  // d[x] := cost "s-x", V := N + 1
            rep(j,sz(g[i])) {
                auto [b,t,_] = g[i][j];
                ll nt = t+d[i]-d[b];
                if (b == V-1) nt = INF; // bが超頂点の場合
                g[i][j] = {b,t,nt};
            }
        }
        return 1;
    }

    vector<ll> run(int start) { // dijkstra
        vecl dist(V,INF), dist2(V,INF);
        dist[start] = dist2[start] = 0;

        priority_queue<pairl,vector<pairl>,greater<pairl>> PQ;
        PQ.emplace(dist[start],start);

        while (!PQ.empty()) {
            auto [_,current] = PQ.top(); PQ.pop();

            for (auto [to,cost,ncost] : g[current]) {
                if (!chmin(dist[to],dist[current] + ncost)) continue;
                PQ.emplace(dist[to],to);
                dist2[to] = dist2[current] + cost; // 本来のコストで更新したものも記憶
            }
        }
        
        return dist2;
    }

    void bellman_ford() {
        d.assign(V, INF);
        ll start = V-1; // start: 超頂点
        d[start] = 0;
        rep(i,V-1) { // 負の経路が存在しないなら高々V-1回で終了する
            rep(j,V) {
                for (auto [to,cost,_] : g[j]) {
                    if(d[j] == INF) continue;
                    chmin(d[to], d[j] + cost);
                }
            }
        }

        rep(i,V) {
            rep(j,V) {
                for (auto [to,cost,_] : g[j]) {
                    if(d[j] == INF) continue;
                    has_negative_cycle |= d[to] != min(d[to], d[j] + cost); // V + i回目で更新されるなら負の閉路あり
                    if (has_negative_cycle) break;
                }
            }
            if (has_negative_cycle) break;
        }
    }
};

// #PORT_END# 



int solve(ostringstream &cout) { // AOJ: All Pairs Shortest Path
    #ifdef LOCAL
    ifstream in("../../Atcoder/input.txt");
    cin.rdbuf(in.rdbuf());
    #endif

    ll N,M;
    cin>>N>>M;

    vector<vector<pair<int,ll>>> g(N);

    rep(i,M) {
        ll a,b,t;
        cin>>a>>b>>t;
        g[a].emplace_back(b,t);
    }

    Johnson john(N,g);
    if (!john.rebuild()) {
        cout << "NEGATIVE CYCLE" << endl;
        return 0;
    }

    rep(start,N) {
        auto dist = john.run(start);
        rep(i,N) {
            if (dist[i] == INF) cout << "INF" << " \n"[i==N-1];
            else cout << dist[i] << " \n"[i==N-1];
        }
    }

    return 0;
}

int main() {
    ostringstream oss;
    int res = solve(oss);

    cout << oss.str();
    return res;
}
