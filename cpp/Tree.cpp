#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const double EPS = 1e-9;
const ll INF = ((1LL<<62)-(1LL<<31));
typedef vector<ll> vecl;
typedef pair<ll, ll> pairl;
#define REP(i, x, n) for (ll i = x; i < n; i++)
#define rep(i, n) REP(i, 0, n)

// #PORT#
// name: "tree"
// prefix: "tree"
// description: "木のライブラリ"

class Tree {
    public:
    vecl leaves;

    Tree(ll V) {
        LOG = log2(V + 10) + 10;
        parents.resize(LOG + 1, vecl(V + 1));
        depths.resize(V + 1);
        children_count.assign(V + 1, 0LL);
        G.resize(V + 1);
    }

    void add(ll x, ll y) {
        G[x].push_back(y);
        G[y].push_back(x);
    }

    void build(ll root = 1) {
        ll V = G.size() - 1;
        dfs(root, -1, 0);
        rep(k, LOG) {
            rep(v, V) {
                if (parents[k][v] < 0)
                    parents[k + 1][v] = -1;
                else
                    parents[k + 1][v] = parents[k][parents[k][v]];
            }
        }
        set_leaves();
        set_children_count();
    }

    ll lca(ll u, ll v) {
        if (depths[u] > depths[v]) swap(u, v);
        rep(k, LOG) {
            if ((depths[v] - depths[u]) >> k & 1) v = parents[k][v];
        }
        if (u == v) return u;
        for (ll k = LOG; k >= 0; k--) {
            if (parents[k][u] != parents[k][v]) {
                u = parents[k][u];
                v = parents[k][v];
            }
        }
        return parents[0][u];
    }

    ll distance(ll u, ll v) {
        ll t = lca(u, v);
        return depths[u] + depths[v] - 2 * depths[t];
    }

    ll get_depth(ll u) { return depths[u]; }

    ll get_parent(ll u) { return parents[0][u]; }

    ll count_children(ll u) { return children_count[u]; }

    private:
    vector<vecl> G;
    vector<vecl> parents;
    vector<ll> depths;
    vector<ll> children_count;
    ll LOG;

    void dfs(ll v, ll parent, ll depth) {
        parents[0][v] = parent;
        depths[v] = depth;
        for (auto& i : G[v]) {
            if (i != parent) dfs(i, v, depth + 1);
        }
    }

    void set_leaves() {
        ll V = G.size() - 1;
        vecl res;
        rep(i, V) {
            if (G[i].size() == 1) res.push_back(i);
        }
        leaves = res;
    }

    void set_children_count() {
        priority_queue<pairl> nexts;
        for (auto leaf : leaves) nexts.push({depths[leaf], leaf});

        unordered_set<ll> seen;
        while (!nexts.empty()) {
            auto [d, current] = nexts.top(); nexts.pop();
            ll parent = parents[0][current];
            if (parent == -1) continue;

            children_count[parent] += children_count[current] + 1;

            if (seen.count(parent) > 0) continue;
            nexts.push({depths[parent], parent});
            seen.insert(parent);
        }
    }

};

// #PORT_END#


// #PORT#
// name: "rerooting"
// prefix: "reroot"
// description: "全方位DP"

///// 全方位DP: https://algo-logic.info/tree-dp/

template <class DP, DP (*e)(), DP (*merge)(DP, DP), DP (*add_root)(DP)>
struct Rerooting {
    struct Edge {
        int to;
    };
    using Graph = vector<vector<Edge>>;
    vector<vector<DP>> dp;  // dp[v][i]: vから出るi番目の有向辺に対応する部分木のDP
    vector<DP> ans;         // ans[v]: 頂点vを根とする木の答え
    Graph G;
    Rerooting(int N) : G(N) {
        dp.resize(N);
        ans.assign(N, identity);
    }
    void add_edge(int a, int b) {
        G[a].push_back({b});
    }
    void build() {
        dfs(0);            // 普通に木DP
        bfs(0, identity);  // 残りの部分木に対応するDPを計算
    }
    DP dfs(int v, int p = -1) {  // 頂点v, 親p
        DP dp_cum = identity;
        int deg = G[v].size();
        dp[v] = vector<DP>(deg, identity);
        for (int i = 0; i < deg; i++) {
            int u = G[v][i].to;
            if (u == p) continue;
            dp[v][i] = dfs(u, v);
            dp_cum = merge(dp_cum, dp[v][i]);
        }
        return add_root(dp_cum);
    }
    void bfs(int v, const DP& dp_p, int p = -1) { 
        int deg = G[v].size();
        for (int i = 0; i < deg; i++) {  // 前のbfsで計算した有向辺に対応する部分木のDPを保存
            if (G[v][i].to == p) dp[v][i] = dp_p;
        }
        vector<DP> dp_l(deg + 1, identity), dp_r(deg + 1, identity);  // 累積merge
        for (int i = 0; i < deg; i++) {
            dp_l[i + 1] = merge(dp_l[i], dp[v][i]);
        }
        for (int i = deg - 1; i >= 0; i--) {
            dp_r[i] = merge(dp_r[i + 1], dp[v][i]);
        }
        ans[v] = add_root(dp_l[deg]);  // 頂点 v の答え
        for (int i = 0; i < deg; i++) {  // 一つ隣の頂点に対しても同様に計算
            int u = G[v][i].to;
            if (u == p) continue;
            bfs(u, add_root(merge(dp_l[i], dp_r[i + 1])), v);
        }
    }
};


struct DP {  // DP の型 (モノイド)
    ll dp;
    DP(ll dp_) : dp(dp_) {}
};

DP e() { return DP(-1); }  // 単位元(末端の値は add_root(e) になるので注意)

DP merge(DP dp_cum, DP d) { // 累積計算をするための二項演算 (親へとDP情報を「上げる」前に子の情報を全てmergeする)
    return DP(max(dp_cum.dp, d.dp));
}

DP add_root(DP d) { // mergeしたDP情報を親に渡す前に加工する
    return DP(d.dp + 1);
};


// #PORT_END#

