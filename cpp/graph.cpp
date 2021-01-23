#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
#define REP(i, x, n) for (int i = x; i < n; i++)
#define rep(i, n) REP(i, 0, n)

// #PORT#
// name: "Graph"
// prefix: "graph"
// description: "グラフテンプレート"
struct Point { ll i,j; };

ll H,W;

int di[] = {1, 0, -1, 0, 1, -1, -1, 1, 0};
int dj[] = {0, 1, 0, -1, 1, 1, -1, -1, 0};

bool valid_point(ll i, ll j) { return 0 <= i && i < H && 0 <= j && j < W; }
ll transform(ll i, ll j) { return i * W + j; }
Point transform(ll transed) { return {transed / W, transed % W}; }
bool is_alphabet(char c) { return 'a' <= c && c <= 'z'; }


struct Edge {
    ll from, to; // 2次元の場合は一次に変換しておく
    ll cost;
    ll idx;

    Edge() = default;
    Edge(ll from, ll to, ll cost = 1, ll idx = -1) : from(from), to(to), cost(cost), idx(idx) {}
};

struct ShortestPath { // 最短経路
    vecl d;
    vecl update_list; // update_list[i]: 頂点iを更新した辺番号
    bool diverge; //start → goal でコストが負に発散する
};

struct Graph {
    vector<vector<Edge>> g;
    vector<Edge> allE; // allE[es]: 辺番号esの辺
    int es; // 辺番号

    Graph() = default;
    explicit Graph(int n) : g(n), es(0) {}

    size_t size() const {
        return g.size();
    }

    void add_directed_edge(int from, int to, ll cost = 1) {
        Edge from_e = {from, to, cost, es++};
        g[from].push_back(from_e);
        allE.push_back(from_e);
    }

    void add_edge(int from, int to, ll cost = 1) {
        Edge from_e = {from, to, cost, es++};
        Edge to_e = {to, from, cost, es++};   
        // Edge to_e = {to, from, cost, es++}; // 行きと帰りを区別しないならこっち
        g[from].push_back(from_e);
        g[to].push_back(to_e);

        allE.push_back(from_e);
        allE.push_back(to_e);
    }

    vector<Edge> get_edges(int v) {
        return g[v];
    }

    vector<Edge> get_all_edges() {
        return allE;
    }

    // wallを避けて辺情報を追加する
    void auto_make(vector<vector<char>> G, char wall = '#') {
        rep(i,H) {
            rep(j,W) {
                if (G[i][j] == wall) continue;
                rep(k,2) {
                    ll ni = i + di[k];
                    ll nj = j + dj[k];
                    if (!valid_point(ni,nj)) continue;
                    if (G[ni][nj] == wall) continue;
                    add_edge(transform(i,j),transform(ni,nj),1);
                }
            }
        }
    }

    ////////////////////////////////////////////////////

    // BFS: O(E+V)
    void bfs(Point start, vector<vecl> &dp, ll initial_cost = 0) {
        dp[start.i][start.j] = initial_cost;

        queue<ll> S;
        S.push(transform(start.i,start.j));

        while (!S.empty()) {
            ll current = S.front(); S.pop();
            auto [i,j] = transform(current);

            for (auto e : get_edges(current)) {
                auto [ni,nj] = transform(e.to);
                if (!chmin(dp[ni][nj],dp[i][j] + e.cost)) continue;
                S.push(e.to);
            }
        }
    }

    // 01BFS: O(E+V)
    void zero_one_bfs(Point start, vector<vecl> &dp, ll initial_cost = 0) {
        deque<pairl> S;
        S.push_front({transform(start.i, start.j),initial_cost});

        while (!S.empty()) {
            auto [current,cost] = S.front(); S.pop_front();
            auto [i,j] = transform(current);

            if (!chmin(dp[i][j], cost)) continue;
            for (auto e : get_edges(current)) {
                auto next = e.to;
                auto [ni,nj] = transform(next);
                if (e.cost) {
                    S.push_back({next,dp[i][j]+e.cost}); // コスト1の辺は後ろにpush
                }
                else {
                    S.push_front({next,dp[i][j]+e.cost}); // コスト0の辺は前にpush
                }
            }
        }
    }

    // ベルマンフォード: O(V*E), 負の閉路が存在し最短距離が負へと発散するときは ShortestPath(false) 返す
    ShortestPath bellman_ford(ll start, ll goal, ll V) {
        vecl d(V, INF);
        vecl update_list(V+10,-1); // update_list[i]: 頂点iを更新した辺番号, 最初は-1
        d[start] = 0;
        rep(i,V-1) { // 負の経路が存在しないなら高々V-1回で終了する
            for(auto &e : allE) {
                if(d[e.from] == INF) continue;
                
                bool updated = chmin(d[e.to], d[e.from] + e.cost);
                if (updated) {
                    update_list[e.to] = e.idx;
                }
            }
        }

        rep(i,V) {
            for(auto &e : allE) {
                if(d[e.from] == INF) continue;
                
                bool updated = d[e.to] != min(d[e.to], d[e.from] + e.cost); // V + i回目で更新されるなら負の閉路あり

                // goalを含む閉路が存在する場合、追加したV回ループ内で確実にgoalへとリーチするためにコストを-infにしておく
                if (updated) { 
                    d[e.to] = -INF; 
                    if (e.to == goal) return ShortestPath{vecl(),vecl(),true}; // 負の閉路からgoalに到達できるなら負に発散するからダメ
                }
            }
        }
        
        return ShortestPath {d, update_list, false}; // update_list[goal] から辺番号をさかのぼれば経路復元できる
    }
    
    // ダイクストラ: O((E+V)logV)
    ShortestPath dijkstra(ll start, ll V, ll initial_cost = 0) {
        vecl dp(V, INF);
        dp[start] = initial_cost;

        priority_queue<pairl,vector<pairl>,greater<pairl>> PQ;
        PQ.emplace(dp[start],start);

        vecl update_list(V+10,-1); // update_list[i]: 頂点iを更新した辺番号, 最初は-1
        while (!PQ.empty()) {
            auto [_,current] = PQ.top(); PQ.pop();

            for (auto e : get_edges(current)) {
                if (!chmin(dp[e.to],dp[current] + e.cost)) continue;
                PQ.emplace(dp[e.to],e.to);
                update_list[e.to] = e.idx;
            }
        }

        return ShortestPath {dp, update_list, false};
    }

    // プリム法: O(ElogV), 辺番号の格納されたuset<ll>を返す
    uset<ll> prim(ll V) {
        using P = tuple<ll,ll,ll>;
        priority_queue<P,vector<P>,greater<P>> PQ;
        PQ.push({-1,0,-1}); // 頂点0からスタート

        vecl used(V+10,0LL);
        uset<ll> E;
        ll count = 0;
        while (count < V && !PQ.empty()) {
            auto [cost,current,es] = PQ.top(); PQ.pop();

            if (used[current]) continue;
            for (auto e : get_edges(current)) {
                if (e.idx == es) continue;
                PQ.push({e.cost,e.to,e.idx});
            }

            if (es > -1) E.insert(es);
            used[current] = true;
            count++;
        }

        if (count < V) return uset<ll>(); // 使われてないものがあるなら全域木構成できず
        
        return E; 
    }

    // DFSで適当な全域木を取ってくる {root, Graph}
    pair<ll,Graph> spanning(ll s) { 
        Graph g2(size());
        vecl deg(size(),0LL);
        stack<ll> S;
        uset<ll> seen;
        S.push(s), seen.insert(s);

        while (!S.empty()) {
            ll current = S.top(); S.pop();
            ll count = 0;

            for (auto e : get_edges(current)) {
                if (seen.count(e.to) > 0) continue;
                g2.add_edge(current, e.to, e.cost);
                S.push(e.to);
                deg[e.to]++; // from
                count++;
                seen.insert(e.to);
            }

            deg[current] += count; // to
        }

        rep(i,size()) if (deg[i] == 1) return {i,g2};
        return {0,g2};
    }
};



// #PORT_END#