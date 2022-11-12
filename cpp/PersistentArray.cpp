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
#define pb(x) push_back(x)
#define eb(x) emplace_back(x)
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

#ifndef DEBUG_
#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#endif


// #PORT#
// name: "persistent"
// prefix: "persistent"
// description: "永続配列"

template <typename T>
struct PersistentArray {
    struct node_t;
    using node_ptr = node_t*;
    struct node_t {
        T data;
        node_ptr left = nullptr;
        node_ptr right = nullptr;
 
        node_t(const T& _data, const node_ptr& _left, const node_ptr& _right)
            : data(_data), left(_left), right(_right) {}
    };
 
    vector<pair<node_ptr, int>> m_roots;
 
    explicit PersistentArray(int _size, const T& _init = {})
        : m_roots({{new node_t(_init, nullptr, nullptr), _size}}) {
        auto dfs = [&](auto& F, node_ptr& parent, int pidx, int pdep) -> void {
            if (pidx + (1 << pdep) < _size) {
                parent->left = new node_t(_init, nullptr, nullptr);
                F(F, parent->left, pidx + (1 << pdep), pdep + 1);
            }
            if (pidx + (1 << (pdep + 1)) < _size) {
                parent->right = new node_t(_init, nullptr, nullptr);
                F(F, parent->right, pidx + (1 << (pdep + 1)), pdep + 1);
            }
        };
 
        dfs(dfs, m_roots[0].first, 0, 0);
    }
 
    explicit PersistentArray(const vector<T>& _vec)
        : m_roots({{new node_t(_vec.empty() ? T{} : _vec[0], nullptr, nullptr), _vec.size()}}) {
        int _size = static_cast<int>(_vec.size());
 
        auto dfs = [&](auto F, node_ptr& parent, int pidx, int pdep) -> void {
            if (int nindex = pidx + (1 << pdep); nindex < _size) {
                parent->left = new node_t(_vec[nindex], nullptr, nullptr);
                F(F, parent->left, nindex, pdep + 1);
            }
            if (int nindex = pidx + (1 << (pdep + 1)); nindex < _size) {
                parent->right = new node_t(_vec[nindex], nullptr, nullptr);
                F(F, parent->right, nindex, pdep + 1);
            }
        };
 
        dfs(dfs, m_roots[0].first, 0, 0);
    }
 
    int size(int time) const { return m_roots[time].second; }
 
    bool empty(int time) const { return m_roots[time].second == 0; }
 
    int now() const { return static_cast<int>(m_roots.size()) - 1; }
 
    int set(int time, int index, const T& val) { // 時刻timeにおける配列を格納する
        auto dfs = [&](auto& F, node_ptr& node) -> node_ptr {
            if (index == 0) {
                return new node_t(val, node->left, node->right);
            }
 
            if (index & 1) {
                index = (index >> 1);
                return new node_t(node->data, F(F, node->left), node->right);
            } else {
                index = ((index - 1) >> 1);
                return new node_t(node->data, node->left, F(F, node->right));
            }
        };
 
        m_roots.emplace_back(dfs(dfs, m_roots[time].first), m_roots[time].second);
        return now();
    }
 
    T get(int time, int index) const { // 時刻timeにおける配列を返す
        T ret;
 
        auto dfs = [&](const auto& F, const node_ptr& node) -> void {
            if (index == 0) {
                ret = node->data;
                return;
            }
 
            if (index & 1) {
                index = (index >> 1);
                F(F, node->left);
            } else {
                index = ((index - 1) >> 1);
                F(F, node->right);
            }
        };
 
        dfs(dfs, m_roots[time].first);
        return ret;
    }
  
    template <class... Args>
    int emplace_back(int time, Args&&... args) { // argsを格納する・timeを返す
        int index = m_roots[time].second;
 
        auto dfs = [&](auto& F, node_ptr& node) -> node_ptr {
            if (index == 0) return new node_t(T(forward<Args>(args)...), nullptr, nullptr);
            else if (index & 1) {
                index = (index >> 1);
                return new node_t(node->data, F(F, node->left), node->right);
            } else {
                index = ((index - 1) >> 1);
                return new node_t(node->data, node->left, F(F, node->right));
            }
        };
 
        m_roots.emplace_back(dfs(dfs, m_roots[time].first), m_roots[time].second + 1);
        return now();
    }
 
    int pop_back(int time) { // timeを返す
        m_roots.emplace_back(m_roots[time].first, m_roots[time].second - 1);
        return now();
    }

    T tail(int time) const { return get(time, m_roots[time].second - 1); } // 時刻timeの最後尾を返す
 
    T front(int time) const { return m_roots[time].first->data; } // 時刻timeの先頭を返す
};
// #PORT_END#



int solve(ostringstream &cout) {
    #ifdef LOCAL
    ifstream in("../../input.txt");
    cin.rdbuf(in.rdbuf());
    #endif

    ll Q;
    cin>>Q;

    PersistentArray<ll> pa(0);
    umap<ll,ll> note;
    ll time = pa.now();
    rep(q,Q) {
        string S;
        cin>>S;
        char s = S[0];
        if (s == 'A') {
            ll x; cin>>x;
            time = pa.emplace_back(time,x);
        }
        else if (s == 'D') {
            if (!pa.empty(time)) {
                time = pa.pop_back(time);
            }
        }
        else if (s == 'S') {
            ll x; cin>>x;
            note[x] = time;
        }
        else if (s == 'L') {
            ll x; cin>>x;
            time = note[x];
        }

        cout << (pa.empty(time) ? -1 : pa.tail(time)) << " ";
    }
    cout << endl;

    return 0;
}

int main() {
    ostringstream oss;
    int res = solve(oss);

    cout << oss.str();
    return res;
}




