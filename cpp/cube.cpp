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
// name: "cube"
// prefix: "cube"
// description: "6面サイコロの回転"

struct Cube {
    using T = tuple<int,int,int,int,int,int>; // u,d,f,b,l,r
    T ed;
    Cube() : ed({1,2,3,4,5,6}) {}
    Cube(T _ed) : ed(_ed) {}
    Cube(vecl _ed) { ed = {_ed[0],_ed[1],_ed[2],_ed[3],_ed[4],_ed[5]}; }

    void rotate(char type) {
        auto [u,d,f,b,l,r] = ed;
        int buff;
        switch (type) {
            case 'N': buff = d, d = f, f = u, u = b, b = buff;
            case 'S': buff = d, d = b, b = u, u = f, f = buff;
            case 'L': buff = f, f = l, l = b, b = r, r = buff;
            case 'R': buff = f, f = r, r = b, b = l, l = buff;
            case 'E': buff = d, d = l, l = u, u = r, r = buff;
            case 'W': buff = d, d = r, r = u, u = l, l = buff;
            default: break;
        }
        ed = {u,d,f,b,l,r};
    }

    vector<Cube> get_all() {
        char procedures[5][2] = {{'N','X'},{'S','X'},{'S','S'},{'L','X'},{'R','X'}};
        vector<Cube> res(1,Cube(ed));
        rep(i,5) {
            rep(j,2) rotate(procedures[i][j]);
            res.push_back(Cube(ed));
        }
        return res;
    }

    vecl make_vector() {
        auto [u,d,f,b,l,r] = ed;
        return vecl{u,d,f,b,l,r};
    }

    bool operator==(Cube &Q) {
        auto S = Q.get_all();
        for (auto cube : this->get_all()) {
            rep(i,sz(S)) {
                if (cube.ed == S[i].ed) return true;
            }
        }
        return false;
    }

    bool operator!=(Cube &Q) { return !(*this == Q); }
};

// #PORT_END#

int solve(ostringstream &cout) {
    #ifdef LOCAL
    ifstream in("../../../Atcoder/input.txt");
    cin.rdbuf(in.rdbuf());
    #endif

    // 和がNであるようなサイコロ面の数え上げ(愚直解)

    ll N;
    cin>>N;

    vector<Cube> S;
    ll n = N-5;

    rep(i,pow(n,6)) {
        vecl X;
        ll x = i, sum = 0;
        rep(j,6) {
            ll y = x % n + 1;
            X.push_back(y), sum += y, x /= n;
        }

        if (sum != N) continue;
        S.push_back(Cube(X));
    }   

    UnionFind uf(sz(S));
    rep(i,sz(S)) REP(j,i+1,sz(S)) if (S[i] == S[j]) uf.unite(i,j);

    cout << sz(uf.lands()) << endl;

    return 0;
}

int main() {
    ostringstream oss;
    int res = solve(oss);

    cout << oss.str();
    return res;
}
