#include <bits/stdc++.h>
using namespace std;
#define LOCAL // 提出時はコメントアウト
#define DEBUG_

typedef long long ll;
const double EPS = 1e-9;
const ll INF = ((1LL<<62)-(1LL<<31));
typedef vector<ll> vecl;
typedef pair<ll, ll> pairl;
template<typename T, typename U> using mapv = map<T,vector<U>>;

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
// name: "run"
// prefix: "run"
// description: "ランレングス圧縮"

template<class T>
struct Gyu { 
    T item;
    ll count;
};

template<class T>
vector<Gyu<T>> rle(const vector<T> &X) { // ランレングス圧縮
    T pre = X[0];
    ll current_count = 1;
    vector<Gyu<T>> res;
    rep(i,X.size()-1) {
        T x = X[i+1];
        if (x == pre) {
            current_count++;
        }
        else {
            res.push_back({pre,current_count});
            current_count = 1;
        }
        pre = x;
    }

    res.push_back({pre,current_count});
    return res;
}

template<class T>
vector<T> un_rle(const vector<Gyu<T>> &X) { // 展開
    vector<T> res;
    for (auto &x : X) rep(i,x.count) res.push_back(x.item);
    return res;
}

// #PORT_END#


// #PORT#
// name: "zarts"
// prefix: "zarts_pure"
// description: "座圧単体"

ll zarts(const vector<ll> &X, ll x) {
    return lower_bound(ALL(X), x) - X.begin();
}

// #PORT_END#

// #PORT#
// name: "zarts"
// prefix: "zarts"
// description: "座圧"


//// 座標圧縮

template<class T, class U>
struct Zarted {
    unordered_map<T,ll> table; // table[x] → z
    vector<U> Z; 
};

template<class T, class U>
U id_convert(T t) { return t; } 

// 座標圧縮 & T → U 
// Xの位置情報は固定したままXを変換してZに格納する: Z[i] = convert(X[j])
// Xはコピーされる
template<class T = ll, class U = ll, class Compare = less<T>>
Zarted<T,U> zarts(vector<T> X, bool already_sorted, bool already_unique, function<U(T)> convert=id_convert<T,U>) { 
    Compare sort_algorism = Compare();
    umap<T,ll> table;
    
    if (!already_sorted) sort(ALL(X),sort_algorism);
    if (!already_unique) X.erase(unique(ALL(X)), X.end());

    vector<U> Z(sz(X));
    rep(i,sz(X)) {
        T x = X[i];
        ll z = lower_bound(ALL(X), x, sort_algorism) - X.begin();
        table[x] = z; // x → z
        Z[z] = convert(x); // Z[z] = ?
    }
    
    return Zarted<T,U>{table,Z};
}

// #PORT_END#