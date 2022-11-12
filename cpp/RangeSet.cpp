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

// #PORT#
// name: "RangeSet"
// prefix: "rangeset"
// description: "区間管理"

struct RangeSet {
    set<pairl> st;

    RangeSet() {
        st.emplace(-INF, -INF + 1);
        st.emplace(INF, INF + 1);
    }
 
    bool covered(const ll l, const ll r) { // [l,r)
        assert(l < r);
        auto itr = prev(st.lower_bound({l + 1, -INF}));
        return itr->first <= l && r <= itr->second;
    }

    bool covered(const ll x) { return covered(x, x + 1); }
 
    pairl covered_by(const ll l, const ll r) { // 存在しなかったら[-INF, -INF)
        assert(l < r);
        auto itr = prev(st.lower_bound({l + 1, -INF}));
        if(itr->first <= l && r <= itr->second) return *itr;
        return {-INF, -INF};
    }
 
    pairl covered_by(const ll x) { return covered_by(x, x + 1); } // 存在しなかったら[-INF, -INF)
 
    ll insert(ll l, ll r) { // [l,r)
        if(l >= r) return 0;
        auto itr = prev(st.lower_bound({l + 1, -INF}));
        if(itr->first <= l && r <= itr->second) return 0;
        ll sum_erased = 0;
        if(itr->first <= l && l <= itr->second) {
            l = itr->first;
            sum_erased += itr->second - itr->first;
            itr = st.erase(itr);
        } else
            itr = next(itr);
        while(r > itr->second) {
            sum_erased += itr->second - itr->first;
            itr = st.erase(itr);
        }
        if(itr->first <= r) {
            sum_erased += itr->second - itr->first;
            r = itr->second;
            st.erase(itr);
        }
        st.emplace(l, r);
        return r - l - sum_erased;
    }
 
    ll insert(const ll x) { return insert(x, x + 1); }
 
    ll erase(const ll l, const ll r) { // [l,r)
        assert(l < r);
        auto itr = prev(st.lower_bound({l + 1, -INF}));
        if(itr->first <= l && r <= itr->second) {
            if(itr->first < l) st.emplace(itr->first, l);
            if(r < itr->second) st.emplace(r, itr->second);
            st.erase(itr);
            return r - l;
        }
 
        ll ret = 0;
        if(itr->first <= l && l < itr->second) {
            ret += itr->second - l;
            if(itr->first < l) st.emplace(itr->first, l);
            itr = st.erase(itr);
        } else
            itr = next(itr);
        while(itr->second <= r) {
            ret += itr->second - itr->first;
            itr = st.erase(itr);
        }
        if(itr->first < r) {
            ret += r - itr->first;
            st.emplace(r, itr->second);
            st.erase(itr);
        }
        return ret;
    }
 
    ll erase(const ll x) { return erase(x, x + 1); }
 
    int size() { return (int)st.size() - 2; }
 
    int mex(const ll x = 0) {
        auto itr = prev(st.lower_bound({x + 1, -INF}));
        if(itr->first <= x && x < itr->second)
            return itr->second;
        else
            return x;
    }
 
    ll sum() const {
        ll res = 0;
        for(auto &p : st) {
            if(abs(p.first) == INF) continue;
            res += p.second - p.first;
        }
        return res;
    }
 
    set<pairl> get() const {
        set<pairl> res;
        for(auto &p : st) {
            if(abs(p.first) == INF) continue;
            res.emplace(p.first, p.second);
        }
        return res;
    }

    friend ostream &operator<<(ostream &os, const RangeSet &rs) {
        for (auto iter = rs.st.begin(); iter != rs.st.end(); iter++) {
            auto [l,r] = *iter;
            if(abs(l) == INF) continue;
            os << "[" << l << "," << r << "),";
        }
        return os << endl;
    }
};

// #PORT_END#


 
int main() {
    #ifdef LOCAL
    ifstream in("../../Atcoder/input.txt");
    cin.rdbuf(in.rdbuf());
    #endif

    RangeSet rs;

    rs.insert(-1,100);
    rs.insert(-20,2);

    dump(rs);
    return 0;
}
