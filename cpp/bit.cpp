#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
#define REP(i, x, n) for (int i = x; i < n; i++)
#define rep(i, n) REP(i, 0, n)

// #PORT#
// name: "msb"
// prefix: "msb"
// description: "最上位bitの位置を返す"

// 最上位bitの位置を返す
ll msb(ll x) {
    return 63 - __builtin_clzll(x); // (0-index)
}

// #PORT_END#


// #PORT#
// name: "bit"
// prefix: "bit"
// description: "bit全探索"

// 0~n-1でbit全探索
vector<vecl> bit_search(ll n) {
    vector<vecl> res;
    rep(i,pow(2,n)) {
        vecl st;
        rep(j,n) {
            if ((i>>j)&1) 
                st.push_back(j);
        }
        res.push_back(st);
    }

    return res;
}


// #PORT_END#

// #PORT#
// name: "binary_indexed_tree"
// prefix: "bit_binary_indexed_tree"
// description: "Binary Indexed Tree"

template <class T>
struct Bit { // Tによっては演算変えなきゃいけないので注意
    ll size;
    vector<T> v;
    vector<T> X;

    Bit(ll n): size(n), v(n+1,T(0)), X(n+1,T(0)) {}
    Bit(vector<T> vec, ll N = -1) {
        if (N == -1) size = (ll)vec.size() + 1;
        else size = N;

        v.assign(size,T(0)), X.assign(size,T(0));
        rep(i,(ll)vec.size()) add(i,vec[i]);
    }

    T get_sum(ll index) { // [0,index): O(logN)
        T res = T(0);
        for (int i = index; i > 0; i -= (i & -i)) { // (i & -i): LSB
            res += v[i];
        }
        return res;
    }

    T get(ll l, ll r) { // [l,r): O(logN)
        return get_sum(r) - get_sum(l);
    }

    T get(ll index) {
        return get(index,index+1);
    }

    void add(ll index, T x) { // O(logN)
        X[index] += x;
        index++;
        for (int i = index; i < size + 1; i += (i & -i)) { // (i & -i): LSB
            v[i] += x;
        }
    }

    void sub(ll index, T x) {  // O(logN)
        add(index, -x);
    }

    void update(ll index, T x) {  // O(logN)
        T old = get(index);
        add(index, -old + x);
    }

    ll lower_bound(ll x) { // min{i} s.t. ∑a >= x: O(logN)
        if (x <= 0) return 0;
        ll res = 0, r = 1;

        while(r < size) r <<= 1; // sizeに最も近い2の冪乗を探す
        for (int len = r; len > 0; len >>= 1) { // 区間幅ごとに探索
            if (res + len >= size) continue;
            if(v[res+len] < x){
                x -= v[res+len];
                res += len;
            }
        }
        return res;
    }

    friend ostream &operator<<(ostream &os, const Bit<T> &p) {
        return os << p.X;
    }

};



// #PORT_END#