#include <bits/stdc++.h>
using namespace std;

// #PORT#
// name: "Mo"
// prefix: "mo"
// description: "Mo"

// クエリ先読み可能で、要素が更新されないならOK
struct Mo { // O(k(N+Q)√N)
    int n;
    vector<pair<int, int> > lr;

    explicit Mo(int n) : n(n) {}

    // クエリの半開区間[l,r)をaddする
    void add(int l, int r) {
        lr.emplace_back(l, r);
    }

    // sqrt(N)の区間に分割して、クエリを左端の区間場所で昇順にソート
    // 区間を伸縮させて、区間内の点を順になめていく
    // 点の到達をadd(position), erase(position)で伝播していって、クエリ完遂時にはout(query_id)を呼ぶ
    template <typename AL, typename AR, typename EL, typename ER, typename O>
    void build(const AL &add_left, const AR &add_right, const EL &erase_left,
               const ER &erase_right, const O &out) {
        int q = (int)lr.size();
        int bs = n / min<int>(n, sqrt(q));
        vector<int> ord(q);
        iota(begin(ord), end(ord), 0);
        sort(begin(ord), end(ord), [&](int a, int b) {
            int ablock = lr[a].first / bs, bblock = lr[b].first / bs;
            if (ablock != bblock) return ablock < bblock;
            return (ablock & 1) ? lr[a].second > lr[b].second
                                : lr[a].second < lr[b].second;
        });
        int l = 0, r = 0;
        for (auto idx : ord) {
            while (l > lr[idx].first) add_left(--l);
            while (r < lr[idx].second) add_right(r++);
            while (l < lr[idx].first) erase_left(l++);
            while (r > lr[idx].second) erase_right(--r);
            out(idx);
        }
    }

    template <typename A, typename E, typename O>
    void build(const A &add, const E &erase, const O &out) {
        build(add, add, erase, erase, out);
    }
};

// #PORT_END#