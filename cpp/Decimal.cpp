#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
#define REP(i, x, n) for (int i = x; i < n; i++)
#define rep(i, n) REP(i, 0, n)

// #PORT#
// name: "decimal"
// prefix: "decimal"
// description: "文字列のまま数値比較する"
enum Compared { greater, less, equal };

// 文字列のまま数値比較 O(1)
Compared compare(string x, string y) { // greater = x < y
    if (x.length() != y.length()) {
        if (x.length() < y.length()) return Compared::greater;
        else return Compared::less;
    } 

    // 文字列の大小は辞書順で定義されるので、サイズが同じ場合は直接比較すれば良い
    if (x < y) {  
        return Compared::greater;
    }
    else if (x > y) {
        return Compared::less;
    }
    else {
        return Compared::equal;
    }
}
// #PORT_END#