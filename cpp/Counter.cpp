#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
#define REP(i, x, n) for (int i = x; i < n; i++)
#define rep(i, n) REP(i, 0, n)

// #PORT#
// name: "counter"
// prefix: "counter"
// description: "Counter"

template<typename T>
class Counter {
  public:
    unordered_map<T,ll> memo;
    void add(T x) {
      if (memo.count(x) == 0) memo[x] = 0;
      memo[x]++;
    }

    void decrement(T x) {
      memo[x]--;
      if (memo[x]<=0) memo.erase(x);
    }
};

// #PORT_END#