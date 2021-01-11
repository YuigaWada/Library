#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
#define REP(i, x, n) for (int i = x; i < n; i++)
#define rep(i, n) REP(i, 0, n)

// #PORT#
// name: "accum"
// prefix: "accum"
// description: "累積和"
template<typename T>
class Accumulation {
    public: 
        vector<T> v;
        Accumulation(const vector<ll> &X,ll N, ll init_val=0) {
            v.assign(N+1,0);
            v[0] = init_val;
            rep(i,N) {
                v[i+1] = v[i] + X[i];
            }
        }

        T query(ll l,ll r) { // [0,r]
            return v[r+1] - v[l];
        }
};

// #PORT_END#


// FOR DEBUG
void printVector(const vector<ll>& vec) {
    for (int value : vec) {
        cout << value << " ";
    }
    cout << endl;
}

int main() {
    vector<ll> v {1,2,3,4};
    auto acc = Accumulation<4>(v);

    auto a = acc.v;
    rep(i,5) {
        cout << a[i] << endl;
    }
    cout << endl;
    

    return 1;
}