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
// name: "test"
// prefix: "testcase"
// description: "testcase"
int judgecheck = 0, greedycheck = 1;
#ifdef DEBUG_

int solve(ostringstream &cout);
int greedy(ostringstream &cout);

enum casetype { 
    smallcorner = 0,
    bigcorner = 1,
    normal = 2,
};

using Range = pairl; // [m,M]

struct Param {
    using uid = uniform_int_distribution<ll>;
    
    ll m,M;
    mt19937_64 rnd;
    uid dist;
    Param(Range range) : m(range.first), M(range.second) { // m <= x <= M
        random_device seed_gen;
        rnd = mt19937_64(seed_gen()); 
        dist = uid(m,M);
    }

    ll make(casetype type = normal) {
        if (type == casetype::smallcorner) return m;
        else if (type == casetype::bigcorner) return M;
        return dist(rnd);
    }
};

struct ParamGen {
    int N;
    vector<Param> params;
    ParamGen(vector<Range> ranges) : N(sz(ranges)) {
        rep(i,sz(ranges)) params.push_back(Param(ranges[i]));
    }

    vector<vecl> make() { // O(3^N)
        vector<vecl> res;
        rep(i,pow(N,3)) {
            vecl sub;
            int _i = i;
            rep(j,N) {
                sub.push_back(params[j].make(casetype(_i % 3)));
                _i /= 3;
            }
            res.push_back(sub);
        }
        return res;
    }
};

enum arraytype {
    permutation = 0,
    randomness = 1,
};

struct ArrayGen {
    Range range;
    arraytype type;

    ArrayGen(Range range, arraytype type = randomness) : range(range), type(type) {}

    vecl make(int size) {
        vecl X(size,0LL);
        if (type == arraytype::permutation) iota(ALL(X),range.first);
        else rep(i,size) X[i] = Param(range).make();

        // shuffle
        random_device grd;
        mt19937 grm(grd());
        shuffle(ALL(X),grm);

        return X;
    }
};


struct StressTester { // stress_test
    enum checktype {
        judgec = 0,
        greedyc = 1,
    };

    StressTester() = default;

    void print_testcase() {
        ifstream rf("../../Atcoder/input.txt");
        while (!rf.eof()) {
            string str;
            std::getline(rf, str);
            cout << str << endl;
        }
        rf.close();
    }

    pair<string,string> get_answer() {
        ostringstream oss,oss_greedy;
        solve(oss), greedy(oss_greedy);

        string ans = oss.str();
        string ans_greedy = oss_greedy.str();
        return {ans,ans_greedy};
    }

    bool greedy_check(string ans, string ans_greedy) { // solve vs greedy
        return ans == ans_greedy;
    }

    bool judge_check(ll N) { // solve vs judge
        ostringstream oss;
        solve(oss);
        istringstream cin(oss.str());

        return true;
    }

    void run(checktype ctype = checktype::greedyc, bool pausable = false) { // 入力データを書き換えるので注意
        vector<Range> ranges;
        ranges.emplace_back(1,5); // N
        ranges.emplace_back(1,5); // K

        ParamGen pg(ranges);
        int count = 0;
        int pass = 0;
        for (auto params : pg.make()) {
            // gen params
            ll N = params[0];
            ArrayGen ag({0,N-1},arraytype::randomness);
            auto A = ag.make(N);

            // write testcases 
            ofstream wf; wf.open("../../Atcoder/input.txt");
            wf << N << endl;
            rep(i,N) wf << A[i]+1 << " \n"[i==N-1];
            wf.close(); 

            // solve vs greedy
            auto [ans,ans_greedy] = get_answer();
            bool is_ok = ctype == checktype::greedyc ? greedy_check(ans,ans_greedy) : judge_check(N);
            cout << (count++) << "... ";
            if (is_ok) {
                cout << "pass!" << endl;
                pass++;
            }
            else {
                cout << "fail." << endl;
                cout << "TestCase: " << endl;
                print_testcase();
                
                cout << endl << "ans = " << ans << endl;
                if (ctype == checktype::greedyc) cout << "ans_greedy = " << ans_greedy << endl;
                cout << "------------" << endl;
                if (pausable) break;
            }
        }

        if (pausable && count != pass) {
            cout << "fail: 1" << endl;
            return;
        }

        cout << endl << "---- results ----" << endl;
        cout << "pass: " << pass << endl;
        cout << "fail: " << count - pass << endl;
        cout << "---" << endl;
        cout << "all: " << count << endl;
    }
};

void stress_test(int ctype = greedycheck, bool pausable = false) {
    auto st = StressTester();
    auto ct = StressTester::checktype::greedyc;
    if (ctype == judgecheck) ct = StressTester::checktype::judgec;
    st.run(ct, pausable);
}

#endif

#ifndef DEBUG_
void stress_test(int ctype = greedycheck, bool pausable = false) {}
#endif

// #PORT_END#