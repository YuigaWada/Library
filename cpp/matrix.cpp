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
// name: "Matrix"
// prefix: "matrix"
// description: "行列演算"


template <class T>
struct Matrix {
    vector<vector<T>> A;

    Matrix() {}

    Matrix(size_t n, size_t m) : A(n, vector<T>(m, 0)) {}

    Matrix(size_t n) : A(n, vector<T>(n, 0)){};

    size_t height() const { return (A.size()); }

    size_t width() const { return (A[0].size()); }

    inline const vector<T> &operator[](int k) const { return (A.at(k)); }

    inline vector<T> &operator[](int k) { return (A.at(k)); }

    static Matrix I(size_t n) {
        Matrix mat(n);
        for (int i = 0; i < n; i++) mat[i][i] = 1;
        return (mat);
    }

    Matrix &operator+=(const Matrix &B) {
        size_t n = height(), m = width();
        assert(n == B.height() && m == B.width());
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) (*this)[i][j] += B[i][j];
        return (*this);
    }

    Matrix &operator-=(const Matrix &B) {
        size_t n = height(), m = width();
        assert(n == B.height() && m == B.width());
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) (*this)[i][j] -= B[i][j];
        return (*this);
    }

    Matrix &operator*=(const Matrix &B) {
        size_t n = height(), m = B.width(), p = width();
        assert(p == B.height());
        vector<vector<T>> C(n, vector<T>(m, 0));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                for (int k = 0; k < p; k++)
                    C[i][j] = (C[i][j] + (*this)[i][k] * B[k][j]);
        A.swap(C);
        return (*this);
    }

    Matrix &operator^=(long long k) {
        Matrix B = Matrix::I(height());
        while (k > 0) {
            if (k & 1) B *= *this;
            *this *= *this;
            k >>= 1LL;
        }
        A.swap(B.A);
        return (*this);
    }

    Matrix operator+(const Matrix &B) const { return (Matrix(*this) += B); }

    Matrix operator-(const Matrix &B) const { return (Matrix(*this) -= B); }

    Matrix operator*(const Matrix &B) const { return (Matrix(*this) *= B); }

    Matrix operator^(const long long k) const { return (Matrix(*this) ^= k); }

    friend ostream &operator<<(ostream &os, Matrix &p) {
        size_t n = p.height(), m = p.width();
        for (int i = 0; i < n; i++) {
            os << "[";
            for (int j = 0; j < m; j++) {
                os << p[i][j] << (j + 1 == m ? "]\n" : ",");
            }
        }
        return (os);
    }

    Matrix pow(ll n) {
        Matrix<T> B = I(height());
        Matrix<T> C(width(), height());
        C.A = A;

        while (n > 0) {
            if (n & 1) B *= C;
            C *= C;
            n >>= 1;
        }

        return B;
    }

    T determinant() {
        Matrix B(*this);
        assert(width() == height());
        T ret = 1;
        for (int i = 0; i < width(); i++) {
            int idx = -1;
            for (int j = i; j < width(); j++) {
                if (B[j][i] != 0) idx = j;
            }
            if (idx == -1) return (0);
            if (i != idx) {
                ret *= -1;
                swap(B[i], B[idx]);
            }
            ret *= B[i][i];
            T vv = B[i][i];
            for (int j = 0; j < width(); j++) {
                B[i][j] /= vv;
            }
            for (int j = i + 1; j < width(); j++) {
                T a = B[j][i];
                for (int k = 0; k < width(); k++) {
                    B[j][k] -= B[i][k] * a;
                }
            }
        }
        return (ret);
    }
};


// #PORT_END#

// #PORT#
// name: "vec"
// prefix: "vec"
// description: "3次元ベクトル"

double pi = 3.141592653589793238;
struct vec {
    using ld = long double;
    ld x, y, z;
    vec(ld x = 0, ld y = 0, ld z = 0) : x(x), y(y), z(z) {}

    ld dot(const vec &q) { // dot product
        return x * q.x + y * q.y + z * q.z;
    }  

    vec cross(const vec &q) { // cross product = p × q
        vec res;
        res.x = y * q.z - q.y - z;
        res.y = -x * q.z + q.x - z;
        res.z = x * q.y - q.x - y;
        return res;
    }

    vec rot2D(ld rad) { // xy平面で回転
        assert(z == 0);
        ld c = cos(rad), s = sin(rad);
        return vec {c*x-s*y,s*x+c*y};
    }

    ld norm() { // |p|^2
        return x * x + y * y + z * z;
    }

    ld abs() { // |p|
        return sqrt(norm());
    }

    vec operator+(const vec &q) {
        return vec {x+q.x, y+q.y, z+q.z};
    }

    vec operator-(const vec &q) {
        return vec {x-q.x, y-q.y, z-q.z};
    }   

    vec operator*(ld c) {
        return vec {x*c,y*c,z*c};
    }
    
    vec operator/(ld c) {
        return vec {x/c,y/c,z/c};
    }
    
    friend ostream &operator<<(ostream &os, const vec &p) {
        return os << "(" << p.x << "," << p.y << "," << p.z << ")";
    }

};

// #PORT_END#