#include<iostream>
#include<string>
#include<algorithm>
#include<vector>
#include<queue>
#include<map>
#include<math.h>
#include<iomanip>
#include<set>
#include<numeric>
#include<cstring>
#include<cstdio>
#include<functional>
#include<bitset>
#include<limits.h>
#include<cassert>
#include<iterator>
#include<complex>
#include<stack>
#include<unordered_map>
#include<unordered_set>
#include<time.h>
#include<random>
#include<array>
using namespace std;
using ll = long long;
using ull = unsigned long long;
#define rep(i, a, b) for(int i = a; i < b; i++)
#define rrep(i, a, b) for(int i = b - 1; i >= a; i--)
#define repl(i, a, b) for(long long i = a; i < b; i++)
#define rrepl(i, a, b) for(long long i = b - 1; i >= a; i--)
#define ALL(a) a.begin(), a.end()
using pii = pair<int,int>;
using piii = pair<pii,int>;
using pll = pair<long long, long long>;
using plll = pair<pll, long long>;
// #pragma GCC optimize("Ofast")
// #define _GLIBCXX_DEBUG
#define pcnt __builtin_popcount
#define buli(x) __builtin_popcountll(x)
#define pb push_back
#define mp make_pair
#define UNIQUE(v) v.erase( unique(v.begin(), v.end()), v.end() );
#define isSquare(x) (sqrt(x)*sqrt(x) == x)
template<class T>inline bool chmax(T &a, const T &b) {if(a<b){a = b; return 1;} return 0; };
template<class T>inline bool chmin(T &a, const T &b) {if(a>b){a = b; return 1;} return 0; };
inline void in(void){return;}
template <typename First, typename... Rest> void in(First& first, Rest&... rest){cin >> first;in(rest...);return;}
inline void out(void){cout << "\n";return;}
template <typename First, typename... Rest> void out(First first, Rest... rest){cout << first << " ";out(rest...);return;}
const double EPS = 1e-9;
const int mod = 1e9 + 7;
// const int mod = 998244353;
const int INF = 1e9;
const long long INFLL = 1e18;
void iosetup() {
    cin.tie(nullptr);ios::sync_with_stdio(false);
    cout << fixed << setprecision(10);
    cerr << fixed << setprecision(10);
}
template< typename T1, typename T2 >
ostream &operator<<(ostream &os, const pair< T1, T2 >& p) {
    os << p.first << " " << p.second;
    return os;
}
template< typename T1, typename T2 >
istream &operator>>(istream &is, pair< T1, T2 > &p) {
    is >> p.first >> p.second;
    return is;
}
template< typename T >
ostream &operator<<(ostream &os, const vector< T > &v) {
    for(int i = 0; i < (int) v.size(); i++) {
        os << v[i] << (i + 1 != v.size() ? " " : "");
    }
    return os;
}
template< typename T >
istream &operator>>(istream &is, vector< T > &v) {
    for(T &in : v) is >> in;
    return is;
}
template<class T> vector<T> make_vec(size_t a) {return vector<T>(a); }
template<class T, class... Ts> auto make_vec(size_t a, Ts... ts){
    return vector<decltype(make_vec<T>(ts...))>(a, make_vec<T>(ts...));
}
template<class S, class T> pair<S,T> operator+(const pair<S,T> &s, const pair<S, T>& t){return pair<S,T>(s.first+t.first, s.second+t.second);}
template<class S, class T> pair<S,T> operator-(const pair<S,T> &s, const pair<S, T>& t){return pair<S,T>(s.first-t.first, s.second-t.second);}
template<class S, class T> pair<S,T> operator*(const pair<S,T> &s, const S& t){return pair<S,T>(s.first*t, s.second*t);}
template <typename T> void Exit(T first){cout << first << endl;exit(0); };
template< int mod > struct ModInt {
    unsigned x; ModInt() : x(0) {}
    ModInt(int64_t y) : x(y >= 0 ? y % mod : (mod - (-y) % mod) % mod) {}
    ModInt &operator+=(const ModInt &p) {if((x += p.x) >= mod) x -= mod;return *this;}
    ModInt &operator-=(const ModInt &p) {if((x += mod - p.x) >= mod) x -= mod;return *this;}
    ModInt &operator*=(const ModInt &p) {x = (int) (1LL * x * p.x % mod);return *this;}
    ModInt &operator/=(const ModInt &p) {*this *= p.inverse();return *this;}
    ModInt operator-() const { return ModInt(-x); }
    ModInt operator+(const ModInt &p) const { return ModInt(*this) += p; }
    ModInt operator-(const ModInt &p) const { return ModInt(*this) -= p; }
    ModInt operator*(const ModInt &p) const { return ModInt(*this) *= p; }
    ModInt operator/(const ModInt &p) const { return ModInt(*this) /= p; }
    bool operator==(const ModInt &p) const { return x == p.x; }
    bool operator!=(const ModInt &p) const { return x != p.x; }
    ModInt inverse() const {int a = x, b = mod, u = 1, v = 0, t;
    while(b > 0) { t = a / b; swap(a -= t * b, b); swap(u -= t * v, v); }return ModInt(u);}
    ModInt pow(int64_t n) const {ModInt ret(1), mul(x); while(n > 0) {if(n & 1) ret *= mul;mul *= mul;n >>= 1;}return ret;}
    friend ostream &operator<<(ostream &os, const ModInt &p) { return os << p.x;}
    friend istream &operator>>(istream &is, ModInt &a) { int64_t t; is >> t; a = ModInt< mod >(t); return (is); }
    static int get_mod() { return mod; }
}; using modint = ModInt< mod >;
const int dx[4] = {1, 0, -1, 0};
const int dy[4] = {0, 1, 0, -1};
const pii dxy[4] = {pii(1,0), pii(0, 1), pii(-1, 0), pii(0, -1)};
 
//----------------------- edit from here ---------------------------------
 
template< typename T > struct edge {
   int src, to; T cost;
   edge(int to, T cost) : src(-1), to(to), cost(cost) {}
   edge(int src, int to, T cost) : src(src), to(to), cost(cost) {}
   edge &operator=(const int &x) { to = x; return *this;}
   operator int() const { return to; }
};
template< typename T > using Edges = vector< edge< T > >;
template< typename T > using WeightedGraph = vector< Edges< T > >;
using UnWeightedGraph = vector< vector< int > >;
template< typename T > using Matrix = vector< vector< T > >;
void add_edge(UnWeightedGraph& g, int x, int y){
    g[x].push_back(y);
    g[y].push_back(x);
}
template< typename T > 
void add_edge(WeightedGraph< T >& g, int x, int y, T c){
    g[x].push_back({y, c});
    g[y].push_back({x, c});
}
void add_edge_bi(UnWeightedGraph& g, int x, int y){
    g[x].push_back(y);
}
template< typename T > 
void add_edge_bi(WeightedGraph< T >& g, int x, int y, T c){
    g[x].push_back({y, c});
}
 
template< typename T >
void warshall_floyd(Matrix< T > &g, T INF) {
  for(int k = 0; k < g.size(); k++) {
    for(int i = 0; i < g.size(); i++) {
      for(int j = 0; j < g.size(); j++) {
        if(g[i][k] == INF || g[k][j] == INF) continue;
        g[i][j] = min(g[i][j], g[i][k] + g[k][j]);
      }
    }
  }
}
// Matrix< int > mat(V, vector< int >(V, INT_MAX));
// warshall_floyd(mat, INT_MAX);
int main(){
    iosetup();
    int n; cin >> n;
    vector<string> S(n); cin >> S;
    Matrix< int > mat(n, vector< int >(n, INT_MAX));
    rep(i, 0, n) mat[i][i] = 0;
    rep(i, 0, n) rep(j, 0, n) if(S[i][j] == '1') mat[i][j] = 0;
    warshall_floyd(mat, INT_MAX);
    long double ans = 0;
    vector<int> cnt(n, 0);
    rep(i, 0, n){
        rep(j, 0, n) if(mat[i][j] == 0) cnt[j]++;
    }
    // cerr << cnt << endl;
    rep(i, 0, n) ans += 1.0 / cnt[i];
    cout << ans << endl;
 
    return 0;
}