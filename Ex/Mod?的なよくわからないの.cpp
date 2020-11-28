/* #region header */
#ifdef LOCAL
#include "cxx-prettyprint-master/prettyprint.hpp"
#define debug(x) cout << x << endl
#else
#define debug(...) 42
#endif

#pragma GCC optimize("Ofast")
#include <bits/stdc++.h>
using namespace std;
//#include <atcoder/all>
// types
using ll = long long;
using ull = unsigned long long;
using ld = long double;
typedef pair<ll, ll> Pl;
typedef pair<int, int> Pi;
typedef vector<ll> vl;
typedef vector<int> vi;
typedef vector<char> vc;
template <typename T>
using mat = vector<vector<T>>;
typedef vector<vector<int>> vvi;
typedef vector<vector<long long>> vvl;
typedef vector<vector<char>> vvc;
template <int mod>
struct modint {
    int x;

    modint() : x(0) {}

    modint(int64_t y) : x(y >= 0 ? y % mod : (mod - (-y) % mod) % mod) {}

    modint& operator+=(const modint& p) {
        if ((x += p.x) >= mod) x -= mod;
        return *this;
    }

    modint& operator-=(const modint& p) {
        if ((x += mod - p.x) >= mod) x -= mod;
        return *this;
    }

    modint& operator*=(const modint& p) {
        x = (int)(1LL * x * p.x % mod);
        return *this;
    }

    modint& operator/=(const modint& p) {
        *this *= p.inverse();
        return *this;
    }

    modint operator-() const { return modint(-x); }

    modint operator+(const modint& p) const { return modint(*this) += p; }

    modint operator-(const modint& p) const { return modint(*this) -= p; }

    modint operator*(const modint& p) const { return modint(*this) *= p; }

    modint operator/(const modint& p) const { return modint(*this) /= p; }

    bool operator==(const modint& p) const { return x == p.x; }

    bool operator!=(const modint& p) const { return x != p.x; }

    modint inverse() const {
        int a = x, b = mod, u = 1, v = 0, t;
        while (b > 0) {
            t = a / b;
            swap(a -= t * b, b);
            swap(u -= t * v, v);
        }
        return modint(u);
    }

    modint pow(int64_t n) const {
        modint ret(1), mul(x);
        while (n > 0) {
            if (n & 1) ret *= mul;
            mul *= mul;
            n >>= 1;
        }
        return ret;
    }

    friend ostream& operator<<(ostream& os, const modint& p) {
        return os << p.x;
    }

    friend istream& operator>>(istream& is, modint& a) {
        int64_t t;
        is >> t;
        a = modint<mod>(t);
        return (is);
    }

    static int get_mod() { return mod; }
};
// abreviations
#define all(x) (x).begin(), (x).end()
#define rall(x) (x).rbegin(), (x).rend()
#define rep_(i, a_, b_, a, b, ...) for (ll i = (a), max_i = (b); i < max_i; i++)
#define rep(i, ...) rep_(i, __VA_ARGS__, __VA_ARGS__, 0, __VA_ARGS__)
#define rrep_(i, a_, b_, a, b, ...) \
    for (ll i = (b - 1), min_i = (a); i >= min_i; i--)
#define rrep(i, ...) rrep_(i, __VA_ARGS__, __VA_ARGS__, 0, __VA_ARGS__)
#define srep(i, a, b, c) for (ll i = (a), max_i = (b); i < max_i; i += c)
#define SZ(x) ((int)(x).size())
#define pb(x) push_back(x)
#define eb(x) emplace_back(x)
#define mp make_pair
//入出力
#define print(x) cout << x << endl
template <class T>
ostream& operator<<(ostream& os, const vector<T>& v) {
    for (auto& e : v) cout << e << " ";
    cout << endl;
    return os;
}
void scan(int& a) { cin >> a; }
void scan(long long& a) { cin >> a; }
void scan(char& a) { cin >> a; }
void scan(double& a) { cin >> a; }
void scan(string& a) { cin >> a; }
template <class T>
void scan(vector<T>& a) {
    for (auto& i : a) scan(i);
}
#define vsum(x) accumulate(all(x), 0LL)
#define vmax(a) *max_element(all(a))
#define vmin(a) *min_element(all(a))
#define lb(c, x) distance((c).begin(), lower_bound(all(c), (x)))
#define ub(c, x) distance((c).begin(), upper_bound(all(c), (x)))
// functions
// gcd(0, x) fails.
ll gcd(ll a, ll b) { return b ? gcd(b, a % b) : a; }
ll lcm(ll a, ll b) { return a / gcd(a, b) * b; }
template <class T>
bool chmax(T& a, const T& b) {
    if (a < b) {
        a = b;
        return 1;
    }
    return 0;
}
template <class T>
bool chmin(T& a, const T& b) {
    if (b < a) {
        a = b;
        return 1;
    }
    return 0;
}
template <typename T>
T mypow(T x, ll n) {
    T ret = 1;
    while (n > 0) {
        if (n & 1) (ret *= x);
        (x *= x);
        n >>= 1;
    }
    return ret;
}
ll modpow(ll x, ll n, const ll mod) {
    ll ret = 1;
    while (n > 0) {
        if (n & 1) (ret *= x);
        (x *= x);
        n >>= 1;
        x %= mod;
        ret %= mod;
    }
    return ret;
}
template <typename T>
uint64_t my_rand(void) {
    static uint64_t x = 88172645463325252ULL;
    x = x ^ (x << 13);
    x = x ^ (x >> 7);
    return x = x ^ (x << 17);
}
int popcnt(ull x) { return __builtin_popcountll(x); }
// graph template
template <typename T>
struct edge {
    int src, to;
    T cost;

    edge(int to, T cost) : src(-1), to(to), cost(cost) {}

    edge(int src, int to, T cost) : src(src), to(to), cost(cost) {}

    edge& operator=(const int& x) {
        to = x;
        return *this;
    }

    bool operator<(const edge<T>& r) const { return cost < r.cost; }

    operator int() const { return to; }
};
template <typename T>
using Edges = vector<edge<T>>;
template <typename T>
using WeightedGraph = vector<Edges<T>>;
using UnWeightedGraph = vector<vector<int>>;
struct Timer {
    clock_t start_time;
    void start() { start_time = clock(); }
    int lap() {
        // return x ms.
        return (clock() - start_time) * 1000 / CLOCKS_PER_SEC;
    }
};
/* #endregion*/
// constant
#define inf 4000000000ll
#define INF 4000000004000000000LL
#define mod 998244353ll
using mint = modint<mod>;
typedef vector<mint> vmint;
typedef vector<vector<mint>> vvmint;
#define endl '\n'
const long double eps = 0.000000000000001;
const long double PI = 3.141592653589793;
// library
void solve() {
    int n;
    cin >> n;
    vl a(n);
    scan(a);
    if(n%2==1){
        print("Second");
    }else{
        bool flag = true;
        sort(all(a));
        rep(i, n/2)if(a[i*2]!=a[i*2+1])flag = false;
        if(flag){
            print("Second");
        }else{
            print("First");
        }
    }
}
int main() {
    cin.tie(0);
    ios::sync_with_stdio(0);
    cout << setprecision(30) << fixed;
    int t;
    cin >> t;
    while (t--) {
        solve();
    }
}