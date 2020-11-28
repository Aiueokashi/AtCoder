#include<bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace __gnu_pbds;
using namespace std;
#define ll long long
#define endl "\n"
#define f first
#define s second
#define ar array
#define pb push_back
#define eb emplace_back
#define mp make_pair
#define sz(X) ((int)(X).size())
#define rsz resize
#define pcnt __builtin_popcount
#define sort_unique(c) (sort(c.begin(),c.end()), c.resize(distance(c.begin(),unique(c.begin(),c.end()))))
#define get_pos(c, x) (lower_bound(c.begin(),c.end(),x)-c.begin())
#define all(X) (X).begin(), (X).end()
#define rall(X) (X).rbegin(), (X).rend()
#define ms(c, x)  memset(c,x,sizeof c)
#define  No(x) cout<<(x?"NO":"No")<<endl;
#define  Yes(x) cout<<(x?"YES":"Yes")<<endl;
#define  nl cout<<endl;
#define forn(i, n) for(int i = 0; i < int(n); i++)
#define fore(i, l, r) for(int i = l; i <r; i++)
#define fored(i, l, r) for(int i = r-1; i >=l; i--)
#define ford(i, n) for (int i = n - 1; i >= 0; --i)
using ld = long double;
using db = double;
using str = string;
using pi = pair<int, int>;
using pl = pair<ll, ll>;
using pd = pair<db, db>;

using vi = vector<int>;
using vb = vector<bool>;
using vl = vector<ll>;
using vd = vector<db>;
using vs = vector<str>;
using vpi = vector<pi>;
using vpl = vector<pl>;
using vpd = vector<pd>;

template<typename T>
using osl = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
template<typename T>
using osg = tree<T, null_type, greater<T>, rb_tree_tag, tree_order_statistics_node_update>;
// for greater use greater<T>
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());

template<typename T>
ostream &operator+(ostream &out, const vector<T> &vec) {
    for (const auto &x : vec) {
        out << x << " ";
    }
    out << "\n";
    return out;
}

template<typename T>
ostream &operator*(ostream &out, const vector<T> &vec) {
    for (const auto &x : vec) {
        out + x;
    }
    return out;
}

template<typename T>
istream &operator>>(istream &in, vector<T> &vec) {
    for (auto &x : vec) {
        in >> x;
    }
    return in;
}

template<class T>
bool ckmin(T &a, const T &b) { return b < a ? a = b, 1 : 0; }

template<class T>
bool ckmax(T &a, const T &b) { return a < b ? a = b, 1 : 0; }

#ifdef LOCAL

void debug() { cerr << endl; }

template<class T, class... Args>
void debug(const T &t, const Args &... args) {
    cerr << t << " ";
    debug(args...);
}

#define dbg(...) {cerr<<__LINE__<<" [[DEBUG]] ";debug(__VA_ARGS__);};
#else
#define dbg(...) void(0)
#endif
const ll mod = 1e9 + 7;
//const int mod = 998244353;
const ll INF = 1e17 + 6;

inline ll power(ll x, ll y) {
    ll res = 1;
    x = x;
    while (y > 0) {
        if (y & 1)
            res = (res * x);
        y = y >> 1;
        x = (x * x);
    }
    return res;
}

int dx[4] = {0, 1, 0, -1};
int dy[4] = {1, 0, -1, 0};
const int MXN = 2e5 + 5;

void solve() {
    int n;
    cin>>n;
    vector<vi>a(n,vi(n));
    forn(i,n){

        string s;
        cin>>s;
        forn(j,n){
            a[i][j]=s[j]-'0';
        }
        a[i][i]=1;
    }
    forn(k,n)forn(i,n)forn(j,n){
        a[i][j]|=a[i][k]&&a[k][j];
    }
    long double ans=0;
    forn(i,n){
        int d=0;
        forn(j,n){
            d+=a[j][i];
        }
        ans+=(long  double)1/d;
    }
    cout<<ans<<endl;
}

int main() {

    ios_base::sync_with_stdio(0);
    cin.tie(0);
    cout.tie(0);
    cout << fixed << setprecision(10);
#ifdef LOCAL
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
#endif
    int t = 1, tc = 1;
//    cin >> t;
    while (t--) {
//        cout<<"Case #"<<tc<<": " ;
        solve();
        tc++;
    }
#ifdef LOCAL
    cerr << endl << "Time elapsed : " << clock() * 1000.0 / CLOCKS_PER_SEC << " ms" << '\n';
#endif
    return 0;
}
//look if it requires ll