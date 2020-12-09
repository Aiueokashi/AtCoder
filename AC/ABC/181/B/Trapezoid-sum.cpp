
#pragma GCC optimize("Ofast")
//#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")
//浮動小数点数計算時
#include <bits/stdc++.h>
using namespace std;
//脳筋int定義
#define int long long int
const int inf = 1e9 + 7;
const int MOD = 1000000007;
const double pi = 3.141592653589793238;
using pii = pair<int, int>;
#define ll long long
#define IO                                                                     \
	ios_base::sync_with_stdio(false);                                          \
	cin.tie(0);                                                                \
	cout.tie(0)
// vector定義
using vi = vector<int>;
#define vec1(type, name, ...) vector<type> name(__VA_ARGS__)
#define VEC1(type, ...) vector<type>(__VA_ARGS__)
#define vec2(type, name, a, ...)                                               \
	vector<vector<type>> name(a, VEC1(type, __VA_ARGS__))
#define VEC2(type, a, ...) vector<vector<type>>(a, VEC1(type, __VA_ARGS__))
#define vec3(type, name, a, b, ...)                                            \
	vector<vector<vector<type>>> name(a, VEC2(type, b, __VA_ARGS__))
#define VEC3(type, a, b, ...)                                                  \
	vector<vector<vector<type>>>(a, VEC2(type, b, __VA_ARGS__))
#define vec4(type, name, a, b, c, ...)                                         \
	vector<vector<vector<vector<type>>>> name(a, VEC3(type, b, c, __VA_ARGS__))
#define VEC4(type, a, b, c, ...)                                               \
	vector<vector<vector<vector<type>>>>(a, VEC3(type, b, c, __VA_ARGS__))
#define vec5(type, name, a, b, c, d, ...)                                      \
	vector<vector<vector<vector<vector<type>>>>> name(                         \
		a, VEC4(type, b, c, d, __VA_ARGS__))
#define VEC5(type, a, b, c, d, ...)                                            \
	vector<vector<vector<vector<vector<type>>>>>(                              \
		a, VEC4(type, b, c, d, __VA_ARAGS__))
//stack定義
#define stk(type, name) stack<type> name;
//queue定義
#define que(type, name) queue<type> name;
//点、線分、円
using DD = double;
const DD INF = 1LL<<60;      // 最適化
const DD EPS = 1e-10;        // 最適化
const DD PI = acosl(-1.0);
DD torad(int deg) {return (DD)(deg) * PI / 180;}
DD todeg(DD ang) {return ang * 180 / PI;}
#define INF 2000000000
#define MAX_V 10000
 
//point operator
struct Point {
    DD x, y;
    Point(DD x = 0.0, DD y = 0.0) : x(x), y(y) {}
    friend ostream& operator << (ostream &s, const Point &p) {return s << '(' << p.x << ", " << p.y << ')';}
};
inline Point operator + (const Point &p, const Point &q) {return Point(p.x + q.x, p.y + q.y);}
inline Point operator - (const Point &p, const Point &q) {return Point(p.x - q.x, p.y - q.y);}
inline Point operator * (const Point &p, DD a) {return Point(p.x * a, p.y * a);}
inline Point operator * (DD a, const Point &p) {return Point(a * p.x, a * p.y);}
inline Point operator * (const Point &p, const Point &q) {return Point(p.x * q.x - p.y * q.y, p.x * q.y + p.y * q.x);}
inline Point operator / (const Point &p, DD a) {return Point(p.x / a, p.y / a);}
inline Point conj(const Point &p) {return Point(p.x, -p.y);}
inline Point rot(const Point &p, DD ang) {return Point(cos(ang) * p.x - sin(ang) * p.y, sin(ang) * p.x + cos(ang) * p.y);}
inline Point rot90(const Point &p) {return Point(-p.y, p.x);}
inline DD cross(const Point &p, const Point &q) {return p.x * q.y - p.y * q.x;}
inline DD dot(const Point &p, const Point &q) {return p.x * q.x + p.y * q.y;}
inline DD norm(const Point &p) {return dot(p, p);}
inline DD abs(const Point &p) {return sqrt(dot(p, p));}
inline DD amp(const Point &p) {DD res = atan2(p.y, p.x); if (res < 0) res += PI*2; return res;}
inline bool eq(const Point &p, const Point &q) {return abs(p - q) < EPS;}
inline bool operator < (const Point &p, const Point &q) {return (abs(p.x - q.x) > EPS ? p.x < q.x : p.y < q.y);}
inline bool operator > (const Point &p, const Point &q) {return (abs(p.x - q.x) > EPS ? p.x > q.x : p.y > q.y);}
inline Point operator / (const Point &p, const Point &q) {return p * conj(q) / norm(q);}
//line operator
struct Line : vector<Point> {
    Line(Point a = Point(0.0, 0.0), Point b = Point(0.0, 0.0)) {
        this->push_back(a);
        this->push_back(b);
    }
    friend ostream& operator << (ostream &s, const Line &l) {return s << '{' << l[0] << ", " << l[1] << '}';}
};
//circle operator
struct Circle : Point {
    DD r;
    Circle(Point p = Point(0.0, 0.0), DD r = 0.0) : Point(p), r(r) {}
    friend ostream& operator << (ostream &s, const Circle &c) {return s << '(' << c.x << ", " << c.y << ", " << c.r << ')';}
};
// rep系
#define rep(i, n) for (int i = 0; i < n; i++)
#define rep2(i, x, n) for (int i = x; i <= n; i++)
#define rep3(i, x, n) for (int i = x; i >= n; i--)
#define each(e, v) for (auto &e : v)
//便利系
#define ff first
#define ss second
#define len(x) (x.size())
#define all(x) x.begin(), x.end()
#define rall(x) x.rbegin(), x.rend()
#define pb push_back
#define eb emplace_back
#define sz(x) (int)x.size()
#define acm(x) accumulate(all(x), 0)
#define maximum(...) max({__VA_ARGS__})
#define minimum(...) min({__VA_ARGS__})
//二分探索系
#define LB lower_bound
#define UB upper_bound
// quick_binary_search
#define lb(c, x) distance((c).begin(), lower_bound(all(c), (x)))
#define ub(c, x) distance((c).begin(), upper_bound(all(c), (x)))
// type定義
#define pq(type, name) priority_queue<type> name
#define iq(type, name) priority_queue<type, vector<type>, greater<type>> name
// bit演算系
#define bitsearch(x, process)                                                  \
	for (int tmp = 0; tmp < (1 << x); tmp++) {                                 \
		bitset<x> s(tmp);                                                      \
		process                                                                \
	}
#define get_pos(c, x) (LB(c.begin(), c.end(), x) - c.begin())
//入出力
#define absprint(x) cout << abs(x) << endl;
#define nullE cout << endl;
#define No(x) cout << (x ? "NO" : "No") << endl
#define Yes(x) cout << (x ? "YES" : "Yes") << endl
#define TF(x) cout << (x ? "true" : "false") << endl
// function list
ll Factorial(ll k);
ll modpower(ll a, ll n, ll mod);
ll vecmax(vector<ll> v);
ll vecmin(vector<ll> v);
vector<pair<ll, ll>> primeFunc(ll N);
ll dgitget(ll num);
vector<ll> divisor(ll n);
ll gcd(ll a, ll b);
ll lcm(ll a, ll b);
bool ntimes_of(ll a, ll b);
bool nplus_of(ll a, ll b);
bool boolswitch(bool x); // bool反転
// function list
// function

// to_string func
template <class T> string to_string(T s);
template <class S, class T> string to_string(pair<S, T> p);
string to_string(char c) {
	return string(1, c);
}
string to_string(string s) {
	return s;
}
string to_string(const char s[]) {
	return string(s);
}
template <class T> string to_string(T v) {
	if (v.empty())
		return "{}";
	string ret = "{";
	for (auto x : v)
		ret += to_string(x) + ",";
	ret.back() = '}';
	return ret;
}
template <class S, class T> string to_string(pair<S, T> p) {
	return "{" + to_string(p.first) + ":" + to_string(p.second) + "}";
}

// in out func
void in() {
}
template <typename Head, typename... Tail>
void in(Head &&head, Tail &&... tail) {
	cin >> head;
	in(forward<Tail>(tail)...);
}
void out() {
	cout << '\n';
}
template <typename Head, typename... Tail>
void out(Head &&head, Tail &&... tail) {
	cout << head << ' ';
	out(forward<Tail>(tail)...);
}
void outn() {
}
template <typename Head, typename... Tail>
void outn(Head &&head, Tail &&... tail) {
	cout << head << '\n';
	outn(forward<Tail>(tail)...);
}
template <typename T, typename U> void in(pair<T, U> &p) {
	cin >> p.first >> p.second;
}
template <typename T, typename U> void out(pair<T, U> p) {
	cout << p.first << ' ' << p.second << '\n';
}

// vin vout func
template <typename T> void vin(vector<T> &a) {
	rep(i, sz(a)) cin >> a[i];
}
template <typename T> void vout(const vector<T> &a) {
	for (auto &e : a)
		cout << e << ' ';
	cout << '\n';
}
template <typename T> void voutn(const vector<T> &a) {
	for (auto &e : a)
		cout << e << '\n';
}
template <typename T, typename U> void vin(vector<pair<T, U>> &p) {
	rep(i, sz(p)) cin >> p[i].first >> p[i].second;
}
template <typename T, typename U> void vout(vector<pair<T, U>> &p) {
	rep(i, sz(p)) cout << p[i].first << ' ' << p[i].second;
}
template <typename T, typename U> void voutn(vector<pair<T, U>> &p) {
	rep(i, sz(p)) cout << p[i].first << ' ' << p[i].second << endl;
}

template <typename T> void unique(vector<T> &a) {
	sort(all(a)), a.erase(unique(all(a)), a.end());
}
int vector_finder(vector<int> vec, int number) {
	auto itr = find(vec.begin(), vec.end(), number);
	size_t index = distance(vec.begin(), itr);
	if (index != vec.size()) {
		return itr - vec.begin();
	} else {
		return -1;
	}
}
vector<int> iota(int n) {
	vector<int> ret(n);
	iota(all(ret), 0);
	return ret;
}
template <typename T>
vector<int> iota(const vector<T> &a, bool greater = false) {
	vector<int> ret = iota(sz(a));
	sort(all(ret), [&](int i, int j) { return (a[i] < a[j]) ^ greater; });
	return ret;
}
struct io_setup {
	io_setup() {
		ios_base::sync_with_stdio(false);
		cin.tie(NULL);
		cout << fixed << setprecision(15);
	}
} io_setup;

//演算系
template <typename T> bool chmax(T &x, const T &y) {
	if (x < y) {
		x = y;
		return true;
	}
	return false;
}
template <typename T> bool chmin(T &x, const T &y) {
	if (x > y) {
		x = y;
		return true;
	}
	return false;
}
template <typename T> bool judgePrime(T n) {
	for (T i = 2; i * i <= n; i++) {
		if (n % i == 0)
			return false;
	}
	return n != 1;
}
ll Factorial(ll k) {
	ll sum = 1;
	for (ll i = 1; i <= k; ++i) {
		sum *= i;
	}
	return sum;
}
ll modpower(ll a, ll n, ll mod) {
	ll res = 1;
	while (n > 0) {
		if (n & 1)
			res = res * a % mod;
		a = a * a % mod;
		n >>= 1;
	}
	return res;
}
ll vecmax(vector<ll> v) {
	sort(all(v));
	return v.at(v.size() - 1);
}
ll vecmin(vector<ll> v) {
	sort(all(v));
	return v.at(0);
}
vector<ll> simplePF(ll N){
  vector<ll> res;
 	for (ll z = 2; z * z <= N; ++z) {
		if (N % z != 0)
			continue;
		ll exnum = 0;

		while (N % z == 0) {
			++exnum;
			N /= z;
		}

		res.pb(z);
	}

	if (N != 1)
		res.pb(N);
	return res; 
}
vector<pair<ll, ll>> primeFunc(ll N) {
	vector<pair<ll, ll>> res;
	for (ll z = 2; z * z <= N; ++z) {
		if (N % z != 0)
			continue;
		ll exnum = 0;

		while (N % z == 0) {
			++exnum;
			N /= z;
		}

		res.pb({z, exnum});
	}

	if (N != 1)
		res.pb({N, 1});
	return res;
}
ll dgitget(ll num) {
	ll ans = 0;
	while (num != 0) {
		num /= 10;
		++ans;
	}

	return ans;
}
vector<ll> divisor(ll n) {
	vector<ll> ret;
	for (ll i = 1; i * i <= n; i++) {
		if (n % i == 0) {
			ret.pb(i);
			if (i * i != n)
				ret.pb(n / i);
		}
	}
	sort(ret.begin(), ret.end());
	return ret;
}

ll gcd(ll a, ll b) {
	if (a % b == 0) {
		return b;
	} else {
		return gcd(b, a % b);
	}
}

ll lcm(ll a, ll b) {
	return (a / gcd(a, b)) * b;
}

bool ntimes_of(ll a, ll b) {
	ll prod = a * b;

	return (prod / b == a);
}

bool nplus_of(ll a, ll b) {
	ll prod = a + b;
	return (prod > 0);
}

bool boolswitch(bool x) {
	if (x == true) {
		return false;
	} else {
		return true;
	}
}

//文字列処理
//split(文字列,カットに使う文字)
vector<string> split(string s, char delim) {
	int startIdx = 0;
	int len = 0;
	vector<string> ret;
	for (int i = 0; i < (int)s.size(); i++) {
		if (s[i] == delim) {
			if (len != 0)
				ret.push_back(s.substr(startIdx, len));
			startIdx = i + 1;
			len = 0;
		} else {
			len++;
		}
	}
	if (len != 0)
		ret.push_back(s.substr(startIdx, len));
	return ret;
}
string toLower(string s) {
	transform(s.cbegin(), s.cend(), s.begin(), ::tolower);
	return s;
}
string toUpper(string s) {
	transform(s.cbegin(), s.cend(), s.begin(), ::toupper);
	return s;
}
//最長回文
//vector< int > manacher(const string &s) {//iを中心とした回文半径を出力
//  vector< int > radius(s.size());
//  int i = 0, j = 0;
//  while(i < s.size()) {
//    while(i - j >= 0 && i + j < s.size() && s[i - j] == s[i + j]) {
//      ++j;
//    }
//   radius[i] = j;
//    int k = 1;
//    while(i - k >= 0 && i + k < s.size() && k + radius[i - k] < j) {
//      radius[i + k] = radius[i - k];
//    ++k;
//    }
//    i += k;
//    j -= k;
//  }
//  return radius;
//}
//一次不定方程式 /ax + by = gcd(a,b) /返り値 => (x,y)
ll diophantine_Equation(ll a, ll b, ll &x, ll &y) {
    if (b == 0) { x = 1; y = 0; return a; }
    long long d = diophantine_Equation(b, a%b, y, x);
    y -= a/b * x;
    return d;
}
//二数値間の数列

//点と線分の位置
// 1：a-bから見てcは左側(反時計回り)、-1：a-bから見てcは右側(時計回り)、0：一直線上
int simple_ccw(const Point &a, const Point &b, const Point &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    return 0;
}
// 1：a-bから見てcは左側(反時計回り)、-1：a-bから見てcは右側(時計回り)
// 2：c-a-bの順に一直線上、-2：a-b-cの順に一直線上、0：a-c-bの順に一直線上
int ccw(const Point &a, const Point &b, const Point &c) {
    if (cross(b-a, c-a) > EPS) return 1;
    if (cross(b-a, c-a) < -EPS) return -1;
    if (dot(b-a, c-a) < -EPS) return 2;
    if (norm(b-a) < norm(c-a) - EPS) return -2;
    return 0;
}
// 点と三角形の包含関係(辺上については判定していない)
bool is_contain(const Point &p, const Point &a, const Point &b, const Point &c) {
    int r1 = simple_ccw(p, b, c), r2 = simple_ccw(p, c, a), r3 = simple_ccw(p, a, b);
    if (r1 == 1 && r2 == 1 && r3 == 1) return true;
    if (r1 == -1 && r2 == -1 && r3 == -1) return true;
    return false;
}
// 2点の比率a:bのアポロニウスの円
Circle Apporonius(const Point &p, const Point &q, DD a, DD b) {
    if ( abs(a-b) < EPS ) return Circle(Point(0,0),0);
    Point c1 = (p * b + q * a) / (a + b);
    Point c2 = (p * b - q * a) / (b - a);
    Point c = (c1 + c2) / 2;
    DD r = abs(c - c1);
    return Circle(c, r);
}
//即席function
// main
int32_t main() {
int n;
  in(n);
vec1(pii,vec,n);
vin(vec);
vec1(int,vecc,n);
rep(i,len(vec)){
  vecc[i] = (vec[i].ff + vec[i].ss) * (vec[i].ss - vec[i].ff + 1) / 2;
}
int ans=0;
rep(i,len(vecc)){
  ans = ans + vecc[i];
}
out(ans);
}
//https://atcoder.jp/contests/abc181/tasks/abc181_b
