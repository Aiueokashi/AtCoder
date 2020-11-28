#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable: 4244) // 最悪をします
#include <atcoder/all>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric>
#include <queue>
#include <stack>
#include <map> 
#include <set>
#include <string>
#include <functional>
#include <list>
#include <random>
#include <time.h>
#include <iomanip>
#include <assert.h>
#include <numeric>
#include <sstream>
#define BIT(nr) (1UL << (nr))
#define int long long
//#define ll long long
#define double long double
#define mod 1000000007
#define MAXN (int)1e+5 * 2+1
#define LL_MAX 9223372036854775807	//ない環境用
#define LL_HALFMAX 9223372036854775807 / 2	//ない環境用
#define MIN -(9223372036854775807 / 2)
#define REP(i,a,n) for(int i=(a); i<(int)(n); i++)
#define rep(i,n) REP(i,0,n)
#define FOR(it,c) for(__typeof((c).begin()) it=(c).begin(); it!=(c).end(); ++it)
#define ALLOF(c) (c).begin(), (c).end()
#define REPS(i,x) for(int i=1;i<=(int)(x);i++)
#define RREP(i,x) for(int i=((int)(x)-1);i>=0;i--)
#define RREPS(i,x) for(int i=((int)(x));i>0;i--)
#define repl(i,a,b) for(int i=(int)(a);i<(int)(b);i++)
#define mp make_pair
template<typename T1, typename T2> inline void chmin(T1 & a, T2 b) { if (a > b) a = b; }
template<typename T1, typename T2> inline void chmax(T1& a, T2 b) { if (a < b) a = b; }
 
 
using namespace std;
 
//デバッグ用カッコの有無
#ifdef DEBUG
template <class T>ostream &operator<<(ostream &o, const vector<T>&v)
{
	o << "{"; for (int i = 0; i<(int)v.size(); i++)o << (i>0 ? ", " : "") << v[i]; o << "}"; return o;
}
#endif // DEBUG
 
template <class T>ostream &operator<<(ostream &o, const vector<T>&v)
{
	for (int i = 0; i<(int)v.size(); i++)o << (i>0 ? " " : "") << v[i]; return o;
}
 
int dx[4] = { 0, 1, 0, -1 }; // x軸方向への変位
int dy[4] = { 1, 0, -1, 0 }; // y軸方向への変位
 
int dxp[4] = { 0, 1 }; // x軸方向への変位(正のみ)
int dyp[4] = { 1, 0 }; // y軸方向への変位(負のみ)
 
using Weight = int;
using Flow = int;
struct Edge {
	int src, dst;
 
	// libalgo のものに追加、メンバを追加するだけなので互換性は崩さないはず、逆辺のG[e.dstの]インデックスを保持
	int rev;
	Weight weight;
	Flow cap;
	Edge() : src(0), dst(0), weight(0) {}
	Edge(int s, int d, Weight w) : src(s), dst(d), weight(w) {}
};
using Edges = std::vector<Edge>;
using Graph = std::vector<Edges>;
using Array = std::vector<Weight>;
using Matrix = std::vector<Array>;
 
void add_edge(Graph& g, int a, int b, Weight w = 1) {
	g[a].emplace_back(a, b, w);
	g[b].emplace_back(b, a, w);
}
void add_arc(Graph& g, int a, int b, Weight w = 1) { g[a].emplace_back(a, b, w); }
 
// 辺がメンバに src と dst を持つ隣接リスト表記のグラフをダンプ(https://hello-world-494ec.firebaseapp.com/) に投げることを想定
template <typename T>
void dump_graph(T G) {
	int V = G.size();
	int E = 0;
	ostringstream os;
 
	for (auto es : G) {
		for (auto e : es) {
			E++;
			os << e.src << " " << e.dst << "\n";
		}
	}
	cout << V << " " << E << "\n";
	cout << os.str() << "\n";
}
 
// グリッドからグラフを構築
// @pre: gはノード数H*Wのグラフ
void create_from_grid(Graph& g, int h, int w, vector<string>& mapData, char wall) {
	//グラフ構築 O(HW)
	rep(y, h) {
		rep(x, w) {
			if (mapData[y][x] == wall) {
				continue;
			}
 
			int id = y * w + x;
			//右と下(変位が正)のみ見る(辺の重複を回避するため)
			rep(i, 2) {
				int nx = x + dxp[i];
				int ny = y + dyp[i];
				int nid = ny * w + nx;
				if (nx < 0 || nx >= w) {
					continue;
				}
				if (ny < 0 || ny >= h) {
					continue;
				}
				if (mapData[ny][nx] != wall) {
					add_edge(g, id, nid);
				}
			}
		}
	}
}
 
// マスに重みが定義されるグリッドから重み付きグラフを構築、ダイクストラなどをするとき始点のこすとは入らないことに注意
// @pre: gはノード数H*Wのグラフ
void create_weighted_from_grid(Graph& g, int h, int w, vector<vector<int>>& mapData) {
	//グラフ構築 O(HW)
	rep(y, h) {
		rep(x, w) {
			int id = y * w + x;
			// こんどは全方向見る(行きと帰りで重みが違うはず)
			rep(i, 4) {
				int nx = x + dx[i];
				int ny = y + dy[i];
				int nid = ny * w + nx;
				if (nx < 0 || nx >= w) {
					continue;
				}
				if (ny < 0 || ny >= h) {
					continue;
				}
 
				// 移動先のコストを足す
				add_arc(g, id, nid, mapData[ny][nx]);
			}
		}
	}
}
 
// グリッドにおいて座標をグラフのノード番号に変換する
int point_to_node_num(int x, int y, int W) {
	return y * W + x;
}
 
struct uf_tree {
	std::vector<int> parent;
	int __size;
	uf_tree(int size_) : parent(size_, -1), __size(size_) {}
	void unite(int x, int y) {
		if ((x = find(x)) != (y = find(y))) {
			if (parent[y] < parent[x]) std::swap(x, y);
			parent[x] += parent[y];
			parent[y] = x;
			__size--;
		}
	}
	bool is_same(int x, int y) { return find(x) == find(y); }
	int find(int x) { return parent[x] < 0 ? x : parent[x] = find(parent[x]); }
	int size(int x) { return -parent[find(x)]; }
	int size() { return __size; }
};
 
 
 
//!!!問題をちゃんと読む!!!
//!!!問題をちゃんと読め!!!
//!!!問題は読みましたか？!!!
 
template <signed M, unsigned T>
struct mod_int {
	constexpr static signed MODULO = M;
	constexpr static unsigned TABLE_SIZE = T;
 
	signed x;
 
	mod_int() : x(0) {}
 
	mod_int(long long y) : x(static_cast<signed>(y >= 0 ? y % MODULO : MODULO - (-y) % MODULO)) {}
 
	mod_int(signed y) : x(y >= 0 ? y % MODULO : MODULO - (-y) % MODULO) {}
 
	mod_int& operator+=(const mod_int& rhs) {
		if ((x += rhs.x) >= MODULO) x -= MODULO;
		return *this;
	}
 
	mod_int& operator-=(const mod_int& rhs) {
		if ((x += MODULO - rhs.x) >= MODULO) x -= MODULO;
		return *this;
	}
 
	mod_int& operator*=(const mod_int& rhs) {
		x = static_cast<signed>(1LL * x * rhs.x % MODULO);
		return *this;
	}
 
	mod_int& operator/=(const mod_int& rhs) {
		x = static_cast<signed>((1LL * x * rhs.inv().x) % MODULO);
		return *this;
	}
 
	mod_int operator-() const { return mod_int(-x); }
 
	mod_int operator+(const mod_int& rhs) const { return mod_int(*this) += rhs; }
 
	mod_int operator-(const mod_int& rhs) const { return mod_int(*this) -= rhs; }
 
	mod_int operator*(const mod_int& rhs) const { return mod_int(*this) *= rhs; }
 
	mod_int operator/(const mod_int& rhs) const { return mod_int(*this) /= rhs; }
 
	bool operator<(const mod_int& rhs) const { return x < rhs.x; }
 
	mod_int inv() const {
		assert(x != 0);
		if (x <= static_cast<signed>(TABLE_SIZE)) {
			if (_inv[1].x == 0) prepare();
			return _inv[x];
		}
		else {
			signed a = x, b = MODULO, u = 1, v = 0, t;
			while (b) {
				t = a / b;
				a -= t * b;
				std::swap(a, b);
				u -= t * v;
				std::swap(u, v);
			}
			return mod_int(u);
		}
	}
 
	mod_int pow(long long t) const {
		assert(!(x == 0 && t == 0));
		mod_int e = *this, res = mod_int(1);
		for (; t; e *= e, t >>= 1)
			if (t & 1) res *= e;
		return res;
	}
 
	mod_int fact() {
		if (_fact[0].x == 0) prepare();
		return _fact[x];
	}
 
	mod_int inv_fact() {
		if (_fact[0].x == 0) prepare();
		return _inv_fact[x];
	}
 
	mod_int choose(mod_int y) {
		assert(y.x <= x);
		return this->fact() * y.inv_fact() * mod_int(x - y.x).inv_fact();
	}
 
	static mod_int _inv[TABLE_SIZE + 1];
 
	static mod_int _fact[TABLE_SIZE + 1];
 
	static mod_int _inv_fact[TABLE_SIZE + 1];
 
	static void prepare() {
		_inv[1] = 1;
		for (int i = 2; i <= (int)TABLE_SIZE; ++i) {
			_inv[i] = 1LL * _inv[MODULO % i].x * (MODULO - MODULO / i) % MODULO;
		}
		_fact[0] = 1;
		for (unsigned i = 1; i <= TABLE_SIZE; ++i) {
			_fact[i] = _fact[i - 1] * signed(i);
		}
		_inv_fact[TABLE_SIZE] = _fact[TABLE_SIZE].inv();
		for (int i = (int)TABLE_SIZE - 1; i >= 0; --i) {
			_inv_fact[i] = _inv_fact[i + 1] * (i + 1);
		}
	}
};
 
template <signed M, unsigned F>
std::ostream& operator<<(std::ostream& os, const mod_int<M, F>& rhs) {
	return os << rhs.x;
}
 
template <signed M, unsigned F>
std::istream& operator >> (std::istream& is, mod_int<M, F>& rhs) {
	long long s;
	is >> s;
	rhs = mod_int<M, F>(s);
	return is;
}
 
template <signed M, unsigned F>
mod_int<M, F> mod_int<M, F>::_inv[TABLE_SIZE + 1];
 
template <signed M, unsigned F>
mod_int<M, F> mod_int<M, F>::_fact[TABLE_SIZE + 1];
 
template <signed M, unsigned F>
mod_int<M, F> mod_int<M, F>::_inv_fact[TABLE_SIZE + 1];
 
template <signed M, unsigned F>
bool operator==(const mod_int<M, F>& lhs, const mod_int<M, F>& rhs) {
	return lhs.x == rhs.x;
}
 
template <int M, unsigned F>
bool operator!=(const mod_int<M, F>& lhs, const mod_int<M, F>& rhs) {
	return !(lhs == rhs);
}
 
const signed MF = 1000010;
const signed MOD = 1000000007;
 
using mint = mod_int<MOD, MF>;
 
mint binom(int n, int r) { return (r < 0 || r > n || n < 0) ? 0 : mint(n).choose(r); }
 
mint fact(int n) { return mint(n).fact(); }
 
mint inv_fact(int n) { return mint(n).inv_fact(); }
 
//出典 http://beet-aizu.hatenablog.com/entry/2017/12/01/225955
/*
コンストラクタ引数説明
int n_
要素数。
f
2つの要素Tをマージするための関数。
区間MAX区間更新の時: max
区間Sum区間Addの時: +
g
1つの要素Tに作用素Eを適用するための関数。
区間MAX区間更新の時: =
区間Sum区間Addの時: +
h
2つの作用素Eをマージするための関数。
区間MAX区間更新の時: =
区間Sum区間Addの時: +
T d1
演算fの単位元。
区間MAX区間更新の時: -INF　
区間Sum区間Addの時: 0
E d0,
g, hの単位元。
区間MAX区間更新の時: 定義域外のどこか
区間Sum区間Addの時: 0
vector<T> v = vector<T>()
セグ木を構成するときのvector
P p = [](E a, int b) {return a; }
区間の長さbを引数に取り、区間の長さによって変化する作用素E'を返す関数。
例えば、区間MAX区間Addの時なんかは区間長によって足すべき数が変化するので必要
区間Sum区間Addの時: *
 
//具体例
//区間chmin, 区間min
auto myMin = [](int a, int b) {return min(a, b); };
SegmentTree<int, int> seg(n, myMin, myMin, myMin, LL_HALFMAX, LL_HALFMAX);
//区間update、区間min
SegmentTree<int, int> seg(n, myMin, myMin, myMin, LL_HALFMAX, LL_HALFMAX);
//区間Add、区間Sum
vector<int> v(0, N + 1);
SegmentTree<int, int> segtree(N + 1, plus<int>(), plus<int>(), plus<int>(), 0, 0, v, [](int a, int b) {return a * b; });
//区間Add、区間Min
vector<int> v(0, N + 1);
SegmentTree<int, int> segtree(N + 1, myMin, plus<int>(), plus<int>(), LL_HALFMAX, 0, v, [](int a, int b) {return a; });
*/
 
template <typename T, typename E>
struct SegmentTree {
	typedef function<T(T, T)> F;
	typedef function<T(T, E)> G;
	typedef function<E(E, E)> H;
	typedef function<E(E, int)> P;
	int n;
	F f;
	G g;
	H h;
	P p;
	T d1;
	E d0;
	vector<T> dat;
	vector<E> laz;
	SegmentTree(int n_, F f, G g, H h, T d1, E d0,
		vector<T> v = vector<T>(), P p = [](E a, int b) {return a; }) :
		f(f), g(g), h(h), d1(d1), d0(d0), p(p) {
		init(n_);
		if (n_ == (int)v.size()) build(n_, v);
	}
	//初期化。要素配列と遅延配列を2*n-1個にする
	void init(int n_) {
		n = 1;
		while (n < n_) n *= 2;
		dat.clear();
		dat.resize(2 * n - 1, d1);
		laz.clear();
		laz.resize(2 * n - 1, d0);
	}
	//既存のvectorからセグ木を構築
	void build(int n_, vector<T> v) {
		for (int i = 0; i < n_; i++) dat[i + n - 1] = v[i];
		for (int i = n - 2; i >= 0; i--)
			dat[i] = f(dat[i * 2 + 1], dat[i * 2 + 2]);
	}
	//ノードを評価する。
	inline void eval(int len, int k) {
		//遅延配列に単位元が入ってたら評価済みなのでおしまい
		if (laz[k] == d0) return;
		//葉ノードでないなら遅延伝播する
		if (k * 2 + 1 < n * 2 - 1) {
			//h: 2つの作用素を引数に取り合成した作用素を返す関数
			laz[k * 2 + 1] = h(laz[k * 2 + 1], laz[k]);
			laz[k * 2 + 2] = h(laz[k * 2 + 2], laz[k]);
		}
		//p: このノードに対応する区間長と作用素を引数に取り、区間長に対応する作用素を返す関数
		//dat[k] にlaz に溜めていた作用素を適用(g: 要素型と作用素型を引数に取り、要素に作用素を作用させた結果を返す関数、ここでの作用素とは区間Sum区間Addなら (+ 3) とか)
		dat[k] = g(dat[k], p(laz[k], len));
		//適用し終わったので遅延配列をクリア
		laz[k] = d0;
	}
	//[l,r)の区間を再帰的に見ながら0-indexedの[a, b)を更新する
	T update(int a, int b, E x, int k, int l, int r) {
		//先に評価
		eval(r - l, k);
		//範囲外ならなにもしないでそのノードが持つ値を返す
		if (r <= a || b <= l) return dat[k];
		//完全被覆なら既に遅延配列に入っている作用素と追加したい作用素をマージした後にそれを要素に作用させた結果を返す、pは区間長に対応する作用素を得るための（ｒｙ
		if (a <= l && r <= b) {
			laz[k] = h(laz[k], x);
			return g(dat[k], p(laz[k], r - l));
		}
		//完全被覆でも範囲外でもないなら(中途半端にかぶっているなら)完全被覆と範囲外の境界が見えるまで木を潜って変化後の値を得る
		return dat[k] = f(update(a, b, x, k * 2 + 1, l, (l + r) / 2),
			update(a, b, x, k * 2 + 2, (l + r) / 2, r));
	}
	T update(int a, int b, E x) {
		return update(a, b, x, 0, 0, n);
	}
 
	T update(int a, E x) {
		return update(a, a + 1, x);
	}
 
	T query(int a, int b, int k, int l, int r) {
		eval(r - l, k);
		//範囲外なら単位元を返す
		if (r <= a || b <= l) return d1;
		//完全被覆ならそのまま返す
		if (a <= l && r <= b) return dat[k];
		//一部被覆なら完全被覆と範囲外に分かれるまで木を潜る
		T vl = query(a, b, k * 2 + 1, l, (l + r) / 2);
		T vr = query(a, b, k * 2 + 2, (l + r) / 2, r);
		return f(vl, vr);
	}
	//0-indexedで[a, b)の区間*を求める
	T query(int a, int b) {
		return query(a, b, 0, 0, n);
	}
 
	T query(int a) {
		return query(a, a + 1, 0, 0, n);
	}
 
	void debug_print(int num) {
		vector<T> v;
		rep(i, num) {
			v.push_back(query(i));
		}
		cout << "{" << v << "}\n";
	}
};
 
//座標圧縮
 
class compress {
public:
	map<int, int> zip;
	vector<int> unzip;
 
	compress(vector<int> x)
	{
		sort(x.begin(), x.end());
		x.erase(unique(x.begin(), x.end()), x.end());
		for (int i = 0; i < x.size(); i++) {
			zip[x[i]] = i;
			unzip.push_back(i);
		}
	}
};
 
 
int euclidean_gcd(int a, int b) {
	while (1) {
		if (a < b) swap(a, b);
		if (!b) break;
		a %= b;
	}
	return a;
}
 
//https://ei1333.github.io/luzhiled/snippets/dp/cumulative-sum-2d.html
template< class T >
struct CumulativeSum2D {
	vector< vector< T > > data;
 
	CumulativeSum2D(int W, int H) : data(W + 1, vector< int >(H + 1, 0)) {}
 
	void add(int x, int y, T z) {
		++x, ++y;
		if (x >= data.size() || y >= data[0].size()) return;
		data[x][y] += z;
	}
 
	void build() {
		for (int i = 1; i < data.size(); i++) {
			for (int j = 1; j < data[i].size(); j++) {
				data[i][j] += data[i][j - 1] + data[i - 1][j] - data[i - 1][j - 1];
			}
		}
	}
 
	T query(int sx, int sy, int gx, int gy) {
		return (data[gx][gy] - data[sx][gy] - data[gx][sy] + data[sx][sy]);
	}
};
 
//lib
int nC2(int n) {
	return n * (n - 1) / 2;
}
 
class node {
public:
	int depth;
	int num;
 
	node(int d, int n) {
		depth = d;
		num = n;
	}
};
 
template< class T >
struct CumulativeSum {
	vector< T > data;
 
	CumulativeSum(int sz) : data(sz, 0) {};
 
	void add(int k, T x) {
		data[k] += x;
	}
 
	void build() {
		for (int i = 1; i < data.size(); i++) {
			data[i] += data[i - 1];
		}
	}
 
	T query(int k) {
		if (k < 0) return (0);
		return (data[min(k, (int)data.size() - 1)]);
	}
	//[left, right]の和
	T query(int left, int right) {
		return query(right) - query(left - 1);
	}
};
 
std::vector<int> eratosthenes_sieve(int n) {
	std::vector<int> ps(n + 1);
	std::iota(ps.begin() + 2, ps.end(), 2);
	for (int i = 2; i * i <= n; ++i)
		if (ps[i])
			for (int j = i * i; j <= n; j += i) ps[j] = 0;
	return ps;
}
 
std::vector<int> make_primes(int n) {
	std::vector<int> ps = eratosthenes_sieve(n);
	ps.erase(std::remove(ps.begin(), ps.end(), 0), ps.end());
	return ps;
}
 
// 区間[a, b)の素数判定をする、is_prime[i]: a + i が素数 or not つまり is_prime[i-a] が true: iが素数
std::vector<bool> segment_eratosthenes_sieve(int a, int b) {
	vector<bool> is_prime(b - a, true);
	vector<bool> is_prime_small;
	for (int i = 0; i*i < b; i++)is_prime_small.push_back(true);
 
	for (int i = 2; i*i < b; i++) {
		if (is_prime_small[i]) {
			for (int j = 2 * i; j*j < b; j += i) {
				is_prime_small[j] = false;	// [2, sqrt(b))のふるい
			}
			// (a + i - 1LL) / i * i a以上の最小のiの倍数
			for (int j = max(2LL, (a + i - 1LL) / i) * i; j < b; j += i) {
				is_prime[j - a] = false;	// [a, b)のふるい
			}
		}
	}
	return is_prime;
}
 
vector< int64_t > divisor(int64_t n) {
	vector< int64_t > ret;
	for (int64_t i = 1; i * i <= n; i++) {
		if (n % i == 0) {
			ret.push_back(i);
			if (i * i != n) ret.push_back(n / i);
		}
	}
	sort(begin(ret), end(ret));
	return (ret);
}
 
 
 
 
// 汎用的な二分探索のテンプレ(めぐる式)
int binary_search(function<bool(int)> isOk, int ng, int ok) {
 
	/* ok と ng のどちらが大きいかわからないことを考慮 */
	while (abs(ok - ng) > 1) {
		int mid = (ok + ng) / 2;
 
		if (isOk(mid)) ok = mid;
		else ng = mid;
	}
	return ok;
}
 
std::pair<std::vector<Weight>, bool> bellmanFord(const Graph& g, int s) {
	int n = g.size();
	const Weight inf = std::numeric_limits<Weight>::max() / 8;
	Edges es;
	for (int i = 0; i < n; i++)
		for (auto& e : g[i]) es.emplace_back(e);
 
	//初期化、スタート地点以外の距離は無限大
	std::vector<Weight> dist(n, inf);
	dist[s] = 0;
	bool negCycle = false;
	for (int i = 0;; i++) {
		bool update = false;
		//すべての辺について、その辺をとおった場合に最短経路が更新できる場合は更新する
		for (auto& e : es) {
			if (dist[e.src] != inf && dist[e.dst] > dist[e.src] + e.weight) {
				dist[e.dst] = dist[e.src] + e.weight;
				update = true;
			}
		}
		//更新がなくなったらおはり
		if (!update) break;
		//n回以上更新されてたら負閉路がある
		if (i > n) {
			negCycle = true;
			break;
		}
	}
	return std::make_pair(dist, !negCycle);
}
 
//ゴールを指定して、それまでのパスに負閉路がなかったらOK(嘘修正済)
std::pair<std::vector<Weight>, bool> bellmanFord(const Graph& g, int s, int d) {
	int n = g.size();
	const Weight inf = std::numeric_limits<Weight>::max() / 8;
	Edges es;
	for (int i = 0; i < n; i++)
		for (auto& e : g[i]) es.emplace_back(e);
 
	//初期化、スタート地点以外の距離は無限大
	std::vector<Weight> dist(n, inf);
	dist[s] = 0;
	bool negCycle = false;
	for (int i = 0; i < n * 2; i++) {
		bool update = false;
		//すべての辺について、その辺をとおった場合に最短経路が更新できる場合は更新する
		for (auto& e : es) {
			if (dist[e.src] != inf && dist[e.dst] > dist[e.src] + e.weight) {
				// n回目の更新で d が更新されてたら問答無用で負閉路ありとしてNG
				if (i >= n - 1 && e.dst == d) {
					negCycle = true;
				}
				// 終点以外に負閉路がある場合はそこの距離を十分小さい値に置き換える
				else if (i >= n - 1) {
					dist[e.dst] = -inf;
					update = true;
				}
				else {
					dist[e.dst] = dist[e.src] + e.weight;
					update = true;
				}
			}
		}
		//更新がなくなったらおはり
		if (!update) break;
	}
	return std::make_pair(dist, !negCycle);
}
 
//R[i] == S[i] を中心とした極大回文長 なるvector Rを返す
vector<int> Manachar(string S) {
	int len = S.length();
	vector<int> R(len);
 
	int i = 0, j = 0;
	while (i < S.size()) {
		while (i - j >= 0 && i + j < S.size() && S[i - j] == S[i + j]) ++j;
		R[i] = j;
		int k = 1;
		while (i - k >= 0 && i + k < S.size() && k + R[i - k] < j) R[i + k] = R[i - k], ++k;
		i += k; j -= k;
	}
	return R;
}
 
std::vector<int> tsort(const Graph &g) {
	int n = g.size(), k = 0;
	std::vector<int> ord(n), in(n);
	for (auto &es : g)
		for (auto &e : es) in[e.dst]++;
	std::queue<int> q;
	//入次数0の点をキューに追加
	for (int i = 0; i < n; ++i)
		if (in[i] == 0) q.push(i);
	while (q.size()) {
		int v = q.front();
		//Sから node n を削除する
		q.pop();
		//L に n を追加する
		ord[k++] = v;
		for (auto &e : g[v]) {
			//選択した点から出てる辺を削除、0になったらキューに追加
			if (--in[e.dst] == 0) {
				q.push(e.dst);
			}
		}
 
	}
	return *std::max_element(in.begin(), in.end()) == 0 ? ord : std::vector<int>();
}
 
std::vector<Weight> dijkstra(const Graph &g, int s) {
	const Weight INF = std::numeric_limits<Weight>::max() / 8;
	using state = std::tuple<Weight, int>;
	std::priority_queue<state> q;
	std::vector<Weight> dist(g.size(), INF);
	dist[s] = 0;
	q.emplace(0, s);
	while (q.size()) {
		Weight d;
		int v;
		std::tie(d, v) = q.top();
		q.pop();
		d *= -1;
		/* if(v == t) return d; */
		if (dist[v] < d) continue;
		for (auto &e : g[v]) {
			if (dist[e.dst] > dist[v] + e.weight) {
				dist[e.dst] = dist[v] + e.weight;
				q.emplace(-dist[e.dst], e.dst);
			}
		}
	}
	return dist;
}
 
Matrix WarshallFloyd(const Graph &g) {
	auto const INF = std::numeric_limits<Weight>::max() / 8;
	int n = g.size();
	Matrix d(n, Array(n, INF));
	rep(i, n) d[i][i] = 0;
	rep(i, n) for (auto &e : g[i]) d[e.src][e.dst] = std::min(d[e.src][e.dst], e.weight);
	rep(k, n) rep(i, n) rep(j, n) {
		if (d[i][k] != INF && d[k][j] != INF) d[i][j] = std::min(d[i][j], d[i][k] + d[k][j]);
	}
	return d;
}
 
std::pair<std::vector<int>, std::vector<int>> prime_factor_decomp(int n) {
	std::vector<int> p, e;
	int m = n;
	for (int i = 2; i * i <= n; i++) {
		if (m % i != 0) continue;
		int c = 0;
		while (m % i == 0) c++, m /= i;
		p.push_back(i);
		e.push_back(c);
	}
	if (m > 1) {
		p.push_back(m);
		e.push_back(1);
	}
	return std::make_pair(p, e);
}
 
int extgcd(int a, int b, int &x, int &y) {
	int g = a;
	x = 1;
	y = 0;
	if (b != 0) g = extgcd(b, a % b, y, x), y -= (a / b) * x;
	return g;
}
 
// 不定方程式 ax + by = c の一般整数解(pt + q, rt + s)を求める
/*
* exist: 解が存在するか否か
* p, q, r, s: 存在するならば不定方程式の一般解(pt + q, rt + s)
* ここで、式変形から、p > 0、 q < 0 となることに注意する。(解の条件を絞るときなどに必要になる)
*/
void IndeterminateEq(int a, int b, int c, bool& exist, int& p, int& q, int& r, int& s) {
	int X, Y;
 
	int g = euclidean_gcd(a, b);
 
	// c が最大公約数の整数倍でないならNG
	if (c % g != 0) {
		exist = false;
		return;
	}
	exist = true;
 
	// 拡張ユークリッドの互除法で ax + by = gcd(a, b) なる (X, Y) を求める
	extgcd(a, b, X, Y);
	int m = c / g;
 
	// ax + by = c の解にする
	X *= m;
	Y *= m;
 
	int a2 = a / g;
	int b2 = b / g;
 
	p = b2;
	q = X;
	r = -a2;
	s = Y;
}
 
// x^n mod modulo を繰り返し二乗法で計算する 
// n を 2^k の和で表す -> n を二進表記したとき、kbit目(0-indexed)が立っているときだけx^kをかける
int mod_pow(int x, int n, int modulo) {
	int res = 1;
	while (n > 0) {
		if (n & 1) {
			res = res * x % modulo;
		}
		x = x * x % modulo;
		n >>= 1;
	}
	return res;
}
 
int64_t popcnt(int64_t n)
{
	int64_t c = 0;
	c = (n & 0x5555555555555555) + ((n >> 1) & 0x5555555555555555);
	c = (c & 0x3333333333333333) + ((c >> 2) & 0x3333333333333333);
	c = (c & 0x0f0f0f0f0f0f0f0f) + ((c >> 4) & 0x0f0f0f0f0f0f0f0f);
	c = (c & 0x00ff00ff00ff00ff) + ((c >> 8) & 0x00ff00ff00ff00ff);
	c = (c & 0x0000ffff0000ffff) + ((c >> 16) & 0x0000ffff0000ffff);
	c = (c & 0x00000000ffffffff) + ((c >> 32) & 0x00000000ffffffff);
	return(c);
}
 
/*
行列積と行列累乗
行列積
vector<vector<T>> matrixMultiplies(vector<vector<T>> l, vector<vector<T>> r, F plus = plus<T>(), G multiple = multiplies<T>(), T eplus = 0LL)
行列累乗
vector<vector<T>> matrixPower(vector<vector<T>> m, int n, F plus = std::plus<T>(), G multiple = multiplies<T>(), T eplus = 0LL, T emultiple = 1LL)
T:			考える集合(競プロにおいてはたぶんほぼ整数)
l:			左からかける行列
r:			右からかける行列
plus:		加法演算
multiple:	乗法演算
eplus:		加法の単位元
emultiple:	乗法の単位元
*/
template<typename T = long long, typename F = decltype(std::plus<T>()), typename G = decltype(multiplies<T>())>
vector<vector<T>> matrixMultiplies(vector<vector<T>> l, vector<vector<T>> r, F plus = plus<T>(), G multiple = multiplies<T>(), T eplus = 0LL) {
	int rx = r[0].size();
	int ry = r.size();
	vector<vector<T> > ret;
 
	for (int y = 0; y < ry; y++) {
		vector<T> add;
		for (int x = 0; x < rx; x++) {
			T cell = eplus;
			for (int i = 0; i < ry; i++) {
				T mul = multiple(l[y][i], r[i][x]);
				cell = plus(cell, mul);
			}
			add.push_back(cell);
		}
		ret.push_back(add);
	}
	return ret;
}
 
template<typename T = long long, typename F = decltype(std::plus<T>()), typename G = decltype(multiplies<T>())>
vector<vector<T>> matrixPower(vector<vector<T>> m, int n, F plus = std::plus<T>(), G multiple = multiplies<T>(), T eplus = 0LL, T emultiple = 1LL) {
	int k = m.size();
	if (n == 0) {
		vector<vector<T> > E;
		for (int i = 0; i < k; i++) {
			// 単位行列は対角成分を乗法単位元、非対角成分をゼロ元で埋める
			vector<T> v(k, eplus);
			v[i] = emultiple;
			E.push_back(v);
		}
		return E;
	}
	vector<vector<T>> ret = matrixPower(matrixMultiplies(m, m, plus, multiple, eplus), n / 2, plus, multiple, eplus, emultiple);
	if (n % 2 == 1) {
		ret = matrixMultiplies(m, ret, plus, multiple);
	}
	return ret;
}
 
// フロー系のアルゴリズム
// 最大流
/*
Ford-Fulkerson法(蟻本) O(F|E|)
F: 最大流量
E: 辺数
コンストラクタ引数でノード数nを受け取り初期化し、add_edge で辺と逆辺を追加していく
*/
class Ford_Fulkerson {
private:
	struct Edge {
		int src, dst;
 
		// libalgo のものに追加、メンバを追加するだけなので互換性は崩さないはず、逆辺のG[e.dstの]インデックスを保持
		int rev;
		int cap;
		Edge(int s, int d, int c, int r) : src(s), dst(d), cap(c), rev(r) {}
	};
	vector<vector<Edge> > G;
	vector<bool> used;
public:
	Ford_Fulkerson(int n) :
		G(n),
		used(n, false)
	{}
 
	void add_edge(int s, int d, int cap) {
		G[s].emplace_back(s, d, cap, G[d].size());
		G[d].emplace_back(d, s, 0, G[s].size() - 1);
	}
 
	int dfs(int v, int t, int f) {
		if (v == t) {
			return f;
		}
		used[v] = true;
		for (Edge& e : G[v]) {
			if (!used[e.dst] && e.cap > 0) {
				// 流せる辺があったら流す
				int d = dfs(e.dst, t, min(f, e.cap));
				if (d > 0) {
					// 辺の残り容量を減らす
					e.cap -= d;
					// 逆辺の容量を増やす
					G[e.dst][e.rev].cap += d;
					return d;
				}
			}
		}
		// t にたどり着けなかったら0
		return 0;
	}
	int max_flow(int s, int t) {
		int flow = 0;
 
		while (1) {
			for (int i = 0; i < used.size(); i++) {
				used[i] = false;
			}
			int f = dfs(s, t, LL_HALFMAX);
			if (f == 0) {
				return flow;
			}
			flow += f;
		}
	}
 
};
 
/*
Dinic法 From libalgo O(V^2 * E)
dinic::solve(s, t) : s -> t の最大流を求める
dinic;;flow[u][v] : 辺(u, v)の流量
*/
struct dinic {
	int n, s, t;
	std::vector<int> level, prog, que;
	std::vector<std::vector<Flow>> cap, flow;
	std::vector<std::vector<int>> g;
	Flow inf;
	dinic(const Graph &graph)
		: n(graph.size()),
		cap(n, std::vector<Flow>(n)),
		flow(n, std::vector<Flow>(n)),
		g(n, std::vector<int>()),
		inf(std::numeric_limits<Flow>::max() / 8) {
		for (int i = 0; i < n; i++) {
			for (auto &e : graph[i]) {
				int u = e.src, v = e.dst;
				Flow c = e.cap;
				cap[u][v] += c;
				cap[v][u] += c;
				flow[v][u] += c;
				g[u].push_back(v);
				g[v].push_back(u);
			}
		}
	}
	// 残りを求める
	inline Flow residue(int u, int v) { return cap[u][v] - flow[u][v]; }
 
	// 実際に最大流問題を解く
	Flow solve(int s_, int t_) {
		this->t = t_, this->s = s_;
		que.resize(n + 1);
		Flow res = 0;
		// levelize() == false: bfs で s から t に到達できなかった
		while (levelize()) {
			prog.assign(n, 0);
			res += augment(s, inf);
		}
		return res;
	}
	// bfs でレベルグラフをつくる
	bool levelize() {
		int l = 0, r = 0;
		level.assign(n, -1);
		level[s] = 0;
		que[r++] = s;
		while (l != r) {
			int v = que[l++];
			if (v == t) break;
			for (const int &d : g[v]) {
				// まだレベルが決まっておらず、v -> dの辺に流せるならlevel[d] = level[v] + 1
				if (level[d] == -1 && residue(v, d) != 0) {
					level[d] = level[v] + 1;
					que[r++] = d;
				}
			}
		}
		// t に到達できるなら true を返す
		return level[t] != -1;
	}
	// dfs で実際に流してみる
	Flow augment(int v, Flow lim) {
		Flow res = 0;
		if (v == t) return lim;
		// prog[v]: dfs において、vを展開する際、vの子の何番目まで展開したかを覚えておく
		for (int &i = prog[v]; i < (int)g[v].size(); i++) {
			const int &d = g[v][i];
			// v -> d に流せない or v(流す側) の方がレベルが大きい(=深い)場合NG
			if (residue(v, d) == 0 || level[v] >= level[d]) continue;
			// 流せるなら、流せるだけ流す
			const Flow aug = augment(d, std::min(lim, residue(v, d)));
			flow[v][d] += aug;
			flow[d][v] -= aug;
			res += aug;
			lim -= aug;
			// ノードvに来ている流量を使い切ったら終わり
			if (lim == 0) break;
		}
		return res;
	}
};
 
/*
Primal-Dual法(蟻本版 / ベルマンフォード)
*/
 
class Primal_Dual_BellmanFord {
	using Cost = int;
	struct Edge {
		int src, dst;
 
		// libalgo のものに追加、メンバを追加するだけなので互換性は崩さないはず、逆辺のG[e.dstの]インデックスを保持
		int rev;
		Cost cost;
		Flow cap;
		Edge(int s, int d, int aRev, Cost aCost, Flow aCap) : src(s), dst(d), rev(aRev), cost(aCost), cap(aCap) {}
	};
 
	int V;							//頂点数
	vector<vector<Edge>> G;			// 隣接リスト
	vector<int> dist;				// 最短距離
	vector<int> prevv;				// 直前の頂点
	vector<int> preve;				// 直前の辺
	const int INF;
 
public:
	// 頂点数 n を引数にとって初期化
	Primal_Dual_BellmanFord(int n) :
		V(n),
		G(n),
		dist(n, 0),
		prevv(n, 0),
		preve(n, 0),
		INF(std::numeric_limits<int>::max() / 8) {}
	void add_edge(int src, int dst, int cap, int cost) {
		// cost は weight に入れる
		G[src].emplace_back(src, dst, G[dst].size(), cost, cap);
		G[dst].emplace_back(dst, src, G[src].size() - 1, -cost, 0);
	}
 
	int min_cost_flow(int s, int t, int f) {
		int res = 0;
		while (f > 0) {
			// ベルマンフォードによりs-t最短路をもとめる
			dist.assign(V, INF);
			dist[s] = 0;
			bool update = true;
			while (update) {
				update = false;
				for (int v = 0; v < V; v++) {
					if (dist[v] == INF) continue;
					for (int i = 0; i < G[v].size(); i++) {
						Edge& e = G[v][i];
						if (e.cap > 0 && dist[e.dst] > dist[v] + e.cost) {
							dist[e.dst] = dist[v] + e.cost;
							prevv[e.dst] = v;
							preve[e.dst] = i;
							update = true;
						}
					}
				}
			}
 
			// これ以上流せない
			if (dist[t] == INF) {
				return -1;
			}
 
			// 復元したs-t最短路に沿って流せるだけ流す
			int d = f;
			// 尻からprevvを辿っていき、流せる量を求める
			for (int v = t; v != s; v = prevv[v]) {
				// 一つ手前に戻るための辺
				Edge &e = G[prevv[v]][preve[v]];
				chmin(d, e.cap);
			}
 
			f -= d;
 
			// ここでの dist はコスト和なので、それに流す量をかけると今回見つけた最短パスに流すコストとなる。
			res += d * dist[t];
			for (int v = t; v != s; v = prevv[v]) {
				Edge &e = G[prevv[v]][preve[v]];
				e.cap -= d;
				G[v][e.rev].cap += d;
			}
		}
		return res;
	}
};
 
/*
ダイクストラ版 Primal-Dual
出典: https://ei1333.github.io/luzhiled/snippets/graph/primal-dual.html
*/
 
template< typename flow_t, typename cost_t >
struct PrimalDual {
	const cost_t INF;
 
	struct edge {
		int to;
		flow_t cap;
		cost_t cost;
		int rev;
		bool isrev;
		edge(int aTo, flow_t aCap, cost_t aCost, int aRev, bool aIsRev) : to(aTo), cap(aCap), cost(aCost), rev(aRev), isrev(aIsRev) {}
	};
	vector< vector< edge > > graph;
	vector< cost_t > potential, min_cost;
	vector< int > prevv, preve;
 
	PrimalDual(int V) : graph(V), INF(numeric_limits< cost_t >::max()) {}
 
	void add_edge(int from, int to, flow_t cap, cost_t cost) {
		graph[from].emplace_back(to, cap, cost, (int)graph[to].size(), false);
		graph[to].emplace_back(from, 0, -cost, (int)graph[from].size() - 1, true);
	}
 
	cost_t min_cost_flow(int s, int t, flow_t f) {
		int V = (int)graph.size();
		cost_t ret = 0;
		using Pi = pair< cost_t, int >;
		priority_queue< Pi, vector< Pi >, greater< Pi > > que;
		potential.assign(V, 0);
		preve.assign(V, -1);
		prevv.assign(V, -1);
 
		while (f > 0) {
			min_cost.assign(V, INF);
			que.emplace(0, s);
			min_cost[s] = 0;
			while (!que.empty()) {
				Pi p = que.top();
				que.pop();
				if (min_cost[p.second] < p.first) continue;
				for (int i = 0; i < graph[p.second].size(); i++) {
					edge &e = graph[p.second][i];
					cost_t nextCost = min_cost[p.second] + e.cost + potential[p.second] - potential[e.to];
					if (e.cap > 0 && min_cost[e.to] > nextCost) {
						min_cost[e.to] = nextCost;
						prevv[e.to] = p.second, preve[e.to] = i;
						que.emplace(min_cost[e.to], e.to);
					}
				}
			}
			if (min_cost[t] == INF) return -1;
			for (int v = 0; v < V; v++) potential[v] += min_cost[v];
			flow_t addflow = f;
			for (int v = t; v != s; v = prevv[v]) {
				addflow = min(addflow, graph[prevv[v]][preve[v]].cap);
			}
			f -= addflow;
			ret += addflow * potential[t];
			for (int v = t; v != s; v = prevv[v]) {
				edge &e = graph[prevv[v]][preve[v]];
				e.cap -= addflow;
				graph[v][e.rev].cap += addflow;
			}
		}
		return ret;
	}
 
	void output() {
		for (int i = 0; i < graph.size(); i++) {
			for (auto &e : graph[i]) {
				if (e.isrev) continue;
				auto &rev_e = graph[e.to][e.rev];
				cout << i << "->" << e.to << " (flow: " << rev_e.cap << "/" << rev_e.cap + e.cap << ")" << endl;
			}
		}
	}
};
 
class lca {
public:
	int n, segn;
	vector<int> path;		// 蟻本の vs、オイラーツアーを保持
	vector<int> depth;		// 蟻本の depth、path[i] であるノードの深さを保持
	vector<int> in_order;	// 蟻本の id、ノードiがオイラーツアーで最初に出てくるインデックスを保持
	vector<pair<int, int>> dat;
	const std::pair<int, int> INF = std::make_pair(1000000000, 1000000000);
 
	lca(const Graph& g, int root) : n(g.size()), path(n * 2 - 1), depth(n * 2 - 1), in_order(n) {
		int k = 0;
		dfs(g, root, -1, 0, k);
 
		// セグ木を構築、持つのはpair(depth, index) => depth が最小となる index がわかる 
		for (segn = 1; segn < n * 2 - 1; segn <<= 1);
 
		dat.assign(segn * 2, INF);
		for (int i = 0; i < (int)depth.size(); ++i) dat[segn + i] = std::make_pair(depth[i], i);
		for (int i = segn - 1; i >= 1; --i) dat[i] = min(dat[i * 2], dat[i * 2 + 1]);
	}
 
	int get(int u, int v) const {
		int l = std::min(in_order[u], in_order[v]);
		int r = std::max(in_order[u], in_order[v]) + 1;
		return path[range_min(1, segn, l, r).second];
	}
 
	void dfs(const Graph& g, int v, int p, int d, int& k) {
		// k: オイラーツアーの何番目かを保持する変数
		in_order[v] = k;
		path[k] = v;
		depth[k++] = d;
		for (auto &e : g[v]) {
			if (e.dst != p) {
				dfs(g, e.dst, v, d + 1, k);
				// ここに来た時はノードvの子であるe.dstを展開し終わってvに戻ってきたときなので、再度 path と depth に記録する
				path[k] = v;
				depth[k++] = d;
			}
		}
	}
 
	// v : いまみてるノード、w: 今見てるノードに対応する区間長 l: ? r: ?
	pair<int, int> range_min(int v, int w, int l, int r) const {
		if (r <= l || w == 0) return INF;
		if (r - l == w)
			return dat[v];
 
		int m = w / 2;
		auto rmin = range_min(v * 2, m, l, std::min(r, m));
		auto lmin = range_min(v * 2 + 1, m, std::max(0LL, l - m), r - m);
 
		return min(rmin, lmin);
	}
};
 
// int における ceil と floor、負数対応(a / b の ceil, floor)
int64_t intceil(int64_t a, int64_t b) {
	int sign_a = (a > 0) - (a < 0);
	int sign_b = (b > 0) - (b < 0);
 
	if (sign_a == sign_b) {
		return (a + b - sign_b) / b;
	}
	else {
		return a / b;
	}
}
int64_t intfloor(int64_t a, int64_t b) {
	int sign_a = (a > 0) - (a < 0);
	int sign_b = (b > 0) - (b < 0);
 
	if (sign_a == sign_b) {
		return a / b;
	}
	else {
		return (a - b + sign_b) / b;
	}
}
 
class Point {
public:
	int y, x;
	Point() { y = x = 0; }
	Point(int y0, int x0) {
		y = y0;
		x = x0;
	}
	Point operator+(const Point& p) const { return Point(y + p.y, x + p.x); }
	Point operator-(const Point& p) const { return Point(y - p.y, x - p.x); }
	Point operator*(int a) const { return Point(y * a, x * a); }
	long long length2() const { return y * (long long)y + x * (long long)x; }
	long long dist2(const Point& p) const {
		return (y - p.y) * (long long)(y - p.y) + (x - p.x) * (long long)(x - p.x);
	}
	long long dot(const Point& p) const {
		return y * (long long)p.y + x * (long long)p.x;  // |a|*|b|*cosθ
	}
	long long cross(const Point& p) const {
		return x * (long long)p.y - y * (long long)p.x;  // |a|*|b|*sinθ
	}
 
	static bool Sorter(const Point& p1, const Point& p2) {
		bool a = p1.y > 0 || (p1.y == 0 && p1.x >= 0);
		bool b = p2.y > 0 || (p2.y == 0 && p2.x >= 0);
		if (a != b) return a;
		long long c = p2.x * (long long)p1.y;
		long long d = p1.x * (long long)p2.y;
		if (c != d) return c < d;
		return p1.length2() < p2.length2();
	}
};
 
 
int solve(ostringstream& aout, long long S);
void solve_TLE(ostringstream& aout, long long S);
 
class StressTest {
private:
	mt19937 m_RandEngine;
	bool judge_case(long long S) {
		ostringstream fast, tle;
		solve(fast, S);
		solve_TLE(tle, S);
		if (fast.str() == tle.str()) {
			return true;
		}
		else {
			return false;
		}
	}
	// [l, l+1, ... r] の数列を生成し、シャッフルする
	vector<int> create_range_permutation(int l, int r) {
		vector<int> ret;
		for (int i = l; i <= r; i++) {
			ret.push_back(i);
		}
		shuffle(ret.begin(), ret.end(), m_RandEngine);
		return ret;
	}
	// [1, n] の順列を生成する
	vector<int> create_permutation(int n) {
		create_range_permutation(1, n);
	}
	// 範囲が[l, r] でサイズが n の数列を生成する
	vector<int> create_random_sequence(int l, int r, int n) {
		uniform_int_distribution<> randLR(l, r);
		vector<int> ret;
		for (int i = 0; i < n; i++) {
			ret.push_back(randLR(m_RandEngine));
		}
		return ret;
	}
 
	/*
	* 頂点数 n, 辺数 m で自己ループと多重辺のない無向グラフを生成
	* 慣習的に頂点番号が1-indexed な AtCoder で 1-n の頂点が使えるようにするため n+1 頂点のグラフを生成し、0番を無視することとする
	* weighted を true にすると重み付き、maxWeight で最大重みを指定
	* 連結でないグラフが出力される可能性があることに注意する
	*/
	Graph create_undirected_graph(int n, int m, bool weighted = false, int maxWeight = 10) {
		Graph ret(n + 1);
		set<pair<int, int>> used;
		uniform_int_distribution<> randNode(1, n);
		uniform_int_distribution<> randWeight(1, maxWeight);
		while (used.size() < m * 2) {
			int src = randNode(m_RandEngine);
			int dst = randNode(m_RandEngine);
 
			// 自己ループ、多重辺判定
			if (used.count(make_pair(src, dst)) == 0 && used.count(make_pair(dst, src)) == 0 && src != dst) {
				used.insert(make_pair(src, dst));
				used.insert(make_pair(dst, src));
				add_edge(ret, src, dst, weighted ? randWeight(m_RandEngine) : 1);
			}
		}
		return ret;
	}
 
	/*
	* 頂点数 n, 辺数 m で自己ループと多重辺のない有向グラフを生成
	* 慣習的に頂点番号が1-indexed な AtCoder で 1-n の頂点が使えるようにするため n+1 頂点のグラフを生成し、0番を無視することとする
	* weighted を true にすると重み付き、maxWeight で最大重みを指定
	* 連結でないグラフが出力される可能性があることに注意する
	*/
	Graph create_directed_graph(int n, int m, bool weighted = false, int maxWeight = 10) {
		Graph ret(n + 1);
		set<pair<int, int>> used;
		uniform_int_distribution<> randNode(1, n);
		uniform_int_distribution<> randWeight(1, maxWeight);
		while (used.size() < m) {
			int src = randNode(m_RandEngine);
			int dst = randNode(m_RandEngine);
 
			// 自己ループ、多重辺判定
			if (used.count(make_pair(src, dst)) == 0 && src != dst) {
				used.insert(make_pair(src, dst));
				add_arc(ret, src, dst, weighted ? randWeight(m_RandEngine) : 1);
			}
		}
		return ret;
	}
 
	/*
	* 頂点数nの木(無向)を生成します。
	*/
	Graph create_tree(int n, bool weighted = false, int maxWeight = 10) {
		Graph ret(n + 1);
		uf_tree uf(n + 1);
		int cnt = 0;
 
		uniform_int_distribution<> randNode(1, n);
		uniform_int_distribution<> randWeight(1, maxWeight);
 
		while (cnt < n - 1) {
			int n1 = randNode(m_RandEngine);
			int n2 = randNode(m_RandEngine);
			if (n1 != n2 && !uf.is_same(n1, n2)) {
				cnt++;
				add_edge(ret, n1, n2, weighted ? randWeight(m_RandEngine) : 1);
			}
		}
	}
public:
	StressTest(int seed) :
		m_RandEngine(seed) {}
	void test() {
		while (1) {
			// TODO: generate random case
			//if (!judge_case(S)) {
			// TODO: output case
			//break;
			//}
		}
	}
};
 
int solve(ostringstream& aout, long long S) {
	// 区間sum, 区間chg
	auto myAdd = [](int a, int b) -> int {
		return (a + b) % mod;
	};
	auto myChg = [](int a, int b) -> int {
		return b;
	};
 
	// 区間更新はしない
	vector<SegmentTree<int, int>> dp(S + 10, SegmentTree<int, int>(S + 10, myAdd, myChg, myChg, 0, LL_HALFMAX));
 
	// dp[0][0] = 1;
	dp[0].update(0, 1);
 
	// もらう
	REPS(i, S) {
		REP(j, 3, S + 1) {
			int chg = dp[i - 1].query(0, j - 3 + 1);
			dp[i].update(j, chg);
		}
	}
 
	int ans = 0;
	REPS(i, S) {
		int add = dp[i].query(S);
		ans += add;
		ans %= mod;
	}
 
	return ans;
}
 
void solve_TLE(ostringstream& aout, long long S) {
 
}
 
const int modu = 998244353;
 
int op(int a, int b) {
	return (a + b) % modu;
}
 
int e() {
	return 0;
}
 
int mapping(pair<int, int> f, int x) {
	int ret = x * f.first;
	ret %= modu;
	ret += f.second;
	ret %= modu;
	
	return ret;
}
 
int g(int x, pair<int, int> f) {
	int ret = x * f.first;
	ret %= modu;
	ret += f.second;
	ret %= modu;
 
	return ret;
}
 
pair<int, int> composition(pair<int, int> f, pair<int, int> g) {
	// (b, c) と (d, e) をマージ
	int b = f.first, c = f.second, d = g.first, e = g.second;
 
	return mp((b * d) % modu, (c*d + e) % modu);
}
 
pair<int, int> id() {
	return mp(1, 0);
}
 
signed main() {
	int N, Q;
	cin >> N >> Q;
	vector<int> a(N);
	rep(i, N) {
		cin >> a[i];
	}
 
	vector<vector<int>> query;
	rep(i, Q) {
		int qtype;
		cin >> qtype;
		if (qtype == 0) {
			int l, r, b, c;
			cin >> l >> r >> b >> c;
			query.push_back({ 0, l, r, b, c });
		}
		else {
			int l, r;
			cin >> l >> r;
			query.push_back({ 1, l, r });
		}
	}
	// seg2 でいうところの p わたせなくない？
	atcoder::lazy_segtree<int, op, e, pair<int, int>, mapping, composition, id> seg(a);
	auto p = [](pair<int, int> a, int d) -> pair<int, int> {
		return mp(a.first, a.second*d);
	};
 
	SegmentTree<int, pair<int, int>> seg2(N, op, g, composition, 0, mp(1, 0), a, p);
 
	auto print = [&]() {
		rep(i, N) {
			cout << seg.get(i) << ", ";
		}
		cout << "\n";
	};
 
	rep(i, Q) {
		if (query[i][0] == 1) {
			//int ans = seg.prod(query[i][1], query[i][2]);
			int ans = seg2.query(query[i][1], query[i][2]);
			cout << ans << "\n";
		}
		else {
			int l = query[i][1], r = query[i][2], b = query[i][3], c = query[i][4];
			//seg.apply(l, r, mp(b, c));
			seg2.update(l, r, mp(b, c));
		}
		//print();
	}
 
	return 0;
}