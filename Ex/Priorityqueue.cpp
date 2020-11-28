#include "bits/stdc++.h"
#include <ext/pb_ds/assoc_container.hpp>
using namespace __gnu_pbds;
using namespace std;

#define ub              upper_bound
#define lb              lower_bound
#define isrt            insert
#define clr             clear
#define rsz             resize
#define ff              first
#define ss              second
#define lli             long long int
#define pb              push_back
#define pf              push_front
#define mkp              make_pair
#define pii             pair<lli,lli>
#define vi              vector<int>
#define mii             map<lli,lli>
#define pqb             priority_queue<int>
#define pqs             priority_queue<int,vi,greater<int> >
#define setbits(x)      __builtin_popcountll(x)
#define zrobits(x)      __builtin_ctzll(x)
#define mod             1000000007
#define INF             1e9
#define ps(x,y)         fixed<<setprecision(y)<<x
#define mk(arr,n,type)  type *arr=new type[n];
#define w(x)            int x; cin>>x; while(x--)
#define all(v)          v.begin(),v.end()
#define endl            "\n"
const double PI = acos(-1);

typedef tree<int, null_type, less<int>, rb_tree_tag, tree_order_statistics_node_update> pbds;

lli maxm(lli a, lli b)
{
  return (a >= b ? a : b);
}


lli minm(lli a, lli b)
{
  return (a <= b ? a : b);
}

lli power(lli x, lli y, lli p)
{
  lli res = 1;

  x = x % p;
  if (x == 0)
    return 0;

  while (y > 0)
  {

    if (y & 1)
      res = (res * x) % p;


    y = y >> 1;
    x = (x * x) % p;
  }
  return res;
}

lli nCrmodp(lli n, lli r, lli p)
{
  if (r == 0)
    return 1;
  if (r == 1)
    return n;


  lli pro = 1, k = n;
  while (k > r)
  {
    pro *= k;
    pro %= p;
    --k;
  }

  k = 2;
  while (k <= (n - r))
  {
    pro *= power(k, (p - 2), p);
    pro %= p;
    ++k;
  }

  return pro;
}


lli first_big(lli n, lli k, vector<lli> arr)
{
  lli l, r, ans = n;
  l = 0;
  r = n - 1;

  while (l <= r) {
    lli m = l + (r - l) / 2;

    if (arr[m] > k) {
      ans = m;
      r = m - 1;
    }

    else
      l = m + 1;
  }

  return ans;
}

lli last_small(lli n, lli k, vector<lli> arr)
{
  lli l, r, ans = -1;
  l = 0;
  r = n - 1;

  while (l <= r)
  {
    lli mid = (l + r) / 2;

    if (arr[mid] <= k)
    {
      ans = mid;
      l = mid + 1;
    }
    else
      r = mid - 1;
  }

  return ans;
}

void c_p_c()
{
  ios_base::sync_with_stdio(0); cin.tie(0); cout.tie(0);
#ifndef ONLINE_JUDGE
  freopen("input.txt", "r", stdin);
  freopen("output.txt", "w", stdout);
#endif
}
#define MAXN 110
vector<double>cnt(MAXN);
vector<int>adj[MAXN];
vector<bool>vis;

void dfs(int s)
{
  cnt[s] += 1;
  vis[s] = 1;
  for (auto it : adj[s])
    if (!vis[it])
      dfs(it);
}

int main()
{
  c_p_c();

  int n, x;
  cin >> n;
  string str;
  for (int i = 0; i < n; ++i)
  {
    cin >> str;
    for (int j = 0; j < n; ++j)
    {
      if (str[j] == '1')
        adj[i].pb(j);
    }
  }

  for (int i = 0; i < n; ++i)
  {
    vis.clr();
    vis.rsz(n + 1, false);
    dfs(i);
  }

  double ans = 0.0000000000;
  for (int i = 0; i < n; ++i)
    ans += (1 / cnt[i]);
  cout << ps(ans, 15) << endl;
  return 0;
}