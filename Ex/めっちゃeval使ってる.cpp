#pragma GCC optimize("Ofast")
#include<bits/stdc++.h>
using namespace std;
//#include<boost/multiprecision/cpp_int.hpp>
//#include<boost/multiprecision/cpp_dec_float.hpp>
//namespace mp=boost::multiprecision;
//#define mulint mp::cpp_int
//#define mulfloat mp::cpp_dec_float_100
struct __INIT{__INIT(){cin.tie(0);ios::sync_with_stdio(false);cout<<fixed<<setprecision(15);}} __init;
#define max3(a,b,c) max(a,max(b,c))
#define min3(a,b,c) min(a,min(b,c))
constexpr int MOD=1000000007;
//constexpr int MOD=998244353;
#define INF (1<<30)
#define LINF (lint)(1LL<<56)
#define endl "\n"
#define rep(i,n) for(lint (i)=0;(i)<(n);(i)++)
#define reprev(i,n) for(lint (i)=(n-1);(i)>=0;(i)--)
#define Flag(x) (1<<(x))
#define Flagcount(x) __builtin_popcountll(x)
#define pint pair<int,int>
#define pdouble pair<double,double>
#define plint pair<lint,lint>
#define fi first
#define se second
typedef long long lint;
int dx[8]={1,1,0,-1,-1,-1,0,1};
int dy[8]={0,1,1,1,0,-1,-1,-1};
const int MAX_N=2e5+5;
//struct edge{lint to,num;};
//vector<int> bucket[MAX_N/1000];

lint segdata[(1<<19)-1];
lint lazy[(1<<19)-1];

lint neutral(){
    return 0; //単位元
}
lint calc(lint a,lint b){
    return max(a,b); //演算 デフォルトでは最小値
}
void seginit(){
    rep(i,(1<<19)-1) segdata[i]=neutral();
}
void eval(lint k,lint l,lint r){
    if(lazy[k]){
        segdata[k]+=lazy[k]/(r-l);
        if(r-l>1) lazy[2*k+1]+=lazy[k]/2,lazy[2*k+2]+=lazy[k]/2;
        lazy[k]=0;
    }
}
void add(lint a,lint b,lint x,lint k,lint l,lint r){
    eval(k,l,r);
    if(r<=a || b<=l) return;
    if(a<=l && r<=b){
        lazy[k]+=(r-l)*x;
        eval(k,l,r);
    }
    else{
        add(a,b,x,2*k+1,l,(l+r)/2);
        add(a,b,x,2*k+2,(l+r)/2,r);
        segdata[k]=calc(segdata[2*k+1],segdata[2*k+2]);
    }
}
lint seg(lint a,lint b,lint k,lint l,lint r){
    eval(k,l,r);
    if(r<=a || b<=l) return neutral();
    if(a<=l && r<=b) return segdata[k];
    else{
        lint left=seg(a,b,k*2+1,l,(l+r)/2);
        lint right=seg(a,b,k*2+2,(l+r)/2,r);
        return calc(left,right);
    }
}
lint query(lint a,lint b){
    return seg(a,b+1,0,0,1<<18); //a~b(閉区間）の範囲での演算
}
void rangeadd(lint a,lint b,lint x){ //a~b(閉区間)にxを加算
    add(a,b+1,x,0,0,1<<18); 
    return;
}

int main(void){
    lint N,W;
    cin >> N >> W;
    seginit();
    rep(i,N){
        lint s,t,p;
        cin >> s >> t >> p;
        rangeadd(s,t-1,p);
    }
    lint maxi=query(0,200010);
    if(maxi<=W) cout << "Yes" << endl;
    else cout << "No" << endl;
}