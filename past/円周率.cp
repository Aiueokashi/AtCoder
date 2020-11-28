/*********************************************
 * 円周率計算 by Klingenstierna の公式
 *********************************************/
#include <iostream>  // for cout
#include <math.h>    // for pow()
#include <stdio.h>   // for printf()
#include <time.h>    // for clock()
#define FNAME "pi_klingenstierna.txt"

using namespace std;


class Calc
{
    // 各種変数
    int l, l1, n;    // 計算桁数、配列サイズ、計算項数
    int i, k;        // LOOP インデックス
    int cr, br;      // 繰り上がり、借り
    long w;          // 被乗数、被除数ワーク
    long r;          // 剰余ワーク
    clock_t t1, t2;  // 計算開始・終了CPU時刻
    double tt;       // 計算時間
    FILE *out_file;  // 結果出力ファイル

    public:
        // コンストラクタ
        Calc(int);
        // 計算
        void calc();
        // ロング + ロング
        void ladd(int *, int *, int *);
        // ロング - ロング
        void lsub(int *, int *, int *);
        // ロング / ショート
        void ldiv(int *, int, int *);
        // 結果出力
        void display(double, int *);
};


Calc::Calc(int x)
{
    l  = x;                // 計算桁数
    l1 = (l / 8) + 1;      // 配列サイズ
    n  = (l + 1) / 2 + 1;  // 計算項数
}


void Calc::calc()
{
    // 計算開始時刻
    t1 = clock();

    // 配列宣言
    int s[l1 + 1], a[l1 + 1], b[l1 + 1], c[l1 + 1], q[l1 + 1];

    // 配列初期化
    for (k = 0; k <= l1 + 1; k++)
        s[k] = a[k] = b[k] = c[k] = q[k] = 0;

    // Klingenstierna の公式
    a[0] = 32 *  10;
    b[0] =  4 * 239;
    c[0] = 16 * 515;
    for (k = 1; k <= n; k++) {
        ldiv(a,  10 *  10, a);
        ldiv(b, 239 * 239, b);
        ldiv(c, 515 * 515, c);
        lsub(a, b, q);
        lsub(q, c, q);
        ldiv(q, 2 * k - 1, q);
        if ((k % 2) != 0)
            ladd(s, q, s);
        else
            lsub(s, q, s);
    }

    // 計算終了時刻
    t2 = clock();

    // 計算時間
    tt = (double)(t2 - t1) / CLOCKS_PER_SEC;

    // 結果出力
    display(tt, s);
}


void Calc::ladd(int a[], int b[], int c[])
{
    cr = 0;
    for (i = l1 + 1; i >=0; i--) {
        c[i] = a[i] + b[i] + cr;
        if (c[i] < 100000000) {
            cr = 0;
        } else {
            c[i] -= 100000000;
            cr = 1;
        }
    }
}


void Calc::lsub(int a[], int b[], int c[])
{
    br = 0;
    for (i = l1 + 1; i >=0; i--) {
        c[i] = a[i] - b[i] - br;
        if (c[i] >= 0) {
            br = 0;
        } else {
            c[i] += 100000000;
            br = 1;
        }
    }
}


void Calc::ldiv(int d[], int e, int f[])
{
    r = 0;
    for (i = 0; i < l1 + 1; i++) {
        w = d[i];
        f[i] = (w + r) / e;
        r = ((w + r) % e) * 100000000;
    }
}


void Calc::display(double tt, int s[])
{
    printf("** Pi Computation with the Klingenstierna formula method **\n") ;
    printf("   Digits = %d.\n", l);
    printf("   Time   = %f seconds\n", tt);

    // ファイル出力
    out_file = fopen(FNAME,"w");
    fprintf(out_file, "** Pi Computation with the Klingenstierna formula method **\n") ;
    fprintf(out_file, "   Digits = %d.\n", l);
    fprintf(out_file, "   Time   = %f seconds.\n\n", tt);
    fprintf(out_file, "          %d.\n", s[0]);
    for (i = 1; i < l1; i++) {
        if (i % 10 == 1) fprintf(out_file, "%08d:", (i - 1) * 8 + 1);
        fprintf(out_file, " %08d", s[i]);
        if (i % 10 == 0) fprintf(out_file, "\n");
    }
    fprintf(out_file, "\n\n");
}


int main()
{
    int n;  // 計算桁数

    try
    {
        // 計算桁数入力
        printf("Please input number of Pi Decimal-Digits : ");
        scanf("%d", &n);

        // 計算クラスインスタンス化
        Calc objCalc(n);

        // 円周率計算
        objCalc.calc();
    }
    catch (...) {
        cout << "例外発生！" << endl;
        return -1;
    }

    // 正常終了
    return 0;
}
