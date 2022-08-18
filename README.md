# MSYS2 MinGW64 環境の C++ 開発で LAPACK を使ったエルミート行列対角化をする

## はじめに

Windows における C++ 開発で LAPACK を使う場合のための覚書として記す。
MSYS2 を使い、MinGW64 の gcc コンパイラ (g++) でコンパイルすることを想定する。
ソースコードとコンパイルオプションを明示して、最低限動かすために必要な情報を残す。
また、REFERENCE BLAS と OpenBLAS の速度比較も行う。

## 環境

- Windows 10
- gcc 12.1.0 (Rev3, Built by MSYS2 project)
- Intel Core i7-8500Y 1.50 GHz
- LAPACK v3.10.1
- REFERENCE BLAS v3.10.0
- OpenBlas v0.3.21

LAPACK、REFERENCE BLAS、OpenBLAS のインストールは以下のコマンドでできる。REFERENCE BLAS は LAPACK のパッケージに同梱されている。

```
pacman -S mingw-w64-x86_64-lapack mingw-w64-x86_64-openblas
```

REFERENCE BLAS か OpenBLAS かについては基本的に OpenBLAS を使えばよい。
コンパイルオプションとして `-lopenblas` を指定すると OpenBLAS がリンクされる。
あまり、使う理由はないが REFERENCE BLAS を使いたい場合は `-lblas` を代わりに使う。

## BLAS と LAPACK

[BLAS](https://netlib.org/blas/) は Basic Linear Algebra Subprograms の略で線形代数に必要な基礎的な演算を高速に行う数値計算ライブラリである。
これには内積やノルムの計算、ベクトルや行列の乗算などが含まれている。

BLAS にはいくつかのバリエーションがあり、[Netlib](https://netlib.org/blas/) で公開されているものが公式の BLAS となる。
これは特に REFERENCE BLAS とも呼ばれている。
つまり、すべての BLAS の「お手本」である。
ドキュメントが整っており、どのバリエーションを使う場合も Netlib のドキュメントを確認すればよい。
REFERENCE BLAS は高速な計算アルゴリズムを示すためのもので、コードは可読性に重点を置かれている。
これに対して、より速度に重点をおいて種々の最適化を施された BLAS を最適化 BLAS (Optimized BLAS) という。
最適化 BLAS には現在有名なものに MKL (Intel Math Kernel Library)、[OpenBLAS](https://github.com/xianyi/OpenBLAS) などがある。
（より正確に言えば MKL は BLAS 以外にも LAPACK 相当の機能も含んでおり、更に多くの数値計算ルーチンを包含している。）

[LAPACK](https://netlib.org/lapack/) は Linear Algebra PACKage の略でより一般的な線形代数の問題を解くための数値計算ライブラリである。
これには連立1次方程式、固有値問題、特異値分解をはじめとして、それらを解くために必要となる多くの補助ルーチンが含まれる。

LAPACK は BLAS のルーチンを使って計算を行うために必ずこれが必要になる。
もちろん REFERENCE BLAS でもよいが、最適化 BLAS を使うことで LAPACK の計算も高速になる。

## エルミート行列対角化のためのLAPACKドライバー

LAPACK は 1000 を超えるルーチンからなり、それらを組み合わせて線形代数の諸問題を数値的に計算するためのライブラリである。
多くのルーチンを組み合わせて一般的な問題（連立1次方程式、固有値問題、特異値分解など）を解くためのルーチンは特別にドライバー(driver routine)と呼ばれている。

エルミート行列の対角化には `zheev` 系のドライバーを用いる。
「系」というのには理由があり、同じ問題でも実装によっていくつかの選択肢が存在するためである。

- `zheev` 
    - Simple driver
    - 最も単純で初期から存在するドライバー
- `zheevx`
    - Expert driver
    - Simple driver とアルゴリズムは同じだが、固有値を指定の範囲・個数まで求めて計算を打ち切るための仕組みが提供されている
        - 固有値が全部必要ない場合には計算が早く終了することが期待できる
    - 基本的には Simple driver の上位互換と考えて良い
- `zheevd`
    - Divide and conquer アルゴリズムによって固有値問題を解くドライバー
    - 上の2つのドライバーに比べて高速
- `zheevr`
    - Relatively Robust Representation アルゴリズムによって固有値問題を解くドライバー
    - 比較的新しく、v3.10.1 時点ではまだ実装は完全ではないようである
        - 具体的に言うと引数 ABSTOL のドキュメントに少し書いてある
        - しかし、十分に使えるようではある
    - 高速かつ Expert driver と同じく一部の固有値のみを計算することができる

ここでは例として 1000x1000 のエルミート行列をタブ区切りファイルから読み込み対角化するプログラムを示す。
行列の実部と虚部はそれぞれ別のファイルに格納されており、2つのファイルを読み込んでエルミート行列を作る。
LAPACK はもともと FORTRAN 言語でかかれており、その呼び出し規則に従わなければならない。
つまり、行列は **Column major** で連続したメモリ上に配置されていなければならない。
C/C++ のプログラムでは1000\*1000=1000000要素の1次元配列としてメモリを確保し、Column majorで数値を格納する。

---

コンピュータのメモリ空間は1次元的な配置になっているので、
2次元の行列をメモリ上に展開するには何らかの規則が必要である。
Column majorとは以下のような行列

$$
\begin{matrix}
a_{11} & a_{12} & a_{13}
a_{21} & a_{22} & a_{23}
a_{31} & a_{32} & a_{33}
\begin{matrix}
$$

を配列Aに次のように配置することである。

$$
A[0] = a_{11}
A[1] = a_{21}
A[2] = a_{31}
A[3] = a_{12}
A[4] = a_{22}
A[5] = a_{32}
A[6] = a_{13}
A[7] = a_{23}
A[8] = a_{33}
$$

Cにおける配列（C++における生配列）はメモリ上に連続的にデータを配置する。
つまり、各列 (column) をメモリ上に連続的に配置している。

Column major と対になる言葉として Row major というものがある。
これは同じ行列を次のように配置することである。

$$
A[0] = a_{11}
A[1] = a_{12}
A[2] = a_{13}
A[3] = a_{21}
A[4] = a_{22}
A[5] = a_{23}
A[6] = a_{31}
A[7] = a_{32}
A[8] = a_{33}
$$

つまり、各行 (row) をメモリ上に連続的に配置している。

---

C/C++ から FORTRAN でかかれたサブルーチンを呼び出すためには、FORTRAN でのサブルーチン名の最後にアンダースコア(`_`)をつけなければいけない。
つまり、`zheev` ドライバーを呼び出すためには `zheev_` と記述する。

下記のプログラムでは手抜きのためにファイルにおける行を列として読み込んでいる。
本来であればファイルに書いてあるとおりの行列を対角化したければ、転置しなければならない。
（とはいえ、エルミート行列なので転地しても固有値は変わらないはずである。）

また、計算には上三角行列あるいは下三角行列のどちらかしか使わないので、対角化する行列を作る際は全部の要素を埋める必要はない。
つまり、エルミート対称性から対角要素の値は決まるので、原理的に計算に必要ないのである。
少しでもプログラムを高速化するためには重要である。

### zheevドライバー

ドキュメント: https://netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_gaf23fb5b3ae38072ef4890ba43d5cfea2.html

引数 `A` に対角化する行列を与えると `W` に与えた配列に計算された固有値が格納される。
配列 `A` は対角化計算の過程で破壊的変更を受ける。
つまり計算が終了した時点で `A` は与えた配列と同一の内容ではないので注意。
`JOBZ` に `'V'` を指定した場合は `A` の第*i*列が `W` の*i*番目の要素に対応した固有ベクトルになっている。

`WORK`, `RWORK` は計算のために必要な作業用のメモリ領域であり、`LWORK` は `WORK` の配列サイズである。
`RWORK` は NxN の行列を対角化する場合、3N-2 のサイズ決め打ちでよい。
`WORK` の最適なサイズは `LWORK` に -1 を与えて `zheev` を呼び出すと計算でき、`WORK` の最初の要素に格納される。
つまり、固有値を得るためには 1. `LWORK` の決定 → 2. 固有値の計算というように2回 `zheev` を呼び出す。
それぞれ最適サイズより大きくても問題ない。
最適サイズより小さい場合、計算速度が悪化することが予想される。

何度も固有値問題を解く場合、`WORK`, `RWORK` はサイズが十分な限り使い回せる。
すなわち、ブロック対角化のように何回か対角化をする場合、最も大きい行列で `WORK`, `RWORK` を用意しておけば何度も確保し直す必要はない。

計算終了後に `INFO` に計算の終了状態に関する情報が返される。
`INFO = 0` なら正常終了。
`INFO < 0` なら `INFO` 番目の引数に問題があり、異常終了。
`INFO > 0` なら計算が収束しなかったことを意味する。

- ソースコード `test_LAPACK_zheev.cpp`

```c++
#include <fstream>
#include <sstream>
#include <string>
#include <complex.h>

// 読み取る行列のサイズ
#define SIZE 1000

// im を虚数単位とする
#define im _Complex_I
#undef I

extern "C" {
  void zheev_(const char& JOBZ, const char& UPLO, const int& N, double
        _Complex* A, const int& LDA, double* W, double _Complex* WORK,
        const int& LWORK, double* RWORK, int& INFO);
};


// ２つのファイルを読み込んで Hermitian matrix (a) を作る
void readtsv(double _Complex* a){
    const auto real_file = "real.tsv";
    const auto imag_file = "imag.tsv";

    std::ifstream ifs_real(real_file);
    std::string line;
    int i = 0;
    while (getline(ifs_real, line)) {
        std::istringstream stream(line);
        std::string realstr;
        while (getline(stream, realstr, '\t')) {
            a[i] = std::stod(realstr);
            i++;
        }
    }

    std::ifstream ifs_imag(imag_file);
    i = 0;
    while (getline(ifs_imag, line)) {
        std::istringstream stream(line);
        std::string imagstr;
        while (getline(stream, imagstr, '\t')) {
            a[i] += std::stof(imagstr) * im;
            i++;
        }
    }
}


int main(){
    char jobz = 'V';
    char uplo = 'U';
    int n = SIZE;
    int lda = SIZE;
    double w[SIZE];
    double rwork[3*SIZE - 2];
    int info;

    // 行列はサイズが大きくなりうるので new を使ってヒープに確保する
    auto a = new double _Complex[SIZE*SIZE];
    readtsv(a);

    // 最適な LWORK の計算
    double _Complex work_temp[1];
    zheev_(jobz, uplo, n, a, lda, w, work_temp, -1, rwork, info);
    int lwork = (int) lround(creal(work_temp[0]));
    auto work = new double _Complex[lwork];

    zheev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);

    printf("Exit with info = %d\n", info);
    printf("w[0] = %f\n", w[0]);
    printf("w[1] = %f\n", w[1]);
    printf("...\n");
    printf("w[%d] = %f\n", SIZE-2, w[SIZE-2]);
    printf("w[%d] = %f\n", SIZE-1, w[SIZE-1]);

    delete[] a;
    delete[] work;
    return 0;
}
```

- コンパイルコマンド

```sh
g++ test_LAPACK_zheev.cpp -o a.out -llapack -lopenblas
```


### zheevxドライバー

ドキュメント: https://netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_gaabef68a9c7b10df7aef8f4fec89fddbe.html

引数 `A` に対角化する行列を与えると `W` に与えた配列に計算された固有値が格納される。
配列 `A` は対角化計算の過程で破壊的変更を受ける。
つまり計算が終了した時点で `A` は与えた配列と同一の内容ではないので注意。
`JOBZ` に `'V'` を指定した場合は `Z` の第*i*列が `W` の*i*番目の要素に対応した固有ベクトルになっている。

`zheevx` ドライバーは与えられた行列 `A` の固有値の内の一部のみを計算する機能を持つ。
すべての固有値が必要ではない場合には、計算を打ち切って高速化することができる。
実際に見つかった固有値の個数は `M` に格納される。

- `RANGE` に `'A'` を指定した場合、すべての固有値を計算する。
- `RANGE` に `'V'` を指定した場合、(`VL`, `VU`] の半開区間に入る固有値のみを求める。
- `RANGE` に `'I'` を指定した場合 `IL` 番目から `IU` 番目の固有値のみを求める。

`ABSTOL` は固有値の許容計算誤差に関する引数だが、特に理由がなければ最良の精度を得るために 2\*DLMCH('S') を与えるのがよい。

`WORK`, `RWORK`, `IWORK` は計算のために必要な作業用のメモリ領域であり、`LWORK` は `WORK` の配列サイズである。
`RWORK` は NxN の行列を対角化する場合、7\*N のサイズ決め打ちでよい。
`IWORK` は NxN の行列を対角化する場合、5\*N のサイズ決め打ちでよい。
`WORK` の最適なサイズは `LWORK` に -1 を与えて `zheev` を呼び出すと計算でき、`WORK` の最初の要素に格納される。
つまり、固有値を得るためには 1. `LWORK` の決定 → 2. 固有値の計算というように2回 `zheevx` を呼び出す。
それぞれ最適サイズより大きくても問題ない。
最適サイズより小さい場合、計算速度が悪化することが予想される。

何度も固有値問題を解く場合、`WORK`, `RWORK`, `IWORK` はサイズが十分な限り使い回せる。
すなわち、ブロック対角化のように何回か対角化をする場合、最も大きい行列で `WORK`, `RWORK`, `IWORK` を用意しておけば何度も確保し直す必要はない。

計算終了後に `INFO` に計算の終了状態に関する情報が返される。
`INFO = 0` なら正常終了。
`INFO < 0` なら `INFO` 番目の引数に問題があり、異常終了。
`INFO > 0` なら計算が収束しなかったことを意味する。

- ソースコード `test_LAPACK_zheevx.cpp`

```c++
#include <fstream>
#include <sstream>
#include <string>
#include <complex.h>

// 読み取る行列のサイズ
#define SIZE 1000

// im を虚数単位とする
#define im _Complex_I
#undef I

extern "C" {
  double dlamch_(const char& CMACH);

  void zheevx_(const char& JOBZ, const char& RANGE, const char& UPLO, const int& N,
        double _Complex* A, const int& LDA, const double& VL, const double& VU,
        const int& IL, const int& IU, const double& ABSTOL, int& M, double* W,
        double _Complex* Z, const int& LDZ, double _Complex* WORK,
        const int& LWORK, double* RWORK, int* IWORK, int* IFAIL, int& INFO);
};


// ２つのファイルを読み込んで Hermitian matrix (a) を作る
void readtsv(double _Complex* a){
    const auto real_file = "real.tsv";
    const auto imag_file = "imag.tsv";

    std::ifstream ifs_real(real_file);
    std::string line;
    int i = 0;
    while (getline(ifs_real, line)) {
        std::istringstream stream(line);
        std::string realstr;
        while (getline(stream, realstr, '\t')) {
            a[i] = std::stod(realstr);
            i++;
        }
    }

    std::ifstream ifs_imag(imag_file);
    i = 0;
    while (getline(ifs_imag, line)) {
        std::istringstream stream(line);
        std::string imagstr;
        while (getline(stream, imagstr, '\t')) {
            a[i] += std::stof(imagstr) * im;
            i++;
        }
    }
}


int main(){
    char jobz = 'V';
    char range = 'A';
    char uplo = 'U';
    int n = SIZE;
    int lda = SIZE;
    int vl = -1;
    int vu = -1;
    int il = -1;
    int iu = -1;
    double abstol = dlamch_('S');
    int m = SIZE;
    double w[SIZE];
    int ldz = SIZE;
    double _Complex work_temp[1];
    int lwork = -1;
    double rwork[7*SIZE];
    int iwork[5*SIZE];
    int ifail[2*SIZE];
    int info;

    // 行列はサイズが大きくなりうるので new を使ってヒープに確保する
    auto a = new double _Complex[SIZE*SIZE];
    auto z = new double _Complex[SIZE*SIZE];
    readtsv(a);

    // 最適な LWORK の計算
    zheevx_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
            work_temp, lwork, rwork, iwork, ifail, info);
    if (info != 0){return info;}
    lwork = (int) lround(creal(work_temp[0]));
    auto work = new double _Complex[lwork];

    // 対角化処理
    zheevx_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
            work, lwork, rwork, iwork, ifail, info);

    printf("Exit with info = %d\n", info);
    printf("w[0] = %f\n", w[0]);
    printf("w[1] = %f\n", w[1]);
    printf("...\n");
    printf("w[%d] = %f\n", SIZE-2, w[SIZE-2]);
    printf("w[%d] = %f\n", SIZE-1, w[SIZE-1]);

    delete[] a;
    delete[] z;
    delete[] work;
    return info;
}
```

- コンパイルコマンド

```sh
g++ test_LAPACK_zheevx.cpp -o a.out -llapack -lopenblas
```


### zheevdドライバー

ドキュメント: https://netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_ga9b3e110476166e66f2f62fa1fba6344a.html

引数 `A` に対角化する行列を与えると `W` に与えた配列に計算された固有値が格納される。
配列 `A` は対角化計算の過程で破壊的変更を受ける。
つまり計算が終了した時点で `A` は与えた配列と同一の内容ではないので注意。
`JOBZ` に `'V'` を指定した場合は `A` の第*i*列が `W` の*i*番目の要素に対応した固有ベクトルになっている。

`WORK`, `RWORK`, `IWORK` は計算のために必要な作業用のメモリ領域であり、`LWORK` は `WORK` の、`LRWORK` は `RWORK` の、`LIWORK` は `IWORK` の配列サイズである。
`WORK`, `RWORK`, `IWORK` の最適なサイズは `LWORK`, `LRWORK`, `LIWORK` のどれかに -1 を与えて `zheev` を呼び出すと計算でき、`WORK`, `RWORK`, `IWORK` それぞれの最初の要素に格納される。
つまり、固有値を得るためには 1. `LWORK`, `LRWORK`, `LIWORK` の決定 → 2. 固有値の計算というように2回 `zheevd` を呼び出す。
それぞれ最適サイズより大きくても問題ない。
最適サイズより小さい場合、計算速度が悪化することが予想される。

何度も固有値問題を解く場合、`WORK`, `RWORK`, `IWORK` はサイズが十分な限り使い回せる。
すなわち、ブロック対角化のように何回か対角化をする場合、最も大きい行列で `WORK`, `RWORK`, `IWORK` を用意しておけば何度も確保し直す必要はない。

計算終了後に `INFO` に計算の終了状態に関する情報が返される。
`INFO = 0` なら正常終了。
`INFO < 0` なら `INFO` 番目の引数に問題があり、異常終了。
`INFO > 0` なら計算が収束しなかったことを意味する。

- ソースコード `test_LAPACK_zheevd.cpp`

```c++
#include <fstream>
#include <sstream>
#include <string>
#include <complex.h>

// 読み取る行列のサイズ
#define SIZE 1000

// im を虚数単位とする
#define im _Complex_I
#undef I

extern "C" {
  void zheevd_(const char& JOBZ, const char& UPLO, const int& N,
        double _Complex* A, const int& LDA, double* W, double _Complex* WORK,
        const int& LWORK, double* RWORK, const int& LRWORK, int* IWORK,
        const int& LIWORK, int& INFO);
};


// ２つのファイルを読み込んで Hermitian matrix (a) を作る
void readtsv(double _Complex* a){
    const auto real_file = "real.tsv";
    const auto imag_file = "imag.tsv";

    std::ifstream ifs_real(real_file);
    std::string line;
    int i = 0;
    while (getline(ifs_real, line)) {
        std::istringstream stream(line);
        std::string realstr;
        while (getline(stream, realstr, '\t')) {
            a[i] = std::stod(realstr);
            i++;
        }
    }

    std::ifstream ifs_imag(imag_file);
    i = 0;
    while (getline(ifs_imag, line)) {
        std::istringstream stream(line);
        std::string imagstr;
        while (getline(stream, imagstr, '\t')) {
            a[i] += std::stof(imagstr) * im;
            i++;
        }
    }
}


int main(){
    char jobz = 'V';
    char uplo = 'U';
    int n = SIZE;
    int lda = SIZE;
    double w[SIZE];
    int info;
    // 行列はサイズが大きくなりうるので new を使ってヒープに確保する
    auto a = new double _Complex[SIZE*SIZE];
    readtsv(a);

    // 最適な LWORK の計算
    double _Complex work_temp[1];
    double rwork_temp[1];
    int iwork_temp[1];
    zheevd_(jobz, uplo, n, a, lda, w, work_temp, -1, rwork_temp, -1, iwork_temp, -1, info);
    int lwork = (int) lround(creal(work_temp[0]));
    int lrwork = (int) lround(creal(rwork_temp[0]));
    int liwork = (int) lround(creal(iwork_temp[0]));
    auto work = new double _Complex[lwork];
    auto rwork = new double[lrwork];
    auto iwork = new int[liwork];

    // 対角化処理
    zheevd_(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);

    printf("Exit with info = %d\n", info);
    printf("w[0] = %f\n", w[0]);
    printf("w[1] = %f\n", w[1]);
    printf("...\n");
    printf("w[%d] = %f\n", SIZE-2, w[SIZE-2]);
    printf("w[%d] = %f\n", SIZE-1, w[SIZE-1]);

    delete[] a;
    delete[] work;
    delete[] rwork;
    delete[] iwork;
    return 0;
}

```

- コンパイルコマンド

```sh
g++ test_LAPACK_zheevd.cpp -o a.out -llapack -lopenblas
```


### zheevrドライバー

ドキュメント: https://netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_ga60dd605c63d7183a4c289a4ab3df6df6.html

引数 `A` に対角化する行列を与えると `W` に与えた配列に計算された固有値が格納される。
配列 `A` は対角化計算の過程で破壊的変更を受ける。
つまり計算が終了した時点で `A` は与えた配列と同一の内容ではないので注意。
`JOBZ` に `'V'` を指定した場合は `Z` の第*i*列が `W` の*i*番目の要素に対応した固有ベクトルになっている。

`zheevx` ドライバーは与えられた行列 `A` の固有値の内の一部のみを計算する機能を持つ。
すべての固有値が必要ではない場合には、計算を打ち切って高速化することができる。
実際に見つかった固有値の個数は `M` に格納される。

- `RANGE` に `'A'` を指定した場合、すべての固有値を計算する。
- `RANGE` に `'V'` を指定した場合、(`VL`, `VU`] の半開区間に入る固有値のみを求める。
- `RANGE` に `'I'` を指定した場合 `IL` 番目から `IU` 番目の固有値のみを求める。

`ABSTOL` は固有値の許容計算誤差に関する引数である。特に理由がなければ DLMCH('S') を与えるのがよい。
v3.10.1時点ではまだそのような実装になっていないが、将来のバージョンではこの設定で最良の精度を得るようになる。

`WORK`, `RWORK`, `IWORK` は計算のために必要な作業用のメモリ領域であり、`LWORK` は `WORK` の、`LRWORK` は `RWORK` の、`LIWORK` は `IWORK` の配列サイズである。
`WORK`, `RWORK`, `IWORK` の最適なサイズは `LWORK`, `LRWORK`, `LIWORK` のどれかに -1 を与えて `zheev` を呼び出すと計算でき、`WORK`, `RWORK`, `IWORK` それぞれの最初の要素に格納される。
つまり、固有値を得るためには 1. `LWORK`, `LRWORK`, `LIWORK` の決定 → 2. 固有値の計算というように2回 `zheevd` を呼び出す。
それぞれ最適サイズより大きくても問題ない。
最適サイズより小さい場合、計算速度が悪化することが予想される。

何度も固有値問題を解く場合、`WORK`, `RWORK`, `IWORK` はサイズが十分な限り使い回せる。
すなわち、ブロック対角化のように何回か対角化をする場合、最も大きい行列で `WORK`, `RWORK`, `IWORK` を用意しておけば何度も確保し直す必要はない。

計算終了後に `INFO` に計算の終了状態に関する情報が返される。
`INFO = 0` なら正常終了。
`INFO < 0` なら `INFO` 番目の引数に問題があり、異常終了。
`INFO > 0` なら計算が収束しなかったことを意味する。

- ソースコード `test_LAPACK_zheevr.cpp`

```c++
#include <fstream>
#include <sstream>
#include <string>
#include <complex.h>

// 読み取る行列のサイズ
#define SIZE 1000

// im を虚数単位とする
#define im _Complex_I
#undef I

extern "C" {
  double dlamch_(const char& CMACH);

  void zheevr_(const char& JOBZ, const char& RANGE, const char& UPLO, const int& N,
        double _Complex* A, const int& LDA, const double& VL, const double& VU,
        const int& IL, const int& IU, const double& ABSTOL, int& M, double* W,
        double _Complex* Z, const int& LDZ, int* ISUPPZ, double _Complex* WORK,
        const int& LWORK, double* RWORK, const int& LRWORK, int* IWORK,
        const int& LIWORK, int& INFO);
};


// ２つのファイルを読み込んで Hermitian matrix (a) を作る
void readtsv(double _Complex* a){
    const auto real_file = "real.tsv";
    const auto imag_file = "imag.tsv";

    std::ifstream ifs_real(real_file);
    std::string line;
    int i = 0;
    while (getline(ifs_real, line)) {
        std::istringstream stream(line);
        std::string realstr;
        while (getline(stream, realstr, '\t')) {
            a[i] = std::stod(realstr);
            i++;
        }
    }

    std::ifstream ifs_imag(imag_file);
    i = 0;
    while (getline(ifs_imag, line)) {
        std::istringstream stream(line);
        std::string imagstr;
        while (getline(stream, imagstr, '\t')) {
            a[i] += std::stof(imagstr) * im;
            i++;
        }
    }
}


int main(){
    char jobz = 'V';
    char range = 'A';
    char uplo = 'U';
    int n = SIZE;
    int lda = SIZE;
    int vl = -1;
    int vu = -1;
    int il = -1;
    int iu = -1;
    double abstol = dlamch_('S');
    int m = SIZE;
    double w[SIZE];
    int ldz = SIZE;
    int isuppz[2*SIZE];
    int info;

    // 行列はサイズが大きくなりうるので new を使ってヒープに確保する
    auto a = new double _Complex[SIZE*SIZE];
    auto z = new double _Complex[SIZE*SIZE];
    readtsv(a);

    // 最適な LWORK の計算
    double _Complex work_temp[1];
    double rwork_temp[1];
    int iwork_temp[1];
    zheevr_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz,
            work_temp, -1, rwork_temp, -1, iwork_temp, -1, info);
    if (info != 0){return info;}
    int lwork = (int) lround(creal(work_temp[0]));
    int lrwork = (int) lround(creal(rwork_temp[0]));
    int liwork = (int) lround(creal(iwork_temp[0]));
    auto work = new double _Complex[lwork];
    auto rwork = new double[lrwork];
    auto iwork = new int[liwork];

    // 対角化処理
    zheevr_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
            isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);

    printf("Exit with info = %d\n", info);
    printf("w[0] = %f\n", w[0]);
    printf("w[1] = %f\n", w[1]);
    printf("...\n");
    printf("w[%d] = %f\n", SIZE-2, w[SIZE-2]);
    printf("w[%d] = %f\n", SIZE-1, w[SIZE-1]);

    delete[] a;
    delete[] z;
    delete[] work;
    delete[] rwork;
    delete[] iwork;
    return info;
}

```

- コンパイルコマンド

```sh
g++ test_LAPACK_zheevr.cpp -o a.out -llapack -lopenblas
```


## Driverのパフォーマンス比較

## BLASのパフォーマンス比較
