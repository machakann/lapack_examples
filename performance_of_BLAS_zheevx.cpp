#include <chrono>
#include <complex.h>
#include <iostream>
#include <stdio.h>
#include <string>

// 試行回数
#define NTEST 10

// 読み取る行列のサイズ
#define SIZE 1000

int readtsv(std::string, double _Complex*, int);

extern "C" {
  double dlamch_(const char& CMACH);

  void zheevx_(const char& JOBZ, const char& RANGE, const char& UPLO, const int& N,
        double _Complex* A, const int& LDA, const double& VL, const double& VU,
        const int& IL, const int& IU, const double& ABSTOL, int& M, double* W,
        double _Complex* Z, const int& LDZ, double _Complex* WORK,
        const int& LWORK, double* RWORK, int* IWORK, int* IFAIL, int& INFO);
};


// 配列 a を配列 b の値で上書きする
void fill(double _Complex* a, double _Complex* b, int n){
    for (int i; i<n; i++){
        a[i] = b[i];
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
    double abstol = 2*dlamch_('S');
    int m = SIZE;
    double w[SIZE];
    int ldz = SIZE;
    double rwork[7*SIZE];
    int iwork[5*SIZE];
    int ifail[2*SIZE];
    int info;

    // 行列はサイズが大きくなりうるので new を使ってヒープに確保する
    auto a = new double _Complex[SIZE*SIZE];
    auto a_orig = new double _Complex[SIZE*SIZE];
    auto z = new double _Complex[SIZE*SIZE];
    info = readtsv("matrix.tsv", a_orig, SIZE*SIZE);
    if (info != 0){
        std::cerr << "Exit readtsv with info = " << info << std::endl;
        return info;
    }

    // 最適な LWORK の計算
    double _Complex work_temp[1];
    zheevx_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
            work_temp, -1, rwork, iwork, ifail, info);
    if (info != 0){return info;}
    int lwork = (int) lround(creal(work_temp[0]));
    auto work = new double _Complex[lwork];

    double elapsed;
    std::chrono::system_clock::time_point start, end;
    for (int i; i<NTEST; i++){
        fill(a, a_orig, SIZE*SIZE);
        start = std::chrono::system_clock::now();
        zheevx_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
                work, lwork, rwork, iwork, ifail, info);
        end = std::chrono::system_clock::now();
        elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        printf("%18.15lf\n", elapsed);
    }

    delete[] a;
    delete[] work;
    return 0;
}
