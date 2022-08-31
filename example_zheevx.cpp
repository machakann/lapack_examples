#include <complex>
#include <random>
#include <stdio.h>
#include <vector>
#include "random_Hermitian.h"

// 読み取る行列のサイズ
#define SIZE 1000

using complex = std::complex<double>;
using std::vector;

extern "C" {
  double dlamch_(const char& CMACH);
  void zheevx_(const char& JOBZ, const char& RANGE, const char& UPLO,
               const int& N, complex* A, const int& LDA, const double& VL,
               const double& VU, const int& IL, const int& IU,
               const double& ABSTOL, int& M, double* W, complex* Z,
               const int& LDZ, complex* WORK, const int& LWORK, double* RWORK,
               int* IWORK, int* IFAIL, int& INFO);
};


int main() {
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
    vector<complex> z(SIZE*SIZE);
    int ldz = SIZE;
    int lwork = -1;
    double rwork[7*SIZE];
    int iwork[5*SIZE];
    int ifail[2*SIZE];
    int info;

    vector<complex> a(SIZE*SIZE);
    int seed = 123456;
    std::mt19937 mt(seed);
    random_Hermitian(a, SIZE, mt);

    // 最適な LWORK の計算
    complex work_temp[1];
    zheevx_(jobz, range, uplo, n, a.data(), lda, vl, vu, il, iu, abstol, m, w,
            z.data(), ldz, work_temp, lwork, rwork, iwork, ifail, info);
    lwork = (int) work_temp[0].real();
    vector<complex> work(lwork);

    // 対角化処理
    zheevx_(jobz, range, uplo, n, a.data(), lda, vl, vu, il, iu, abstol, m, w,
            z.data(), ldz, work.data(), lwork, rwork, iwork, ifail, info);

    printf("Exit with info = %d\n", info);
    printf("w[0] = %f\n", w[0]);
    printf("w[1] = %f\n", w[1]);
    printf("...\n");
    printf("w[%d] = %f\n", SIZE-2, w[SIZE-2]);
    printf("w[%d] = %f\n", SIZE-1, w[SIZE-1]);
    return info;
}
