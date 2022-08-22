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
  void zheev_(const char& JOBZ, const char& UPLO, const int& N, complex* A,
              const int& LDA, double* W, complex* WORK, const int& LWORK,
              double* RWORK, int& INFO);
};


int main() {
    char jobz = 'V';
    char uplo = 'U';
    int n = SIZE;
    int lda = SIZE;
    double w[SIZE];
    double rwork[3*SIZE - 2];
    int info;

    vector<complex> a(SIZE*SIZE);
    int seed = 123456;
    std::mt19937 mt(seed);
    random_Hermitian(a, SIZE, mt);

    // 最適な LWORK の計算
    complex work_temp[1];
    zheev_(jobz, uplo, n, a.data(), lda, w, work_temp, -1, rwork, info);
    int lwork = (int) work_temp[0].real();
    vector<complex> work(lwork);

    // 対角化処理
    zheev_(jobz, uplo, n, a.data(), lda, w, work.data(), lwork, rwork, info);

    printf("Exit with info = %d\n", info);
    printf("w[0] = %f\n", w[0]);
    printf("w[1] = %f\n", w[1]);
    printf("...\n");
    printf("w[%d] = %f\n", SIZE-2, w[SIZE-2]);
    printf("w[%d] = %f\n", SIZE-1, w[SIZE-1]);
    return info;
}
