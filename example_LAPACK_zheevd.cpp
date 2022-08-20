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
  void zheevd_(const char& JOBZ, const char& UPLO, const int& N,
               complex* A, const int& LDA, double* W, complex* WORK,
               const int& LWORK, double* RWORK, const int& LRWORK, int* IWORK,
               const int& LIWORK, int& INFO);
};


int main(){
    char jobz = 'V';
    char uplo = 'U';
    int n = SIZE;
    int lda = SIZE;
    double w[SIZE];
    int info;

    vector<complex> a(SIZE*SIZE);
    int seed = 123456;
    std::mt19937 mt(seed);
    random_Hermitian(a, SIZE, mt);

    // 最適な LWORK の計算
    complex work_temp[1];
    double rwork_temp[1];
    int iwork_temp[1];
    zheevd_(jobz, uplo, n, a.data(), lda, w, work_temp, -1,
            rwork_temp, -1, iwork_temp, -1, info);
    int lwork = (int) work_temp[0].real();
    int lrwork = (int) rwork_temp[0];
    int liwork = iwork_temp[0];
    vector<complex> work(lwork);
    vector<double> rwork(lrwork);
    vector<int> iwork(liwork);

    // 対角化処理
    zheevd_(jobz, uplo, n, a.data(), lda, w, work.data(), lwork,
            rwork.data(), lrwork, iwork.data(), liwork, info);

    printf("Exit with info = %d\n", info);
    printf("w[0] = %f\n", w[0]);
    printf("w[1] = %f\n", w[1]);
    printf("...\n");
    printf("w[%d] = %f\n", SIZE-2, w[SIZE-2]);
    printf("w[%d] = %f\n", SIZE-1, w[SIZE-1]);
    return info;
}

