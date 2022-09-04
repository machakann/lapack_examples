#include <complex>
#include <iostream>
#include <random>
#include <vector>
#include "random_Hermitian.h"

using complex = std::complex<double>;
using std::vector;

extern "C" {
  void zheevd_(const char& JOBZ, const char& UPLO, const int& N,
               complex* A, const int& LDA, double* W, complex* WORK,
               const int& LWORK, double* RWORK, const int& LRWORK, int* IWORK,
               const int& LIWORK, int& INFO);
};


int main() {
    constexpr int SIZE = 1000;
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

    std::cout << "Exit with info = " << info << std::endl;
    std::cout << "w[0] = " << w[0] << std::endl;
    std::cout << "w[1] = " << w[1] << std::endl;
    std::cout << "..." << std::endl;
    std::cout << "w[" << SIZE-2 << "] = " << w[SIZE-2] << std::endl;
    std::cout << "w[" << SIZE-1 << "] = " << w[SIZE-1] << std::endl;
    return info;
}

