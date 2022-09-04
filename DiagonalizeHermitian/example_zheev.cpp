#include <complex>
#include <iostream>
#include <random>
#include <vector>
#include "random_Hermitian.h"

using complex = std::complex<double>;

extern "C" {
  void zheev_(const char& JOBZ, const char& UPLO, const int& N, complex* A,
              const int& LDA, double* W, complex* WORK, const int& LWORK,
              double* RWORK, int& INFO);
};


int main() {
    using std::vector;
    constexpr int SIZE = 1000;
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

    std::cout << "Exit with info = " << info << std::endl;
    std::cout << "w[0] = " << w[0] << std::endl;
    std::cout << "w[1] = " << w[1] << std::endl;
    std::cout << "..." << std::endl;
    std::cout << "w[" << SIZE-2 << "] = " << w[SIZE-2] << std::endl;
    std::cout << "w[" << SIZE-1 << "] = " << w[SIZE-1] << std::endl;
    return info;
}
