#include <chrono>
#include <complex>
#include <iostream>
#include <random>
#include <vector>
#include "random_Hermitian.h"

using complex = std::complex<double>;
using std::vector;

extern "C" {
  void zheev_(const char& JOBZ, const char& UPLO, const int& N,
              complex* A, const int& LDA, double* W, complex* WORK, 
              const int& LWORK, double* RWORK, int& INFO);
};


int main() {
    constexpr int NTEST = 10;
    constexpr int SIZE = 1000;
    char jobz = 'V';
    char uplo = 'U';
    int n = SIZE;
    vector<complex> a(SIZE*SIZE);
    int lda = SIZE;
    double w[SIZE];
    double rwork[3*SIZE - 2];
    int info;

    complex work_temp[1];
    zheev_(jobz, uplo, n, a.data(), lda, w, work_temp, -1, rwork, info);
    int lwork = (int) work_temp[0].real();
    vector<complex> work(lwork);

    int seed = 123456;
    std::mt19937 mt(seed);
    std::chrono::system_clock::time_point start, end;
    double elapsed;
    std::cout.precision(8);
    for (int i = 0; i<NTEST; i++) {
        random_Hermitian(a, SIZE, mt);
        start = std::chrono::system_clock::now();
        zheev_(jobz, uplo, n, a.data(), lda, w, work.data(), lwork, rwork, info);
        end = std::chrono::system_clock::now();
        elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        std::cout << elapsed << std::endl;
    }
    return 0;
}
