#include <chrono>
#include <complex>
#include <iostream>
#include <random>
#include <vector>
#include "random_Hermitian.h"

using complex = std::complex<double>;

extern "C" {
  double dlamch_(const char& CMACH);
  void zheevr_(const char& JOBZ, const char& RANGE, const char& UPLO,
               const int& N, complex* A, const int& LDA, const double& VL,
               const double& VU, const int& IL, const int& IU,
               const double& ABSTOL, int& M, double* W, complex* Z,
               const int& LDZ, int* ISUPPZ, complex* WORK, const int& LWORK,
               double* RWORK, const int& LRWORK, int* IWORK, const int& LIWORK,
               int& INFO);
};


int main() {
    using std::vector;
    constexpr int NTEST = 50;
    constexpr int SIZE = 1000;
    char jobz = 'V';
    char range = 'A';
    char uplo = 'U';
    int n = SIZE;
    vector<complex> a(SIZE*SIZE);
    int lda = SIZE;
    int vl = -1;
    int vu = -1;
    int il = -1;
    int iu = -1;
    double abstol = dlamch_('S');
    int m = SIZE;
    double w[SIZE];
    vector<complex> z(SIZE*SIZE);
    int ldz = SIZE;
    int isuppz[2*SIZE];
    int info;

    complex work_temp[1];
    double rwork_temp[1];
    int iwork_temp[1];
    zheevr_(jobz, range, uplo, n, a.data(), lda, vl, vu, il, iu, abstol, m,
            w, z.data(), ldz, isuppz, work_temp, -1, rwork_temp,
            -1, iwork_temp, -1, info);
    int lwork = (int) work_temp[0].real();
    int lrwork = (int) rwork_temp[0];
    int liwork = iwork_temp[0];
    vector<complex> work(lwork);
    vector<double> rwork(lrwork);
    vector<int> iwork(liwork);

    int seed = 123456;
    std::mt19937 mt(seed);
    std::chrono::system_clock::time_point start, end;
    double elapsed;
    std::cout.precision(8);
    for (int i = 0; i<NTEST; i++) {
        random_Hermitian(a, SIZE, mt);
        start = std::chrono::system_clock::now();
        zheevr_(jobz, range, uplo, n, a.data(), lda, vl, vu, il, iu, abstol, m,
                w, z.data(), ldz, isuppz, work.data(), lwork, rwork.data(),
                lrwork, iwork.data(), liwork, info);
        end = std::chrono::system_clock::now();
        elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        std::cout << elapsed << std::endl;
    }
    return 0;
}
