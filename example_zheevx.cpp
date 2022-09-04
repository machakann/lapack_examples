#include <complex>
#include <iostream>
#include <random>
#include <vector>
#include "random_Hermitian.h"

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
    constexpr int SIZE = 1000;
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

    std::cout << "Exit with info = " << info << std::endl;
    std::cout << "w[0] = " << w[0] << std::endl;
    std::cout << "w[1] = " << w[1] << std::endl;
    std::cout << "..." << std::endl;
    std::cout << "w[" << SIZE-2 << "] = " << w[SIZE-2] << std::endl;
    std::cout << "w[" << SIZE-1 << "] = " << w[SIZE-1] << std::endl;
    return info;
}
