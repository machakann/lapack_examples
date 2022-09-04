#include <complex>
#include <iostream>
#include <random>
#include <vector>
#include "random_Hermitian.h"

using complex = std::complex<double>;
using std::vector;

extern "C" {
  double dlamch_(const char& CMACH);
  void zheevr_(const char& JOBZ, const char& RANGE, const char& UPLO, const int& N,
               complex* A, const int& LDA, const double& VL, const double& VU,
               const int& IL, const int& IU, const double& ABSTOL, int& M, double* W,
               complex* Z, const int& LDZ, int* ISUPPZ, complex* WORK,
               const int& LWORK, double* RWORK, const int& LRWORK, int* IWORK,
               const int& LIWORK, int& INFO);
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
    double abstol = dlamch_('S');
    int m = SIZE;
    double w[SIZE];
    vector<complex> z(SIZE*SIZE);
    int ldz = SIZE;
    int isuppz[2*SIZE];
    int info;

    vector<complex> a(SIZE*SIZE);
    int seed = 123456;
    std::mt19937 mt(seed);
    random_Hermitian(a, SIZE, mt);

    // 最適な LWORK の計算
    complex work_temp[1];
    double rwork_temp[1];
    int iwork_temp[1];
    zheevr_(jobz, range, uplo, n, a.data(), lda, vl, vu, il, iu, abstol, m, w,
            z.data(), ldz, isuppz, work_temp, -1, rwork_temp, -1, iwork_temp, -1, info);
    int lwork = (int) work_temp[0].real();
    int lrwork = (int) rwork_temp[0];
    int liwork = iwork_temp[0];
    vector<complex> work(lwork);
    vector<double> rwork(lrwork);
    vector<int> iwork(liwork);

    // 対角化処理
    zheevr_(jobz, range, uplo, n, a.data(), lda, vl, vu, il, iu, abstol, m, w,
            z.data(), ldz, isuppz, work.data(), lwork, rwork.data(), lrwork,
            iwork.data(), liwork, info);

    std::cout << "Exit with info = " << info << std::endl;
    std::cout << "w[0] = " << w[0] << std::endl;
    std::cout << "w[1] = " << w[1] << std::endl;
    std::cout << "..." << std::endl;
    std::cout << "w[" << SIZE-2 << "] = " << w[SIZE-2] << std::endl;
    std::cout << "w[" << SIZE-1 << "] = " << w[SIZE-1] << std::endl;
    return info;
}

