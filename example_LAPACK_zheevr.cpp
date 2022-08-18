#include <complex.h>
#include <iostream>
#include <stdio.h>
#include <string>

// 読み取る行列のサイズ
#define SIZE 1000

int readtsv(std::string, double _Complex*, int);

extern "C" {
  double dlamch_(const char& CMACH);

  void zheevr_(const char& JOBZ, const char& RANGE, const char& UPLO, const int& N,
        double _Complex* A, const int& LDA, const double& VL, const double& VU,
        const int& IL, const int& IU, const double& ABSTOL, int& M, double* W,
        double _Complex* Z, const int& LDZ, int* ISUPPZ, double _Complex* WORK,
        const int& LWORK, double* RWORK, const int& LRWORK, int* IWORK,
        const int& LIWORK, int& INFO);
};


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
    double abstol = dlamch_('S');
    int m = SIZE;
    double w[SIZE];
    int ldz = SIZE;
    int isuppz[2*SIZE];
    int info;

    // 行列はサイズが大きくなりうるので new を使ってヒープに確保する
    auto a = new double _Complex[SIZE*SIZE];
    auto z = new double _Complex[SIZE*SIZE];
    info = readtsv("matrix.tsv", a, SIZE*SIZE);
    if (info != 0){
        std::cerr << "Exit readtsv with info = " << info << std::endl;
        return info;
    }

    // 最適な LWORK の計算
    double _Complex work_temp[1];
    double rwork_temp[1];
    int iwork_temp[1];
    zheevr_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz,
            work_temp, -1, rwork_temp, -1, iwork_temp, -1, info);
    if (info != 0){return info;}
    int lwork = (int) lround(creal(work_temp[0]));
    int lrwork = (int) lround(creal(rwork_temp[0]));
    int liwork = (int) lround(creal(iwork_temp[0]));
    auto work = new double _Complex[lwork];
    auto rwork = new double[lrwork];
    auto iwork = new int[liwork];

    // 対角化処理
    zheevr_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
            isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);

    printf("Exit with info = %d\n", info);
    printf("w[0] = %f\n", w[0]);
    printf("w[1] = %f\n", w[1]);
    printf("...\n");
    printf("w[%d] = %f\n", SIZE-2, w[SIZE-2]);
    printf("w[%d] = %f\n", SIZE-1, w[SIZE-1]);

    delete[] a;
    delete[] z;
    delete[] work;
    delete[] rwork;
    delete[] iwork;
    return info;
}

