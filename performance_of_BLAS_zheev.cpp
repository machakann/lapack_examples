#include <chrono>
#include <complex.h>
#include <iostream>
#include <stdio.h>
#include <string>

// 試行回数
#define NTEST 10

// 読み取る行列のサイズ
#define SIZE 1000

int readtsv(std::string, double _Complex*, int);

extern "C" {
  void zheev_(const char& JOBZ, const char& UPLO, const int& N,
              double _Complex* A, const int& LDA, double* W,
              double _Complex* WORK, const int& LWORK, double* RWORK,
              int& INFO);
};


// 配列 a を配列 b の値で上書きする
void fill(double _Complex* a, double _Complex* b, int n){
    for (int i; i<n; i++){ a[i] = b[i]; }
}


int main(){
    char jobz = 'V';
    char uplo = 'U';
    int n = SIZE;
    int lda = SIZE;
    double w[SIZE];
    double rwork[3*SIZE - 2];
    int info;

    // 行列はサイズが大きくなりうるので new を使ってヒープに確保する
    auto a = new double _Complex[SIZE*SIZE];
    auto a_orig = new double _Complex[SIZE*SIZE];
    info = readtsv("matrix.tsv", a_orig, SIZE*SIZE);
    if (info != 0){
        std::cerr << "Exit readtsv with info = " << info << std::endl;
        return info;
    }

    // 最適な LWORK の計算
    double _Complex work_temp[1];
    zheev_(jobz, uplo, n, a, lda, w, work_temp, -1, rwork, info);
    int lwork = (int) lround(creal(work_temp[0]));
    auto work = new double _Complex[lwork];

    double elapsed;
    std::chrono::system_clock::time_point start, end;
    for (int i; i<NTEST; i++){
        fill(a, a_orig, SIZE*SIZE);
        start = std::chrono::system_clock::now();
        zheev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
        end = std::chrono::system_clock::now();
        elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        printf("%18.15lf\n", elapsed);
    }

    delete[] a;
    delete[] work;
    return 0;
}
