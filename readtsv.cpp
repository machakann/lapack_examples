#include <complex.h>
#include <fstream>
#include <stdlib.h>
#include <string>

// im を虚数単位とする
#define im _Complex_I
#undef I

enum ErrNo {
    Err_Noerror,
    Err_FileNotOpened,
    Err_ParseError,
};


// ファイルを読み込んで配列(a)を埋める
//   std::string        file    読み込むファイルパス
//   double _Complex*   a       読み込んだ数値を入れる配列
//   int                la      a の要素数
int readtsv(std::string file, double _Complex* a, int la) {
    char c;
    char last_c = '\0';
    int i = 0;
    int j = 0;
    char buf[25];
    bool realflg = true;
    double realp = 0.0;
    double imagp = 0.0;
    std::ifstream ifstr(file);
    if (!ifstr) { return Err_FileNotOpened; }

    while (ifstr.get(c)) {
        if (c == '0' || c == '1' || c == '2' || c == '3' || c == '4' ||
            c == '5' || c == '6' || c == '7' || c == '8' || c == '9' ||
            c == '.' || c == 'E' || c == 'e') {
            buf[i] = c;
            i++;
        } else if (c == '+' || c == '-') {
            if (last_c == 'E' || last_c == 'e') {
                buf[i] = c;
                i++;
            } else {
                if (realflg) {
                    buf[i] = '\0';
                    realp = atof(buf);
                    realflg = false;
                    i = 0;
                }
            }
        } else if (c == '\t' || c == '\n') {
            if (!realflg) {
                buf[i] = '\0';
                imagp = atof(buf);
                if (j >= la) {break;}
                a[j] = realp + imagp * im;
                j++;
                realflg = true;
                i = 0;
            } else {
                return Err_ParseError;
            }
        }
        last_c = c;
    }
    return 0;
}
