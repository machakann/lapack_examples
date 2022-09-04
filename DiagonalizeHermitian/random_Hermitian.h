#include <algorithm>
#include <complex>
#include <random>
#include <vector>


template <typename T>
void random_Hermitian(std::vector<std::complex<double>>& v, int n, T& rng) {
    using std::conj;
    using complex = std::complex<double>;
    std::uniform_real_distribution<> udist(-1, 1);
    v.resize(n*n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j < i) {
                v[n*i + j] = conj(v[n*j + i]);
            } else {
                if (i == j) {
                    v[n*i + j] = udist(rng);
                } else {
                    v[n*i + j] = complex(udist(rng), udist(rng));
                }
            }
        }
    }
}


template <typename T>
void random_band_Hermitian(std::vector<std::complex<double>>& v, int n, int w, T& rng) {
    using std::conj;
    using complex = std::complex<double>;
    using namespace std::literals::complex_literals;
    std::uniform_real_distribution<> udist(-1, 1);
    int lo, up;
    v.resize(n*n);
    for (int i = 0; i < n; i++) {
        lo = std::max(i - w, 0);
        up = std::min(i + w, n*n);
        for (int j = 0; j < n; j++) {
            if (lo <= j && j <= up) {
                if (j < i) {
                    v[n*i + j] = conj(v[n*j + i]);
                } else {
                    if (i == j) {
                        v[n*i + j] = udist(rng);
                    } else {
                        v[n*i + j] = complex(udist(rng), udist(rng));
                    }
                }
            } else {
                v[n*i + j] = 0.0i;
            }
        }
    }
}
