#pragma once
// Minimal self-contained FFT (no external deps).
// Forward : X[k] = sum_m x[m] e^{-2πi m k / n}
// Inverse : x[m] = (1/n) sum_k X[k] e^{+2πi m k / n}
// Uses an iterative radix-2 Cooley-Tukey when n is a power of two,
// otherwise Bluestein's algorithm (which reduces the arbitrary-n DFT to
// a power-of-two convolution). Matches numpy.fft results in double precision.
#include <complex>
#include <vector>
#include <cstddef>

namespace simplefft {

inline void fftPow2Raw(std::vector<std::complex<double>> &a, int sign)
{
    const int n = static_cast<int>(a.size());
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (int len = 2; len <= n; len <<= 1) {
        const double ang = sign * 2.0 * M_PI / len;
        const std::complex<double> wlen{std::cos(ang), std::sin(ang)};
        for (int i = 0; i < n; i += len) {
            std::complex<double> w(1.0);
            const int half = len >> 1;
            for (int k = 0; k < half; ++k) {
                const std::complex<double> u = a[i + k];
                const std::complex<double> v = a[i + k + half] * w;
                a[i + k] = u + v;
                a[i + k + half] = u - v;
                w *= wlen;
            }
        }
    }
}

// In-place DFT of arbitrary length. `inv` selects the direction;
// the inverse transform is scaled by 1/n.
inline void fft(std::vector<std::complex<double>> &z, bool inv)
{
    const int n = static_cast<int>(z.size());
    if (n <= 1) return;
    const int sign = inv ? +1 : -1;

    if ((n & (n - 1)) == 0) {
        fftPow2Raw(z, sign);
        if (inv)
            for (auto &c : z) c /= n;
        return;
    }

    // Bluestein: X[y] = w[y] * (a * b)[y] with
    //   w[k] = e^{sign * π i k² / n}, a[k] = z[k] w[k], b[k] = conj(w[k]).
    // Valid because k·y = (k² + y² − (k−y)²)/2.
    int m = 1;
    while (m < 2 * n - 1) m <<= 1;

    std::vector<std::complex<double>> A(m, 0.0), B(m, 0.0);
    auto chirp = [&](int k) -> std::complex<double> {
        const long long k2 = (static_cast<long long>(k) * k) % (2 * n);
        const double ang = sign * M_PI * k2 / n;
        return {std::cos(ang), std::sin(ang)};
    };
    for (int k = 0; k < n; ++k) {
        const auto w = chirp(k);
        A[k] = z[k] * w;
        B[k] = std::conj(w);
    }
    // The kernel must also appear at the tail of the circle so that the
    // circular convolution of length m is the linear convolution of two
    // length-n sequences: b[j] for j = -(n-1) .. n-1 (b has period 2n
    // and is even, so the tail is a mirror copy of the head).
    for (int k = 1; k < n; ++k)
        B[m - k] = B[k];
    fftPow2Raw(A, -1);
    fftPow2Raw(B, -1);
    for (int k = 0; k < m; ++k) A[k] *= B[k];
    fftPow2Raw(A, +1);
    for (auto &c : A) c /= m;              // unscaled core -> exact linear convolution
    for (int k = 0; k < n; ++k) z[k] = A[k] * chirp(k);
    if (inv)
        for (auto &c : z) c /= n;
}

} // namespace simplefft
