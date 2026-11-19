// IQClassify - DSP-Basisfunktionen (FFT via pocketfft, Statistik).

#include "iqclassify/Dsp.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "pocketfft_hdronly.h"

namespace iqclassify {
namespace detail {

std::size_t pow2ceil(std::size_t x) {
    std::size_t n = 1;
    while (n < x) n <<= 1;
    return n;
}

void fftInPlace(std::vector<cd>& a, bool invert) {
    const std::size_t n = a.size();
    if (n == 0) return;
    // pocketfft unterstuetzt beliebige Laengen; die Aufrufer nutzen
    // Zweierpotenzen (Leistungstransformation der Praefixplaene).
    pocketfft::shape_t shape{n};
    pocketfft::stride_t stride(sizeof(cd));
    pocketfft::shape_t axes{0};
    const double fct = invert ? 1.0 / static_cast<double>(n) : 1.0;
    pocketfft::c2c<double>(shape, stride, stride, axes, !invert,
                          reinterpret_cast<const std::complex<double>*>(a.data()),
                          reinterpret_cast<std::complex<double>*>(a.data()), fct,
                          1);
}

std::vector<cd> fftPadded(const cd* x, std::size_t len, std::size_t n,
                          bool invert) {
    if (n < len) n = len;
    std::vector<cd> buf(x, x + len);
    buf.resize(n, cd(0.0, 0.0));
    fftInPlace(buf, invert);
    return buf;
}

double db10(double v) {
    if (v < 1e-30) v = 1e-30;
    return 10.0 * std::log10(v);
}

double median(std::vector<double>& v) {
    // Haertung: nicht-endliche Werte aussortieren. std::sort mit NaN
    // verletzt die strikte schwache Ordnung und waere undefiniertes
    // Verhalten (bis hin zu Abstuerzen). Endliche Stichproben bleiben
    // unberuehrt.
    std::size_t w = 0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (std::isfinite(v[i])) {
            if (w != i) v[w] = v[i];
            ++w;
        }
    }
    v.resize(w);
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const std::size_t m = v.size();
    if (m % 2) return v[m / 2];
    return 0.5 * (v[m / 2 - 1] + v[m / 2]);
}

double medianOf(const std::vector<double>& v) {
    std::vector<double> w = v;
    return median(w);
}

std::size_t protoPad(std::size_t n) { return pow2ceil(2 * n); }

std::vector<double> unwrapPhase(const std::vector<cd>& y) {
    const std::size_t n = y.size();
    std::vector<double> up(n, 0.0);
    if (n == 0) return up;
    const double pi = 3.14159265358979323846;
    double acc = 0.0;
    double prev = std::atan2(y[0].imag(), y[0].real());
    up[0] = prev;
    for (std::size_t i = 1; i < n; ++i) {
        const double p = std::atan2(y[i].imag(), y[i].real());
        double d = p - prev;
        if (d > pi) d -= 2.0 * pi;
        else if (d < -pi) d += 2.0 * pi;
        acc += d;
        prev = p;
        up[i] = acc;
    }
    return up;
}

std::vector<double> rrcTaps(long sps, double alpha, long span) {
    const double pi = 3.14159265358979323846;
    const long ntaps = span * sps + 1;
    std::vector<double> taps(static_cast<std::size_t>(ntaps), 0.0);
    const double fsps = static_cast<double>(sps);
    for (long i = 0; i < ntaps; ++i) {
        const double t =
            (static_cast<double>(i) - (ntaps - 1) / 2.0) / fsps;
        double v;
        if (std::abs(t) < 1e-9) {
            v = 1.0 - alpha + 4.0 * alpha / pi;
        } else if (std::abs(std::abs(t) - 1.0 / (4.0 * alpha)) < 1e-9) {
            v = (alpha / std::sqrt(2.0)) *
                ((1.0 + 2.0 / pi) * std::sin(pi / (4.0 * alpha)) -
                 (1.0 - 2.0 / pi) * std::cos(pi / (4.0 * alpha)));
        } else {
            const double num = std::sin(pi * t * (1.0 - alpha)) +
                4.0 * alpha * t * std::cos(pi * t * (1.0 + alpha));
            const double den = pi * t * (1.0 - (4.0 * alpha * t) * (4.0 * alpha * t));
            v = num / den;
        }
        taps[static_cast<std::size_t>(i)] = v;
    }
    double s = 0.0;
    for (double v : taps) s += v;
    if (std::abs(s) > 1e-30) {
        for (double& v : taps) v /= s;
    }
    return taps;
}

}  // namespace detail
}  // namespace iqclassify
