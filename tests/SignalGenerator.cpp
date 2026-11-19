// IQClassify - Implementierung der Testsignale.

#include "SignalGenerator.hpp"

#include <algorithm>
#include <cmath>

#include "iqclassify/Dsp.hpp"

namespace iqtest {

namespace {

const double kPi = 3.14159265358979323846;

// Kausale Faltungsfilterung (wie die Referenz: y[i] = sum x[i-k] taps[k]).
std::vector<cd> fir(const std::vector<cd>& x,
                    const std::vector<double>& taps) {
    const std::size_t nt = taps.size();
    std::vector<cd> y(x.size(), cd(0.0, 0.0));
    for (std::size_t i = 0; i < x.size(); ++i) {
        cd acc(0.0, 0.0);
        const std::size_t hi = std::min(nt, i + 1);
        for (std::size_t k = 0; k < hi; ++k)
            acc += x[i - k] * taps[k];
        y[i] = acc;
    }
    return y;
}

double meanPower(const std::vector<cd>& v) {
    double p = 0.0;
    for (const cd& s : v) p += std::norm(s);
    return p / static_cast<double>(v.size());
}

}  // namespace

std::vector<cd> genPsk(std::mt19937& rng, int m, std::size_t nsym, int sps,
                       double alpha) {
    std::uniform_int_distribution<int> pick(0, m - 1);
    const std::size_t guard = static_cast<std::size_t>(8 * sps);
    std::vector<cd> up(nsym * static_cast<std::size_t>(sps) + guard,
                       cd(0.0, 0.0));
    for (std::size_t i = 0; i < nsym; ++i) {
        const double a = 2.0 * kPi * pick(rng) / m;
        up[i * static_cast<std::size_t>(sps) +
           static_cast<std::size_t>(4 * sps)] = cd(std::cos(a), std::sin(a));
    }
    const std::vector<double> taps =
        iqclassify::detail::rrcTaps(sps, alpha, 8);
    std::vector<cd> y = fir(up, taps);
    y.assign(y.begin() + static_cast<std::ptrdiff_t>(4 * sps),
             y.begin() + static_cast<std::ptrdiff_t>(4 * sps) +
                 static_cast<std::ptrdiff_t>(nsym * sps));
    const double p = meanPower(y);
    const double sc = 1.0 / std::sqrt(p > 1e-30 ? p : 1e-30);
    for (cd& v : y) v *= sc;
    return y;
}

std::vector<cd> genQam(std::mt19937& rng, int m, std::size_t nsym, int sps,
                       double alpha) {
    const int L = static_cast<int>(std::lround(std::sqrt(static_cast<double>(m))));
    std::vector<int> lv(static_cast<std::size_t>(L));
    for (int i = 0; i < L; ++i) lv[static_cast<std::size_t>(i)] = 2 * i - (L - 1);
    std::uniform_int_distribution<int> pick(0, L - 1);
    std::vector<cd> syms(nsym);
    for (std::size_t i = 0; i < nsym; ++i) {
        const double re = lv[static_cast<std::size_t>(pick(rng))];
        const double im = lv[static_cast<std::size_t>(pick(rng))];
        syms[i] = cd(re, im);
    }
    const double p = meanPower(syms);
    const double sc = 1.0 / std::sqrt(p > 1e-30 ? p : 1e-30);
    for (cd& s : syms) s *= sc;
    const std::size_t guard = static_cast<std::size_t>(8 * sps);
    std::vector<cd> up(nsym * static_cast<std::size_t>(sps) + guard,
                       cd(0.0, 0.0));
    for (std::size_t i = 0; i < nsym; ++i) {
        up[i * static_cast<std::size_t>(sps) +
           static_cast<std::size_t>(4 * sps)] = syms[i];
    }
    const std::vector<double> taps =
        iqclassify::detail::rrcTaps(sps, alpha, 8);
    std::vector<cd> y = fir(up, taps);
    y.assign(y.begin() + static_cast<std::ptrdiff_t>(4 * sps),
             y.begin() + static_cast<std::ptrdiff_t>(4 * sps) +
                 static_cast<std::ptrdiff_t>(nsym * sps));
    const double p2 = meanPower(y);
    const double sc2 = 1.0 / std::sqrt(p2 > 1e-30 ? p2 : 1e-30);
    for (cd& v : y) v *= sc2;
    return y;
}

std::vector<cd> genFsk2(std::mt19937& rng, std::size_t nsym, int sps,
                        double sep) {
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    double ph = u01(rng) * 2.0 * kPi;
    const double d = sep / 2.0;
    std::vector<cd> out;
    out.reserve(nsym * static_cast<std::size_t>(sps));
    for (std::size_t i = 0; i < nsym; ++i) {
        const double f = u01(rng) < 0.5 ? d : -d;
        for (int k = 0; k < sps; ++k) {
            out.push_back(cd(std::cos(ph), std::sin(ph)));
            ph += 2.0 * kPi * f;
        }
    }
    return out;
}

std::vector<cd> genFrank(std::mt19937& rng, int m, int sps, double alpha) {
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    const std::size_t nsym = static_cast<std::size_t>(m) * m;
    const std::size_t guard = static_cast<std::size_t>(8 * sps);
    std::vector<cd> up(nsym * static_cast<std::size_t>(sps) + guard,
                       cd(0.0, 0.0));
    const double ph0 = u01(rng) * 2.0 * kPi;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            const double a = ph0 + 2.0 * kPi * i * j / m;
            up[(static_cast<std::size_t>(i) * m + j) *
                   static_cast<std::size_t>(sps) +
               static_cast<std::size_t>(4 * sps)] =
                cd(std::cos(a), std::sin(a));
        }
    }
    const std::vector<double> taps =
        iqclassify::detail::rrcTaps(sps, alpha, 8);
    std::vector<cd> y = fir(up, taps);
    y.assign(y.begin() + static_cast<std::ptrdiff_t>(4 * sps),
             y.begin() + static_cast<std::ptrdiff_t>(4 * sps) +
                 static_cast<std::ptrdiff_t>(nsym * sps));
    const double p = meanPower(y);
    const double sc = 1.0 / std::sqrt(p > 1e-30 ? p : 1e-30);
    for (cd& v : y) v *= sc;
    return y;
}

std::vector<cd> genChirp(std::size_t nsamp, double bw, bool up) {
    const double rate =
        bw / static_cast<double>(nsamp) * (up ? 1.0 : -1.0);
    std::vector<cd> out(nsamp);
    for (std::size_t i = 0; i < nsamp; ++i) {
        const double t = static_cast<double>(i) -
                          static_cast<double>(nsamp) / 2.0;
        const double ph = kPi * rate * t * t;
        out[i] = cd(std::cos(ph), std::sin(ph));
    }
    const double p = meanPower(out);
    const double sc = 1.0 / std::sqrt(p > 1e-30 ? p : 1e-30);
    for (cd& v : out) v *= sc;
    return out;
}

std::vector<cd> genCw(std::mt19937& rng, std::size_t nsamp) {
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    const double ph = u01(rng) * 2.0 * kPi;
    std::vector<cd> out(nsamp, cd(std::cos(ph), std::sin(ph)));
    return out;
}

std::vector<cd> applyEnv(const std::vector<cd>& sig, int mode, double depth,
                         std::mt19937& rng) {
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    const std::size_t n = sig.size();
    std::vector<cd> out(n);
    for (std::size_t i = 0; i < n; ++i) {
        double e;
        if (mode == 0) {
            e = 1.0 - depth * u01(rng);
        } else if (mode == 1) {
            e = (1.0 - depth) +
                depth * 0.5 * (1.0 + std::sin(2.0 * kPi * i /
                                               static_cast<double>(n)));
        } else {
            const double x = static_cast<double>(i) / static_cast<double>(n);
            e = 1.0 - depth * std::exp(-std::pow((x - 0.5) * 8.0, 2.0));
        }
        out[i] = sig[i] * e;
    }
    const double p = meanPower(out);
    const double sc = 1.0 / std::sqrt(p > 1e-30 ? p : 1e-30);
    for (cd& v : out) v *= sc;
    return out;
}

std::vector<cd> genOfdm(std::mt19937& rng, int nfft, int cp, int nsym,
                        double load) {
    int nact = std::max(16, static_cast<int>(nfft * load));
    if (nact % 2) ++nact;
    const int half = nact / 2;
    std::uniform_int_distribution<int> pick(0, 3);
    std::vector<cd> acc;
    acc.reserve(static_cast<std::size_t>(nsym * (nfft + cp)));
    for (int s = 0; s < nsym; ++s) {
        std::vector<cd> x(static_cast<std::size_t>(nfft), cd(0.0, 0.0));
        for (int k = -half; k <= half; ++k) {
            if (k == 0) continue;
            const double a = kPi / 4.0 +
                             kPi / 2.0 * pick(rng);
            const int kk = ((k % nfft) + nfft) % nfft;
            x[static_cast<std::size_t>(kk)] =
                cd(std::cos(a), std::sin(a));
        }
        // inverse FFT: Zeitsymbol (PocketFFT, 1/n-Skalierung wie die FFT
        // der Referenz)
        iqclassify::detail::fftInPlace(x, true);
        const double p = meanPower(x);
        const double sc = 1.0 / std::sqrt(p > 1e-30 ? p : 1e-30);
        for (cd& v : x) v *= sc;
        acc.insert(acc.end(), x.begin(), x.end());
        acc.insert(acc.end(), x.begin(),
                   x.begin() + static_cast<std::ptrdiff_t>(cp));
    }
    return acc;
}

std::vector<cd> mix(const std::vector<cd>& sig, double fc) {
    std::vector<cd> out(sig.size());
    for (std::size_t i = 0; i < sig.size(); ++i) {
        const double ph = 2.0 * kPi * fc * static_cast<double>(i);
        out[i] = sig[i] * cd(std::cos(ph), std::sin(ph));
    }
    return out;
}

std::vector<cd> noised(std::mt19937& rng, const std::vector<cd>& sig,
                       double snrDb) {
    const double npow = std::pow(10.0, -snrDb / 10.0);
    const double sigP = meanPower(sig);
    const double sigma = std::sqrt(npow * sigP / 2.0);
    std::normal_distribution<double> gauss(0.0, sigma);
    std::vector<cd> out(sig.size());
    for (std::size_t i = 0; i < sig.size(); ++i) {
        out[i] = sig[i] + cd(gauss(rng), gauss(rng));
    }
    return out;
}

void Recorder::silence(std::size_t nsamp) {
    std::mt19937 rng(static_cast<unsigned>(1234 + silenceCount_));
    ++silenceCount_;
    std::normal_distribution<double> gauss(0.0, sigma_);
    for (std::size_t i = 0; i < nsamp; ++i) {
        buf_.push_back(cd(gauss(rng), gauss(rng)));
    }
}

void Recorder::add(const std::vector<cd>& sig, double fc, double snrDb,
                   std::size_t at) {
    std::size_t pos = at;
    if (pos == kAtEnd) {
        pos = buf_.size();
    }
    while (buf_.size() < pos + sig.size())
        buf_.push_back(cd(0.0, 0.0));
    std::mt19937 rng(static_cast<unsigned>(999 + pulseCount_));
    ++pulseCount_;
    std::vector<cd> y = noised(rng, mix(sig, fc), snrDb);
    for (std::size_t i = 0; i < y.size(); ++i) {
        buf_[pos + i] += y[i];
    }
}

}  // namespace iqtest
