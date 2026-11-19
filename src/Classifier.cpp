// IQClassify - Stufe 1 (grobe Klasse), Stufe 2 (Modulationsart) und
// Klassifikator-Fassade.
//
// Alle Schwellwerte und Entscheidungsequenzen entsprechen eins zu eins der
// gegen simulierten Pulsen (SNR 5..20 dB, sps 4/8/16) validierten
// Referenzimplementierung; siehe docs/METHODOLOGY.md.

#include "iqclassify/Classifier.hpp"
#include "iqclassify/PulseSegmenter.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <utility>
#include <vector>

namespace iqclassify {

const char* toString(ModulationType t) {
    switch (t) {
        case ModulationType::Unknown: return "Unknown";
        case ModulationType::Unmodulated: return "Unmodulated";
        case ModulationType::BPSK: return "BPSK";
        case ModulationType::QPSK: return "QPSK";
        case ModulationType::FSK2: return "FSK2";
        case ModulationType::Chirp: return "Chirp";
        case ModulationType::OFDM: return "OFDM";
        case ModulationType::QAM16: return "QAM16";
        case ModulationType::QAM64: return "QAM64";
        case ModulationType::QAM256: return "QAM256";
    }
    return "Unknown";
}

const char* toString(CoarseClass c) {
    switch (c) {
        case CoarseClass::Unknown: return "Unknown";
        case CoarseClass::Tone: return "Tone";
        case CoarseClass::FrequencyModulated: return "FrequencyModulated";
        case CoarseClass::PhaseModulated: return "PhaseModulated";
        case CoarseClass::NoiseLike: return "NoiseLike";
    }
    return "Unknown";
}

namespace {

const double kPi = 3.14159265358979323846;

using cd = detail::cd;

// Rundung fuer Ausgabewerte (Paritaet mit der Referenzausgabe).
double roundTo(double v, int decimals) {
    const double f = std::pow(10.0, decimals);
    return std::round(v * f) / f;
}

// Zyklische Normierung einer fraktionalen Frequenz nach [0, 1).
double wrap01(double f) {
    f = std::fmod(f, 1.0);
    if (f < 0.0) f += 1.0;
    return f;
}

// --------------------------------------------------------------------------
// zyklische Praefix-Evidenz (OFDM)
// --------------------------------------------------------------------------
struct CpEvidence {
    double rho = 0.0;
    long tau = 0;
    double contrast = 0.0;
    double singularity = 0.0;
    double clean = 0.0;
};

CpEvidence cpAutocorrPeak(const std::vector<cd>& y,
                          const ClassifierConfig& cfg) {
    const std::size_t n = y.size();
    CpEvidence e;
    if (n < 256) return e;
    const std::size_t nft = detail::protoPad(n);
    std::vector<cd> buf = y;
    buf.resize(nft, cd(0.0, 0.0));
    detail::fftInPlace(buf, false);
    for (auto& v : buf) {
        const double r = v.real();
        const double im = v.imag();
        v = cd(r * r + im * im, 0.0);
    }
    detail::fftInPlace(buf, true);
    const double r0 = buf[0].real();
    if (r0 <= 0.0) return e;
    auto rv = [&](long t) {
        const double denom =
            r0 * (1.0 - static_cast<double>(t) / static_cast<double>(n));
        return std::abs(buf[static_cast<std::size_t>(t)].real()) /
               (denom > 1e-30 ? denom : 1e-30);
    };
    const long tauMax =
        std::min<long>(static_cast<long>(n / 4), cfg.ofdmTauMax);
    const long tauMin = cfg.ofdmTauMin;
    double best = 0.0;
    double second = 0.0;
    long bt = 0;
    const long tauEnd = std::min<long>(tauMax, static_cast<long>(n) - 1);
    for (long tau = tauMin; tau < tauEnd; ++tau) {
        const double v = rv(tau);
        if (v > best) {
            if (std::abs(tau - bt) > cfg.cpIsolation) second = best;
            best = v;
            bt = tau;
        } else if (v > second) {
            if (std::abs(tau - bt) > cfg.cpIsolation) second = v;
        }
    }
    if (bt == 0) return e;
    const long lo = std::max(tauMin, bt - 64);
    const long hi = std::min(tauMax - 1, bt + 64);
    double s = 0.0;
    long c = 0;
    for (long t = lo; t <= hi; ++t) {
        if (std::abs(t - bt) > 3) {
            s += rv(t);
            ++c;
        }
    }
    const double mean = s / static_cast<double>(c > 0 ? c : 1);
    const double contrast = best / std::max(mean, 1e-9);
    const double sing =
        second > 0.0 ? best / std::max(second, 1e-9) : 99.0;
    const long W = cfg.cpCleanWin;
    long deep = 0;
    long tot = 0;
    const long clo = std::max(tauMin, bt - W);
    const long chi = std::min(tauMax - 1, bt + W);
    for (long t = clo; t <= chi; ++t) {
        if (std::abs(t - bt) <= 3) continue;
        ++tot;
        if (rv(t) < cfg.tCpDeepFrac * best) ++deep;
    }
    const double clean = static_cast<double>(deep) /
                         static_cast<double>(tot > 0 ? tot : 1);
    e.rho = best;
    e.tau = bt;
    e.contrast = contrast;
    e.singularity = sing;
    e.clean = clean;
    return e;
}

bool cpHalfConsistency(const std::vector<cd>& y, const ClassifierConfig& cfg,
                       double rho, long /*tau*/) {
    const std::size_t n = y.size();
    std::pair<long, double> res[2];
    for (int h = 0; h < 2; ++h) {
        const std::size_t lo = h == 0 ? 0 : n / 2;
        const std::size_t hi = h == 0 ? n / 2 : n;
        const std::size_t m = hi - lo;
        if (m < 256) return false;
        const std::size_t nft = detail::protoPad(m);
        std::vector<cd> buf(y.begin() + static_cast<std::ptrdiff_t>(lo),
                            y.begin() + static_cast<std::ptrdiff_t>(hi));
        buf.resize(nft, cd(0.0, 0.0));
        detail::fftInPlace(buf, false);
        for (auto& v : buf) {
            const double r = v.real();
            const double im = v.imag();
            v = cd(r * r + im * im, 0.0);
        }
        detail::fftInPlace(buf, true);
        const double r0 = buf[0].real();
        if (r0 <= 0.0) return false;
        const long tmax =
            std::min<long>(static_cast<long>(m / 4), cfg.ofdmTauMax);
        double best = 0.0;
        long bt = 0;
        for (long t = cfg.ofdmTauMin; t < tmax; ++t) {
            const double num =
                std::abs(buf[static_cast<std::size_t>(t)].real());
            if (num <= 0.0) continue;
            const double v = num / (r0 * (1.0 - static_cast<double>(t) /
                                                static_cast<double>(m)));
            if (v > best) {
                best = v;
                bt = t;
            }
        }
        res[h] = {bt, best};
    }
    if (res[0].first == 0 || res[1].first == 0) return false;
    if (std::abs(res[0].first - res[1].first) > cfg.cpHalfTol) return false;
    if (std::min(res[0].second, res[1].second) < cfg.cpHalfRhoRel * rho)
        return false;
    return true;
}

// --------------------------------------------------------------------------
// Symboltiming
// --------------------------------------------------------------------------
std::pair<double, double> rsFromCyclic(const std::vector<double>& env2,
                                       std::size_t n, double bw) {
    const std::size_t nft = detail::pow2ceil(2 * n);
    double mean = 0.0;
    for (double v : env2) mean += v;
    mean /= static_cast<double>(n);
    std::vector<cd> buf(nft, cd(0.0, 0.0));
    for (std::size_t i = 0; i < n; ++i)
        buf[i] = cd(env2[i] - mean, 0.0);
    detail::fftInPlace(buf, false);
    const long kLo = std::max(1L, static_cast<long>(
                                     std::floor(bw / 1.50 *
                                                static_cast<double>(nft))));
    const long kHi = std::min<long>(
        static_cast<long>(nft / 4),
        static_cast<long>(std::ceil(bw / 1.05 * static_cast<double>(nft))));
    if (kHi <= kLo + 2) return {0.0, 0.0};
    std::vector<double> mags;
    mags.reserve(static_cast<std::size_t>(kHi - kLo));
    long bestK = 0;
    double bestV = -1.0;
    for (long k = kLo; k < kHi; ++k) {
        const cd& f = buf[static_cast<std::size_t>(k)];
        const double v = f.real() * f.real() + f.imag() * f.imag();
        mags.push_back(v);
        if (v > bestV) {
            bestV = v;
            bestK = k;
        }
    }
    std::sort(mags.begin(), mags.end());
    const double med = !mags.empty() ? mags[mags.size() / 2] : 0.0;
    const double prom = bestV / std::max(med, 1e-30);
    if (bestK <= 0 || bestK >= static_cast<long>(nft) - 1)
        return {0.0, prom};
    // parabolische Peakinterpolation
    const double a0 = std::norm(buf[static_cast<std::size_t>(bestK - 1)]);
    const double a1 = std::norm(buf[static_cast<std::size_t>(bestK)]);
    const double a2 = std::norm(buf[static_cast<std::size_t>(bestK + 1)]);
    const double den = a0 - 2.0 * a1 + a2;
    double delta = 0.0;
    if (std::abs(den) > 1e-30) {
        delta = 0.5 * (a0 - a2) / den;
        if (std::abs(delta) > 1.0) delta = 0.0;
    }
    double fr = (static_cast<double>(bestK) + delta) /
                static_cast<double>(nft);
    // Sub-Bin-Verfeinerung: Maximum der zyklischen Kohaerenz
    // |sum env2[i] exp(-i 2 pi fr i)|. Der parabolische Fit des
    // geleakten Dirichlet-Nebenzipfels ist um bis zu einen Bin
    // vorgeneigt; ein Symbolratenfehler von 1e-4/Sample driftet ueber
    // einige tausend Symbole um ein volles Sample.
    auto coh = [&](double f) {
        cd acc(0.0, 0.0);
        double ph = 0.0;
        const double stp = 2.0 * kPi * f;
        for (const cd& v : buf) {
            ph += stp;
            acc += v * cd(std::cos(ph), -std::sin(ph));
        }
        return std::abs(acc);
    };
    double d = 3.0 / static_cast<double>(nft);
    while (d > 1e-7) {
        const double flo = fr - d;
        const double fhi = fr + d;
        fr = coh(flo) > coh(fhi) ? flo : fhi;
        d *= 0.5;
    }
    return {fr, prom};
}

double eyeOffset(const std::vector<double>& env, std::size_t n, double sps) {
    const double step = 0.125;
    std::vector<double> prof;
    double off = 0.0;
    while (off < sps) {
        double acc = 0.0;
        long c = 0;
        long k = 0;
        while (off + static_cast<double>(k) * sps <
               static_cast<double>(n) - 1.0) {
            const double p = off + static_cast<double>(k) * sps;
            const std::size_t i = static_cast<std::size_t>(p);
            const double frac = p - static_cast<double>(i);
            if (i + 1 >= n) break;
            const double v =
                env[i] * (1.0 - frac) + env[i + 1] * frac;
            acc += v * v;
            ++c;
            ++k;
        }
        prof.push_back(acc / static_cast<double>(c > 0 ? c : 1));
        off += step;
    }
    std::size_t imax = 0;
    for (std::size_t i = 1; i < prof.size(); ++i) {
        if (prof[i] > prof[imax]) imax = i;
    }
    double d = 0.0;
    if (imax > 0 && imax + 1 < prof.size()) {
        const double den = prof[imax - 1] - 2.0 * prof[imax] +
                           prof[imax + 1];
        if (den < -1e-30) {
            d = 0.5 * (prof[imax - 1] - prof[imax + 1]) / den;
        }
    }
    const double res = (static_cast<double>(imax) + d) * step;
    return res > 0.0 ? res : 0.0;
}

// Symbolpositionen (Sampleindizes) auf dem geschaetzten Gitter.
std::vector<std::size_t> estimateSymbolsRobust(const Features& f,
                                                const ClassifierConfig& cfg) {
    const std::size_t n = f.n;
    if (n < 256) return {};
    std::vector<double> env2(n);
    for (std::size_t i = 0; i < n; ++i) env2[i] = f.env[i] * f.env[i];
    const auto rp = rsFromCyclic(env2, n, f.bw);
    if (rp.first <= 0.0 || rp.second < cfg.tRsProm) return {};
    const double sps = 1.0 / rp.first;
    const double spsMax =
        std::max(2.0, static_cast<double>(n) / 12.0);
    if (sps < 1.6 || sps > spsMax) return {};
    const double off = eyeOffset(f.env, n, sps);
    std::vector<std::size_t> pos;
    long k = 0;
    while (off + static_cast<double>(k) * sps <
           static_cast<double>(n) - 2.0) {
        const double p = off + static_cast<double>(k) * sps;
        const long ip = std::lround(p);
        if (ip >= 1) pos.push_back(static_cast<std::size_t>(ip));
        ++k;
    }
    if (pos.size() < cfg.timingMinSyms) return {};
    return pos;
}

// --------------------------------------------------------------------------
// Konstellationsstatistik
// --------------------------------------------------------------------------
double foldConc(const std::vector<double>& ph, int m) {
    cd acc(0.0, 0.0);
    for (double p : ph) {
        const double a = static_cast<double>(m) * p;
        acc += cd(std::cos(a), std::sin(a));
    }
    return std::abs(acc) / static_cast<double>(ph.size());
}

double kappa(const std::vector<cd>& s) {
    const double n = static_cast<double>(s.size());
    double m2 = 0.0;
    double m4 = 0.0;
    cd m22(0.0, 0.0);
    for (const cd& v : s) {
        const double a = std::abs(v);
        m2 += a * a;
        m4 += a * a * a * a;
        m22 += v * v;
    }
    m2 /= n;
    m4 /= n;
    m22 /= n;
    const double c42 = m4 - std::norm(m22) - 2.0 * m2 * m2;
    return -c42 / (m2 * m2);
}

double snrCorrectK42(double k42Obs, double snrDb) {
    const double r = std::pow(10.0, snrDb / 10.0);
    return k42Obs * (1.0 + 1.0 / r) * (1.0 + 1.0 / r);
}

std::pair<double, double> carrierFromSymbols(const std::vector<cd>& s) {
    const std::size_t len = s.size();
    if (len < 16) return {0.0, 1.0};
    std::vector<cd> z4(len);
    for (std::size_t i = 0; i < len; ++i) {
        const cd s2 = s[i] * s[i];
        z4[i] = s2 * s2;
    }
    const std::size_t n = detail::pow2ceil(4 * len);
    std::vector<cd> buf = z4;
    buf.resize(n, cd(0.0, 0.0));
    detail::fftInPlace(buf, false);
    std::size_t k = 0;
    double bv = -1.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double v = std::abs(buf[i]);
        if (v > bv) {
            bv = v;
            k = i;
        }
    }
    const double a = std::abs(buf[(k + n - 1) % n]);
    const double b = std::abs(buf[k]);
    const double c = std::abs(buf[(k + 1) % n]);
    const double den = a - 2.0 * b + c;
    double kf = static_cast<double>(k) +
               (std::abs(den) > 1e-30 ? 0.5 * (a - c) / den : 0.0);
    if (kf > static_cast<double>(n) / 2.0) kf -= static_cast<double>(n);
    const double f0 = kf / (4.0 * static_cast<double>(n));
    // Kohaerenz-Verfeinerung |sum z4[k] exp(-i 8 pi f k)|
    auto coh = [&](double f) {
        cd acc(0.0, 0.0);
        for (std::size_t kk = 0; kk < z4.size(); ++kk) {
            const double ph = -8.0 * kPi * f * static_cast<double>(kk);
            acc += z4[kk] * cd(std::cos(ph), std::sin(ph));
        }
        return std::abs(acc);
    };
    double fr = f0;
    double d = 0.5 / (2.0 * static_cast<double>(n));
    while (d > 1e-7) {
        const double c1 = fr - d;
        const double c2 = fr + d;
        fr = coh(c1) > coh(c2) ? c1 : c2;
        d *= 0.5;
    }
    double sAbs = 0.0;
    for (const cd& v : z4) sAbs += std::abs(v);
    const double prom = coh(fr) / std::max(sAbs, 1e-30);
    return {fr, prom};
}

std::vector<cd> mfSymbols(const std::vector<cd>& y,
                          const std::vector<std::size_t>& idx,
                          const ClassifierConfig& cfg) {
    std::vector<cd> out;
    const std::size_t n = y.size();
    if (idx.size() < 2) return out;
    std::vector<long> ds;
    for (std::size_t i = 0; i + 1 < idx.size(); ++i)
        ds.push_back(static_cast<long>(idx[i + 1]) -
                     static_cast<long>(idx[i]));
    std::sort(ds.begin(), ds.end());
    const long sps = ds[ds.size() / 2];
    if (sps < 2) return out;
    std::vector<double> taps =
        detail::rrcTaps(sps, cfg.mfAlpha, static_cast<long>(cfg.mfSpan));
    double e = 0.0;
    for (double t : taps) e += t * t;
    e = std::sqrt(e);
    if (e <= 1e-12) return out;
    for (double& t : taps) t /= e;
    const std::size_t m = taps.size();
    const std::size_t half = (m - 1) / 2;
    const std::size_t tail = m - half;
    for (std::size_t i : idx) {
        if (i < half) continue;
        if (i + tail >= n) continue;
        const std::size_t base = i - half;
        cd acc(0.0, 0.0);
        for (std::size_t j = 0; j < m; ++j)
            acc += y[base + j] * taps[j];
        out.push_back(acc);
    }
    // Rampe am Crop-Rand: nahezu verschwindende Pseudosymbole
    // verfaelschen jede Konstellationsstatistik.
    if (out.size() >= 8) {
        std::vector<double> am;
        am.reserve(out.size());
        for (const cd& v : out) am.push_back(std::abs(v));
        std::sort(am.begin(), am.end());
        const double gate = cfg.mfMinRel * am[am.size() / 2];
        std::vector<cd> kept;
        kept.reserve(out.size());
        for (const cd& v : out) {
            if (std::abs(v) > gate) kept.push_back(v);
        }
        out.swap(kept);
    }
    return out;
}

// --------------------------------------------------------------------------
// Stufe 1
// --------------------------------------------------------------------------
bool stage1(const Features& f, const ClassifierConfig& cfg,
            CoarseClass* out) {
    if (f.snrDb < cfg.minSnrDb) return false;
    // Ton-Test zuerst, unabhaengig von der Bandbreite: starke Toene
    // behalten schmale Baender, doch der Splatter benachbarter Pulse kann
    // den Lauf legitim verbreitern; das darf keinen CW-Puls in die
    // Phasenanalyse schicken.
    std::size_t inside = 0;
    const std::size_t n1 = f.n >= 1 ? f.n - 1 : 0;
    for (double v : f.finst) {
        if (std::abs(v) < cfg.tCwFthr) ++inside;
    }
    const double frac =
        static_cast<double>(inside) /
        static_cast<double>(n1 > 0 ? n1 : 1);
    if (frac > cfg.tCwFrac && f.carrCoh > cfg.tCwCoh) {
        *out = CoarseClass::Tone;
        return true;
    }
    if (f.bandBins <= cfg.tBwCwBins) {
        *out = CoarseClass::FrequencyModulated;
        return true;
    }
    if (f.chirpR2 > cfg.tChirpR2 && f.chirpSweep > cfg.tChirpSweep) {
        *out = CoarseClass::FrequencyModulated;
        return true;
    }
    if (f.nfreqModes >= 2 && f.modeConc > cfg.tModeConc) {
        *out = CoarseClass::FrequencyModulated;
        return true;
    }
    if (f.paprWDb > cfg.tPaprOfdmDb && f.cv > cfg.tCvOfdm) {
        *out = CoarseClass::NoiseLike;
        return true;
    }
    *out = CoarseClass::PhaseModulated;
    return true;
}

// --------------------------------------------------------------------------
// QAM-Ordnung
// --------------------------------------------------------------------------
ModulationType qamTypeFromOrder(int m) {
    switch (m) {
        case 16: return ModulationType::QAM16;
        case 64: return ModulationType::QAM64;
        default: return ModulationType::QAM256;
    }
}

std::pair<ModulationType, double> qamOrder(const std::vector<cd>& s,
                                           const ClassifierConfig& cfg,
                                           double snrDb) {
    const double sig2 = 0.5;
    const double base = 1.0 / 3.0;
    const int orders[3] = {16, 64, 256};
    double scores[3];
    for (int oi = 0; oi < 3; ++oi) {
        const int m = orders[oi];
        const double l = std::round(std::sqrt(static_cast<double>(m)));
        const double d0 = std::sqrt(3.0 * sig2 / (l * l - 1.0));
        double eb = 0.0;
        bool have = false;
        for (int deg = -45; deg <= 45; deg += 2) {
            const double rot = kPi * static_cast<double>(deg) / 180.0;
            const double cs = std::cos(rot);
            const double sn = std::sin(rot);
            double e = 0.0;
            for (const cd& v : s) {
                const double p = (v.real() * cs - v.imag() * sn) / d0;
                const double kk = std::nearbyint((p - 1.0) * 0.5);
                const double dx = p - (2.0 * kk + 1.0);
                e += dx * dx;
            }
            e /= static_cast<double>(s.size());
            if (!have || e < eb) {
                eb = e;
                have = true;
            }
        }
        scores[oi] = base - eb;
    }
    // Python-Seitige Referenz: Sortierung der (Score, Ordnung)-Paare
    // absteigend; bei Score-Gleichheit dominiert die hoehere Ordnung.
    std::vector<std::pair<double, int>> ranked;
    for (int oi = 0; oi < 3; ++oi)
        ranked.emplace_back(scores[oi], orders[oi]);
    std::sort(ranked.begin(), ranked.end(),
              [](const std::pair<double, int>& a,
                 const std::pair<double, int>& b) {
                  if (a.first != b.first) return a.first > b.first;
                  return a.second > b.second;
              });
    const int bestOrder = ranked[0].second;
    const double sc0 = ranked[0].first;
    const double sc1 = ranked[1].first;
    const double lead = sc0 - std::max(sc1, 0.0);
    if (sc0 > cfg.qamScoreStrong && lead > cfg.qamScoreMargin) {
        double conf =
            std::min(0.9, 0.45 + 0.45 * std::sqrt(std::max(0.0, lead)));
        conf = std::max(0.5, conf);
        return {qamTypeFromOrder(bestOrder), conf};
    }
    // keine aufloesbare Levelstruktur: Dichte-Fallback. Die QAM256-Level-
    // Abstaende brauchen jedes dB; unterhalb der Schwelle ist die
    // Konstellation nicht mehr identifizierbar.
    if (snrDb >= cfg.qamDenseSnr)
        return {ModulationType::QAM256, 0.35};
    return {ModulationType::QAM64, 0.4};
}

// --------------------------------------------------------------------------
// Stufe 2
// --------------------------------------------------------------------------
std::pair<ModulationType, double> classifyPhase(const Features& f,
                                                const ClassifierConfig& cfg,
                                                const std::vector<cd>& y) {
    // OFDM via zyklischem Praefix (starke, orthogonale Evidenz)
    if (y.size() >= 256 && f.bandBins > cfg.ofdmMinBandBins) {
        const CpEvidence e = cpAutocorrPeak(y, cfg);
        if (e.rho > cfg.tCpPeak && e.tau >= cfg.tCpMinTau &&
            e.contrast > cfg.tCpContrast && e.clean >= cfg.tCpClean) {
            const bool sure = e.clean >= cfg.tCpCleanSure;
            if (sure || cpHalfConsistency(y, cfg, e.rho, e.tau)) {
                const double conf =
                    std::min(0.97,
                             0.55 + (e.rho - cfg.tCpPeak) * 1.5 +
                                 std::min(0.25, e.clean - cfg.tCpClean));
                return {ModulationType::OFDM, conf};
            }
        }
    }
    // OFDM-Fallback: Hoellkurven-Spielraum + glaettendes Spektrum
    if (f.paprWDb > cfg.tPaprOfdmDb &&
        f.flatDb < cfg.tFlatDb + 1.5) {
        return {ModulationType::OFDM, 0.5};
    }
    const std::vector<std::size_t> idx = estimateSymbolsRobust(f, cfg);
    if (idx.size() < cfg.timingMinSyms) {
        return {ModulationType::Unknown, 0.1};
    }
    std::vector<cd> s0 = mfSymbols(y, idx, cfg);
    if (s0.size() < cfg.timingMinSyms) {
        s0.clear();
        for (std::size_t i : idx) {
            const double a = f.env[i];
            const double p = f.ph[i];
            s0.push_back(cd(a * std::cos(p), a * std::sin(p)));
        }
    }
    const auto fp = carrierFromSymbols(s0);
    const double fres = fp.first;
    double sp = 0.0;
    for (const cd& v : s0) sp += std::norm(v);
    sp /= static_cast<double>(s0.size());
    const double sc = 1.0 / std::sqrt(std::max(sp, 1e-30));
    const std::size_t len = s0.size();
    std::vector<cd> s(len);
    for (std::size_t k = 0; k < len; ++k) {
        const double ph = -2.0 * kPi * fres * static_cast<double>(k);
        s[k] = s0[k] * sc * cd(std::cos(ph), std::sin(ph));
    }
    const double k42 = snrCorrectK42(kappa(s), f.snrDb);
    double am = 0.0;
    {
        double acc = 0.0;
        for (const cd& v : s) acc += std::abs(v);
        am = acc / static_cast<double>(len);
    }
    double vsig = 0.0;
    {
        double acc = 0.0;
        for (const cd& v : s) {
            const double d = std::abs(v) - am;
            acc += d * d;
        }
        acc /= static_cast<double>(len);
        vsig = std::max(0.0, acc - std::pow(10.0, -f.snrDb / 10.0));
    }
    const double cvSym = std::sqrt(vsig) / std::max(am, 1e-9);
    std::vector<double> ph(len);
    for (std::size_t k = 0; k < len; ++k)
        ph[k] = std::atan2(s[k].imag(), s[k].real());
    const double f2 = foldConc(ph, 2);
    const double f4 = foldConc(ph, 4);
    // PSK: starke Phasenfaltung + Kumulantenpositionierung
    if (f2 > cfg.tFold2 && k42 > 1.5) {
        return {ModulationType::BPSK, std::min(0.98, 0.5 + f2 * 0.4)};
    }
    if (f4 > cfg.tFold4 && 0.75 < k42 && k42 < 1.35) {
        return {ModulationType::QPSK, std::min(0.98, 0.5 + f4 * 0.4)};
    }
    if (k42 > 1.5 && cvSym < cfg.tCvPskHi) {
        return {ModulationType::BPSK, 0.55};
    }
    if (0.8 < k42 && k42 < 1.35 && cvSym < cfg.tCvPskHi && f4 > 0.55) {
        return {ModulationType::QPSK, 0.55};
    }
    if (cvSym < cfg.tCvConstMod) {
        // konstanter Modulus ohne erkennbare Faltstruktur
        return {ModulationType::Unknown, 0.2};
    }
    // QAM-Pfad
    auto res = qamOrder(s, cfg, f.snrDb);
    if (res.first == ModulationType::QAM64 ||
        res.first == ModulationType::QAM256) {
        const double d64 = std::abs(k42 - cfg.k42Qam64);
        const double d256 = std::abs(k42 - cfg.k42Qam256);
        if (d64 < d256) {
            res.second = std::min(res.second, 0.5 + d64 * 0.5);
        } else {
            res.second = std::min(res.second, 0.5 + d256 * 0.5);
        }
    }
    return res;
}

std::pair<ModulationType, double> stage2(const Features& f,
                                         const ClassifierConfig& cfg,
                                         CoarseClass coarse,
                                         const std::vector<cd>& y) {
    if (coarse == CoarseClass::Tone) {
        return {ModulationType::Unmodulated, 0.9};
    }
    if (coarse == CoarseClass::FrequencyModulated) {
        if (f.chirpR2 > cfg.tChirpR2 && f.chirpSweep > cfg.tChirpSweep) {
            return {ModulationType::Chirp,
                    std::min(0.98, 0.5 + 0.45 * f.chirpR2)};
        }
        if (f.nfreqModes >= 2) {
            return {ModulationType::FSK2, 0.8};
        }
        if (f.nfreqModes == 1) {
            return {ModulationType::Unmodulated, 0.6};
        }
        return {ModulationType::Unknown, 0.2};
    }
    if (coarse == CoarseClass::NoiseLike) {
        const double conf = 0.55 +
            std::min(0.4, 0.05 * (f.paprWDb - cfg.tPaprOfdmDb));
        if (f.flatDb < cfg.tFlatDb + 1.5) {
            return {ModulationType::OFDM, std::min(0.97, conf)};
        }
        return {ModulationType::OFDM, std::min(0.8, conf)};
    }
    return classifyPhase(f, cfg, y);
}

}  // namespace

// --------------------------------------------------------------------------
// Klassifikator-Fassade
// --------------------------------------------------------------------------
std::vector<PulseResult> ModulationClassifier::classify(
    const std::complex<double>* iq, std::size_t n) const {
    const ClassifierConfig& cfg = cfg_;
    std::vector<PulseResult> results;
    const PulseSegmenter seg(cfg);
    const BandExtractor bxs(cfg);
    for (const auto& ab : seg.detect(iq, n)) {
        const std::size_t a = ab.first;
        const std::size_t b = ab.second;
        BandExtractor::BandDetection det = bxs.bands(iq, a, b);
        std::vector<Band> bands = det.bands;
        if (bands.empty()) {
            bands.push_back(Band{0, static_cast<long>(det.nft)});
        }
        for (const Band& band : bands) {
            const BandExtractor::Slice slice =
                bxs.extract(iq, a, b, band, det.nft);
            for (const auto& crop : timeCrops(slice.y, cfg)) {
                const std::size_t c0 = crop.first;
                const std::size_t c1 = crop.second;
                const std::size_t len = c1 - c0;
                if (len < cfg.minPulseLen) continue;
                std::vector<cd> ys(slice.y.begin() +
                                       static_cast<std::ptrdiff_t>(c0),
                                   slice.y.begin() +
                                       static_cast<std::ptrdiff_t>(c1));
                double snrDb = 0.0;
                double flatDb = 0.0;
                estimateSnrFlat(det.psdDb, band, cfg, bands, &snrDb,
                                &flatDb);
                const Features f =
                    extractFeatures(ys, slice.fc, slice.bw,
                                    static_cast<std::size_t>(band.width()),
                                    snrDb, flatDb, cfg);
                CoarseClass coarse = CoarseClass::Unknown;
                if (!stage1(f, cfg, &coarse)) continue;
                const auto res = stage2(f, cfg, coarse, ys);
                PulseResult r;
                r.startSample = a + c0;
                r.numSamples = len;
                r.centerFreqNorm = static_cast<float>(
                    roundTo(wrap01(slice.fc), 5));
                r.bandwidthNorm =
                    static_cast<float>(roundTo(slice.bw, 5));
                r.snrDb = static_cast<float>(roundTo(f.snrDb, 1));
                r.coarse = coarse;
                r.modulation = res.first;
                r.confidence =
                    static_cast<float>(roundTo(res.second, 2));
                r.envelopeCv = static_cast<float>(roundTo(f.cv, 3));
                r.paprWindowedDb =
                    static_cast<float>(roundTo(f.paprWDb, 1));
                r.chirpR2 = static_cast<float>(roundTo(f.chirpR2, 2));
                r.chirpSweep = static_cast<float>(f.chirpSweep);
                r.dwellFraction =
                    static_cast<float>(roundTo(f.dwell, 2));
                r.freqModes = f.nfreqModes;
                r.spectralFlatnessDb =
                    static_cast<float>(roundTo(f.flatDb, 2));
                r.bandBins = static_cast<std::size_t>(band.width());
                results.push_back(r);
            }
        }
    }
    return results;
}

std::vector<PulseResult> ModulationClassifier::classifyF(
    const std::complex<float>* iq, std::size_t n) const {
    std::vector<std::complex<double>> buf(n);
    for (std::size_t i = 0; i < n; ++i) {
        buf[i] = std::complex<double>(static_cast<double>(iq[i].real()),
                                      static_cast<double>(iq[i].imag()));
    }
    return classify(buf.data(), buf.size());
}

}  // namespace iqclassify
