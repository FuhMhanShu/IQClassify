// IQClassify - Stufe 1 (grobe Klasse), Stufe 2 (Modulationsart) und
// Klassifikator-Fassade.
//
// Alle Schwellwerte und Entscheidungsequenzen entsprechen eins zu eins der
// gegen simulierten Pulsen (SNR 5..20 dB, sps 4/8/16) validierten
// Referenzimplementierung; siehe docs/METHODOLOGY.md.

#include "iqclassify/Classifier.hpp"


#include <atomic>
#include <thread>
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
        case ModulationType::FrankCode: return "FrankCode";
        case ModulationType::PCode: return "PCode";
        case ModulationType::BarkerCode: return "BarkerCode";
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

// Skaliert die sample-basierten Schwellwerte fuer ein um U hochgesampletes
// Signal (Folgestufen sehen die Rate U); symbol- und bin-relative Groessen
// bleiben unveraendert (Bin-Schwellen skalieren mit nft automatisch).
ClassifierConfig scaleForUpsample(const ClassifierConfig& c, std::size_t U) {
    if (U <= 1) return c;
    ClassifierConfig r = c;
    const auto m = [U](std::size_t v) { return v * U; };
    const auto ml = [U](long v) { return v * static_cast<long>(U); };
    // Sample-basierte Groessen skalieren mit U. segFloorGroup zaehlt Bloecke
    // (n/segWin ist rateninvariant) und bleibt unveraendert; ebenso Symbol-
    // und Histogramm-Anzahlen.
    r.segWin = m(c.segWin);
    r.segMargin = m(c.segMargin);
    r.minPulseLen = m(c.minPulseLen);
    r.minAnalyzeLen = m(c.minAnalyzeLen);
    r.detectWin = m(c.detectWin);
    // Bin-basierte Groessen skalieren mit U (nft waechst mit detectWin).
    r.tBwCwBins = m(c.tBwCwBins);
    r.ofdmMinBandBins = m(c.ofdmMinBandBins);
    r.satMergeMin = ml(c.satMergeMin);
    r.ofdmTauMin = ml(c.ofdmTauMin);
    r.ofdmTauMax = ml(c.ofdmTauMax);
    r.tCpMinTau = ml(c.tCpMinTau);
    r.cpCleanWin = ml(c.cpCleanWin);
    r.cpHalfTol = ml(c.cpHalfTol);
    r.cpIsolation = ml(c.cpIsolation);
    return r;
}

// Skaliert die sample-basierten Schwellwerte fuer einen um D dezimierten
// Kanal (die Folgestufen sehen den Kanal mit Rate 1/D); symbol- und
// verhaeltnisbasierte Groessen bleiben unveraendert.
ClassifierConfig scaleForDecimation(const ClassifierConfig& c,
                                    std::size_t D) {
    if (D <= 1) return c;
    ClassifierConfig r = c;
    const auto sd = [D](std::size_t v) {
        const std::size_t q = v / D;
        return q > 0 ? q : 1;
    };
    const auto sl = [D](long v) {
        const long q = v / static_cast<long>(D);
        return q > 0 ? q : 1;
    };
    // segFloorGroup zaehlt Bloecke (n/segWin ist rateninvariant) -> unveraendert.
    r.segWin = sd(c.segWin);
    r.segMargin = sd(c.segMargin);
    r.minPulseLen = sd(c.minPulseLen);
    r.detectWin = sd(c.detectWin);
    r.tBwCwBins = sd(c.tBwCwBins);
    r.ofdmTauMin = sl(c.ofdmTauMin);
    r.ofdmTauMax = sl(c.ofdmTauMax);
    r.ofdmMinBandBins = sd(c.ofdmMinBandBins);
    r.tCpMinTau = sl(c.tCpMinTau);
    r.cpCleanWin = sl(c.cpCleanWin);
    r.cpHalfTol = sl(c.cpHalfTol);
    r.cpIsolation = sl(c.cpIsolation);
    r.satMergeMin = sl(c.satMergeMin);
    return r;
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
    if (n < cfg.minAnalyzeLen) return e;
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
        // Komplexer Betrag: die Autokorrelation bei Lag t traegt die Phase
        // des Resttraegers (e^{-j*2pi*f_res*t}); nur der Realteil wuerde
        // bei sub-Bin-Restoffset den CP-Peak zufaellig ausloeschen.
        return std::abs(buf[static_cast<std::size_t>(t)]) /
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
        if (m < cfg.minAnalyzeLen) return false;
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
                std::abs(buf[static_cast<std::size_t>(t)]);
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
                                       std::size_t n, double bw,
                                       std::size_t nftOverride = 0,
                                       double bwDivLo = 1.5,
                                       double bwDivHi = 1.05,
                                       double promGate = 0.0) {
    const std::size_t nft =
        nftOverride > 0 ? nftOverride : detail::pow2ceil(2 * n);
    double mean = 0.0;
    for (double v : env2) mean += v;
    mean /= static_cast<double>(n);
    std::vector<cd> buf(nft, cd(0.0, 0.0));
    for (std::size_t i = 0; i < n; ++i)
        buf[i] = cd(env2[i] - mean, 0.0);
    // Zeitsignal vor der Transformation sichern: Die sub-Bin-Verfeinerung
    // maximiert die zyklische Kohaerenz UEBER DEM ZEITSIGNAL; buf selbst
    // wird durch die in-place-FFT zum Spektrum.
    const std::vector<cd> tsig = buf;
    detail::fftInPlace(buf, false);
    const long kLo = std::max(1L, static_cast<long>(
                                     std::floor(bw / bwDivLo *
                                                static_cast<double>(nft))));
    const long kHi = std::min<long>(
        static_cast<long>(nft / 4),
        static_cast<long>(std::ceil(bw / bwDivHi * static_cast<double>(nft))));
    if (kHi <= kLo + 2) return {0.0, 0.0};
    // Betragsquadrate ueber einen gegenueber dem Suchbereich erweiterten
    // K-Bereich; die Nachbarschafts-Median der Linien-Auswahl braucht auch
    // am Rand des Suchbereichs echte Nachbarn.
    const long kLoW = kLo / 2 > 1 ? kLo / 2 : 1;
    const long kHiW = std::min<long>(static_cast<long>(nft / 4),
                                     kHi + (kHi - kLo) / 2);
    std::vector<double> magsW;
    magsW.reserve(static_cast<std::size_t>(kHiW - kLoW));
    for (long k = kLoW; k < kHiW; ++k) {
        const cd& f = buf[static_cast<std::size_t>(k)];
        magsW.push_back(f.real() * f.real() + f.imag() * f.imag());
    }
    const auto magW = [&](long k) -> double {
        return magsW[static_cast<std::size_t>(k - kLoW)];
    };
    std::vector<double> mags(magsW.begin() +
                                 static_cast<std::ptrdiff_t>(kLo - kLoW),
                             magsW.begin() +
                                 static_cast<std::ptrdiff_t>(kHi - kLoW));
    std::sort(mags.begin(), mags.end());
    const double med = !mags.empty() ? mags[mags.size() / 2] : 0.0;
    // Linien-Auswahl: Bei kurzen Pulsen (wenige Symbole) liegt unter der
    // Symbolratenlinie ein breiter Buckel der Pulsflanken; die eigentliche
    // Ratenlinie ist dagegen schmal (1-2 Bins). Ist ein Prominenz-Gate
    // gesetzt, wird unter den Kandidaten oberhalb des Gates die schmalste
    // Linie (Staerke zur lokalen Nachbarschafts-Median) gewaehlt, nicht der
    // Buckel. Ohne Gate bleibt es beim reinen Maximum.
    const long win = (kHi - kLo) / 6 > 4 ? (kHi - kLo) / 6 : 4;
    long bestK = 0;
    double bestV = -1.0;
    double bestProm = -1.0;
    {
        long maxK = 0;
        double maxV = -1.0;
        for (long k = kLo; k < kHi; ++k) {
            const double v = magW(k);
            if (v > maxV) {
                maxV = v;
                maxK = k;
            }
        }
        bestK = maxK;
        bestV = maxV;
        bestProm = maxV / std::max(med, 1e-30);
    }
    if (promGate > 0.0 && med > 0.0 && bestK > kLo && bestK < kHi - 1) {
        // Tie-Breaker: Das reine Maximum kann bei kurzen Pulsen auf dem
        // breiten Flanken-Buckel liegen (wenige Symbole -> schwache Linie).
        // Gibt es ein deutlich schaerferes lokales Maximum (lokale Prominenz
        // >= Gate und >= 3x die des Maximums), wird auf dieses gewechselt.
        auto localScore = [&](long k) {
            std::vector<double> nb;
            nb.reserve(static_cast<std::size_t>(2 * win));
            for (long j = k - win; j <= k + win; ++j) {
                if (std::abs(j - k) <= 1) continue;
                nb.push_back(magW(j));
            }
            if (nb.empty()) return 0.0;
            std::sort(nb.begin(), nb.end());
            return magW(k) / std::max(nb[nb.size() / 2], 1e-30);
        };
        const double s0 = localScore(bestK);
        double s1 = -1.0;
        long k1 = bestK;
        for (long k = kLo + 1; k < kHi - 1; ++k) {
            const double v = magW(k);
            if (v <= magW(k - 1) || v <= magW(k + 1)) continue;
            const double sc = localScore(k);
            if (sc >= promGate && sc > s1) {
                s1 = sc;
                k1 = k;
            }
        }
        if (s1 >= 3.0 * s0) {
            bestK = k1;
            bestV = magW(k1);
            bestProm = s1;
        }
    }
    const double prom = bestProm > 0.0 ? bestProm : bestV / std::max(med, 1e-30);
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
    // Schnelle zyklische Kohaerenz: Drehung als inkrementelle Rekurrenz
    // (4 Mults pro Sample statt sin/cos — bei der Sub-Bin-Suche der
    // dominierende Kostenblock). Rotation m <- m * e^{-j*2pi*f}.
    auto coh = [&](double f) {
        cd acc(0.0, 0.0);
        const double cs = std::cos(2.0 * kPi * f);
        const double sn = -std::sin(2.0 * kPi * f);
        double cr = 1.0;
        double ci = 0.0;
        for (const cd& v : tsig) {
            acc += v * cd(cr, ci);
            const double nr = cr * cs - ci * sn;
            ci = cr * sn + ci * cs;
            cr = nr;
        }
        return std::abs(acc);
    };
    {
        const double fr0 = fr;
        const double lim = 1.5 / static_cast<double>(nft);
        double d = 3.0 / static_cast<double>(nft);
        while (d > 1e-7) {
            const double flo = fr - d;
            const double fhi = fr + d;
            fr = coh(flo) > coh(fhi) ? flo : fhi;
            if (fr < fr0 - lim) fr = fr0 - lim;
            if (fr > fr0 + lim) fr = fr0 + lim;
            d *= 0.5;
        }
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
    if (n < cfg.minAnalyzeLen) return {};
    std::vector<double> env2(n);
    for (std::size_t i = 0; i < n; ++i) env2[i] = f.env[i] * f.env[i];
    // Kurzpulse (< minAnalyzeLen): feineres Zero-Padding fuer eine scharfere
    // Ratenlinie (wenige Symbole -> grobes K-Raster sonst).
    const std::size_t nftOv =
        cfg.minAnalyzeLen < 256 && n < 512 ? detail::pow2ceil(8 * n) : 0;
    const auto rp = rsFromCyclic(env2, n, f.bw, nftOv, cfg.rsBwDivLo,
                                 cfg.rsBwDivHi, cfg.tRsProm);
    if (rp.first <= 0.0 || rp.second < cfg.tRsProm) return {};
    const double spsMax =
        std::max(2.0, static_cast<double>(n) / 12.0);
    double sps = 1.0 / rp.first;
    std::vector<std::size_t> pos;
    const auto grid = [&](double s) {
        std::vector<std::size_t> p;
        const double o = eyeOffset(f.env, n, s);
        long k = 0;
        while (o + static_cast<double>(k) * s <
               static_cast<double>(n) - 2.0) {
            const long ip = std::lround(o + static_cast<double>(k) * s);
            if (ip >= 1) p.push_back(static_cast<std::size_t>(ip));
            ++k;
        }
        return p;
    };
    if (sps >= 1.6 && sps <= spsMax) pos = grid(sps);
    if (pos.size() < cfg.timingMinSyms && f.bw > 0.0) {
        // Fallback fuer sehr kurze Pulse (wenige Symbole): die zyklische
        // Linie ist dann im breiten Flanken-Buckel vergraben. Die Bandbreite
        // eines PSK-Pulses ist ~1.2-1.5x die Symbolrate -> sps ~ 1.25/bw.
        const double spsNom = 1.25 / f.bw;
        if (spsNom >= 1.6 && spsNom <= spsMax) {
            auto pos2 = grid(spsNom);
            if (pos2.size() > pos.size()) {
                sps = spsNom;
                pos = std::move(pos2);
            }
        }
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

// 2-PSK-Faltkonzentration auf den energie-tragenden Samples (env-gefiltert):
// robust gegen Rausch-Guards am Crop-Rand (BPSK/Barker vs. OFDM/CP).
double fold2Gated(const Features& f, const ClassifierConfig& cfg) {
    double me = 0.0;
    for (double e : f.env) me += e;
    me /= static_cast<double>(f.env.size());
    std::vector<double> ph;
    ph.reserve(f.ph.size());
    for (std::size_t i = 0; i < f.ph.size(); ++i)
        if (f.env[i] >= cfg.histEnvRel * me) ph.push_back(f.ph[i]);
    return ph.size() >= 8 ? foldConc(ph, 2) : 0.0;
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
    // Auch kurze Pulse (ab 8 Symbolen) bekommen eine Traeger-Schaetzung;
    // ohne sie bleibt ein Resttraeger und die Phasenfaltung kollabiert.
    if (len < 8) return {0.0, 1.0};
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
// Frank-Code-Test (polyphase Radar-Code)
// --------------------------------------------------------------------------
// Die Momentanfrequenz eines Frank-Codes der Ordnung M ist eine Treppe mit
// M Stufen: jede Stufe dauert M*sps Samples, die Stufen-Mitten liegen
// gleichabstaendig auf dem Frequenzkreis und werden in der Zeit monoton
// durchlaufen. Rueckgabe: Ordnung M (>= tFrankMinModes) oder 0; r2Out
// erhaelt das R^2 der Stufen-Mitten gegen den Stufen-Index.
int detectFrankOrder(const Features& f, const ClassifierConfig& cfg,
                     double* r2Out = nullptr, double* platLenOut = nullptr,
                     double* dirOut = nullptr) {
    if (r2Out) *r2Out = 0.0;
    if (platLenOut) *platLenOut = 0.0;
    if (dirOut) *dirOut = 1.0;
    const std::size_t n = f.n;
    if (n < 96) return 0;
    const std::size_t n1 = n - 1;
    double me = 0.0;
    for (double e : f.env) me += e;
    me /= static_cast<double>(n);
    const double gate = cfg.histEnvRel * me;
    // Median-Filter ueber finst (killt die Spruenge an den Chipgrenzen).
    const std::size_t W = n1 / 256 > 3 ? n1 / 256 : 3;
    std::vector<double> fm(n1);
    for (std::size_t i = 0; i < n1; ++i) {
        std::vector<double> w;
        w.reserve(2 * W + 1);
        const std::size_t lo = i >= W / 2 ? i - W / 2 : 0;
        const std::size_t hi = std::min(n1 - 1, i + W / 2);
        for (std::size_t j = lo; j <= hi; ++j)
            if (f.env[j] >= gate) w.push_back(f.finst[j]);
        if (w.empty()) {
            fm[i] = f.finst[i];
            continue;
        }
        std::sort(w.begin(), w.end());
        fm[i] = w[w.size() / 2];
    }
    // Histogramm der gefilterten Werte ueber den beobachteten Bereich
    // [fmin, fmax] (der Sweep eines Frank-Codes belegt nur einen Teil des
    // Frequenzkreises).
    double fmin = 1e300, fmax = -1e300;
    std::size_t tot = 0;
    for (std::size_t i = 0; i < n1; ++i) {
        if (f.env[i] < gate) continue;
        fmin = std::min(fmin, fm[i]);
        fmax = std::max(fmax, fm[i]);
        ++tot;
    }
    if (tot < 48 || !(fmax - fmin > 1e-9)) return 0;
    const int Q = 64;
    const double scale = Q / (fmax - fmin);
    std::vector<double> hist(static_cast<std::size_t>(Q), 0.0);
    for (std::size_t i = 0; i < n1; ++i) {
        if (f.env[i] < gate) continue;
        int b = static_cast<int>((fm[i] - fmin) * scale);
        if (b >= Q) b = Q - 1;
        hist[static_cast<std::size_t>(b)] += 1.0;
    }
    // geglaettete Gruppen (wie das Moden-Histogramm in extractFeatures)
    std::vector<double> sm(static_cast<std::size_t>(Q), 0.0);
    for (int i = 0; i < Q; ++i) {
        double acc = hist[static_cast<std::size_t>(i)];
        if (i > 0) acc += hist[static_cast<std::size_t>(i - 1)];
        if (i < Q - 1) acc += hist[static_cast<std::size_t>(i + 1)];
        sm[static_cast<std::size_t>(i)] = acc / 3.0;
    }
    const double hmax = *std::max_element(sm.begin(), sm.end());
    if (hmax <= 0.0) return 0;
    struct Grp {
        double mass;
        double center;
    };
    std::vector<Grp> groups;
    {
        int i = 0;
        while (i < Q) {
            if (sm[static_cast<std::size_t>(i)] > cfg.modeBinThr * hmax) {
                int j = i;
                while (j + 1 < Q &&
                       sm[static_cast<std::size_t>(j + 1)] > cfg.modeBinThr * hmax)
                    ++j;
                double mass = 0.0, csum = 0.0;
                for (int k = i; k <= j; ++k) {
                    mass += hist[static_cast<std::size_t>(k)];
                    csum += hist[static_cast<std::size_t>(k)] *
                            (fmin + (k + 0.5) / scale);
                }
                groups.push_back({mass, mass > 0 ? csum / mass : 0.0});
                i = j + 1;
            } else {
                ++i;
            }
        }
    }
    if (groups.size() < 2) return 0;
    // Sweep ueber die +-0.5-Grenze: ins lineare Frame aufdrehen
    if (fmax - fmin > 0.4)
        for (Grp& g : groups)
            if (g.center < 0.0) g.center += 1.0;
    std::sort(groups.begin(), groups.end(),
              [](const Grp& a, const Grp& b) { return a.center < b.center; });
    // fein aufgeteilte Stufen (Rausch-Splits) verschmelzen
    {
        bool merged = true;
        while (merged && groups.size() > 2) {
            merged = false;
            std::vector<double> g0;
            for (std::size_t k = 0; k + 1 < groups.size(); ++k)
                g0.push_back(groups[k + 1].center - groups[k].center);
            std::vector<double> gs = g0;
            std::sort(gs.begin(), gs.end());
            const double med = gs[gs.size() / 2];
            for (std::size_t k = 0; k + 1 < groups.size(); ++k) {
                if (g0[k] < 0.45 * med) {
                    Grp& a = groups[k];
                    const Grp b = groups[k + 1];
                    const double m = a.mass + b.mass;
                    a.center = (a.center * a.mass + b.center * b.mass) / m;
                    a.mass = m;
                    groups.erase(groups.begin() + static_cast<std::ptrdiff_t>(k + 1));
                    merged = true;
                    break;
                }
            }
        }
    }
    // mediane Stufendauer (vor dem Verschmelzen fein geteilter Stufen):
    // eine Stufe dauert ~ M Chips -> Chip-Laenge = platLen/M
    if (platLenOut)
        *platLenOut = static_cast<double>(tot) /
                      static_cast<double>(groups.size());
    const int P = static_cast<int>(groups.size());
    if (P < cfg.tFrankMinModes || P > cfg.tFrankMaxModes) return 0;
    // Coverage: Anteil der gated Samples in Stufen
    double msum = 0.0;
    for (const Grp& g : groups) msum += g.mass;
    if (msum < cfg.tFrankCover * static_cast<double>(tot)) return 0;
    // interne Stufen-Abstaende gleichmaessig?
    std::vector<double> gaps;
    for (int k = 0; k < P - 1; ++k)
        gaps.push_back(groups[static_cast<std::size_t>(k) + 1].center -
                       groups[static_cast<std::size_t>(k)].center);
    const double gmin = *std::min_element(gaps.begin(), gaps.end());
    const double gmax = *std::max_element(gaps.begin(), gaps.end());
    if (gmin <= 0.0 || gmax > 2.2 * gmin) return 0;
    // R^2 der Stufen-Mitten gegen den Index
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    for (int k = 0; k < P; ++k) {
        const double x = k, y = groups[static_cast<std::size_t>(k)].center;
        sx += x; sy += y; sxx += x * x; sxy += x * y;
    }
    const double den = P * sxx - sx * sx;
    double r2 = 0.0;
    if (std::abs(den) > 1e-30) {
        const double b = (P * sxy - sx * sy) / den;
        const double a = (sy - b * sx) / P;
        double ssr = 0.0, sst = 0.0;
        for (int k = 0; k < P; ++k) {
            const double x = k, y = groups[static_cast<std::size_t>(k)].center;
            const double yh = a + b * x;
            ssr += (y - yh) * (y - yh);
            sst += (y - sy / P) * (y - sy / P);
        }
        r2 = sst > 1e-30 ? 1.0 - ssr / sst : 0.0;
    }
    if (r2 < cfg.tFrankR2) return 0;
    if (r2Out) *r2Out = r2;
    // Stufendauer in der Zeit: jeder Stufe ein langer Run (beim Frank-Code
    // M*sps Samples; beim PSK-Rauschen nur kurze Runs).
    {
        std::vector<std::size_t> maxRun(static_cast<std::size_t>(P), 0);
        int prev = -1;
        std::size_t run = 0;
        // Zeitfolge der Stufen: Run-Mittel der gefilterten finst je
        // zusammenhaengendem Run (liefert die Richtung der Treppe).
        std::vector<double> runMean;
        std::vector<double> runPos;
        std::size_t runStart = 0;
        double runSum = 0.0;
        for (std::size_t i = 0; i < n1; ++i) {
            if (f.env[i] < gate) continue;
            int best = 0;
            double bd = 1e300;
            for (int k = 0; k < P; ++k) {
                double d = fm[i] - groups[static_cast<std::size_t>(k)].center;
                d = std::fmod(d + 1.5, 1.0) - 0.5;
                d = std::abs(d);
                if (d < bd) { bd = d; best = k; }
            }
            if (prev >= 0 && best != prev) {
                maxRun[static_cast<std::size_t>(prev)] =
                    std::max(maxRun[static_cast<std::size_t>(prev)], run);
                runMean.push_back(runSum / static_cast<double>(run));
                runPos.push_back(static_cast<double>(runStart));
                run = 0;
                runStart = i;
                runSum = 0.0;
            }
            ++run;
            runSum += fm[i];
            prev = best;
        }
        if (prev >= 0) {
            maxRun[static_cast<std::size_t>(prev)] =
                std::max(maxRun[static_cast<std::size_t>(prev)], run);
            runMean.push_back(runSum / static_cast<double>(run));
            runPos.push_back(static_cast<double>(runStart));
        }
        const std::size_t lMin = tot / (3 * static_cast<std::size_t>(P)) > 8
                                     ? tot / (3 * static_cast<std::size_t>(P))
                                     : 8;
        for (int k = 0; k < P; ++k)
            if (maxRun[static_cast<std::size_t>(k)] < lMin) return 0;
        // Richtung der Treppe in der Zeit: die Run-Mittel werden mit der
        // Chip-Periode p = 1/sps entwrappt (sonst kippt das Vorzeichen am
        // Faltungsrand +/-0.5/sps), dann linear gegen die Run-Position
        // gefittet. sps kommt aus der Run-Struktur: eine Stufe (Zeile des
        // Codes) ist ein langer Run der Laenge M*sps, also sps ~ rowLen^2 /
        // span. Frank steigt (> 0), P1/P2 fallen (< 0).
        if (dirOut && runMean.size() >= 3) {
            std::vector<double> lm, lp;
            for (std::size_t k = 0; k < runMean.size(); ++k) {
                const double len =
                    (k + 1 < runMean.size() ? runPos[k + 1]
                                            : static_cast<double>(n1)) -
                    runPos[k];
                if (len >= static_cast<double>(lMin)) {
                    lm.push_back(runMean[k]);
                    lp.push_back(runPos[k]);
                }
            }
            if (lm.size() >= 3) {
                std::vector<double> ls;
                for (std::size_t k = 0; k + 1 < lm.size(); ++k)
                    ls.push_back(lp[k + 1] - lp[k]);
                std::sort(ls.begin(), ls.end());
                const double rowLen = ls[ls.size() / 2];
                const double span = lp.back() - lp.front() + rowLen;
                double p = span / (rowLen * rowLen);
                std::vector<double> steps;
                for (std::size_t k = 0; k + 1 < lm.size(); ++k) {
                    const double d =
                        std::fmod(lm[k + 1] - lm[k] + 1.5, 1.0) - 0.5;
                    steps.push_back(std::abs(d));
                }
                std::sort(steps.begin(), steps.end());
                const double medStep = steps[steps.size() / 2];
                if (p < 2.0 * medStep || p > 0.45) p = 1.0;
                std::vector<double> u(lm.size());
                u[0] = lm[0];
                for (std::size_t k = 1; k < lm.size(); ++k) {
                    const double d = lm[k] - u[k - 1];
                    double r = std::fmod(d + p / 2.0, p);
                    if (r < 0.0) r += p;
                    r -= p / 2.0;
                    u[k] = u[k - 1] + r;
                }
                double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
                const double R = static_cast<double>(lm.size());
                for (std::size_t k = 0; k < lm.size(); ++k) {
                    sx += lp[k];
                    sy += u[k];
                    sxx += lp[k] * lp[k];
                    sxy += lp[k] * u[k];
                }
                const double den = R * sxx - sx * sx;
                if (std::abs(den) > 1e-30)
                    *dirOut = (R * sxy - sx * sy) / den;
            }
        }
    }
    return P;
}

// --------------------------------------------------------------------------
// Phasen-Gesetze der M x M polyphase Codes (Chip k = i*M + j)
// --------------------------------------------------------------------------
double polyLawPhase(int law, int M, int i, int j) {
    switch (law) {
        case 0:
            return 2.0 * kPi * i * j / M;               // Frank
        case 1:
            return -(kPi / M) * (M - 2 * i) * (M - 2 * j);  // P1
        default:
            return -(kPi / M) * (M - 2 * i) * (M - 2 * j - 1);  // P2
    }
}

// RMS-Phasenfehler der an den Chip-Mitten rekonstruierten Phasen gegen ein
// Phasen-Gesetz (Offset-frei ueber den zirkularen Mittelwert).
double lawErr(const std::vector<cd>& y, const Features& f, double chipLen,
              int M, int law) {
    const std::size_t N = static_cast<std::size_t>(M) * M;
    double me = 0.0;
    for (double e : f.env) me += e;
    me /= static_cast<double>(f.env.size());
    std::vector<cd> z;
    z.reserve(N);
    for (std::size_t k = 0; k < N; ++k) {
        std::size_t ci =
            static_cast<std::size_t>((static_cast<double>(k) + 0.5) * chipLen);
        if (ci >= y.size()) ci = y.size() - 1;
        if (ci >= f.env.size() || f.env[ci] < 0.4 * me) continue;  // Rand-Chips
        double ph = std::atan2(y[ci].imag(), y[ci].real()) -
                    2.0 * kPi * f.residF * static_cast<double>(ci);
        ph -= polyLawPhase(law, M, static_cast<int>(k / M),
                           static_cast<int>(k % M));
        z.push_back(cd(std::cos(ph), std::sin(ph)));
    }
    if (z.size() < 9) return 1e300;
    cd acc(0.0, 0.0);
    for (const cd& v : z) acc += v;
    const double meanPh = std::atan2(acc.imag(), acc.real());
    double err = 0.0;
    for (const cd& v : z) {
        double d = std::atan2(v.imag(), v.real()) - meanPh;
        d = std::fmod(d + 3.0 * kPi, 2.0 * kPi) - kPi;
        err += d * d;
    }
    return std::sqrt(err / static_cast<double>(z.size()));
}

// Bestes Phasen-Gesetz fuer die Code-Phasen des Slices: probiert Frank/P1/P2
// mit Ordnungen um Mapprox und gibt 0/1/2 zurueck (0 = Frank).
int polyLawMatch(const std::vector<cd>& y, const Features& f, int Mapprox,
                 double platLen, int* Mbest = nullptr,
                 double* errBest = nullptr) {
    (void)platLen;  // Chip-Laenge wird aus y.size()/M^2 geschaetzt
    int bestLaw = 0;
    int bestM = Mapprox > 3 ? Mapprox : 3;
    double bestErr = 1e300;
    double bestFrankErr = 1e300;
    int bestFrankM = bestM;
    const int lo = Mapprox > 4 ? Mapprox : 4;  // Merging unterschaetzt M nie
    const int hi = 2 * Mapprox + 2 < 48 ? 2 * Mapprox + 2 : 48;
    for (int M = lo; M <= hi; ++M) {
        const std::size_t N = static_cast<std::size_t>(M) * M;
        if (N > y.size() || N < 9) continue;
        // Chip-Laenge: der Crop enthaelt ca. M^2 Chips; um n/N fein suchen
        // (1%-Schritte, damit die Drift ueber M^2 Chips klein bleibt).
        const double chip0 = static_cast<double>(y.size()) / N;
        for (int law = 0; law < 3; ++law) {
            double bestE = 1e300;
            for (int o = -8; o <= 8; ++o) {
                double chipLen = chip0 * (1.0 + 0.01 * o);
                if (chipLen < 1.0) chipLen = 1.0;
                bestE = std::min(bestE, lawErr(y, f, chipLen, M, law));
            }
            const double err = bestE;
            if (err < bestErr) {
                bestErr = err;
                bestLaw = law;
                bestM = M;
            }
            if (law == 0 && err < bestFrankErr) {
                bestFrankErr = err;
                bestFrankM = M;
            }
        }
    }
    if (Mbest) *Mbest = bestM;
    if (errBest) *errBest = bestErr;
    // Nur bei klarem Vorsprung gegen das beste Frank-Gesetz (sonst FrankCode)
    if (bestLaw != 0 && bestFrankErr < bestErr + 0.25) {
        if (Mbest) *Mbest = bestFrankM;
        if (errBest) *errBest = bestFrankErr;
        return 0;
    }
    return bestLaw;
}

// --------------------------------------------------------------------------
// Matched-Filter-Vergleich der Phasen-Gesetze (robust gegen Impulsformung):
// gleitende Korrelation des Slices mit dem ungefilterten Chip-Template.
// --------------------------------------------------------------------------
double mfLawScore(const std::vector<cd>& y, int M, int law, double chipLen,
                  double offDiv = 4.0) {
    const std::size_t N = static_cast<std::size_t>(M) * M;
    if (N < 9 || chipLen < 1.0) return 0.0;
    const std::size_t T =
        static_cast<std::size_t>(static_cast<double>(N) * chipLen + chipLen);
    // Template darf laenger als der Slice sein (Rand-Guards): Korrelation
    // laeuft ueber den verfuegbaren Ueberlapp, solange der groesste Teil des
    // Templates abgedeckt ist.
    if (static_cast<double>(T) > 2.0 * static_cast<double>(y.size())) return 0.0;
    std::vector<cd> t(T, cd(0.0, 0.0));
    for (std::size_t k = 0; k < N; ++k) {
        const double ph =
            polyLawPhase(law, M, static_cast<int>(k / M),
                         static_cast<int>(k % M));
        const std::size_t c = static_cast<std::size_t>(k * chipLen);
        if (c < T) t[c] = cd(std::cos(ph), std::sin(ph));
    }
    double ey = 0.0;
    for (const cd& v : y) ey += std::norm(v);
    if (ey < 1e-30) return 0.0;
    double step = chipLen / offDiv;
    if (step < 1.0) step = 1.0;
    double best = 0.0;
    for (double off = 0.0; off < static_cast<double>(y.size()); off += step) {
        const std::size_t o0 = static_cast<std::size_t>(off);
        const std::size_t W = std::min(T, y.size() - o0);
        if (static_cast<double>(W) < 0.6 * static_cast<double>(T)) break;
        cd acc(0.0, 0.0);
        for (std::size_t i = 0; i < W; ++i)
            acc += y[o0 + i] * std::conj(t[i]);
        const double s =
            std::abs(acc) /
            std::sqrt(ey * static_cast<double>(N)) *
            std::sqrt(static_cast<double>(W) / static_cast<double>(y.size()));
        if (s > best) best = s;
    }
    return best;
}

// Bestes Phasen-Gesetz per Matched Filter ueber Ordnungen um Mapprox.
// Rueckgabe: 0 = Frank, 1 = P1, 2 = P2. P gewinnt nur, wenn sein Score bei
// derselben Ordnung M klar ueber dem Frank-Score liegt (Verhaeltnis 1.15).
// scoreBest/frankBest liefern die Scores der Best-Ordnung.
int polyLawMf(const std::vector<cd>& y, int Mapprox, double platLen,
              int* Mbest = nullptr, double* scoreBest = nullptr,
              double* frankBest = nullptr, bool fineGrid = true) {
    (void)platLen;  // Chip-Laenge wird aus y.size()/M^2 geschaetzt
    int bestLaw = 0;
    int bestM = Mapprox > 3 ? Mapprox : 3;
    double bestP = 0.0;
    double maxF = 0.0;
    const int lo = Mapprox > 4 ? Mapprox : 4;
    const int hi = 2 * Mapprox + 2 < 20 ? 2 * Mapprox + 2 : 20;
    for (int M = lo; M <= hi; ++M) {
        const std::size_t N = static_cast<std::size_t>(M) * M;
        if (N + 8 > y.size()) break;
        const double chip0 = static_cast<double>(y.size()) / N;
        double scF = 0.0;
        double scP[2] = {0.0, 0.0};
        for (int law = 0; law < 3; ++law) {
            double sc = 0.0;
            if (fineGrid) {
                for (int o = -8; o <= 8; ++o) {
                    double cl = chip0 * (1.0 + 0.01 * o);
                    if (cl < 1.0) cl = 1.0;
                    sc = std::max(sc, mfLawScore(y, M, law, cl, 4.0));
                }
            } else {
                for (int o = -4; o <= 4; ++o) {
                    double cl = chip0 * (1.0 + 0.02 * o);
                    if (cl < 1.0) cl = 1.0;
                    sc = std::max(sc, mfLawScore(y, M, law, cl, 2.0));
                }
            }
            if (law == 0) scF = sc;
            else scP[law - 1] = sc;
        }
        maxF = std::max(maxF, scF);
        const double p = std::max(scP[0], scP[1]);
        if (p > bestP) {
            bestP = p;
            bestM = M;
            bestLaw = scP[0] >= scP[1] ? 1 : 2;
        }
    }
    if (Mbest) *Mbest = bestM;
    if (scoreBest) *scoreBest = bestP;
    if (frankBest) *frankBest = maxF;
    // P gewinnt nur bei klarem Vorsprung des besten P-Scores gegen den
    // besten Frank-Score im ganzen Ordnungsbereich (sonst Frank).
    if (bestLaw != 0 && bestP < 1.15 * maxF) return 0;
    return bestLaw;
}

// --------------------------------------------------------------------------
// Barker-Match direkt auf dem Slice: Symboltakt um die Schaetzung fein suchen
// (die Barker-Periodizitaet stoert die Symbolratenschaetzung), Zeichen an den
// Chip-Mitten (env-gefiltert, auf den Resttraeger derotiert via f.ph), gegen
// alle Rotationen und beide Vorzeichen der Barker-Folgen testen.
// Rueckgabe: Code-Laenge (7/11/13) oder 0; bestReps/bestScore/bestSps optional.
// --------------------------------------------------------------------------
int barkerMatchSlice(const std::vector<cd>& y, const Features& f,
                     const ClassifierConfig& cfg, double spsEst,
                     int* bestReps = nullptr, double* bestScore = nullptr,
                     double* bestSps = nullptr) {
    static const int b13[] = {1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1};
    static const int b11[] = {1, 1, 1, -1, -1, -1, 1, -1, -1, 1, -1};
    static const int b7[] = {1, 1, 1, -1, -1, 1, -1};
    static const int* seqs[] = {b13, b11, b7};
    static const int lens[] = {13, 11, 7};
    if (spsEst < 2.0 || y.size() < 96) return 0;
    double me = 0.0;
    for (double e : f.env) me += e;
    me /= static_cast<double>(f.env.size());
    const double gate = cfg.histEnvRel * me;
    int bestLen = 0;
    int bestRep = 0;
    double bestSc = 0.0;
    double bestSp = spsEst;
    const double lo = spsEst * 0.65;
    const double hi = spsEst * 1.45;
    for (int t = 0; t < 3; ++t) {
        const int l = lens[t];
        if (l < cfg.tBarkerMinLen) continue;
        for (double sps = lo; sps <= hi; sps += 0.1) {
            const int nsym =
                static_cast<int>(static_cast<double>(y.size()) / sps);
            if (nsym < 2 * l) continue;
            for (double off = 0.0; off < sps; off += sps / 8.0) {
                std::vector<int> sg(static_cast<std::size_t>(nsym), 0);
                int valid = 0;
                for (int k = 0; k < nsym; ++k) {
                    const std::size_t ci = static_cast<std::size_t>(
                        off + (static_cast<double>(k) + 0.5) * sps);
                    if (ci >= y.size() || ci >= f.env.size()) continue;
                    if (f.env[ci] < gate) continue;
                    sg[static_cast<std::size_t>(k)] =
                        std::cos(f.ph[ci]) >= 0.0 ? 1 : -1;
                    ++valid;
                }
                if (valid < 2 * l) continue;
                for (int rot = 0; rot < l; ++rot) {
                    for (int sign = 0; sign < 2; ++sign) {
                        int match = 0;
                        for (int k = 0; k < nsym; ++k) {
                            if (sg[static_cast<std::size_t>(k)] == 0) continue;
                            const int expect = (sign ? -1 : 1) *
                                               seqs[t][(k + rot) % l];
                            if (sg[static_cast<std::size_t>(k)] == expect)
                                ++match;
                        }
                        const double sc = static_cast<double>(match) /
                                          static_cast<double>(valid);
                        if (sc > bestSc) {
                            bestSc = sc;
                            bestLen = l;
                            bestRep = valid / l;
                            bestSp = sps;
                        }
                    }
                }
            }
        }
    }
    if (bestSc >= 0.95 && bestLen >= cfg.tBarkerMinLen) {
        if (bestReps) *bestReps = bestRep;
        if (bestScore) *bestScore = bestSc;
        if (bestSps) *bestSps = bestSp;
        return bestLen;
    }
    return 0;
}

// --------------------------------------------------------------------------
// Billiger Vorab-Test fuer den P-Code-Fallback: lange Runs einzelner
// f_inst-Level (Stufen-Signatur eines Polyphase-Codes). FSK2 und Chirp
// haben keine langen Runs und sparen so den teuren Matched-Filter.
// --------------------------------------------------------------------------
bool hasLongLevelRuns(const Features& f, const ClassifierConfig& cfg) {
    const std::size_t n = f.n;
    if (n < 128) return false;
    const std::size_t n1 = n - 1;
    double me = 0.0;
    for (double e : f.env) me += e;
    me /= static_cast<double>(n);
    const double gate = cfg.histEnvRel * me;
    const std::size_t W = n1 / 256 > 3 ? n1 / 256 : 3;
    std::vector<double> fm(n1);
    for (std::size_t i = 0; i < n1; ++i) {
        std::vector<double> w;
        w.reserve(2 * W + 1);
        const std::size_t lo = i >= W / 2 ? i - W / 2 : 0;
        const std::size_t hi = std::min(n1 - 1, i + W / 2);
        for (std::size_t j = lo; j <= hi; ++j)
            if (f.env[j] >= gate) w.push_back(f.finst[j]);
        if (w.empty()) {
            fm[i] = f.finst[i];
            continue;
        }
        std::sort(w.begin(), w.end());
        fm[i] = w[w.size() / 2];
    }
    std::vector<double> g;
    for (std::size_t i = 0; i < n1; ++i)
        if (f.env[i] >= gate) g.push_back(fm[i]);
    if (g.size() < 96) return false;
    double fmin = 1e300, fmax = -1e300;
    for (double v : g) {
        fmin = std::min(fmin, v);
        fmax = std::max(fmax, v);
    }
    if (fmax - fmin < 1e-9) return true;  // ein Level = ein langer Run
    const int Q = 16;
    const double scale = Q / (fmax - fmin);
    std::vector<int> bin(n1, -1);
    std::vector<double> h(static_cast<std::size_t>(Q), 0.0);
    std::size_t tot = 0;
    for (std::size_t i = 0; i < n1; ++i) {
        if (f.env[i] < gate) continue;
        int b = static_cast<int>((fm[i] - fmin) * scale);
        if (b >= Q) b = Q - 1;
        bin[i] = b;
        h[static_cast<std::size_t>(b)] += 1.0;
        ++tot;
    }
    std::vector<std::size_t> run(static_cast<std::size_t>(Q), 0);
    std::vector<std::size_t> cur(static_cast<std::size_t>(Q), 0);
    for (std::size_t i = 0; i < n1; ++i) {
        if (bin[i] < 0) {
            std::fill(cur.begin(), cur.end(), 0);
            continue;
        }
        const int b = bin[i];
        ++cur[static_cast<std::size_t>(b)];
        run[static_cast<std::size_t>(b)] =
            std::max(run[static_cast<std::size_t>(b)], cur[static_cast<std::size_t>(b)]);
        for (int k = 0; k < Q; ++k)
            if (k != b) cur[static_cast<std::size_t>(k)] = 0;
    }
    const std::size_t lMin = tot / 32 > 8 ? tot / 32 : 8;
    for (int b = 0; b < Q; ++b)
        if (h[static_cast<std::size_t>(b)] > 0.02 * static_cast<double>(tot) &&
            run[static_cast<std::size_t>(b)] >= lMin)
            return true;
    return false;
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
    if (detectFrankOrder(f, cfg) >= cfg.tFrankMinModes) {
        *out = CoarseClass::FrequencyModulated;
        return true;
    }
    if ((f.chirpR2 > cfg.tChirpR2 &&
         (f.chirpSweep > cfg.tChirpSweep ||
          f.chirpSweepAbs > cfg.tChirpSweepAbs)) ||
        (f.chirpFoldR2 > cfg.tChirpFoldR2 &&
         (f.chirpFoldSweep > cfg.tChirpSweep ||
          f.chirpFoldSweepAbs > cfg.tChirpSweepAbs))) {
        *out = CoarseClass::FrequencyModulated;
        return true;
    }
    if (f.nfreqModes >= 2 && f.modeConc > cfg.tModeConc) {
        *out = CoarseClass::FrequencyModulated;
        return true;
    }
    // NoiseLike (OFDM-Verdacht) nur ohne klare 2-PSK-Faltung: BPSK/Barker
    // haben ein Linienspektrum und erhoehten Huellkurven-Spielraum.
    if (f.paprWDb > cfg.tPaprOfdmDb && f.cv > cfg.tCvOfdm &&
        fold2Gated(f, cfg) < cfg.tFold2) {
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
    // OFDM via zyklischem Praefix (starke, orthogonale Evidenz). Nicht bei
    // klarer 2-PSK-Faltung: ein repetitiver Barker-Code hat dieselbe
    // Periodizitaet wie ein zyklischer Praefix.
    if (y.size() >= cfg.minAnalyzeLen && f.bandBins > cfg.ofdmMinBandBins &&
        fold2Gated(f, cfg) < cfg.tFold2) {
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
    // OFDM-Fallback: Hoellkurven-Spielraum + glaettendes Spektrum. Nicht bei
    // klarer 2-PSK-Faltung: BPSK/Barker haben ein Linienspektrum (niedriges
    // flatDb) und duerfen nicht als OFDM fehlgehen.
    if (f.paprWDb > cfg.tPaprOfdmDb &&
        f.flatDb < cfg.tFlatDb + 1.5 &&
        fold2Gated(f, cfg) < cfg.tFold2) {
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
    const double spsEst = idx.size() >= 2
        ? static_cast<double>(idx.back() - idx.front()) /
              static_cast<double>(idx.size() - 1)
        : 8.0;
    if (f2 > cfg.tFold2 && k42 > 1.5) {
        // Feinklasse: BPSK-Puls mit Barker-Code?
        const int bl = barkerMatchSlice(y, f, cfg, spsEst);
        if (bl > 0)
            return {ModulationType::BarkerCode,
                    std::min(0.98, 0.65 + 0.02 * bl)};
        return {ModulationType::BPSK, std::min(0.98, 0.5 + f2 * 0.4)};
    }
    if (f4 > cfg.tFold4 && 0.75 < k42 && k42 < 1.35) {
        return {ModulationType::QPSK, std::min(0.98, 0.5 + f4 * 0.4)};
    }
    if (k42 > 1.5 && cvSym < cfg.tCvPskHi) {
        const int bl = barkerMatchSlice(y, f, cfg, spsEst);
        if (bl > 0) return {ModulationType::BarkerCode, 0.8};
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
        // Frank-Code zuerst: seine Treppe hat zwar eine hohe lineare
        // Komponente, aber die Stufen-Struktur ist spezifischer als der
        // Chirp-Fit.
        // enger Sub-Slice wie beim Fein-Timing (0.25*max): der Matched
        // Filter arbeitet auf dem Puls ohne Rand-Guards.
        std::vector<cd> yt = y;
        {
            double emax = 0.0;
            for (const cd& v : y) emax = std::max(emax, std::abs(v));
            if (emax > 0.0 && y.size() > 2) {
                const double thr = 0.25 * emax;
                std::size_t a0 = 0;
                std::size_t a1 = y.size();
                while (a0 < y.size() && std::abs(y[a0]) <= thr) ++a0;
                while (a1 > a0 && std::abs(y[a1 - 1]) <= thr) --a1;
                if (a1 > a0 && a1 - a0 >= 96)
                    yt.assign(y.begin() + static_cast<std::ptrdiff_t>(a0),
                              y.begin() + static_cast<std::ptrdiff_t>(a1));
            }
        }
        double frankR2 = 0.0;
        double platLen = 0.0;
        const int frankM = detectFrankOrder(f, cfg, &frankR2, &platLen);
        if (frankM >= cfg.tFrankMinModes) {
            // Frank vs. P1/P2 per Matched Filter der Phasen-Gesetze
            // (robust gegen Impulsformung): 0 = Frank, 1/2 = P1/P2.
            int Mb = frankM;
            double sc = 0.0, fsc = 0.0;
            const int law = polyLawMf(yt, frankM, platLen, &Mb, &sc, &fsc);
            const double conf = std::min(0.95, 0.55 + 0.4 * frankR2);
            if (law != 0) return {ModulationType::PCode, conf};
            return {ModulationType::FrankCode, conf};
        }
        if ((f.chirpR2 > cfg.tChirpR2 &&
         (f.chirpSweep > cfg.tChirpSweep ||
          f.chirpSweepAbs > cfg.tChirpSweepAbs)) ||
        (f.chirpFoldR2 > cfg.tChirpFoldR2 &&
         (f.chirpFoldSweep > cfg.tChirpSweep ||
          f.chirpFoldSweepAbs > cfg.tChirpSweepAbs))) {
            return {ModulationType::Chirp,
                    std::min(0.98, 0.5 + 0.45 * f.chirpR2)};
        }
        // P1/P2 kleiner Ordnung: die Stufen sind durch die Impulsformung so
        // verschmiert, dass das Histogramm keine Treppe findet — die
        // Phasen-Gesetze per Matched Filter direkt testen, bevor es FSK2 wird.
        {
            // ohne lange f_inst-Runs (FSK2/Chirp) lohnt der teure
            // Matched-Filter nicht — billiger Vorab-Test first.
            if (hasLongLevelRuns(f, cfg)) {
                int Mb = 0;
                double sc = 0.0, fsc = 0.0;
                const int law = polyLawMf(yt, 8, 0.0, &Mb, &sc, &fsc, false);
                // deutlich staerkeres Kriterium als im Treppen-Pfad: nur bei
                // klarem Score-Vorsprung (echte P-Codes: sc > 0.18,
                // ratio > 2.7; FSK2/PSK-Rauschen: sc < 0.09, ratio < 1.5).
                if (law != 0 && sc > 0.10 && sc > 1.8 * fsc) {
                    return {ModulationType::PCode,
                            std::min(0.9, 0.5 + 2.0 * sc)};
                }
            }
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

// Vollsuche der Symbolratenlinie im zyklischen Spektrum eines reellen
// Signals (z. B. |Delta f_inst| bei FSK, wo die Huellkurve konstant ist):
// K in [12*nft/n, 0.625*nft], also sps in [1.6, n/12].
std::pair<double, double> rsLineFull(const std::vector<double>& v) {
    const std::size_t n = v.size();
    if (n < 8) return {0.0, 0.0};
    const std::size_t nft = detail::pow2ceil(2 * n);
    double mean = 0.0;
    for (double x : v) mean += x;
    mean /= static_cast<double>(n);
    std::vector<cd> buf(nft, cd(0.0, 0.0));
    for (std::size_t i = 0; i < n; ++i) buf[i] = cd(v[i] - mean, 0.0);
    const std::vector<cd> tsig = buf;
    detail::fftInPlace(buf, false);
    const long kLo = std::max(1L, static_cast<long>(12.0 * static_cast<double>(nft) /
                                                   static_cast<double>(n)));
    const long kHi = std::min<long>(static_cast<long>(0.625 * static_cast<double>(nft)),
                                    static_cast<long>(nft) - 2);
    if (kHi <= kLo + 2) return {0.0, 0.0};
    std::vector<double> mags;
    mags.reserve(static_cast<std::size_t>(kHi - kLo));
    long bestK = 0;
    double bestV = -1.0;
    for (long k = kLo; k < kHi; ++k) {
        const cd& f = buf[static_cast<std::size_t>(k)];
        const double vv = f.real() * f.real() + f.imag() * f.imag();
        mags.push_back(vv);
        if (vv > bestV) {
            bestV = vv;
            bestK = k;
        }
    }
    std::sort(mags.begin(), mags.end());
    const double med = mags[mags.size() / 2];
    const double prom = bestV / std::max(med, 1e-30);
    if (bestK <= kLo || bestK >= kHi - 1) return {0.0, prom};
    const double a0 = std::norm(buf[static_cast<std::size_t>(bestK - 1)]);
    const double a1 = std::norm(buf[static_cast<std::size_t>(bestK)]);
    const double a2 = std::norm(buf[static_cast<std::size_t>(bestK + 1)]);
    const double den = a0 - 2.0 * a1 + a2;
    double delta = 0.0;
    if (std::abs(den) > 1e-30) {
        delta = 0.5 * (a0 - a2) / den;
        if (std::abs(delta) > 1.0) delta = 0.0;
    }
    double fr = (static_cast<double>(bestK) + delta) / static_cast<double>(nft);
    // Schnelle zyklische Kohaerenz via inkrementeller Drehung (s. o.).
    auto coh = [&](double f) {
        cd acc(0.0, 0.0);
        const double cs = std::cos(2.0 * kPi * f);
        const double sn = -std::sin(2.0 * kPi * f);
        double cr = 1.0;
        double ci = 0.0;
        for (const cd& vv : tsig) {
            acc += vv * cd(cr, ci);
            const double nr = cr * cs - ci * sn;
            ci = cr * sn + ci * cs;
            cr = nr;
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

// Haertung: nicht-endliche Samples eines Slices als Null behandeln.
void sanitizeSlice(std::vector<cd>& y) {
    for (cd& v : y) {
        if (!std::isfinite(v.real()) || !std::isfinite(v.imag()))
            v = cd(0.0, 0.0);
    }
}

}  // namespace

// --------------------------------------------------------------------------
// Klassifikator-Fassade
// --------------------------------------------------------------------------
std::vector<PulseResult> ModulationClassifier::classify(
    const std::complex<double>* iq, std::size_t n) const {
    // Haertung: leere Eingabe und Nullzeiger sind definiert leer.
    if (iq == nullptr || n == 0) return {};
    // Erster Pipelineschritt: spektrale Kanalisierung (Kanaele um die
    // detektierten Energiebereiche), sofern in der Konfiguration aktiv.
    // Das Upsampling fuer sehr kurze Pulse sitzt nicht hier am Eingang,
    // sondern in processSlice auf dem extrahierten Basisband-Slice (dort
    // verliert der Anti-Image-Filter kein Nutzsignal, weil der Slice schmal
    // ist, und die Banddetektion bleibt ratenunabhaengig).
    if (cfg_.channelize) return classifyChannelized(cfg_, iq, n);
    return classifyLegacy(cfg_, iq, n);
}

std::vector<PulseResult> ModulationClassifier::classifyChannelized(
    const ClassifierConfig& cfg, const std::complex<double>* iq,
    std::size_t n) const {
    std::vector<PulseResult> results;
    // Haertung: leere Eingabe und Nullzeiger sind definiert leer.
    if (iq == nullptr || n == 0) return results;
    const BandExtractor bxs(cfg);
    // Erster Pipelineschritt: spektrale Kanalisierung auf dem ganzen Puffer.
    // Die Kanaele sind die detektierten Energiebereiche des Spektrums (nicht
    // gleichverteilt, sondern um die Stellen mit Energie gelegt).
    const BandExtractor::BandDetection det0 = bxs.bands(iq, 0, n, PulseSegmenter::noiseFloor(iq, n, cfg));
    std::vector<Band> chans = det0.bands;
    if (chans.empty())
        chans.push_back(Band{0, static_cast<long>(det0.nft)});
    // Zeitsegmente global ueber die Gesamtenergie des Puffers.
    const PulseSegmenter seg(cfg);
    const std::vector<std::pair<std::size_t, std::size_t>> segs =
        seg.detect(iq, n);
    // globaler Rauschboden als Gate-Referenz fuer die Detektion
    const double nf = PulseSegmenter::noiseFloor(iq, n, cfg);
    // Pro Segment: Detektion (mit korrektem Rauschboden des Originals) und
    // Zuordnung jedes Segment-Bandes zum ueberlappenden Kanal.
    struct SegDet {
        std::size_t a = 0;
        std::size_t b = 0;
        BandExtractor::BandDetection det;
        std::vector<Band> bands;
        std::vector<std::size_t> owner;
    };
    // Ab dieser Pufferlaenge lohnt sich die Parallelverarbeitung der
    // Segmente/Kanaele (darunter dominiert der Thread-Start).
    const unsigned thr =
        n >= 32768 ? detail::effectiveThreads(cfg.threads) : 1u;
    std::vector<SegDet> sd(segs.size());
    // Pro Segment: Detektion (mit korrektem Rauschboden des Originals) und
    // Zuordnung jedes Segment-Bandes zum ueberlappenden Kanal. Segmente sind
    // unabhaengig -> parallelisierbar.
    detail::parallelFor(
        segs.size(), thr,
        [&](std::size_t i) {
            sd[i].a = segs[i].first;
            sd[i].b = segs[i].second;
            sd[i].det = bxs.bands(iq, sd[i].a, sd[i].b, nf);
            sd[i].bands = sd[i].det.bands;
            if (sd[i].bands.empty())
                sd[i].bands.push_back(
                    Band{0, static_cast<long>(sd[i].det.nft)});
            // Kanalraster (det0) ins Segmentraster (detS) abbilden.
            const long sc = static_cast<long>(det0.nft) /
                            static_cast<long>(sd[i].det.nft);
            for (const Band& bS : sd[i].bands) {
                std::size_t owner = 0;
                long bestOv = -1;
                for (std::size_t j = 0; j < chans.size(); ++j) {
                    const Band cj{chans[j].lo / sc, chans[j].hi / sc};
                    const long ov = std::min(bS.hi, cj.hi) -
                                    std::max(bS.lo, cj.lo);
                    if (ov > bestOv) {
                        bestOv = ov;
                        owner = j;
                    }
                }
                sd[i].owner.push_back(owner);
            }
        });
    // Folgestufen je Kanal: Baender des Kanals werden pro Zeitsegment
    // extrahiert (Bandpass), auf die Bandmitte gemischt und optional
    // resamplet (Dezimierung mit Anti-Alias-Filter); danach Crops,
    // Merkmale und Stufe 1/2 auf dem Kanal-Basisband. Jedes (Segment, Band)
    // ist unabhaengig -> parallelisierbar; Ergebnisse landen pro Thread in
    // eigenen Buckets und werden danach deterministisch sortiert.
    const std::size_t Dreq =
        cfg.channelDecimate < 1
            ? 1
            : static_cast<std::size_t>(cfg.channelDecimate);
    struct Work {
        std::size_t seg;
        std::size_t band;
    };
    std::vector<Work> items;
    for (std::size_t si = 0; si < sd.size(); ++si)
        for (std::size_t bi = 0; bi < sd[si].bands.size(); ++bi)
            items.push_back(Work{si, bi});
    std::vector<std::vector<PulseResult>> buckets(thr);
    {
        std::atomic<std::size_t> next{0};
        const auto worker = [&](std::size_t tid) {
            for (;;) {
                const std::size_t k = next.fetch_add(1);
                if (k >= items.size()) return;
                const Work& w = items[k];
                const SegDet& s = sd[w.seg];
                const Band& bS = s.bands[w.band];
                double snrDb = 0.0;
                double flatDb = 0.0;
                estimateSnrFlat(s.det.psdDb, bS, cfg, s.bands, &snrDb,
                                &flatDb);
                BandExtractor::Slice slice =
                    bxs.extract(iq, s.a, s.b, bS, s.det.nft);
                // Resample-Faktor des Kanals (Nyquist der Bandbreite).
                std::size_t D = 1;
                if (Dreq > 1 && slice.bw > 0.0) {
                    const double maxD = std::floor(1.0 / slice.bw);
                    D = maxD >= 1.0 ? static_cast<std::size_t>(maxD) : 1;
                    if (D > Dreq) D = Dreq;
                }
                std::vector<cd> y = slice.y;
                if (D > 1) y = Channelizer::resample(y, D);
                const ClassifierConfig c2 = scaleForDecimation(cfg, D);
                processSlice(c2, y, s.a / D, slice.fc,
                             slice.bw * static_cast<double>(D),
                             static_cast<std::size_t>(bS.width()), snrDb,
                             flatDb, s.det.psdDb, bS, s.det.nft, D,
                             buckets[tid]);
            }
        };
        std::vector<std::thread> pool;
        pool.reserve(thr > 0 ? thr - 1 : 0);
        for (unsigned t = 1; t < thr; ++t)
            pool.emplace_back(worker, t);
        worker(0);
        for (auto& th : pool) th.join();
    }
    for (auto& b : buckets)
        for (auto& r : b) results.push_back(std::move(r));
    // Deterministische Reihenfolge unabhaengig von der Threadanzahl:
    // feine Zeitlage, dann Frequenz, dann Laenge.
    std::stable_sort(results.begin(), results.end(),
                     [](const PulseResult& a, const PulseResult& b) {
                         if (a.startSampleFine != b.startSampleFine)
                             return a.startSampleFine < b.startSampleFine;
                         if (a.centerFreqRefined != b.centerFreqRefined)
                             return a.centerFreqRefined < b.centerFreqRefined;
                         return a.numSamplesFine < b.numSamplesFine;
                     });
    return results;
}

void ModulationClassifier::processSlice(
    const ClassifierConfig& cfg, const std::vector<cd>& sliceY,
    std::size_t segA, double fcGlobal, double bwLocal, std::size_t bandBins,
    double snrDb, double flatDb, const std::vector<double>& psdDb,
    const Band& band, std::size_t detNft, std::size_t timeMul,
    std::vector<PulseResult>& results) const {
    // Upsampling des Basisband-Slices fuer sehr kurze Pulse (z. B. ~180
    // Samples Puffer mit 13 BPSK-Bit bei sps=8): automatisch, wenn der Slice
    // kuerzer als 1024 Samples ist (U = clamp(ceil(1024/len),1,8)), sonst
    // gemaess cfg.upsample (1 = aus, >1 = fester Faktor). Hier (und nicht am
    // Puffer-Eingang) ist der Anti-Image-Filter unkritisch, weil der Slice
    // schmalbandig ist, und die Banddetektion bleibt in der Originalrate.
    // Sample-basierte Schwellwerte werden mit U skaliert; die normierte
    // Bandbreite in Slice-Rate wird /U gerechnet, Bin-Groessen xU, und die
    // Ergebniszeiten am Ende /U auf die Originalrate zurueckgerechnet.
    std::size_t U = 1;
    if (cfg.upsample > 1) {
        U = static_cast<std::size_t>(cfg.upsample);
    } else if (cfg.upsample == 0 && sliceY.size() < 1024) {
        const double need = 1024.0 / static_cast<double>(sliceY.size());
        const double u = std::ceil(need);
        U = u < 1.0 ? 1 : (u > 8.0 ? 8 : static_cast<std::size_t>(u));
    }
    ClassifierConfig cfgLoc = cfg;
    std::vector<cd> up;
    if (U > 1) {
        cfgLoc = scaleForUpsample(cfg, U);
        up = Channelizer::upsample(sliceY, U);
    }
    const ClassifierConfig& cfgA = cfgLoc;
    const std::vector<cd>& yAll = U > 1 ? up : sliceY;
    const double bwFeat = bwLocal / static_cast<double>(U);
    const std::size_t bandBinsFeat = bandBins * U;
    const double tScale = static_cast<double>(timeMul) / static_cast<double>(U);
    for (const auto& crop : timeCrops(yAll, cfgA)) {
        const std::size_t c0 = crop.first;
        const std::size_t c1 = crop.second;
        const std::size_t len = c1 - c0;
        if (len < cfgA.minPulseLen) continue;
        std::vector<cd> ys(yAll.begin() + static_cast<std::ptrdiff_t>(c0),
                           yAll.begin() + static_cast<std::ptrdiff_t>(c1));
        // Haertung: nicht-endliche Samples als Null behandeln.
        // NaN/Inf duerfen weder Statistiken vergiften noch in
        // Indexumrechnungen einfliessen; ein Dropout verhaelt
        // sich damit wie Stille. (Konsistent mit dem Oracle.)
        double amax = 0.0;
        for (cd& v : ys) {
            if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) {
                v = cd(0.0, 0.0);
            } else {
                const double a = std::abs(v);
                if (a > amax) amax = a;
            }
        }
        // Haertung: Extremamplituden normalisieren (Klassifikation ist
        // skaleninvariant; sonst wuerden Betragsquadrate/-potenzen in den
        // Statistiken ueberlaufen).
        if (amax > 1e100)
            for (cd& v : ys) v /= amax;
        const Features f =
            extractFeatures(ys, 0.0, bwFeat, bandBinsFeat, snrDb, flatDb, cfgA);
        CoarseClass coarse = CoarseClass::Unknown;
        if (!stage1(f, cfgA, &coarse)) continue;
        const auto res = stage2(f, cfgA, coarse, ys);
        // Feine Zeitlage (samplegenau): Huellkurven-Schwelle im
        // Crop, relativ zum Crop-Maximum. Die Raender des Crops
        // (Sicherheitsrand, Blockraster) fallen darunter weg.
        std::size_t i0 = 0;
        std::size_t i1 = ys.size();
        {
            double emax = 0.0;
            for (const cd& v : ys) emax = std::max(emax, std::abs(v));
            if (emax > 0.0 && ys.size() > 2) {
                const double thr = 0.25 * emax;
                std::size_t a0 = 0;
                while (a0 < ys.size() && std::abs(ys[a0]) <= thr) ++a0;
                std::size_t a1 = ys.size();
                while (a1 > a0 && std::abs(ys[a1 - 1]) <= thr) --a1;
                if (a1 > a0) {
                    i0 = a0;
                    i1 = a1;
                }
            }
        }
        // Sub-Bin-Mittenfrequenz: leistungsgewichteter Schwerpunkt des
        // In-Band-PSD (zyklische Bins ab band.lo), globale Koordinaten.
        double fcRef = wrap01(fcGlobal);
        {
            const long nftL = static_cast<long>(detNft);
            if (nftL > 0 && !psdDb.empty() && band.width() > 0) {
                double num = 0.0;
                double den = 0.0;
                for (long k = band.lo; k < band.hi; ++k) {
                    const long kk = ((k % nftL) + nftL) % nftL;
                    const double p = std::pow(
                        10.0, psdDb[static_cast<std::size_t>(kk)] / 10.0);
                    num += static_cast<double>(k - band.lo) * p;
                    den += p;
                }
                if (den > 0.0) {
                    fcRef = wrap01((static_cast<double>(band.lo) + num / den) /
                                   static_cast<double>(detNft));
                }
            }
        }
        PulseResult r;
        const double base = static_cast<double>(segA) * timeMul;
        r.startSample = static_cast<std::size_t>(
            std::lround(base + static_cast<double>(c0) * tScale));
        r.numSamples = static_cast<std::size_t>(std::lround(len * tScale));
        r.startSampleFine = static_cast<std::size_t>(
            std::lround(base + static_cast<double>(c0 + i0) * tScale));
        r.numSamplesFine =
            static_cast<std::size_t>(std::lround((i1 - i0) * tScale));
        r.centerFreqRefined = fcRef;
        r.centerFreqNorm =
            static_cast<float>(roundTo(wrap01(fcGlobal), 5));
        r.bandwidthNorm = static_cast<float>(
            roundTo(bwLocal / static_cast<double>(timeMul), 5));
        r.snrDb = static_cast<float>(roundTo(f.snrDb, 1));
        r.coarse = coarse;
        r.modulation = res.first;
        r.confidence = static_cast<float>(roundTo(res.second, 2));
        r.envelopeCv = static_cast<float>(roundTo(f.cv, 3));
        r.paprWindowedDb = static_cast<float>(roundTo(f.paprWDb, 1));
        r.chirpR2 = static_cast<float>(roundTo(f.chirpR2, 2));
        r.chirpSweep = static_cast<float>(f.chirpSweep);
        r.dwellFraction = static_cast<float>(roundTo(f.dwell, 2));
        r.freqModes = f.nfreqModes;
        r.spectralFlatnessDb = static_cast<float>(roundTo(f.flatDb, 2));
        r.bandBins = bandBins / timeMul;
        results.push_back(r);
    }
}

std::vector<PulseResult> ModulationClassifier::classifyLegacy(
    const ClassifierConfig& cfg, const std::complex<double>* iq,
    std::size_t n) const {
    std::vector<PulseResult> results;
    // Haertung: leere Eingabe und Nullzeiger sind definiert leer.
    if (iq == nullptr || n == 0) return results;
    const PulseSegmenter seg(cfg);
    const BandExtractor bxs(cfg);
    const double nf = PulseSegmenter::noiseFloor(iq, n, cfg);
    for (const auto& ab : seg.detect(iq, n)) {
        const std::size_t a = ab.first;
        const std::size_t b = ab.second;
        BandExtractor::BandDetection det = bxs.bands(iq, a, b, nf);
        std::vector<Band> bands = det.bands;
        if (bands.empty()) {
            bands.push_back(Band{0, static_cast<long>(det.nft)});
        }
        for (const Band& band : bands) {
            const BandExtractor::Slice slice =
                bxs.extract(iq, a, b, band, det.nft);
            double snrDb = 0.0;
            double flatDb = 0.0;
            estimateSnrFlat(det.psdDb, band, cfg, bands, &snrDb, &flatDb);
            processSlice(cfg, slice.y, a, slice.fc, slice.bw,
                         static_cast<std::size_t>(band.width()), snrDb, flatDb,
                         det.psdDb, band, det.nft, 1, results);
        }
    }
    return results;
}

std::vector<PulseResult> ModulationClassifier::classifyF(
    const std::complex<float>* iq, std::size_t n) const {
    // Haertung: Nullzeiger und Laenge 0 sind definiert leer (wie classify).
    if (iq == nullptr || n == 0) return {};
    std::vector<std::complex<double>> buf(n);
    for (std::size_t i = 0; i < n; ++i) {
        buf[i] = std::complex<double>(static_cast<double>(iq[i].real()),
                                      static_cast<double>(iq[i].imag()));
    }
    return classify(buf.data(), buf.size());
}

std::vector<std::complex<double>> ModulationClassifier::extractPulse(
    const std::complex<double>* iq, std::size_t n,
    const PulseResult& pulse) const {
    std::vector<std::complex<double>> out;
    // Haertung: Nullzeiger/leere Eingabe sind definiert leer.
    if (iq == nullptr || n == 0) return out;
    std::size_t s = pulse.startSampleFine;
    std::size_t len = pulse.numSamplesFine;
    if (len == 0 || s >= n) {
        // Fallback auf das grobe Fenster, falls das feine fehlt.
        s = pulse.startSample;
        len = pulse.numSamples;
    }
    if (len == 0 || s >= n) return out;
    len = std::min(len, n - s);
    std::size_t nft = detail::pow2ceil(len);
    if (nft < 16) nft = 16;
    double fc = std::fmod(pulse.centerFreqRefined, 1.0);
    if (fc < 0.0) fc += 1.0;
    double bw = static_cast<double>(pulse.bandwidthNorm);
    if (!(bw > 0.0)) bw = 1.0 / static_cast<double>(nft);
    // Bandmaske um die feine Mittenfrequenz (zyklische Bin-Koordinaten,
    // BandExtractor::extract fuegt den Waechter selbst hinzu).
    const long lo = static_cast<long>(std::lround(
        (fc - 0.5 * bw) * static_cast<double>(nft)));
    long hi = static_cast<long>(std::lround(
        (fc + 0.5 * bw) * static_cast<double>(nft)));
    if (hi - lo < 1) hi = lo + 1;
    const BandExtractor bxs(cfg_);
    const BandExtractor::Slice slice =
        bxs.extract(iq, s, s + len, Band{lo, hi}, nft);
    return slice.y;
}

std::vector<std::complex<float>> ModulationClassifier::extractPulseF(
    const std::complex<float>* iq, std::size_t n,
    const PulseResult& pulse) const {
    std::vector<std::complex<float>> out;
    if (iq == nullptr || n == 0) return out;
    std::vector<std::complex<double>> buf(n);
    for (std::size_t i = 0; i < n; ++i) {
        buf[i] = std::complex<double>(static_cast<double>(iq[i].real()),
                                      static_cast<double>(iq[i].imag()));
    }
    const std::vector<std::complex<double>> s =
        extractPulse(buf.data(), buf.size(), pulse);
    out.reserve(s.size());
    for (const auto& v : s) {
        out.push_back(std::complex<float>(static_cast<float>(v.real()),
                                          static_cast<float>(v.imag())));
    }
    return out;
}

// ---------------------------------------------------------------------------
// PSK-Synchronisation (batchweise, ohne Einschwingvorgang)
// ---------------------------------------------------------------------------
PskSync ModulationClassifier::synchronizePsk(const std::complex<double>* iq,
                                             std::size_t n,
                                             const PulseResult& pulse) const {
    PskSync r;
    if (iq == nullptr || n == 0) return r;
    if (pulse.modulation != ModulationType::BPSK &&
        pulse.modulation != ModulationType::QPSK)
        return r;
    std::vector<cd> y = extractPulse(iq, n, pulse);
    sanitizeSlice(y);
    if (y.size() < cfg_.minAnalyzeLen) return r;
    {   // Haertung: Extremamplituden normalisieren (skaleninvariant).
        double amax = 0.0;
        for (const cd& v : y) amax = std::max(amax, std::abs(v));
        if (amax > 1e100)
            for (cd& v : y) v /= amax;
    }
    const Features f =
        extractFeatures(y, pulse.centerFreqRefined, pulse.bandwidthNorm,
                        pulse.bandBins, pulse.snrDb, pulse.spectralFlatnessDb,
                        cfg_);
    const std::vector<std::size_t> idx = estimateSymbolsRobust(f, cfg_);
    if (idx.size() < cfg_.timingMinSyms) return r;
    std::vector<cd> s0 = mfSymbols(y, idx, cfg_);
    if (s0.size() < cfg_.timingMinSyms) {
        s0.clear();
        for (std::size_t i : idx)
            s0.push_back(cd(f.env[i] * std::cos(f.ph[i]),
                            f.env[i] * std::sin(f.ph[i])));
    }
    // Resttraeger pro Symbol aus dem z^M-Spektrum (batchweise).
    const auto fp = carrierFromSymbols(s0);
    const double fresSym = fp.first;
    // Symboltakt + Auge: batchweise ueber die zyklische Leistungslinie.
    std::vector<double> env2(f.n);
    for (std::size_t i = 0; i < f.n; ++i) env2[i] = f.env[i] * f.env[i];
    const auto rp = rsFromCyclic(env2, f.n, f.bw, 0, cfg_.rsBwDivLo, cfg_.rsBwDivHi, cfg_.tRsProm);
    if (rp.first <= 0.0) return r;
    const double sps = 1.0 / rp.first;
    const double off = eyeOffset(f.env, f.n, sps);
    // Symbole derotieren und normieren.
    double sp = 0.0;
    for (const cd& v : s0) sp += std::norm(v);
    sp /= static_cast<double>(s0.size());
    const double sc = 1.0 / std::sqrt(std::max(sp, 1e-30));
    const std::size_t len = s0.size();
    std::vector<cd> s(len);
    for (std::size_t k = 0; k < len; ++k) {
        const double ph = -2.0 * kPi * fresSym * static_cast<double>(k);
        s[k] = s0[k] * sc * cd(std::cos(ph), std::sin(ph));
    }
    // M-PSK-Restambiguität 2pi/M aufloesen: Kandidat mit minimaler Distanz
    // der Konstellation zum festen Gitter {2*pi*j/M}.
    const int M = (pulse.modulation == ModulationType::BPSK) ? 2 : 4;
    cd acc(0.0, 0.0);
    for (const cd& v : s) {
        cd p = v;
        for (int i = 1; i < M; ++i) p *= v;
        acc += p;
    }
    const double psi = std::atan2(acc.imag(), acc.real());
    double bestScore = -1e300;
    double phiC = 0.0;
    for (int q = 0; q < M; ++q) {
        const double phc = (psi + 2.0 * kPi * q / M) / M;
        double score = 0.0;
        for (const cd& v : s) {
            double a = std::atan2(v.imag(), v.real()) - phc;
            a = std::fmod(a, 2.0 * kPi);
            if (a < 0.0) a += 2.0 * kPi;
            if (a >= kPi) a -= 2.0 * kPi;
            double d = 1e300;
            for (int j = 0; j < M; ++j)
                d = std::min(d, std::abs(a - 2.0 * kPi * j / M));
            if (d > kPi) d = 2.0 * kPi - d;
            score -= d;
        }
        if (score > bestScore) {
            bestScore = score;
            phiC = phc;
        }
    }
    for (cd& v : s) v *= cd(std::cos(-phiC), std::sin(-phiC));
    // Traeger pro Sample und absolute Traegerfrequenz; Phase am Pulsanfang:
    // y[0] (Slice-Index 0 == startSampleFine) traegt phiC + 2*pi*fRes*off.
    const double fRes = fresSym / sps;
    r.valid = true;
    r.sps = sps;
    r.timingOffset = off;
    r.carrierFreqNorm = wrap01(pulse.centerFreqRefined + fRes);
    double phiStart = phiC - 2.0 * kPi * fRes * off;
    phiStart = std::fmod(phiStart, 2.0 * kPi);
    if (phiStart > kPi) phiStart -= 2.0 * kPi;
    if (phiStart < -kPi) phiStart += 2.0 * kPi;
    r.carrierPhaseAtStart = phiStart;
    r.snrDb = pulse.snrDb;
    r.symbols = std::move(s);
    return r;
}

PskSync ModulationClassifier::synchronizePskF(const std::complex<float>* iq,
                                              std::size_t n,
                                              const PulseResult& pulse) const {
    if (iq == nullptr || n == 0) return PskSync{};
    std::vector<std::complex<double>> buf(n);
    for (std::size_t i = 0; i < n; ++i)
        buf[i] = std::complex<double>(static_cast<double>(iq[i].real()),
                                      static_cast<double>(iq[i].imag()));
    return synchronizePsk(buf.data(), buf.size(), pulse);
}

// ---------------------------------------------------------------------------
// Chirp-Parameter
// ---------------------------------------------------------------------------
ChirpParams ModulationClassifier::analyzeChirp(
    const std::complex<double>* iq, std::size_t n,
    const PulseResult& pulse) const {
    ChirpParams r;
    if (iq == nullptr || n == 0) return r;
    if (pulse.modulation != ModulationType::Chirp) return r;
    std::vector<cd> y = extractPulse(iq, n, pulse);
    sanitizeSlice(y);
    if (y.size() < cfg_.minAnalyzeLen) return r;
    {   // Haertung: Extremamplituden normalisieren (skaleninvariant).
        double amax = 0.0;
        for (const cd& v : y) amax = std::max(amax, std::abs(v));
        if (amax > 1e100)
            for (cd& v : y) v /= amax;
    }
    const Features f =
        extractFeatures(y, pulse.centerFreqRefined, pulse.bandwidthNorm,
                        pulse.bandBins, pulse.snrDb, pulse.spectralFlatnessDb,
                        cfg_);
    // Gewichtete Regression auf f_inst (Resttraeger eingerechnet, da der
    // Achsenabschnitt die absolute Frequenz am Pulsanfang liefern soll),
    // mittlere 80 % des Pulses - analog zum Chirp-Fit in extractFeatures.
    const std::size_t n1 = f.n >= 1 ? f.n - 1 : 0;
    const std::size_t iLo = n1 / 10;
    const std::size_t iHi = n1 - n1 / 10;
    const std::size_t iHiEff = std::min(iHi, n1);
    double sw = 0.0, xm = 0.0, ym = 0.0;
    for (std::size_t i = iLo; i < iHiEff; ++i) {
        const double w = f.env[i] * f.env[i];
        sw += w;
        xm += w * static_cast<double>(i);
        ym += w * (f.finst[i] + f.residF);
    }
    xm /= sw > 1e-30 ? sw : 1e-30;
    ym /= sw > 1e-30 ? sw : 1e-30;
    double num = 0.0, den = 0.0;
    for (std::size_t i = iLo; i < iHiEff; ++i) {
        const double w = f.env[i] * f.env[i];
        const double dx = static_cast<double>(i) - xm;
        num += w * dx * (f.finst[i] + f.residF - ym);
        den += w * dx * dx;
    }
    const double slope = num / (den > 1e-30 ? den : 1e-30);
    const double intercept = ym - slope * xm;  // f_inst am Slice-Anfang
    const double fStartRel = intercept;
    const double fEndRel = intercept + slope * static_cast<double>(n1 - 1);
    r.valid = f.chirpR2 > cfg_.tChirpR2;
    r.startSample = pulse.startSampleFine;
    r.numSamples = pulse.numSamplesFine;
    r.fStartNorm = wrap01(pulse.centerFreqRefined + fStartRel);
    r.fEndNorm = wrap01(pulse.centerFreqRefined + fEndRel);
    r.centerFreqNorm = pulse.centerFreqRefined;
    r.bandwidthNorm = std::abs(fEndRel - fStartRel);
    r.slopeNorm = slope;
    r.snrDb = pulse.snrDb;
    r.r2 = f.chirpR2;
    r.up = slope > 0.0;
    return r;
}

ChirpParams ModulationClassifier::analyzeChirpF(
    const std::complex<float>* iq, std::size_t n,
    const PulseResult& pulse) const {
    if (iq == nullptr || n == 0) return ChirpParams{};
    std::vector<std::complex<double>> buf(n);
    for (std::size_t i = 0; i < n; ++i)
        buf[i] = std::complex<double>(static_cast<double>(iq[i].real()),
                                      static_cast<double>(iq[i].imag()));
    return analyzeChirp(buf.data(), buf.size(), pulse);
}

FrankParams ModulationClassifier::analyzeFrank(
    const std::complex<double>* iq, std::size_t n,
    const PulseResult& pulse) const {
    FrankParams r;
    if (iq == nullptr || n == 0) return r;
    if (pulse.modulation != ModulationType::FrankCode) return r;
    std::vector<cd> y = extractPulse(iq, n, pulse);
    sanitizeSlice(y);
    if (y.size() < cfg_.minAnalyzeLen) return r;
    {   // Haertung: Extremamplituden normalisieren (skaleninvariant).
        double amax = 0.0;
        for (const cd& v : y) amax = std::max(amax, std::abs(v));
        if (amax > 1e100)
            for (cd& v : y) v /= amax;
    }
    // gleiches Slice-Upsampling wie in processSlice (kurze Pulse)
    std::size_t U = 1;
    if (cfg_.upsample > 1) {
        U = static_cast<std::size_t>(cfg_.upsample);
    } else if (cfg_.upsample == 0 && y.size() < 1024) {
        const double need = 1024.0 / static_cast<double>(y.size());
        const double u = std::ceil(need);
        U = u < 1.0 ? 1 : (u > 8.0 ? 8 : static_cast<std::size_t>(u));
    }
    ClassifierConfig cfgA = cfg_;
    if (U > 1) {
        cfgA = scaleForUpsample(cfg_, U);
        y = Channelizer::upsample(y, U);
    }
    const double bwFeat = pulse.bandwidthNorm / static_cast<double>(U);
    const Features f =
        extractFeatures(y, pulse.centerFreqRefined, bwFeat,
                        pulse.bandBins * U, pulse.snrDb, pulse.spectralFlatnessDb,
                        cfgA);
    double r2 = 0.0;
    const int M = detectFrankOrder(f, cfgA, &r2);
    if (M < cfg_.tFrankMinModes) return r;
    r.valid = true;
    r.order = M;
    r.codeLength = static_cast<std::size_t>(M) * M;
    r.startSample = pulse.startSampleFine;
    r.numSamples = pulse.numSamplesFine;
    r.centerFreqNorm = pulse.centerFreqRefined;
    r.bandwidthNorm = pulse.bandwidthNorm;
    r.snrDb = pulse.snrDb;
    r.r2 = r2;
    // Chip-Rate aus der Pulslaenge: N = M^2 Chips im Slice (Originalrate).
    const double chipLen = static_cast<double>(y.size()) /
                           static_cast<double>(r.codeLength) /
                           static_cast<double>(U);
    r.chipRateNorm = chipLen > 1e-9 ? 1.0 / chipLen : 0.0;
    // Code-Phasen an den Chip-Mitten abtasten (auf den Resttraeger
    // f.residF derotiert, zyklisch [0, 2pi)).
    r.codePhases.resize(r.codeLength);
    const double phStep = 2.0 * kPi * f.residF;
    const double chipLenUp = chipLen * static_cast<double>(U);
    for (std::size_t k = 0; k < r.codeLength; ++k) {
        const double c = (static_cast<double>(k) + 0.5) * chipLenUp;
        const std::size_t ci = static_cast<std::size_t>(c);
        if (ci >= y.size()) {
            r.codePhases[k] = 0.0;
            continue;
        }
        const double ph = std::atan2(y[ci].imag(), y[ci].real()) -
                          phStep * static_cast<double>(ci);
        double w = std::fmod(ph, 2.0 * kPi);
        if (w < 0.0) w += 2.0 * kPi;
        r.codePhases[k] = w;
    }
    return r;
}

FrankParams ModulationClassifier::analyzeFrankF(
    const std::complex<float>* iq, std::size_t n,
    const PulseResult& pulse) const {
    if (iq == nullptr || n == 0) return FrankParams{};
    std::vector<std::complex<double>> buf(n);
    for (std::size_t i = 0; i < n; ++i)
        buf[i] = std::complex<double>(static_cast<double>(iq[i].real()),
                                      static_cast<double>(iq[i].imag()));
    return analyzeFrank(buf.data(), buf.size(), pulse);
}


// ---------------------------------------------------------------------------
// P-Codes (P1..P4, polyphase Radar-Codes anknuepfend an Frank)
// ---------------------------------------------------------------------------
// P1/P2 haben M x M Chips mit den Phasen-Gesetzen in polyLawPhase (ihre
// f_inst-Treppe ist symmetrisch und faellt im entfalteten Frame); P3/P4
// sind lineare Chirps der Laenge M (werden als Chirp klassifiziert).
PcodeParams ModulationClassifier::analyzePcode(
    const std::complex<double>* iq, std::size_t n,
    const PulseResult& pulse) const {
    PcodeParams r;
    if (iq == nullptr || n == 0) return r;
    if (pulse.modulation != ModulationType::PCode) return r;
    std::vector<cd> y = extractPulse(iq, n, pulse);
    sanitizeSlice(y);
    if (y.size() < cfg_.minAnalyzeLen) return r;
    {
        double amax = 0.0;
        for (const cd& v : y) amax = std::max(amax, std::abs(v));
        if (amax > 1e100)
            for (cd& v : y) v /= amax;
    }
    // gleiches Slice-Upsampling wie in processSlice
    std::size_t U = 1;
    if (cfg_.upsample > 1) {
        U = static_cast<std::size_t>(cfg_.upsample);
    } else if (cfg_.upsample == 0 && y.size() < 1024) {
        const double need = 1024.0 / static_cast<double>(y.size());
        const double u = std::ceil(need);
        U = u < 1.0 ? 1 : (u > 8.0 ? 8 : static_cast<std::size_t>(u));
    }
    ClassifierConfig cfgA = cfg_;
    if (U > 1) {
        cfgA = scaleForUpsample(cfg_, U);
        y = Channelizer::upsample(y, U);
    }
    const double bwFeat = pulse.bandwidthNorm / static_cast<double>(U);
    const Features f =
        extractFeatures(y, pulse.centerFreqRefined, bwFeat,
                        pulse.bandBins * U, pulse.snrDb, pulse.spectralFlatnessDb,
                        cfgA);
    double r2 = 0.0;
    double platLen = 0.0;
    const int M0 = detectFrankOrder(f, cfgA, &r2, &platLen);
    int Mb = M0 > 0 ? M0 : 8;
    double sc = 0.0, fsc = 0.0;
    const int law = polyLawMf(y, Mb, platLen, &Mb, &sc, &fsc);
    if (law == 0 || Mb < cfg_.tFrankMinModes) return r;  // ist Frank, nicht P
    int MbL = Mb;
    double lerr = 0.0;
    polyLawMatch(y, f, Mb, platLen, &MbL, &lerr);
    r.valid = true;
    r.subtype = law;  // 1 = P1, 2 = P2
    r.order = Mb;
    r.codeLength = static_cast<std::size_t>(Mb) * Mb;
    r.startSample = pulse.startSampleFine;
    r.numSamples = pulse.numSamplesFine;
    r.centerFreqNorm = pulse.centerFreqRefined;
    r.bandwidthNorm = pulse.bandwidthNorm;
    r.snrDb = pulse.snrDb;
    r.r2 = r2;
    r.lawErr = lerr;
    const double chipLen = static_cast<double>(y.size()) /
                           static_cast<double>(r.codeLength) /
                           static_cast<double>(U);
    r.chipRateNorm = chipLen > 1e-9 ? 1.0 / chipLen : 0.0;
    // Code-Phasen an den Chip-Mitten (Resttraeger derotiert)
    r.codePhases.resize(r.codeLength);
    const double phStep = 2.0 * kPi * f.residF;
    const double chipLenUp = chipLen * static_cast<double>(U);
    for (std::size_t k = 0; k < r.codeLength; ++k) {
        const std::size_t ci = static_cast<std::size_t>(
            (static_cast<double>(k) + 0.5) * chipLenUp);
        if (ci >= y.size()) {
            r.codePhases[k] = 0.0;
            continue;
        }
        double w = std::atan2(y[ci].imag(), y[ci].real()) -
                   phStep * static_cast<double>(ci);
        w = std::fmod(w, 2.0 * kPi);
        if (w < 0.0) w += 2.0 * kPi;
        r.codePhases[k] = w;
    }
    return r;
}

PcodeParams ModulationClassifier::analyzePcodeF(
    const std::complex<float>* iq, std::size_t n,
    const PulseResult& pulse) const {
    if (iq == nullptr || n == 0) return PcodeParams{};
    std::vector<std::complex<double>> buf(n);
    for (std::size_t i = 0; i < n; ++i)
        buf[i] = std::complex<double>(static_cast<double>(iq[i].real()),
                                      static_cast<double>(iq[i].imag()));
    return analyzePcode(buf.data(), buf.size(), pulse);
}

// ---------------------------------------------------------------------------
// Barker-Code (BPSK-Puls mit Barker-Folge)
// ---------------------------------------------------------------------------
BarkerParams ModulationClassifier::analyzeBarker(
    const std::complex<double>* iq, std::size_t n,
    const PulseResult& pulse) const {
    BarkerParams r;
    if (iq == nullptr || n == 0) return r;
    if (pulse.modulation != ModulationType::BarkerCode) return r;
    std::vector<cd> y = extractPulse(iq, n, pulse);
    sanitizeSlice(y);
    if (y.size() < cfg_.minAnalyzeLen) return r;
    {
        double amax = 0.0;
        for (const cd& v : y) amax = std::max(amax, std::abs(v));
        if (amax > 1e100)
            for (cd& v : y) v /= amax;
    }
    // gleiches Slice-Upsampling wie in processSlice
    std::size_t U = 1;
    if (cfg_.upsample > 1) {
        U = static_cast<std::size_t>(cfg_.upsample);
    } else if (cfg_.upsample == 0 && y.size() < 1024) {
        const double need = 1024.0 / static_cast<double>(y.size());
        const double u = std::ceil(need);
        U = u < 1.0 ? 1 : (u > 8.0 ? 8 : static_cast<std::size_t>(u));
    }
    ClassifierConfig cfgA = cfg_;
    if (U > 1) {
        cfgA = scaleForUpsample(cfg_, U);
        y = Channelizer::upsample(y, U);
    }
    const double bwFeat = pulse.bandwidthNorm / static_cast<double>(U);
    const Features f =
        extractFeatures(y, pulse.centerFreqRefined, bwFeat,
                        pulse.bandBins * U, pulse.snrDb, pulse.spectralFlatnessDb,
                        cfgA);
    const std::vector<std::size_t> idx = estimateSymbolsRobust(f, cfgA);
    const double spsEst = idx.size() >= 2
        ? static_cast<double>(idx.back() - idx.front()) /
              static_cast<double>(idx.size() - 1)
        : 8.0;
    int reps = 0;
    double sc = 0.0, spsB = spsEst;
    const int bl = barkerMatchSlice(y, f, cfgA, spsEst, &reps, &sc, &spsB);
    if (bl <= 0) return r;
    static const int b13[] = {1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1};
    static const int b11[] = {1, 1, 1, -1, -1, -1, 1, -1, -1, 1, -1};
    static const int b7[] = {1, 1, 1, -1, -1, 1, -1};
    const int* seq = bl == 13 ? b13 : (bl == 11 ? b11 : b7);
    r.valid = true;
    r.length = bl;
    r.code.assign(seq, seq + bl);
    r.repeats = reps;
    r.sps = spsB / static_cast<double>(U);
    r.carrierFreqNorm = pulse.centerFreqRefined;
    r.carrierPhaseAtStart = 0.0;
    r.snrDb = pulse.snrDb;
    return r;
}

BarkerParams ModulationClassifier::analyzeBarkerF(
    const std::complex<float>* iq, std::size_t n,
    const PulseResult& pulse) const {
    if (iq == nullptr || n == 0) return BarkerParams{};
    std::vector<std::complex<double>> buf(n);
    for (std::size_t i = 0; i < n; ++i)
        buf[i] = std::complex<double>(static_cast<double>(iq[i].real()),
                                      static_cast<double>(iq[i].imag()));
    return analyzeBarker(buf.data(), buf.size(), pulse);
}

// ---------------------------------------------------------------------------
// Symbol-Schnittstelle (PSK + FSK)
// ---------------------------------------------------------------------------
SymbolStream ModulationClassifier::extractSymbols(
    const std::complex<double>* iq, std::size_t n,
    const PulseResult& pulse) const {
    SymbolStream r;
    if (iq == nullptr || n == 0) return r;
    if (pulse.modulation == ModulationType::BPSK ||
        pulse.modulation == ModulationType::QPSK ||
        pulse.modulation == ModulationType::BarkerCode) {
        PskSync s = synchronizePsk(iq, n, pulse);
        if (!s.valid) return r;
        r.valid = true;
        r.type = pulse.modulation;
        r.sps = s.sps;
        r.timingOffset = s.timingOffset;
        r.carrierFreqNorm = s.carrierFreqNorm;
        r.carrierPhaseAtStart = s.carrierPhaseAtStart;
        r.symbols = std::move(s.symbols);
        return r;
    }
    if (pulse.modulation == ModulationType::FSK2) {
        std::vector<cd> y = extractPulse(iq, n, pulse);
        sanitizeSlice(y);
        if (y.size() < cfg_.minAnalyzeLen) return r;
        const Features f =
            extractFeatures(y, pulse.centerFreqRefined, pulse.bandwidthNorm,
                            pulse.bandBins, pulse.snrDb,
                            pulse.spectralFlatnessDb, cfg_);
        // FSK hat konstante Huellkurve: Symboltakt aus den f_inst-
        // Umschaltspitzen (Vollsuche der Ratenlinie in |Delta f_inst|).
        // f.finst hat f.n-1 Eintraege, daher n1 = f.n-2 Difference-Werte.
        const std::size_t n1 = f.n >= 2 ? f.n - 2 : 0;
        if (n1 < 64) return r;
        std::vector<double> tr(n1);
        for (std::size_t i = 0; i < n1; ++i)
            tr[i] = std::abs(f.finst[i + 1] - f.finst[i]);
        const auto rp = rsLineFull(tr);
        if (rp.first <= 0.0 || rp.second < cfg_.tRsProm) return r;
        const double sps = 1.0 / rp.first;
        if (sps < 1.6 || sps > std::max(2.0, static_cast<double>(f.n) / 12.0))
            return r;
        // Timing-Offset: Gitter so legen, dass es die Umschaltspitzen trifft.
        double bestOff = 0.0;
        double bestScore = -1.0;
        for (double off = 0.0; off < sps; off += 0.25) {
            double sc = 0.0;
            long c = 0;
            for (double p = off; p < static_cast<double>(n1) - 1.0;
                 p += sps) {
                const std::size_t i = static_cast<std::size_t>(std::lround(p));
                if (i < n1) {
                    sc += tr[i];
                    ++c;
                }
            }
            if (c > 0 && sc / static_cast<double>(c) > bestScore) {
                bestScore = sc / static_cast<double>(c);
                bestOff = off;
            }
        }
        r.valid = true;
        r.type = ModulationType::FSK2;
        r.sps = sps;
        r.timingOffset = bestOff;
        r.carrierFreqNorm = pulse.centerFreqRefined;
        // Symbolzentren liegen bei Umschaltspitze + sps/2; Frequenz pro
        // Symbol als Median ueber +/-2 Samples um die Mitte (robust gegen
        // die Umschaltflanken).
        for (double p = bestOff + 0.5 * sps; p < static_cast<double>(n1) - 3.0;
             p += sps) {
            const std::size_t c = static_cast<std::size_t>(std::lround(p));
            std::vector<double> w;
            for (long d = -2; d <= 2; ++d) {
                const std::size_t i = c + static_cast<std::size_t>(d);
                if (i < n1) w.push_back(f.finst[i]);
            }
            if (w.empty()) continue;
            std::sort(w.begin(), w.end());
            const double fs = w[w.size() / 2];
            r.symbolFreqs.push_back(fs);
            r.fskTone.push_back(fs > 0.0 ? 1 : 0);
        }
        if (static_cast<int>(r.symbolFreqs.size()) <
            static_cast<int>(cfg_.timingMinSyms))
            r.valid = false;
        return r;
    }
    return r;
}

SymbolStream ModulationClassifier::extractSymbolsF(
    const std::complex<float>* iq, std::size_t n,
    const PulseResult& pulse) const {
    if (iq == nullptr || n == 0) return SymbolStream{};
    std::vector<std::complex<double>> buf(n);
    for (std::size_t i = 0; i < n; ++i)
        buf[i] = std::complex<double>(static_cast<double>(iq[i].real()),
                                      static_cast<double>(iq[i].imag()));
    return extractSymbols(buf.data(), buf.size(), pulse);
}

}  // namespace iqclassify
