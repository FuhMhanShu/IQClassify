// IQClassify - Merkmalsextraktion und zeitliche Crops.

#include "iqclassify/FeatureExtractor.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace iqclassify {

namespace {

const double kPi = 3.14159265358979323846;

}  // namespace

Features extractFeatures(const std::vector<std::complex<double>>& y,
                          double fc, double bw, std::size_t bandBins,
                          double snrDb, double flatDb,
                          const ClassifierConfig& cfg) {
    const std::size_t n = y.size();
    Features f;
    f.n = n;
    f.fc = fc;
    f.bw = bw;
    f.bandBins = bandBins;
    f.snrDb = snrDb;
    f.flatDb = flatDb;

    f.env.resize(n);
    double me = 0.0;
    double envSum = 0.0;
    {   // Hoellkurve
        double acc = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            f.env[i] = std::abs(y[i]);
            acc += f.env[i];
        }
        me = n ? acc / static_cast<double>(n) : 0.0;
        envSum = acc;
    }
    // Momentane Phase/Frequenz und Resttraeger: zuerst berechnen, da die
    // Traegerkohaerenz den Resttraeger (Median von f_inst) entfernen muss.
    std::vector<double> up;
    double med = 0.0;
    {
        up = detail::unwrapPhase(y);
        f.finst.resize(n >= 1 ? n - 1 : 0);
        for (std::size_t i = 1; i < n; ++i)
            f.finst[i - 1] = (up[i] - up[i - 1]) / (2.0 * kPi);
        med = detail::medianOf(f.finst);
        f.residF = med;
        for (double& v : f.finst) v -= med;
        f.ph.resize(n);
        for (std::size_t i = 0; i < n; ++i)
            f.ph[i] = up[i] - 2.0 * kPi * med * static_cast<double>(i);
    }
    // Traegerkohaerenz nach Entfernung des Resttraegers:
    // |mean(z * e^{-j*2pi*med*i})| / mean(|z|). Ohne Entfernung wuerde ein
    // sub-Bin-Restoffset (Bandmitte liegt auf quantisiertem FFT-Raster)
    // den Zeiger ueber lange Crop-Laengen wegmitteln (Sinc-Ausloeschung)
    // und die CW-Erkennung crop-laengen- und positionsabhaengig machen.
    // Trennt weiterhin BPSK/QPSK (Phasenspruenge mitteln den Zeiger weg)
    // vom unmodulierten Traeger.
    {
        std::complex<double> acc(0.0, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const double ang = up[i] - 2.0 * kPi * med *
                static_cast<double>(i);
            acc += std::polar(f.env[i], ang);
        }
        f.carrCoh = std::abs(acc) / (envSum > 1e-12 ? envSum : 1e-12);
    }
    {   // CV / PAPR
        double ve = 0.0;
        double pmax = 0.0;
        double pmAcc = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double v = f.env[i];
            ve += (v - me) * (v - me);
            pmax = std::max(pmax, v * v);
            pmAcc += v * v;
        }
        ve /= n ? static_cast<double>(n) : 1.0;
        f.cv = std::sqrt(ve) / (me > 1e-12 ? me : 1e-12);
        const double pm = n ? pmAcc / static_cast<double>(n) : 1e-30;
        f.paprDb = detail::db10(pmax / (pm > 1e-30 ? pm : 1e-30));
        // fenstergemitteltes PAPR (Fensterlaenge 16)
        const std::size_t W = std::max<std::size_t>(
            2, std::min<std::size_t>(16, n / 32));
        double acc2 = 0.0;
        for (std::size_t i = 0; i < std::min(W, n); ++i)
            acc2 += f.env[i] * f.env[i];
        double wmax = acc2;
        if (n > W) {
            for (std::size_t i = 0; i + W < n; ++i) {
                acc2 += f.env[i + W] * f.env[i + W] - f.env[i] * f.env[i];
                if (acc2 > wmax) wmax = acc2;
            }
        }
        f.paprWDb = detail::db10(wmax / (W * pm > 1e-30 ? W * pm : 1e-30));
    }
    // Chirp-Fit auf f_inst, hoellkurven^2-gewichtet, mittlere 80 %
    const std::size_t n1 = n >= 1 ? n - 1 : 0;
    const std::size_t iLo = n1 / 10;
    const std::size_t iHi = n1 - n1 / 10;
    {
        double sw = 0.0;
        double xm = 0.0;
        double ym = 0.0;
        const std::size_t iHiEff = std::min(iHi, n1);
        for (std::size_t i = iLo; i < iHiEff; ++i) {
            const double w = f.env[i] * f.env[i];
            sw += w;
            xm += w * static_cast<double>(i);
            ym += w * f.finst[i];
        }
        xm /= sw > 1e-30 ? sw : 1e-30;
        ym /= sw > 1e-30 ? sw : 1e-30;
        double num = 0.0;
        double den = 0.0;
        double ssTot = 0.0;
        for (std::size_t i = iLo; i < iHiEff; ++i) {
            const double w = f.env[i] * f.env[i];
            const double dx = static_cast<double>(i) - xm;
            num += w * dx * (f.finst[i] - ym);
            den += w * dx * dx;
            ssTot += w * (f.finst[i] - ym) * (f.finst[i] - ym);
        }
        const double slope = num / (den > 1e-30 ? den : 1e-30);
        double ssRes = 0.0;
        for (std::size_t i = iLo; i < iHiEff; ++i) {
            const double w = f.env[i] * f.env[i];
            const double resid = f.finst[i] -
                (ym + slope * (static_cast<double>(i) - xm));
            ssRes += w * resid * resid;
        }
        f.chirpR2 = 1.0 - ssRes / (ssTot > 1e-30 ? ssTot : 1e-30);
        f.chirpSweep = std::abs(slope) * static_cast<double>(iHi - iLo) /
                       (bw > 1e-9 ? bw : 1e-9);
    }
    // Dwell-Anteil: Schwelle relativ zum 90. Perzentil von |f_inst|
    {
        std::vector<double> absf;
        for (std::size_t i = 0; i < n1; ++i)
            if (f.env[i] > 0.4 * me) absf.push_back(std::abs(f.finst[i]));
        double thr = 0.40 * (bw > 1e-3 ? bw : 1e-3);
        if (absf.size() > 16) {
            std::sort(absf.begin(), absf.end());
            const double q90 = absf[static_cast<std::size_t>(
                0.90 * static_cast<double>(absf.size() - 1))];
            thr = 0.40 * q90;
        }
        f.dwellThr = thr;
        int run = 0;
        int runSum = 0;
        for (double v : f.finst) {
            if (std::abs(v) > thr) {
                ++run;
            } else {
                if (run >= cfg.tDwellRun) runSum += run;
                run = 0;
            }
        }
        if (run >= cfg.tDwellRun) runSum += run;
        f.dwell = runSum / (n1 > 0 ? static_cast<double>(n1) : 1.0);
    }
    // dominante f_inst-Moden (nur Hoellekurve > Gate)
    const int nbins = cfg.histBins;
    const double bwEff = bw > 1e-3 ? bw : 1e-3;
    std::vector<double> hist(static_cast<std::size_t>(nbins), 0.0);
    double tot = 0.0;
    for (std::size_t i = 0; i < n1; ++i) {
        if (f.env[i] < cfg.histEnvRel * me) continue;
        const double pos = (f.finst[i] / bwEff + 0.5) * static_cast<double>(nbins);
        // Haertung: NaN-Schutz; ein Cast von NaN nach long waere
        // undefiniert (Defense in Depth, Input ist bereits sanitisiert).
        if (!std::isfinite(pos)) continue;
        const long bi = static_cast<long>(pos);
        if (bi >= 0 && bi < nbins) {
            hist[static_cast<std::size_t>(bi)] += 1.0;
            tot += 1.0;
        }
    }
    f.nfreqModes = 0;
    f.modeConc = 0.0;
    if (tot >= 16) {
        std::vector<double> sm(static_cast<std::size_t>(nbins), 0.0);
        for (int i = 0; i < nbins; ++i) {
            double acc = hist[static_cast<std::size_t>(i)];
            if (i > 0) acc += hist[static_cast<std::size_t>(i - 1)];
            if (i < nbins - 1) acc += hist[static_cast<std::size_t>(i + 1)];
            // Randbins teilen durch 3 wie das Referenzmodell: die
            // angehaengten Nullnachbarn zaehlen als Masse mit.
            sm[static_cast<std::size_t>(i)] = acc / 3.0;
        }
        const double hmax = *std::max_element(sm.begin(), sm.end());
        if (hmax > 0) {
            std::vector<std::pair<double, std::pair<int, int>>> groups;
            int i = 0;
            while (i < nbins) {
                if (sm[static_cast<std::size_t>(i)] > cfg.modeBinThr * hmax) {
                    int j = i;
                    while (j + 1 < nbins &&
                           sm[static_cast<std::size_t>(j + 1)] > cfg.modeBinThr * hmax)
                        ++j;
                    double mass = 0.0;
                    for (int k = i; k <= j; ++k)
                        mass += hist[static_cast<std::size_t>(k)];
                    groups.push_back({mass, {i, j}});
                    i = j + 1;
                } else {
                    ++i;
                }
            }
            std::sort(groups.begin(), groups.end(),
                      [](const std::pair<double, std::pair<int, int>>& a,
                         const std::pair<double, std::pair<int, int>>& b) {
                          return a.first > b.first;
                      });
            f.nfreqModes = static_cast<int>(groups.size());
            const int reach =
                std::max(1, static_cast<int>(std::lround(
                                static_cast<double>(nbins) / 12.0)));
            double mass = 0.0;
            const std::size_t topCount = groups.size() < 2 ? groups.size() : 2;
            for (int gi = 0; gi < nbins; ++gi) {
                for (std::size_t t = 0; t < topCount; ++t) {
                    const double center =
                        0.5 * static_cast<double>(groups[t].second.first +
                                                  groups[t].second.second);
                    if (std::abs(static_cast<double>(gi) - center) <=
                        reach + 0.5) {
                        mass += hist[static_cast<std::size_t>(gi)];
                        break;
                    }
                }
            }
            f.modeConc = mass / tot;
        }
    }
    return f;
}

std::vector<std::pair<std::size_t, std::size_t>> timeCrops(
    const std::vector<std::complex<double>>& y, const ClassifierConfig& cfg) {
    const std::size_t n = y.size();
    const std::size_t w = cfg.segWin;
    const std::size_t nb = (n + w - 1) / w;
    std::vector<double> energy(nb, 0.0);
    for (std::size_t i = 0; i < nb; ++i) {
        const std::size_t lo = i * w;
        const std::size_t hi = std::min(n, lo + w);
        double s = 0.0;
        for (std::size_t k = lo; k < hi; ++k)
            s += y[k].real() * y[k].real() + y[k].imag() * y[k].imag();
        energy[i] = s / static_cast<double>(hi - lo);
    }
    const std::size_t g = cfg.segFloorGroup;
    const std::size_t ng = (nb + g - 1) / g;
    bool haveFloor = false;
    double floor = 0.0;
    constexpr double kSilent = 1e-20;
    for (std::size_t gi = 0; gi < ng; ++gi) {
        const std::size_t lo = gi * g;
        const std::size_t hi = std::min(nb, lo + g);
        std::vector<double> blk(energy.begin() + static_cast<std::ptrdiff_t>(lo),
                                energy.begin() + static_cast<std::ptrdiff_t>(hi));
        if (blk.empty()) continue;
        const double m = detail::median(blk);
        if (m <= kSilent) continue;
        if (!haveFloor || m < floor) {
            floor = m;
            haveFloor = true;
        }
    }
    const double emax = energy.empty() ? 0.0 : *std::max_element(energy.begin(), energy.end());
    std::vector<std::pair<std::size_t, std::size_t>> crops;
    // Kein eigenes Schweigen im Slice (ein einzelner Puls fuellt das
    // extrahierte Segment): dann ist kein Split moeglich, der ganze Slice
    // gilt als ein Puls.
    if (!haveFloor || floor > 0.1 * emax) {
        crops.emplace_back(0, n);
        return crops;
    }
    if (floor < 1e-30) floor = 1e-30;
    const double onThr = std::pow(10.0, cfg.segOnDb / 10.0);
    const double offThr = std::pow(10.0, cfg.segOffDb / 10.0);
    const double onRel = std::pow(10.0, -cfg.cropOnRelDb / 10.0);
    const double offRel = std::pow(10.0, -cfg.cropOffRelDb / 10.0);
    bool inCrop = false;
    std::size_t start = 0;
    for (std::size_t i = 0; i < nb; ++i) {
        const double e = energy[i];
        const double r = e / floor;
        const double rel = e / (emax > 1e-30 ? emax : 1e-30);
        if (!inCrop) {
            if (r > onThr && rel > onRel) {
                start = i;
                inCrop = true;
            }
        } else {
            if (r < offThr || rel < offRel) {
                const std::size_t i0 =
                    start * w > cfg.segMargin ? start * w - cfg.segMargin : 0;
                const std::size_t i1 = std::min(n, i * w + w + cfg.segMargin);
                if (i1 - i0 >= cfg.minPulseLen) crops.emplace_back(i0, i1);
                inCrop = false;
            }
        }
    }
    if (inCrop) {
        const std::size_t i0 =
            start * w > cfg.segMargin ? start * w - cfg.segMargin : 0;
        const std::size_t i1 = std::min(n, nb * w);
        if (i1 - i0 >= cfg.minPulseLen) crops.emplace_back(i0, i1);
    }
    return crops;
}

}  // namespace iqclassify
