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
        double emax = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            f.env[i] = std::abs(y[i]);
            acc += f.env[i];
            if (f.env[i] > emax) emax = f.env[i];
        }
        // Haertung: Extremamplituden normalisieren (sonst wuerde env^2 in
        // den Statistiken ueberlaufen -> inf/NaN). Alle Nutzer von f.env
        // sind skaleninvariant (CV/PAPR als Verhältnisse, Chirp-Fit und
        // Symbolratensuche mittelwertfrei, Gates relativ zum Median).
        if (emax > 1e100) {
            for (double& v : f.env) v /= emax;
            acc /= emax;
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
    // Chirp-Fit auf geglaetteter f_inst (finstSm), hoellkurven^2-gewichtet,
    // mittlere 80 %. Die Glaettung (gleitendes Mittel ueber ~n/32 Samples,
    // O(n)) unterdrueckt Ausreisser/Entpackungsfehler, die bei stark
    // verrauschtem f_inst sonst das Fit-R^2 ruinieren.
    const std::size_t n1 = n >= 1 ? n - 1 : 0;
    const std::size_t iLo = n1 / 10;
    const std::size_t iHi = n1 - n1 / 10;
    {
        const std::size_t Wsm = n1 / 32 > 8 ? n1 / 32 : 8;
        f.finstSm.assign(n1, 0.0);
        double acc = 0.0;
        for (std::size_t i = 0; i < n1; ++i) {
            acc += f.finst[i];
            if (i >= Wsm) acc -= f.finst[i - Wsm];
            const std::size_t cnt = std::min(i + 1, Wsm);
            f.finstSm[i] = acc / static_cast<double>(cnt);
        }
    }
    {
        double sw = 0.0;
        double xm = 0.0;
        double ym = 0.0;
        const std::size_t iHiEff = std::min(iHi, n1);
        for (std::size_t i = iLo; i < iHiEff; ++i) {
            const double w = f.env[i] * f.env[i];
            sw += w;
            xm += w * static_cast<double>(i);
            ym += w * f.finstSm[i];
        }
        xm /= sw > 1e-30 ? sw : 1e-30;
        ym /= sw > 1e-30 ? sw : 1e-30;
        double num = 0.0;
        double den = 0.0;
        double ssTot = 0.0;
        for (std::size_t i = iLo; i < iHiEff; ++i) {
            const double w = f.env[i] * f.env[i];
            const double dx = static_cast<double>(i) - xm;
            num += w * dx * (f.finstSm[i] - ym);
            den += w * dx * dx;
            ssTot += w * (f.finstSm[i] - ym) * (f.finstSm[i] - ym);
        }
        const double slope = num / (den > 1e-30 ? den : 1e-30);
        double ssRes = 0.0;
        for (std::size_t i = iLo; i < iHiEff; ++i) {
            const double w = f.env[i] * f.env[i];
            const double resid = f.finstSm[i] -
                (ym + slope * (static_cast<double>(i) - xm));
            ssRes += w * resid * resid;
        }
        f.chirpR2 = 1.0 - ssRes / (ssTot > 1e-30 ? ssTot : 1e-30);
        f.chirpSweep = std::abs(slope) * static_cast<double>(iHi - iLo) /
                       (bw > 1e-9 ? bw : 1e-9);
        f.chirpSweepAbs = std::abs(slope) * static_cast<double>(iHi - iLo);
        // Gefalteter Chirp (V-Form): Extremlage von f_inst suchen, beide
        // Haelften getrennt linear fitten. Nur gesetzt, wenn beide Haelften
        // hochkorreliert sind und gegenlaeufige Steigungen haben.
        auto fitHalf = [&](const std::vector<double>& fi, std::size_t a,
                           std::size_t b, double* r2, double* sl) {
            double sw = 0.0, xm = 0.0, ym2 = 0.0;
            for (std::size_t i = a; i < b; ++i) {
                const double w = f.env[i] * f.env[i];
                sw += w;
                xm += w * static_cast<double>(i);
                ym2 += w * fi[i];
            }
            if (sw <= 1e-30) {
                *r2 = 0.0;
                *sl = 0.0;
                return;
            }
            xm /= sw;
            ym2 /= sw;
            double num = 0.0, den = 0.0, ssTot = 0.0;
            for (std::size_t i = a; i < b; ++i) {
                const double w = f.env[i] * f.env[i];
                const double dx = static_cast<double>(i) - xm;
                num += w * dx * (fi[i] - ym2);
                den += w * dx * dx;
                ssTot += w * (fi[i] - ym2) * (fi[i] - ym2);
            }
            *sl = num / (den > 1e-30 ? den : 1e-30);
            double ssRes = 0.0;
            for (std::size_t i = a; i < b; ++i) {
                const double w = f.env[i] * f.env[i];
                const double resid =
                    fi[i] - (ym2 + *sl * (static_cast<double>(i) - xm));
                ssRes += w * resid * resid;
            }
            *r2 = ssTot > 1e-30 ? 1.0 - ssRes / ssTot : 0.0;
        };
        // Saegezahn-Erkennung: groesster Abfall von f_inst an einer inneren
        // Stelle (Faltung am Chip-Rand, z. B. P3-Code); beide Seiten muessen
        // linear mit gleichsinniger Steigung sein (entfalteter Chirp).
        std::size_t iPk = iLo;
        double dropBest = 0.0;
        std::vector<double> gf;
        for (std::size_t i = iLo; i < iHiEff; ++i) {
            if (f.env[i] < 0.4 * me) continue;
            // Auslenkung aus der geglaetteten Kurve (robust gegen Rausch-
            // Ausreisser im rohen f_inst), Sprung aber auf rohem f_inst.
            gf.push_back(f.finstSm[i]);
            if (i + 1 < iHiEff && f.env[i + 1] >= 0.4 * me) {
                const double d = f.finst[i] - f.finst[i + 1];
                if (d > dropBest) {
                    dropBest = d;
                    iPk = i + 1;
                }
            }
        }
        // robuste Auslenkung (Perzentile, gegen Ausreisser an Chip-Raendern)
        double range = 0.0;
        if (gf.size() > 16) {
            std::sort(gf.begin(), gf.end());
            const double q05 = gf[static_cast<std::size_t>(
                0.05 * static_cast<double>(gf.size() - 1))];
            const double q95 = gf[static_cast<std::size_t>(
                0.95 * static_cast<double>(gf.size() - 1))];
            range = q95 - q05;
        }
        // Rand um den Sprung herum aussparen (der Abfall ist durch die
        // Impulsformung ueber einige Samples verschmiert).
        const std::size_t marg = (iHiEff - iLo) / 10;
        if (dropBest > 0.5 * range && range > 0.0 &&
            iPk > iLo + marg + (iHiEff - iLo) / 5 &&
            iPk < iHiEff - marg - (iHiEff - iLo) / 5) {
            double r2L = 0.0, r2R = 0.0, sL = 0.0, sR = 0.0;
            fitHalf(f.finstSm, iLo, iPk - marg, &r2L, &sL);
            fitHalf(f.finstSm, iPk + marg, iHiEff, &r2R, &sR);
            if (r2L > 0.0 && r2R > 0.0 && sL * sR > 0.0) {
                f.chirpFoldR2 = std::min(r2L, r2R);
                const double foldHub =
                    std::abs(sL) * static_cast<double>(iPk - marg - iLo) +
                    std::abs(sR) * static_cast<double>(iHiEff - iPk - marg) +
                    dropBest;
                f.chirpFoldSweep = foldHub / (bw > 1e-9 ? bw : 1e-9);
                f.chirpFoldSweepAbs = foldHub;
            }
        }
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
    // Fade-Merge: ein Puls mit schwankender Amplitude kann durch die relative
    // Aus-Schwelle (cropOffRelDb) in mehrere Crops zerlegt werden. Liegt die
    // Energie in der Luecke zwischen zwei Crops nie nahe dem slice-eigenen
    // Minimum (energetisch kein Schweigen, sondern ein Amplituden-Einbruch),
    // gehoeren beide Crops zum selben Puls und werden verschmolzen.
    if (crops.size() > 1 && !energy.empty()) {
        double emin = *std::min_element(energy.begin(), energy.end());
        if (emin < 1e-30) emin = 1e-30;
        const double mergeThr =
            std::pow(10.0, cfg.cropFadeMergeDb / 10.0) * emin;
        std::vector<std::pair<std::size_t, std::size_t>> merged;
        merged.push_back(crops.front());
        for (std::size_t i = 1; i < crops.size(); ++i) {
            // Luecke auf die rohen Hysterese-Grenzen (ohne angelagerten
            // Sicherheitsrand) erweitern, damit der Fade-Tiefpunkt erfasst wird.
            std::size_t gapLo = 0;
            if (merged.back().second >= w + cfg.segMargin)
                gapLo = (merged.back().second - w - cfg.segMargin) / w;
            const std::size_t gapHi =
                std::min(nb, (crops[i].first + cfg.segMargin) / w);
            double gmin = 1e300;
            for (std::size_t g = gapLo; g < gapHi; ++g)
                gmin = std::min(gmin, energy[g]);
            if (gapHi > gapLo && gmin > mergeThr) {
                merged.back().second = crops[i].second;
            } else {
                merged.push_back(crops[i]);
            }
        }
        crops.swap(merged);
    }
    return crops;
}

}  // namespace iqclassify
