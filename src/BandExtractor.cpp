// IQClassify - Banddetektion, Basisbandextraktion, SNR/Flachheit.

#include "iqclassify/BandExtractor.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace iqclassify {

BandExtractor::BandDetection BandExtractor::bands(const std::complex<double>* x,
                                                  std::size_t a, std::size_t b,
                                                  double floorRef) const {
    const ClassifierConfig& cfg = cfg_;
    BandDetection det;
    const std::size_t L = b - a;
    const std::size_t wlen = std::min(L, cfg.detectWin);
    if (wlen < 16) {
        // Haertung: Auch im Fruehabruch muss psdDb eine definierende
        // Fuellung tragen (nft Bins), damit estimateSnrFlat nie auf ein
        // leeres PSD stoesst (Modulo-Division durch 0).
        det.nft = detail::pow2ceil(L) < 16 ? 16 : detail::pow2ceil(L);
        det.psdDb.assign(det.nft, -300.0);
        return det;
    }
    const std::size_t nft = detail::pow2ceil(wlen);
    const std::size_t hop = wlen / 2;
    std::vector<double> win(wlen);
    for (std::size_t i = 0; i < wlen; ++i)
        win[i] = 0.5 - 0.5 * std::cos(2.0 * 3.14159265358979323846 *
                                      static_cast<double>(i) /
                                      static_cast<double>(wlen));
    // Energie-Gate (Performance): Fenster, deren Energie klar unter dem
    // Rauschboden des Puffers liegt, tragen zum Max-Hold nichts bei (aktive
    // Fenster liefern den Rauschboden ohnehin) und werden uebersprungen.
    // Referenz ist der globale Puffer-Boden (floorRef, z. B. vom Segmenter),
    // nicht der des Segments — ein Segment kann vollstaendig vom Puls
    // ausgefuellt sein. Gate-Faktor 1.5: selbst ein 0-dB-Puls (Gesamtleistung
    // Signal = Rauschen) bleibt liegen, reine Rauschfenster nicht.
    const double gate =
        floorRef > 0.0 ? floorRef * static_cast<double>(wlen) * 1.5 : -1.0;
    // per-Bin-Maximum ueber alle Hann-Teilfenster (50 % Ueberlappung);
    // buf wird wiederverwendet (kein Neu-Allokieren pro Fenster).
    std::vector<double> psd(nft, 0.0);
    std::vector<detail::cd> buf(nft);
    const std::size_t lastStart = L >= wlen ? L - wlen + 1 : 1;
    for (std::size_t s = 0; s < lastStart; s += (hop ? hop : 1)) {
        std::fill(buf.begin(), buf.end(), detail::cd(0.0, 0.0));
        const std::size_t lim = std::min(wlen, L - s);
        double ewin = 0.0;
        double wmax = 0.0;
        for (std::size_t t = 0; t < lim; ++t) {
            const detail::cd v = x[a + s + t];
            // Haertung: NaN/Inf vergiften via FFT das ganze Fenster,
            // als Ausfall (Null) behandeln.
            if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) {
                buf[t] = detail::cd(0.0, 0.0);
            } else {
                buf[t] = v * win[t];
                ewin += v.real() * v.real() + v.imag() * v.imag();
                const double av = std::abs(v);
                if (av > wmax) wmax = av;
            }
        }
        if (gate >= 0.0 && ewin < gate) continue;  // stilles Fenster
        // Haertung: Extremamplituden vor der FFT normalisieren (PSD ist
        // skaleninvariant; sonst wuerde |X|^2 ueberlaufen).
        if (wmax > 1e100)
            for (std::size_t t = 0; t < lim; ++t) buf[t] /= wmax;
        detail::fftInPlace(buf, false);
        for (std::size_t k = 0; k < nft; ++k) {
            const double r = buf[k].real();
            const double im = buf[k].imag();
            // Haertung: Leistung auf 1e300 begrenzen (sonst inf bei
            // Extremamplituden; db10 begrenzt anschliessend auf +300 dB).
            double v = r * r + im * im;
            if (!(v < 1e300)) v = 1e300;
            if (v > psd[k]) psd[k] = v;
        }
    }
    det.nft = nft;
    det.psdDb.resize(nft);
    for (std::size_t k = 0; k < nft; ++k) det.psdDb[k] = detail::db10(psd[k]);
    const std::vector<double>& psdDb = det.psdDb;

    const double floorDb = detail::medianOf(psdDb);
    const double thr = floorDb + cfg.bandGuardDb;
    // stillstes Bin (Start der zyklischen Koordinaten)
    std::size_t qb = 0;
    double qv = psdDb[0];
    for (std::size_t i = 1; i < nft; ++i) {
        if (psdDb[i] < qv) {
            qv = psdDb[i];
            qb = i;
        }
    }
    const long minBins =
        static_cast<long>(nft / 512) > 8 ? static_cast<long>(nft / 512) : 8;
    const long minGap =
        static_cast<long>(nft / 512) > 4 ? static_cast<long>(nft / 512) : 4;

    // Laeufe ueber der Schwelle, ausgehend vom stillsten Bin
    std::vector<std::pair<long, long>> runs;
    {
        std::size_t i = 0;
        while (i < nft) {
            const bool up = psdDb[(qb + 1 + i) % nft] > thr;
            if (up) {
                std::size_t j = i;
                while (j < nft && psdDb[(qb + 1 + j) % nft] > thr) ++j;
                runs.emplace_back(static_cast<long>(i), static_cast<long>(j));
                i = j;
            } else {
                ++i;
            }
        }
    }
    // Laenge von Luecken <= minGap: direkt verschmelzen
    std::vector<std::pair<long, long>> mr;
    for (const auto& r : runs) {
        if (!mr.empty() && r.first - mr.back().second <= minGap)
            mr.back().second = r.second;
        else
            mr.push_back(r);
    }
    // Low-SNR-Bruecke: Liegt das mittlere PSD in der Luecke zwischen zwei
    // Laeufen noch klar ueber dem Rauschboden, enthaelt die Luecke Signal
    // (kein Rauschen) -> verschmelzen. So bleibt ein breitbandiger Puls auch
    // bei ~2 dB SNR ein einziges Band, waehrend echte Rauschluecken (Mittel
    // auf Boden) zwei Pulse getrennt halten.
    {
        bool changed = mr.size() > 1;
        while (changed) {
            changed = false;
            std::size_t i = 0;
            while (i + 1 < mr.size()) {
                const long gapLo = mr[i].second;
                const long gapHi = mr[i + 1].first;
                if (gapHi - gapLo > 0) {
                    double s = 0.0;
                    long c = 0;
                    for (long g = gapLo; g < gapHi; ++g) {
                        s += psdDb[(qb + 1 + static_cast<std::size_t>(g)) % nft];
                        ++c;
                    }
                    if (c > 0 &&
                        s / static_cast<double>(c) > floorDb + cfg.bandBridgeDb) {
                        mr[i].second = mr[i + 1].second;
                        mr.erase(mr.begin() + i + 1);
                        changed = true;
                        continue;
                    }
                }
                ++i;
            }
        }
    }
    std::vector<Band> bands;
    for (const auto& r : mr) {
        long lo = static_cast<long>(qb) + 1 + r.first;
        long hi = static_cast<long>(qb) + 1 + r.second;
        if (hi - lo < 1) continue;
        const double mthr = floorDb + cfg.bandMeanDb;
        while (hi - lo > minBins && psdDb[static_cast<std::size_t>(lo % static_cast<long>(nft))] < mthr)
            ++lo;
        while (hi - lo > minBins &&
               psdDb[static_cast<std::size_t>((hi - 1) % static_cast<long>(nft))] < mthr)
            --hi;
        if (hi - lo < 1) continue;
        double peak = -1e300;
        double sumv = 0.0;
        for (long k = lo; k < hi; ++k) {
            const double v = psdDb[static_cast<std::size_t>(k % static_cast<long>(nft))];
            if (v > peak) peak = v;
            sumv += v;
        }
        const double meanv = sumv / static_cast<double>(hi - lo);
        const bool strong = peak > floorDb + cfg.bandStrongDb;
        const bool normal = (hi - lo) >= minBins &&
                             peak > floorDb + cfg.bandPeakDb && meanv >= mthr;
        if (strong || normal) bands.push_back(Band{lo, hi});
    }

    // Spektrumssatelliten zurueck in die emittierende Band verschmelzen.
    // Phasen-/Frequenz-Keying erzeugt Beulen etwa eine Symbolrate neben
    // dem Hauptzipfel; tiefe Einbrueche zerteilen sie zu eigenen Laeufen.
    // relatives Kriterium: Luecke <= Anteil * max(Breite). Bei wenig SNR
    // tragen Satelliten nur schmale Luecken: zusaetzlich eine kleine
    // Absolutluecke fuer die Anhaengung eines schmalen Laufs an eine
    // mindestens sechsmal breitere Band; zwei schmale Nachbarn bleiben
    // getrennt (dicht beieinanderliegende Toene ueberleben).
    const double mf = cfg.bandMergeFrac;
    const long mabs = static_cast<long>(nft / static_cast<std::size_t>(cfg.satMergeDenom)) > cfg.satMergeMin
                          ? static_cast<long>(nft / static_cast<std::size_t>(cfg.satMergeDenom))
                          : cfg.satMergeMin;
    auto merged = bands;  // Kopie fuer die Verschmelzung
    bool changed = merged.size() > 1;
    while (changed) {
        changed = false;
        std::size_t i = 0;
        while (i + 1 < merged.size()) {
            Band& ai = merged[i];
            const Band& bi = merged[i + 1];
            const long wa = ai.width();
            const long wb = bi.width();
            const long gap = bi.lo - ai.hi;
            const long wmin = wa < wb ? wa : wb;
            const long wmax = wa < wb ? wb : wa;
            const bool rel = static_cast<double>(gap) <= mf * static_cast<double>(wmax);
            const bool sat = gap <= mabs &&
                             static_cast<double>(wmax) >=
                                 cfg.satMergeWidthRatio * static_cast<double>(wmin < 1 ? 1 : wmin);
            if (rel || sat) {
                ai.hi = bi.hi;
                merged.erase(merged.begin() + i + 1);
                changed = true;
            } else {
                ++i;
            }
        }
        if (merged.size() > 1) {
            Band& first = merged.front();
            const Band& last = merged.back();
            const long gap = static_cast<long>(first.lo + static_cast<long>(nft) - last.hi);
            const long wf = first.width();
            const long wl = last.width();
            const long wmin = wf < wl ? wf : wl;
            const long wmax = wf < wl ? wl : wf;
            const bool rel = static_cast<double>(gap) <= mf * static_cast<double>(wmax);
            const bool sat = gap <= mabs &&
                             static_cast<double>(wmax) >=
                                 cfg.satMergeWidthRatio * static_cast<double>(wmin < 1 ? 1 : wmin);
            if (rel || sat) {
                first.lo = last.lo - static_cast<long>(nft);
                merged.pop_back();
                changed = true;
            }
        }
    }
    // Normierung der virtuellen Koordinaten: hi <= nft halten
    det.bands.clear();
    for (const Band& bd : merged) {
        if (bd.hi > static_cast<long>(nft))
            det.bands.push_back(Band{bd.lo - static_cast<long>(nft),
                                     bd.hi - static_cast<long>(nft)});
        else
            det.bands.push_back(bd);
    }
    return det;
}

BandExtractor::Slice BandExtractor::extract(const std::complex<double>* x,
                                           std::size_t a, std::size_t b,
                                           const Band& band,
                                           std::size_t nft) const {
    Slice out;
    const std::size_t L = b - a;
    std::size_t padded = detail::pow2ceil(L);
    if (padded < nft) padded = nft;
    if (padded % nft) padded = detail::pow2ceil(padded);
    const long scale = static_cast<long>(padded / nft);
    std::vector<detail::cd> buf(L);
    // Haertung: nicht-endliche Samples als Ausfall (Null) behandeln,
    // sonst vergiften sie das komplette Spektrum des Puffers.
    double amax = 0.0;
    for (std::size_t i = 0; i < L; ++i) {
        const std::complex<double> v = x[a + i];
        if (std::isfinite(v.real()) && std::isfinite(v.imag())) {
            buf[i] = detail::cd(v.real(), v.imag());
            const double av = std::abs(v);
            if (av > amax) amax = av;
        } else {
            buf[i] = detail::cd(0.0, 0.0);
        }
    }
    // Haertung: Extremamplituden vor der FFT normalisieren (FFT ist linear;
    // sonst wuerden Betragsquadrate stromabwaerts ueberlaufen). Nach der IFFT
    // wird zurueckskaliert, damit Slice.y die Originalamplitude behaelt.
    const double xscale = amax > 1e100 ? amax : 1.0;
    if (xscale != 1.0)
        for (std::size_t i = 0; i < L; ++i) buf[i] /= xscale;
    buf.resize(padded, detail::cd(0.0, 0.0));
    detail::fftInPlace(buf, false);
    // Maske mit Waechter: die vollstaendige Uebergangszone des Pulses
    // mitnehmen; eine exakt auf dem detektierten Band liegende Maske
    // schneidet den spektralen Uebergang ab und injiziert Ringing-ISI.
    long lo = band.lo;
    long hi = band.hi;
    long width = hi - lo;
    const long guard = width / 12 > 3 ? width / 12 : 3;
    lo -= guard;
    hi += guard;
    width = hi - lo;
    const long lox = lo * scale;
    const long widx = width * scale;
    std::vector<detail::cd> masked(padded, detail::cd(0.0, 0.0));
    const long lp = static_cast<long>(padded);
    for (std::size_t k = 0; k < padded; ++k) {
        const long kk = static_cast<long>(k);
        const long rel =
            ((kk - lox) % lp + lp) % lp;  // python-artiges modulo, lox < 0 moeglich
        if (rel < widx) masked[k] = buf[k];
    }
    detail::fftInPlace(masked, true);
    const double center = 0.5 * static_cast<double>(lo + hi);
    out.fc = center / static_cast<double>(nft);
    out.bw = static_cast<double>(width) / static_cast<double>(nft);
    const double phStep = -2.0 * 3.14159265358979323846 * out.fc;
    out.y.resize(L);
    for (std::size_t i = 0; i < L; ++i) {
        const double ang = phStep * static_cast<double>(i);
        const double cs = std::cos(ang);
        const double sn = std::sin(ang);
        const std::complex<double> rot(cs, sn);
        out.y[i] = masked[i] * rot * xscale;
    }
    return out;
}

void estimateSnrFlat(const std::vector<double>& psdDb, const Band& band,
                     const ClassifierConfig& cfg,
                     const std::vector<Band>& allBands, double* snrDb,
                     double* flatDb) {
    // Haertung: leeres PSD oder degenerierte Bandbreite sind definiert zu
    // behandeln statt durch 0 zu teilen.
    if (psdDb.empty() || snrDb == nullptr || flatDb == nullptr) {
        if (snrDb) *snrDb = -99.0;
        if (flatDb) *flatDb = 0.0;
        return;
    }
    const long nftL = static_cast<long>(psdDb.size());
    const long lo = band.lo;
    const long hi = band.hi;
    const long width = hi - lo;
    if (width <= 0 || nftL <= 0) {
        *snrDb = -99.0;
        *flatDb = 0.0;
        return;
    }
    (void)cfg;
    std::vector<double> inb;
    inb.reserve(static_cast<std::size_t>(width));
    for (long k = lo; k < hi; ++k)
        inb.push_back(psdDb[static_cast<std::size_t>(((k % nftL) + nftL) % nftL)]);
    // Aussenreferenz ohne die Bloecke fremder Bander (+- Breite dieses
    // Bandes als Waechter); sonst kollabiert die SNR-Schaetzung im
    // Mehrpuls-Fall. Markierung als Byte-Feld statt std::set: O(nft) statt
    // O(Baender * Breite * log) pro Aufruf.
    std::vector<char> excl(static_cast<std::size_t>(nftL), 0);
    for (const Band& ob : allBands) {
        if (ob.lo == band.lo && ob.hi == band.hi) continue;
        for (long k = ob.lo - width; k < ob.hi + width; ++k)
            excl[static_cast<std::size_t>(((k % nftL) + nftL) % nftL)] = 1;
    }
    std::vector<double> outb;
    outb.reserve(psdDb.size());
    for (long k = 0; k < nftL; ++k) {
        const long rel = ((k - lo) % nftL + nftL) % nftL;
        if (rel < width) continue;
        if (excl[static_cast<std::size_t>(k)]) continue;
        outb.push_back(psdDb[static_cast<std::size_t>(k)]);
    }
    if (outb.empty()) {
        for (long k = 0; k < nftL; ++k) {
            const long rel = ((k - lo) % nftL + nftL) % nftL;
            if (rel < width) continue;
            outb.push_back(psdDb[static_cast<std::size_t>(k)]);
        }
    }
    double pIn = 0.0;
    for (double v : inb) pIn += std::pow(10.0, v / 10.0);
    pIn /= inb.empty() ? 1.0 : static_cast<double>(inb.size());
    double pOut = 0.0;
    for (double v : outb) pOut += std::pow(10.0, v / 10.0);
    pOut /= outb.empty() ? 1.0 : static_cast<double>(outb.size());
    // Spektraldichte-SNR: Signaldichte ueber der Rauschdichte des
    // Segments. Bandbreitenunabhaengig, schmale Toene werden nicht
    // verduennt.
    double ratio = pIn / (pOut > 1e-30 ? pOut : 1e-30) - 1.0;
    *snrDb = ratio > 0 ? detail::db10(ratio > 1e-12 ? ratio : 1e-12) : -99.0;
    // Flachheit: stark geglaettete in-band-dB-Kurve, Std des Inneren
    const std::size_t sm = static_cast<std::size_t>(width / 24) > 3
                               ? static_cast<std::size_t>(width / 24)
                               : 3;
    std::vector<double> curve;
    double acc = 0.0;
    for (std::size_t i = 0; i < inb.size(); ++i) {
        acc += inb[i];
        if (i >= sm) acc -= inb[i - sm];
        if (i >= sm / 2) curve.push_back(acc / static_cast<double>(sm));
    }
    if (curve.size() < 4) {
        *flatDb = 0.0;
        return;
    }
    const std::size_t c0 = curve.size() / 8;
    const std::size_t c1 = curve.size() - c0;
    double m = 0.0;
    double cnt = 0.0;
    for (std::size_t i = c0; i < c1; ++i) {
        m += curve[i];
        cnt += 1.0;
    }
    m /= cnt > 0 ? cnt : 1.0;
    double acc2 = 0.0;
    for (std::size_t i = c0; i < c1; ++i) acc2 += (curve[i] - m) * (curve[i] - m);
    *flatDb = std::sqrt(acc2 / (cnt > 0 ? cnt : 1.0));
}

}  // namespace iqclassify
