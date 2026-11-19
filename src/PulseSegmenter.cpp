// IQClassify - Stufe 0: Pulssegmentierung.

#include "iqclassify/PulseSegmenter.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace iqclassify {

double PulseSegmenter::noiseFloor(const std::complex<double>* x,
                                  std::size_t n,
                                  const ClassifierConfig& cfg) {
    if (x == nullptr || n == 0) return 0.0;
    const std::size_t w = cfg.segWin;
    const std::size_t nb = n / w < 1 ? 1 : n / w;
    std::vector<double> energy(nb, 0.0);
    for (std::size_t i = 0; i < nb; ++i) {
        double s = 0.0;
        const std::size_t lo = i * w;
        const std::size_t hi = std::min(n, lo + w);
        for (std::size_t k = lo; k < hi; ++k) {
            s += x[k].real() * x[k].real() + x[k].imag() * x[k].imag();
        }
        energy[i] = s / static_cast<double>(w);
    }
    // Rauschboden: Minimum-Statistik ueber Gruppen von Bloecken.
    // Gruppen mit exakter digitaler Stille (Median <= 1e-20) werden
    // uebersprungen; sie wuerden den Boden auf null druecken und jede
    // physikalische Rauschgruppe als Puls erscheinen lassen.
    const std::size_t g = cfg.segFloorGroup;
    const std::size_t ng = (nb + g - 1) / g;
    bool haveFloor = false;
    double floor = 0.0;
    for (std::size_t gi = 0; gi < ng; ++gi) {
        const std::size_t lo = gi * g;
        const std::size_t hi = std::min(nb, lo + g);
        std::vector<double> blk(energy.begin() + lo, energy.begin() + hi);
        if (blk.empty()) continue;
        const double m = detail::median(blk);
        if (m <= 1e-20) continue;
        if (!haveFloor || m < floor) {
            floor = m;
            haveFloor = true;
        }
    }
    if (!haveFloor) return 0.0;
    return floor < 1e-30 ? 1e-30 : floor;
}

std::vector<std::pair<std::size_t, std::size_t>> PulseSegmenter::detect(
    const std::complex<double>* x, std::size_t n) const {
    const ClassifierConfig& cfg = cfg_;
    const std::size_t w = cfg.segWin;
    const std::size_t nb = n / w < 1 ? 1 : n / w;
    std::vector<double> energy(nb, 0.0);
    for (std::size_t i = 0; i < nb; ++i) {
        double s = 0.0;
        const std::size_t lo = i * w;
        const std::size_t hi = std::min(n, lo + w);
        for (std::size_t k = lo; k < hi; ++k) {
            s += x[k].real() * x[k].real() + x[k].imag() * x[k].imag();
        }
        energy[i] = s / static_cast<double>(w);
    }
    double floor = noiseFloor(x, n, cfg);
    std::vector<std::pair<std::size_t, std::size_t>> segs;
    if (floor <= 0.0) return segs;
    if (floor < 1e-30) floor = 1e-30;
    const double onThr = std::pow(10.0, cfg.segOnDb / 10.0);
    const double offThr = std::pow(10.0, cfg.segOffDb / 10.0);
    bool inPulse = false;
    std::size_t start = 0;
    for (std::size_t i = 0; i < nb; ++i) {
        const double r = energy[i] / floor;
        if (!inPulse) {
            if (r > onThr) {
                start = i;
                inPulse = true;
            }
        } else {
            if (r < offThr) {
                const std::size_t a =
                    start * w > cfg.segMargin ? start * w - cfg.segMargin : 0;
                const std::size_t bEnd =
                    std::min(n, (i + 1) * w + cfg.segMargin);
                if (bEnd - a >= cfg.minPulseLen)
                    segs.emplace_back(a, bEnd);
                inPulse = false;
            }
        }
    }
    if (inPulse) {
        const std::size_t a =
            start * w > cfg.segMargin ? start * w - cfg.segMargin : 0;
        const std::size_t bEnd = std::min(n, nb * w + cfg.segMargin);
        if (bEnd - a >= cfg.minPulseLen) segs.emplace_back(a, bEnd);
    }
    // direkt aneinanderstossende Segmente verschmelzen (Flacker-Schutz)
    for (std::size_t i = 1; i < segs.size();) {
        const std::ptrdiff_t gap =
            static_cast<std::ptrdiff_t>(segs[i].first) -
            static_cast<std::ptrdiff_t>(segs[i - 1].second);
        if (gap < static_cast<std::ptrdiff_t>(w)) {
            segs[i - 1].second = segs[i].second;
            segs.erase(segs.begin() + i);
        } else {
            ++i;
        }
    }
    return segs;
}

}  // namespace iqclassify
