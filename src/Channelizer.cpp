// IQClassify - Spektrale Kanalisierung: adaptive Kanaele um die detektierten
// Energiebereiche, Extraktion (Bandpass + Mischung) und Resampling.

#include "iqclassify/Channelizer.hpp"

#include <algorithm>
#include <cmath>

namespace iqclassify {

namespace {

const double kPi = 3.14159265358979323846;

// Zyklische Normierung einer fraktionalen Frequenz nach [0, 1).
double wrap01(double f) {
    f = std::fmod(f, 1.0);
    if (f < 0.0) f += 1.0;
    return f;
}

}  // namespace

std::vector<Band> Channelizer::energyBands(
    const BandExtractor::BandDetection& det) const {
    std::vector<Band> bands = det.bands;
    if (bands.empty()) {
        // Wie der Legacy-Pfad: ein Kanal ueber das volle Band (die SNR-
        // Schwelle der Folgestufen entscheidet dann gegen Rauschen).
        bands.push_back(Band{0, static_cast<long>(det.nft)});
    }
    std::vector<Band> out;
    out.reserve(bands.size());
    const long nftL = static_cast<long>(det.nft);
    for (const Band& b : bands) {
        // Kanal zentrieren auf den leistungsgewichteten Schwerpunkt des
        // Bandes (zyklische Aufwicklung um die Bandmitte).
        const double center = 0.5 * static_cast<double>(b.lo + b.hi);
        double num = 0.0;
        double den = 0.0;
        for (long k = b.lo; k < b.hi; ++k) {
            const long kk = ((k % nftL) + nftL) % nftL;
            const double p = std::pow(
                10.0, det.psdDb[static_cast<std::size_t>(kk)] / 10.0);
            double rel = static_cast<double>(k) - center;
            if (rel > 0.5 * nftL) rel -= nftL;
            if (rel < -0.5 * nftL) rel += nftL;
            num += rel * p;
            den += p;
        }
        const double centroid = den > 0.0 ? center + num / den : center;
        const long shift =
            static_cast<long>(std::llround(centroid - center));
        out.push_back(Band{b.lo + shift, b.hi + shift});
    }
    return out;
}

std::vector<Channelizer::Channel> Channelizer::channelize(
    const std::complex<double>* iq, std::size_t n,
    std::size_t decimate) const {
    std::vector<Channel> out;
    if (iq == nullptr || n == 0) return out;
    BandExtractor bxs(cfg_);
    // Spektrale Detektion auf dem ganzen Puffer: die Kanaele legen sich um
    // die Bereiche, in denen Energie detektiert wurde.
    const BandExtractor::BandDetection det = bxs.bands(iq, 0, n);
    const std::vector<Band> ebs = energyBands(det);
    out.reserve(ebs.size());
    for (const Band& cb : ebs)
        out.push_back(extract(iq, n, cb, det.nft, decimate));
    return out;
}

Channelizer::Channel Channelizer::extract(const std::complex<double>* iq,
                                          std::size_t n, const Band& band,
                                          std::size_t nft,
                                          std::size_t decimate) const {
    Channel ch;
    ch.band = band;
    ch.nft = nft;
    if (iq == nullptr || n == 0 || nft == 0) return ch;
    BandExtractor bxs(cfg_);
    // Bandpass + zyklische Verschiebung auf die Kanalmitte + IFFT.
    const BandExtractor::Slice slice = bxs.extract(iq, 0, n, band, nft);
    ch.centerNorm = wrap01(slice.fc);
    ch.widthNorm = slice.bw;
    // Resample-Faktor durch die Kanalbreite begrenzen (Nyquist des Kanals,
    // inkl. Waechter der Extraktionsmaske).
    std::size_t D = decimate < 1 ? 1 : decimate;
    if (ch.widthNorm > 0.0) {
        const double maxD = std::floor(1.0 / ch.widthNorm);
        const std::size_t cap =
            maxD >= 1.0 ? static_cast<std::size_t>(maxD) : 1;
        if (D > cap) D = cap;
    }
    ch.decimate = D;
    ch.y = resample(slice.y, D);
    return ch;
}

std::vector<std::complex<double>> Channelizer::resample(
    const std::vector<std::complex<double>>& y, std::size_t D) {
    if (D <= 1 || y.size() < 2) return y;
    // Anti-Alias-Tiefpass: Hann-gefensterte Sinc, Grenzfrequenz 0.45/D
    // (auf die volle Rate normiert), ungerade symmetrische Laenge -> zero
    // delay, DC-Verstaerkung auf 1 normiert.
    const std::size_t taps = 8 * D + 1;
    const std::size_t M = taps / 2;
    const double fc = 0.45 / static_cast<double>(D);
    std::vector<double> h(taps);
    double sum = 0.0;
    for (std::size_t i = 0; i < taps; ++i) {
        const double t = static_cast<double>(i) - static_cast<double>(M);
        const double s = t == 0.0 ? 2.0 * fc
                                  : std::sin(2.0 * kPi * fc * t) / (kPi * t);
        const double w =
            0.5 - 0.5 * std::cos(2.0 * kPi * static_cast<double>(i) /
                                 static_cast<double>(taps - 1));
        h[i] = s * w;
        sum += h[i];
    }
    if (sum != 0.0)
        for (double& v : h) v /= sum;
    std::vector<std::complex<double>> out((y.size() + D - 1) / D);
    for (std::size_t m = 0; m < out.size(); ++m) {
        const long center = static_cast<long>(m * D);
        std::complex<double> acc(0.0, 0.0);
        for (std::size_t i = 0; i < taps; ++i) {
            const long idx = center + static_cast<long>(i) -
                             static_cast<long>(M);
            if (idx >= 0 && idx < static_cast<long>(y.size()))
                acc += h[i] * y[static_cast<std::size_t>(idx)];
        }
        out[m] = acc;
    }
    return out;
}

}  // namespace iqclassify
