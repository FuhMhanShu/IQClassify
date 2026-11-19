// IQClassify - Merkmalsextraktion und zeitliche Crops.
//
// Merkmale eines Basisband-Slices: Hoellkurvenstatistik (CV, PAPR,
// fenstergemitteltes PAPR), Traegerkohaerenz, momentane Frequenz (gepoolt
// um den Resttraeger), Chirp-Fit (hoellkurven^2-gewichtet, mittlere
// 80 % des Pulses), Dwell-Statistik und die dominanten f_inst-Moden
// (Histogramm der hohen Hoellkurve).

#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "Config.hpp"
#include "Dsp.hpp"
#include "Types.hpp"

namespace iqclassify {

// Merkmalsvektor eines Pulses (Stufe 1/2 arbeiten darauf).
struct Features {
    std::size_t n = 0;
    double fc = 0.0;
    double bw = 0.0;
    std::size_t bandBins = 0;
    double snrDb = 0.0;
    double flatDb = 0.0;

    std::vector<double> env;    // Hoellkurve |y|
    std::vector<double> ph;     // kontinuierliche Phase ohne
                                // Resttraeger (f.ph entspricht der
                                // akkumulierten Phase)
    std::vector<double> finst;  // momentane Frequenz (n-1 Werte),
                                // Resttraeger entfernt
    double residF = 0.0;        // entfernte Resttraegerfrequenz
    double carrCoh = 0.0;        // |mean(z)| / mean(|z|)
    double cv = 0.0;
    double paprDb = 0.0;
    double paprWDb = 0.0;
    double chirpR2 = 0.0;
    double chirpSweep = 0.0;
    double dwell = 0.0;
    double dwellThr = 0.0;
    int nfreqModes = 0;
    double modeConc = 0.0;
};

Features extractFeatures(const std::vector<std::complex<double>>& y,
                         double fc, double bw, std::size_t bandBins,
                         double snrDb, double flatDb,
                         const ClassifierConfig& cfg);

// Zeitliche On/Off-Crops eines Band-Slices (mehrere Pulse koennen sich ein
// Frequenzband teilen). Liefert [i0, i1)-Indizes in y.
std::vector<std::pair<std::size_t, std::size_t>> timeCrops(
    const std::vector<std::complex<double>>& y, const ClassifierConfig& cfg);

}  // namespace iqclassify
