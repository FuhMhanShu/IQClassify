// IQClassify - Stufe 0: Pulssegmentierung (Energie im Block).
//
// Der IQ-Puffer wird in Bloecke fester Laenge zerlegt; der Rauschboden wird
// per Minimum-Statistik ueber Gruppen von Bloecken geschaetzt (Median je
// Gruppe, Minimum ueber alle Gruppen). Hysterese-Schwellwerte ueber dem
// Boden markieren Pulsbeginn und Pulsende; jeder Puls erhaelt einen
// Sicherheitsrand.

#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "Config.hpp"
#include "Dsp.hpp"
#include "Types.hpp"

namespace iqclassify {

class PulseSegmenter {
public:
    explicit PulseSegmenter(const ClassifierConfig& cfg) : cfg_(cfg) {}

    // Liefert [start, end)-Intervalle (in Samples) aller Pulse im Puffer.
    std::vector<std::pair<std::size_t, std::size_t>> detect(
        const std::complex<double>* x, std::size_t n) const;

    // Rauschboden des Puffers (pro Sample, Minimum-Statistik ueber
    // Blockgruppen-Mediane; 0 wenn keiner bestimmbar). Wird z. B. als
    // Referenz fuer das Energie-Gate der Banddetektion genutzt.
    static double noiseFloor(const std::complex<double>* x, std::size_t n,
                             const ClassifierConfig& cfg);

private:
    const ClassifierConfig& cfg_;
};

}  // namespace iqclassify
