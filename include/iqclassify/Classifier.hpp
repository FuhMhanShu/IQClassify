// IQClassify - Stufe 1 (grobe Klasse) und Stufe 2 (Modulationsart).
//
// Stufe 1 ordnet jeden (Band, Zeitcrop)-Abschnitt einer groben Klasse zu:
//   Tone                -> unmodulierter Traeger
//   FrequencyModulated  -> Chirp / FSK
//   PhaseModulated      -> PSK / QAM (konstellationsbasiert)
//   NoiseLike           -> OFDM-artig flaechiges Signal
//
// Stufe 2 entscheidet innerhalb der Klasse: OFDM ueber zyklisches
// Praefix (Autokorrelationsspitze + Talreinheit + Halbseiten-
// Konsistenz), PSK ueber Symboltiming, Resttraeger (z^4-Kohaerenz),
// Phasenfaltung und Kumulanten, QAM-Ordnung ueber einen rotations-
// invarianten Rasterfit der analytischen Levelgitter.

#pragma once

#include <cstddef>
#include <vector>

#include "BandExtractor.hpp"
#include "Config.hpp"
#include "FeatureExtractor.hpp"
#include "Types.hpp"

namespace iqclassify {

class ModulationClassifier {
public:
    explicit ModulationClassifier(const ClassifierConfig& cfg)
        : cfg_(cfg) {}

    // Klassifiziert alle Pulse im IQ-Puffer. Rueckgabe: ein Ergebnis je
    // erkanntem Puls, zeitlich sortiert.
    std::vector<PulseResult> classify(const std::complex<double>* iq,
                                      std::size_t n) const;

    // Komfortvariante fuer float-Eingaenge (interne Rechnung bleibt
    // in Doppelpraezision, wie in der Validierungsreferenz).
    std::vector<PulseResult> classifyF(const std::complex<float>* iq,
                                      std::size_t n) const;

    // Liefert den Basisband-Slice eines bereits klassifizierten Pulses
    // aus dem IQ-Puffer (zum Wiederfinden des Pulses im Signal):
    // Zeitfenster [startSampleFine, startSampleFine + numSamplesFine),
    // Bandpass um centerFreqRefined mit bandwidthNorm. Rueckgabe[k]
    // entspricht iq[startSampleFine + k], frequenzverschoben auf die
    // Mittenfrequenz des Pulses (Rueckverschiebung mit
    // e^{+j*2pi*centerFreqRefined*k}). Nullzeiger/leere Eingabe oder ein
    // Puls ausserhalb des Puffers liefern ein leeres Ergebnis.
    std::vector<std::complex<double>> extractPulse(
        const std::complex<double>* iq, std::size_t n,
        const PulseResult& pulse) const;

    // float-Variante von extractPulse.
    std::vector<std::complex<float>> extractPulseF(
        const std::complex<float>* iq, std::size_t n,
        const PulseResult& pulse) const;

    const ClassifierConfig& config() const { return cfg_; }

private:
    ClassifierConfig cfg_;
};

}  // namespace iqclassify
