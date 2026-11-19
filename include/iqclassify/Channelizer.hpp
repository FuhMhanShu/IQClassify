// IQClassify - Spektrale Kanalisierung (erster Pipelineschritt).
//
// Das Spektrum des Puffers wird nicht gleichverteilt, sondern adaptiv in
// Kanaele zerlegt: Eine Banddetektion auf dem ganzen Puffer findet die
// Bereiche mit Energie, und um jeden solchen Bereich wird ein Kanal gelegt
// (zentriert auf den leistungsgewichteten Schwerpunkt des Bandes, so breit
// wie das Band inkl. Waechter). Jeder Kanal wird anschliessend einzeln
// untersucht: extrahiert (Bandpass), auf seine Mitte heruntergemischt und
// optional resamplet (dezimiert). Danach laeuft je Kanal die weitere
// Pipeline (Zeitsegmentierung, Crops, Merkmale, Klassifikation) auf dem
// Kanal-Basisband.

#pragma once

#include <complex>
#include <cstddef>
#include <vector>

#include "BandExtractor.hpp"
#include "Config.hpp"
#include "Dsp.hpp"
#include "Types.hpp"

namespace iqclassify {

class Channelizer {
public:
    struct Channel {
        std::vector<std::complex<double>> y;  // Kanal-Basisband (extrahiert,
                                              // gemischt, resamplet)
        double centerNorm = 0.0;  // Kanalmitte, zyklisch [0,1) der Originalrate
        double widthNorm = 0.0;   // Kanalbreite inkl. Waechter (Originalrate)
        Band band{0, 0};          // Kanalband in nft-Bin-Koordinaten
        std::size_t nft = 0;      // verwendete FFT-Laenge des Puffers
        std::size_t decimate = 1; // Resample-Faktor des Kanals
    };

    explicit Channelizer(const ClassifierConfig& cfg) : cfg_(cfg) {}

    // Kanalbaender zu einer Detektion: je detektiertem Energieband ein Band,
    // zentriert auf dessen leistungsgewichteten Schwerpunkt (die Kanaele
    // bewegen sich zu den Stellen im Spektrum, wo Energie detektiert wurde,
    // sie sind nicht gleichverteilt). Bei leerer Detektion ein Band ueber
    // das volle Spektrum (wie der Legacy-Pfad).
    std::vector<Band> energyBands(
        const BandExtractor::BandDetection& det) const;

    // Zerlegt den Puffer in Kanaele um die detektierten Energiebereiche
    // (energyBands) und extrahiert jeden Kanal (extract + Mischung +
    // optional Resample). decimate: Resample-Faktor (1 = volle Abtastrate);
    // wird pro Kanal auf <= floor(1/Kanalbreite) begrenzt (Nyquist des
    // Kanals). Nullzeiger/leere Eingabe -> leer.
    std::vector<Channel> channelize(const std::complex<double>* iq,
                                    std::size_t n,
                                    std::size_t decimate = 1) const;

    // Extrahiert einen einzelnen Kanal (Band in nft-Koordinaten) aus dem
    // Puffer: Bandpass + Mischung auf die Kanalmitte + optional Resample.
    // Rueckgabe.y[k] entspricht iq[k] * e^{-j*2pi*centerNorm*k} (bei
    // decimate > 1: y[m] entspricht iq[m*decimate] ... im Kanal-Basisband).
    Channel extract(const std::complex<double>* iq, std::size_t n,
                    const Band& band, std::size_t nft,
                    std::size_t decimate = 1) const;

    // Anti-Alias-Tiefpass (Hann-Sinc, Grenzfrequenz 0.45/D auf die volle
    // Rate normiert) + Dezimierung um D. Symmetrisches FIR -> kein
    // zusaetzlicher Phasengang/Versatz (zero delay).
    static std::vector<std::complex<double>> resample(
        const std::vector<std::complex<double>>& y, std::size_t D);

    // Upsampling um U: Zero-Stuffing + Anti-Image-Tiefpass (Hann-Sinc,
    // Grenzfrequenz 0.45/U auf die Ausgaberate normiert, Verstaerkung U,
    // symmetrisch -> zero delay). out.size() == y.size()*U; Frequenzen sind
    // rateninvariant, Zeiten skalieren mit U.
    static std::vector<std::complex<double>> upsample(
        const std::vector<std::complex<double>>& y, std::size_t U);

private:
    const ClassifierConfig& cfg_;
};

}  // namespace iqclassify
