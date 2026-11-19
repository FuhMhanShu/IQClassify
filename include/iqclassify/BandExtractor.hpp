// IQClassify - Banddetektion und Basisbandextraktion.
//
// Banddetektion: ueber gleitende Hann-gefensterte Teilfenster wird das
// per-Bin-Maximum der Leistungsdichte gebildet; daraus werden zusammen-
// haengende Bereiche ueber dem Rauschboden extrahiert (zyklische Koordina-
// ten, lo darf negativ und hi > nft sein). Spektrale Satelliten der
// Eigenstreuung (Keying-Splatter bei Phasen-/Frequenzmodulation) werden
// zur emittierenden Band zurueckverschmolzen.
//
// Extraktion: ungefensterte FFT des ganzen Segments, Bandmaske mit
// Waechter, zyklische Verschiebung zur DC und IFFT; das Ergebnis ist der
// komplexe Basisband-Slice des Segments (Zeitbereich), aus dem der Aufrufer
// die Pulse cropt.

#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "Config.hpp"
#include "Dsp.hpp"
#include "Types.hpp"

namespace iqclassify {

// Band in virtuellen Bin-Koordinaten: gueltig sind lo < hi mit beliebigen
// ganzzahligen Werten; der Bezug auf die FFT ist (k mod nft). Negative lo
// entstehen zyklisch (Band ueberschreitet die Bin 0).
struct Band {
    long lo = 0;
    long hi = 0;
    long width() const { return hi - lo; }
};

class BandExtractor {
public:
    explicit BandExtractor(const ClassifierConfig& cfg) : cfg_(cfg) {}

    // Detektionsergebnis: Bandaender (zyklische virtuelle Koordinaten,
    // aufsteigend ab dem stillsten Bin), Laenge der Detektions-FFT und das
    // dB-Beugungs-PSD (Max-Hold ueber Hann-Teilfenster).
    struct BandDetection {
        std::vector<Band> bands;
        std::size_t nft = 0;
        std::vector<double> psdDb;  // Laenge nft, in dB
    };

    // floorRef: Rauschboden des Puffers pro Sample (z. B. PulseSegmenter::
    // noiseFloor) als Referenz fuer das Energie-Gate; < 0 deaktiviert das
    // Gate (alle Fenster werden gerechnet).
    BandDetection bands(const std::complex<double>* x, std::size_t a,
                        std::size_t b, double floorRef = -1.0) const;

    struct Slice {
        std::vector<std::complex<double>> y;  // Basisband-Slice des Segments
        double fc = 0.0;                       // Bandmitte, fraktional (kann
                                               // ausserhalb [0,1) liegen)
        double bw = 0.0;                       // Maskenbreite inkl. Waechter
    };

    // Extrahiert den Basisband-Slice eines Bandes aus dem Segment [a, b).
    Slice extract(const std::complex<double>* x, std::size_t a,
                  std::size_t b, const Band& band, std::size_t nft) const;

private:
    const ClassifierConfig& cfg_;
};

// Spektraldichte-SNR und in-band-Flachheit eines Bandes aus der
// Detektions-PSD. allBands schliesst die Bloecke anderer Baender (+- Breite)
// aus der Aussen-Referenz aus; in Mehrpuls-Szenarien ist das notwendig,
// weil fremde Baender keine Rauschreferenz sind.
void estimateSnrFlat(const std::vector<double>& psdDb, const Band& band,
                     const ClassifierConfig& cfg,
                     const std::vector<Band>& allBands, double* snrDb,
                     double* flatDb);

}  // namespace iqclassify
