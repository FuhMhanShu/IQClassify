// IQClassify - Modulationsklassifikation von Pulsen in IQ-Signalen.
//
// Oeffentliche Basistypen: Modulationsart, grobe Klasse und Ergebnisstruktur
// eines klassifizierten Pulses. Alle Frequenzen sind auf die Abtastrate
// normiert (fs = 1); Frequenzen werden als fraktionale Frequenz relativ zum
// Mittelpunkt des erfassten Bandes bzw. als [0,1)-antizyklische Position im
// Spektrum angegeben.

#pragma once

#include <cstddef>
#include <cstdint>

namespace iqclassify {

// Feinklassifikation (Stufe 2)
enum class ModulationType : std::uint8_t {
    Unknown = 0,      // nicht entscheidbar (zu wenig SNR / Struktur)
    Unmodulated = 1,  // unmodulierter Traeger (CW)
    BPSK = 2,
    QPSK = 3,
    FSK2 = 4,
    Chirp = 5,
    OFDM = 6,
    QAM16 = 7,
    QAM64 = 8,
    QAM256 = 9,
};

// Grobklassifikation (Stufe 1)
enum class CoarseClass : std::uint8_t {
    Unknown = 0,           // (interner Fehlwert)
    Tone = 1,              // tonartig -> Unmodulated-Pfad
    FrequencyModulated,    // Chirp / FSK
    PhaseModulated,        // PSK / QAM (konstellationsbasiert)
    NoiseLike,             // flaechig, hohes PAPR -> OFDM
};

// Konvertierung der Aufzaehlung in Text (Diagnose / Logging).
const char* toString(ModulationType t);
const char* toString(CoarseClass c);

// Ergebnis eines erkannten und klassifizierten Pulses.
struct PulseResult {
    std::size_t startSample = 0;   // erstes Sample des Pulses im Puffer
    std::size_t numSamples = 0;   // Pulslaenge in Samples
    float centerFreqNorm = 0.0f;  // Bandmitte, fraktional, in [0,1) (zyklisch)
    float bandwidthNorm = 0.0f;   // Bandbreite inkl. Masken-Waechter, fraktional
    float snrDb = 0.0f;           // geschaetzte Spektraldichte-SNR des Pulses

    // Feine Lage des Pulses im IQ-Puffer (zum Wiederfinden):
    // startSampleFine/numSamplesFine sind samplegenau (Huellkurven-Schwelle
    // im Crop), centerFreqRefined ist die Sub-Bin-Mittenfrequenz
    // (leistungsgewichteter Schwerpunkt des In-Band-PSD), zyklisch [0,1).
    // Umrechnung: t = startSampleFine / fs, f = centerFreqRefined * fs.
    std::size_t startSampleFine = 0;  // samplegenauer Pulsbeginn im Puffer
    std::size_t numSamplesFine = 0;   // samplegenaue Pulslaenge
    double centerFreqRefined = 0.0;   // Sub-Bin-Mittenfrequenz, zyklisch [0,1)
    CoarseClass coarse = CoarseClass::Unknown;   // Stufe-1-Entscheidung
    ModulationType modulation = ModulationType::Unknown;  // Stufe-2-Entscheidung
    float confidence = 0.0f;      // Konfidenz der Stufe-2-Entscheidung [0,1]

    // Diagnosemerkmale (Stufe-1-/Stufe-2-Zwischengroessen, see docs/API.md)
    float envelopeCv = 0.0f;       // Variationskoeffizient der Hoellkurve
    float paprWindowedDb = 0.0f;   // PAPR, Fensterlaenge 16
    float chirpR2 = 0.0f;          // Chirp-Fitbestimmtheitsmass
    float chirpSweep = 0.0f;       // Chirp-Hub relativ zur Bandbreite
    float dwellFraction = 0.0f;   // Anteil der Samples mit hohem |f_inst|
    int freqModes = 0;             // Anzahl dominanter f_inst-Moden
    float spectralFlatnessDb = 0.0f; // Standardabweichung der in-band-PSD
    std::size_t bandBins = 0;      // Breite des detektierten Bandes (FFT-Bins)
};

}  // namespace iqclassify
