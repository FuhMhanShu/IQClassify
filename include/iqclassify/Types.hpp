// IQClassify - Modulationsklassifikation von Pulsen in IQ-Signalen.
//
// Oeffentliche Basistypen: Modulationsart, grobe Klasse und Ergebnisstruktur
// eines klassifizierten Pulses. Alle Frequenzen sind auf die Abtastrate
// normiert (fs = 1); Frequenzen werden als fraktionale Frequenz relativ zum
// Mittelpunkt des erfassten Bandes bzw. als [0,1)-antizyklische Position im
// Spektrum angegeben.

#pragma once

#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>

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
    FrankCode = 10,  // polyphase Frank-Code (Radar): Treppe der
                     // Momentanfrequenz mit M Stufen (Ordnung M, Laenge M^2)
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

// PSK-Synchronisation eines erkannten PSK-Pulses. Alle Schaetzungen sind
// batchweise ueber den ganzen Puls (Spektrallinien + gewichtete Regression),
// es gibt keine Regelschleife (kein PLL/Clock-Synchronizer) und damit keinen
// Einschwingvorgang; die Symbole am Pulsanfang sind sofort korrekt.
struct PskSync {
    bool valid = false;
    double sps = 0.0;                  // Samples pro Symbol
    double timingOffset = 0.0;         // erste Symbolzeit relativ zu
                                       // startSampleFine, in [0, sps)
    double carrierFreqNorm = 0.0;      // Traegerfrequenz des Pulses, zyklisch [0,1)
    double carrierPhaseAtStart = 0.0;  // Traegerphase am Pulsanfang
                                       // (startSampleFine), rad, auf das M-PSK-
                                       // Raster aufgeloest (Restambiguität 2pi/M)
    double snrDb = 0.0;
    std::vector<std::complex<double>> symbols;  // derotierte, normierte Symbole
};

// Alle Parameter eines erkannten Chirp-Pulses.
struct ChirpParams {
    bool valid = false;
    std::size_t startSample = 0;   // Pulsanfang im Puffer (startSampleFine)
    std::size_t numSamples = 0;    // Pulslaenge (numSamplesFine)
    double fStartNorm = 0.0;       // Frequenz am Pulsanfang, zyklisch [0,1)
    double fEndNorm = 0.0;         // Frequenz am Pulsende, zyklisch [0,1)
    double centerFreqNorm = 0.0;   // Mittenfrequenz
    double bandwidthNorm = 0.0;    // |fEnd - fStart| (Hub)
    double slopeNorm = 0.0;        // Frequenzhub pro Sample^2 (Vorzeichen =
                                   // Richtung; in Hz/s: slopeNorm * fs^2)
    double snrDb = 0.0;
    double r2 = 0.0;               // Bestimmtheitsmass des linearen Chirp-Fits
    bool up = true;                // true = aufwaerts, false = abwaerts
};

// Alle Parameter eines erkannten Frank-Code-Pulses (polyphase Radar-Code).
struct FrankParams {
    bool valid = false;
    std::size_t startSample = 0;   // Pulsanfang im Puffer (startSampleFine)
    std::size_t numSamples = 0;    // Pulslaenge (numSamplesFine)
    int order = 0;                 // Ordnung M (Code-Laenge N = M^2)
    std::size_t codeLength = 0;    // N = M^2
    double chipRateNorm = 0.0;     // Chip-Rate fraktional (1/(sps_chip))
    double centerFreqNorm = 0.0;   // Mittenfrequenz des Pulses, zyklisch [0,1)
    double bandwidthNorm = 0.0;    // belegtes Band (Sweep der Treppe)
    double snrDb = 0.0;
    double r2 = 0.0;               // R^2 der Stufen-Mitten gegen Stufen-Index
    std::vector<double> codePhases;  // geschaetzte Code-Phasen, N Werte in
                                     // [0, 2pi), Chip k = i*M+j -> Phase
                                     // ~ 2*pi*i*j/M (mod 2pi)
};

// Symbol-Schnittstelle fuer PSK- und FSK-Pulse (Symbolextraktion).
// PSK (BPSK/QPSK): symbols gefuellt (derotiert auf das M-PSK-Raster,
// carrierFreqNorm/carrierPhaseAtStart synchronisieren den Pulsanfang).
// FSK (FSK2): symbolFreqs gefuellt (Mittenfrequenz pro Symbol relativ zur
// Pulsmitte) und fskTone (0/1-Entscheidung je Symbol).
struct SymbolStream {
    bool valid = false;
    ModulationType type = ModulationType::Unknown;
    double sps = 0.0;                  // Samples pro Symbol
    double timingOffset = 0.0;         // erste Symbolzeit relativ zu startSampleFine
    double carrierFreqNorm = 0.0;      // PSK: Traeger [0,1); FSK: Pulsmitte
    double carrierPhaseAtStart = 0.0;  // PSK: Phase am Pulsanfang (rad)
    std::vector<std::complex<double>> symbols;  // PSK: Symbole
    std::vector<double> symbolFreqs;            // FSK: Frequenz pro Symbol
    std::vector<int> fskTone;                   // FSK2: Ton 0/1 pro Symbol
};

}  // namespace iqclassify
