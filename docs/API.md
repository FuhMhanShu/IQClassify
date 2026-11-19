# Programmierschnittstelle (API)

Namespace `iqclassify`. Schirmheader:

    #include "iqclassify/iqclassify.hpp"

## Aufzählungen

    enum class ModulationType : std::uint8_t {
        Unknown,        // nicht entscheidbar (zu wenig SNR/Struktur)
        Unmodulated,    // unmodulierter Träger (CW)
        BPSK, QPSK, FSK2, Chirp, OFDM,
        QAM16, QAM64, QAM256
    };

    enum class CoarseClass : std::uint8_t {
        Unknown,            // interner Fehlwert
        Tone,               // tonartig
        FrequencyModulated, // Chirp / FSK
        PhaseModulated,     // PSK / QAM
        NoiseLike           // flächig, hohes PAPR (OFDM-Gruppe)
    };

Diagnosetexte: `const char* toString(ModulationType)`,
`const char* toString(CoarseClass)`.

## PulseResult

Ergebnis je erkanntem Puls:

| Feld | Typ | Bedeutung |
|---|---|---|
| startSample | size_t | erstes Sample des Pulses im Puffer (Blockraster) |
| numSamples | size_t | Pulslänge in Samples (Blockraster) |
| startSampleFine | size_t | samplegenauer Pulsbeginn im Puffer (Wiederfinden) |
| numSamplesFine | size_t | samplegenaue Pulslänge (Wiederfinden) |
| centerFreqNorm | float | Bandmitte auf dem FFT-Raster, fraktional, zyklisch in [0,1) |
| centerFreqRefined | double | Sub-Bin-Mittenfrequenz (In-Band-Schwerpunkt), zyklisch in [0,1) |
| bandwidthNorm | float | Bandbreite inkl. Masken-Wächter, fraktional |
| snrDb | float | geschätzte Spektraldichte-SNR des Pulses |
| coarse | CoarseClass | Stufe-1-Entscheidung |
| modulation | ModulationType | Stufe-2-Entscheidung |
| confidence | float | Konfidenz der Stufe-2-Entscheidung in [0,1] |
| envelopeCv | float | Variationskoeffizient der Hüllkurve |
| paprWindowedDb | float | PAPR, Fensterlänge 16 |
| chirpR2 | float | Chirp-Fit-Bestimmtheitsmaß |
| chirpSweep | float | Chirp-Hub relativ zur Bandbreite |
| dwellFraction | float | Anteil der Samples mit hohem |f_inst| |
| freqModes | int | Anzahl dominanter f_inst-Moden |
| spectralFlatnessDb | float | Standardabweichung der in-band-PSD |
| bandBins | size_t | Breite des detektierten Bandes (FFT-Bins) |

## ModulationClassifier

    class ModulationClassifier {
    public:
        explicit ModulationClassifier(const ClassifierConfig& cfg);
        std::vector<PulseResult> classify(const std::complex<double>* iq,
                                          std::size_t n) const;
        std::vector<PulseResult> classifyF(const std::complex<float>* iq,
                                           std::size_t n) const;
        std::vector<std::complex<double>> extractPulse(
            const std::complex<double>* iq, std::size_t n,
            const PulseResult& pulse) const;
        std::vector<std::complex<float>> extractPulseF(
            const std::complex<float>* iq, std::size_t n,
            const PulseResult& pulse) const;
        const ClassifierConfig& config() const;
    };

- `classify` / `classifyF` durchsuchen den gesamten Puffer und liefern ein
  Ergebnis je erkanntem Puls, zeitlich sortiert.
- Härtung: Nullzeiger und Länge 0 liefern ein leeres Ergebnis. Nicht
  endliche Samples (NaN, Inf) werden an jeder Stelle der Pipeline als
  Ausfall (Null) behandelt und vergiften dadurch weder FFT-Fenster noch
  Statistiken; ein einzelner NaN-Dropout verhält sich wie eine kurze
  Schweigephase im Puls (siehe docs/METHODOLOGY.md, Abschnitt Härtung).
- Die Methoden sind `const` und thread-sicher, solange jeder Thread eine
  eigene Instanz nutzt (PocketFFT-Plancache ist global und
  thread-sicher ab modifizierten Versionen; Aufrufe erfolgen
  single-threaded).
- `classifyF` konvertiert die Eingabe nach double; die interne Rechnung
  bleibt doppelpräzise.
- Lokalisierung (Wiederfinden im IQ-Signal): jedes `PulseResult` trägt
  neben dem groben Blockraster-Fenster (`startSample`/`numSamples`) das
  samplegenaue Zeitfenster (`startSampleFine`/`numSamplesFine`) und die
  Sub-Bin-Mittenfrequenz (`centerFreqRefined`). Umrechnung in
  physikalische Einheiten: `t = startSampleFine / fs`,
  `f = centerFreqRefined * fs` (zyklisch; negative Frequenzen erscheinen
  als Werte nahe 1).
- `extractPulse`/`extractPulseF` liefern den Basisband-Slice eines
  klassifizierten Pulses aus dem Puffer: Zeitfenster
  `[startSampleFine, startSampleFine + numSamplesFine)`, Bandpass um
  `centerFreqRefined` mit `bandwidthNorm`. `slice[k]` entspricht
  `iq[startSampleFine + k]`, frequenzverschoben auf die Pulsmitte
  (Rückverschiebung: `slice[k] * e^{+j*2pi*centerFreqRefined*k}`).
  Nullzeiger, leere Eingabe oder ein Puls außerhalb des Puffers liefern
  ein leeres Ergebnis.

## ClassifierConfig

Alle Schwellwerte der Referenzvalidierung, gruppiert wie in
docs/METHODOLOGY.md (Pulssegmentierung, Banddetektion, SNR/Stufe 1,
OFDM-Präfix, Symboltiming, Phasenanalyse, QAM-Ordnung). Jedes Feld ist
in include/Config.hpp individuell kommentiert. Standardwerte entsprechen der
validierten Kalibrierung; Änderungen erfordern eine erneute Validierung
gegen die Testsuite.

## Beispiel

    #include "iqclassify/iqclassify.hpp"
    #include <vector>

    std::vector<std::complex<float>> iq = /* Abtaststrom, fs = 1 */;

    iqclassify::ClassifierConfig cfg;          // Defaults = Kalibrierung
    iqclassify::ModulationClassifier clf(cfg);

    for (const auto& p : clf.classifyF(iq.data(), iq.size())) {
        std::printf("Pulse @%zu: %s (SNR %.1f dB, Konfidenz %.2f)\n",
                    p.startSample,
                    iqclassify::toString(p.modulation),
                    p.snrDb, p.confidence);
        // Puls zeitlich/frequenztechnisch wiederfinden:
        std::printf("  fein: Samples [%zu, %zu), f = %.6f * fs\n",
                    p.startSampleFine,
                    p.startSampleFine + p.numSamplesFine,
                    p.centerFreqRefined);
        auto slice = clf.extractPulseF(iq.data(), iq.size(), p);
        // slice[k] == iq[startSampleFine + k], auf die Pulsmitte verschoben
    }

## Schwellwerte anpassen

- Schwaches SNR-Regime: minSnrDb senken, gleichzeitig tCpPeak/tRsProm
  herabsetzen und die QAM-Dichte-Fallback-Schwelle qamDenseSnr an das
  reale SNR-Regime anpassen.
- Kürzere Pulse als 400 Samples: minPulseLen senken und die Testsuite
  um derartige Szenarien erweitern (die Kalibrierung gilt für die
  getesteten Pulslängen).
- Andere Oversamplingfaktoren: Der Symboltiming-Kandidatenbereich ist
  generisch (1.6..n/12 Samples/Symbol); die Kalibrierung prüfte 4/8/16.
