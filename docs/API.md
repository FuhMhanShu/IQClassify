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
        PskSync synchronizePsk(const std::complex<double>* iq, std::size_t n,
                               const PulseResult& pulse) const;
        PskSync synchronizePskF(const std::complex<float>* iq, std::size_t n,
                                const PulseResult& pulse) const;
        ChirpParams analyzeChirp(const std::complex<double>* iq, std::size_t n,
                                 const PulseResult& pulse) const;
        ChirpParams analyzeChirpF(const std::complex<float>* iq, std::size_t n,
                                  const PulseResult& pulse) const;
        SymbolStream extractSymbols(const std::complex<double>* iq,
                                    std::size_t n,
                                    const PulseResult& pulse) const;
        SymbolStream extractSymbolsF(const std::complex<float>* iq,
                                     std::size_t n,
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

## Analyse-API (PSK-Synchronisation, Chirp-Parameter, Symbole)

Drei Methoden verfeinern einen bereits klassifizierten Puls. Alle nehmen
den Originalpuffer plus das `PulseResult` und liefern bei Erfolg
`valid == true`; bei falscher Modulationsart, Nullzeiger oder zu kurzem
Puls ist das Ergebnis definiert `invalid`.

`PskSync synchronizePsk(F)` — nur für BPSK/QPSK. Synchronisiert den Puls
phasenrichtig auf seinen Anfang: Batch-Schätzung über den gesamten Puls
(Symbolratenlinie + RRC-Matched-Filter + Träger-Faltung z^M), daher kein
Einschwingvorgang ("Einschwingheizen") eines Clock-Synchronizers.

| Feld | Bedeutung |
|---|---|
| sps | Samples pro Symbol (fein, z. B. 8.000) |
| timingOffset | erste Symbolzeit relativ zu `startSampleFine`, in [0, sps) |
| carrierFreqNorm | Trägerfrequenz des Pulses, zyklisch in [0,1) |
| carrierPhaseAtStart | Trägerphase am Pulsanfang (`startSampleFine`), rad, auf das M-PSK-Raster aufgelöst (Restambiguität 2π/M) |
| snrDb | geschätzte SNR des Pulses |
| symbols | derotierte, normierte Symbole (Betrag ~1, auf das M-PSK-Raster gedreht) |

`ChirpParams analyzeChirp(F)` — nur für Chirp. Lineare Regression der
Momentanfrequenz über den Puls (hüllkurven²-gewichtet).

| Feld | Bedeutung |
|---|---|
| startSample / numSamples | Pulsfenster (`startSampleFine`/`numSamplesFine`) |
| fStartNorm / fEndNorm | Frequenz am Pulsanfang/-ende, zyklisch in [0,1) |
| centerFreqNorm | Mittenfrequenz (fStart+fEnd)/2 |
| bandwidthNorm | Hub |fEnd-fStart| |
| slopeNorm | Frequenzhub pro Sample (Vorzeichen = Richtung; in Hz/s: slopeNorm * fs²) |
| snrDb, r2 | SNR und Bestimmtheitsmaß des Chirp-Fits |
| up | true = aufwärts, false = abwärts |

`SymbolStream extractSymbols(F)` — Symbol-Schnittstelle für die spätere
Symbolgewinnung; PSK-Signale liefern `symbols` (wie `synchronizePsk`),
FSK2 liefert `symbolFreqs` (Momentanfrequenz pro Symbolzentrum) und
`fskTone` (Ton-Entscheidung 0/1). `type` ist BPSK/QPSK bzw. FSK2.

| Feld | Bedeutung |
|---|---|
| sps / timingOffset | Symboltakt und erste Symbolzeit relativ zu `startSampleFine` |
| carrierFreqNorm | PSK: Träger [0,1); FSK2: Pulsmitte |
| carrierPhaseAtStart | PSK: Trägerphase am Pulsanfang (rad); FSK2: 0 |
| symbols | PSK: derotierte Symbole |
| symbolFreqs / fskTone | FSK2: Frequenz bzw. Ton (0/1) pro Symbol |

## Channelizer (spektrale Kanalisierung)

Erster Pipelineschritt: das Spektrum wird in Kanäle zerlegt, die sich **nicht
gleichverteilt, sondern um die detektierten Energiebereiche** legen — eine
Banddetektion auf dem ganzen Puffer findet die Energiebereiche, und um jeden
Bereich wird ein Kanal gelegt (zentriert auf dessen leistungsgewichteten
Schwerpunkt, so breit wie der Bereich inkl. Waechter).

    class Channelizer {
    public:
        struct Channel {
            std::vector<std::complex<double>> y;  // Kanal-Basisband
            double centerNorm;   // Kanalmitte, zyklisch [0,1) der Originalrate
            double widthNorm;    // Kanalbreite inkl. Waechter (Originalrate)
            Band band;           // Kanalband in nft-Bin-Koordinaten
            std::size_t nft;     // FFT-Laenge des Puffers
            std::size_t decimate;// Resample-Faktor des Kanals
        };
        std::vector<Channel> channelize(const std::complex<double>* iq,
                                        std::size_t n,
                                        std::size_t decimate = 1) const;
        Channel extract(const std::complex<double>* iq, std::size_t n,
                        const Band& band, std::size_t nft,
                        std::size_t decimate = 1) const;
        static std::vector<std::complex<double>> resample(
            const std::vector<std::complex<double>>& y, std::size_t D);
    };

- `channelize` liefert alle Kanäle: jeder Kanal ist extrahiert (Bandpass per
  FFT-Maske + IFFT), auf seine Mitte heruntergemischt und optional resamplet.
  `decimate` wird pro Kanal auf `floor(1/Kanalbreite)` begrenzt (Nyquist).
- `extract` liefert einen einzelnen Kanal; `resample` ist der Anti-Alias-
  Tiefpass (Hann-Sinc) + Dezimierung (symmetrisch, zero delay).
- Bei `decimate = D` hat `ch.y` die Länge `ceil(n/D)`; `ch.y[m]` entspricht
  `iq[m*D]`, frequenzverschoben auf `ch.centerNorm`. Umkehrung:
  `iq[m*D] ≈ ch.y[m] * e^{+j*2pi*centerNorm*m*D}`.
- `ModulationClassifier::classify` nutzt dieselbe Kanalisierung intern als
  ersten Schritt (`ClassifierConfig::channelize`, Default an) und rechnet
  alle Ergebnisse in die Koordinaten des Originalsignals zurück
  (Zeit ×D, Frequenz = Kanalmitte + lokal/D zyklisch). Mit
  `channelDecimate > 1` werden die Kanäle in `classify` vor der
  Merkmalsextraktion dezimiert (sample-basierte Schwellwerte werden skaliert,
  Ergebnisse zurückgerechnet).

- Sehr kurze Pulse (~180 Samples, z. B. 13 BPSK-Bit bei sps=8): Das
  Upsampling sitzt auf dem extrahierten Basisband-Slice in `processSlice`
  (`ClassifierConfig::upsample`: 0 = automatisch fuer Slice < 1024 Samples,
  1 = aus, >1 = fester Faktor; Zeiten der Ergebnisse werden /U gerechnet,
  Frequenzen bleiben in der Originalrate). Kombiniert mit
  `ClassifierConfig::shortPulseProfile()` (feinere Segmentierung,
  `minPulseLen=96`, `timingMinSyms=8`, angepasste Symbolratensuche) werden
  solche Pulse zuverlaessig detektiert und lokalisiert (Test 20).

## ClassifierConfig

Alle Schwellwerte der Referenzvalidierung, gruppiert wie in
docs/METHODOLOGY.md (Kanalisierung, Parallelverarbeitung (`threads`: 0 =
automatisch, 1 = seriell), Pulssegmentierung, Banddetektion, SNR/Stufe 1,
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

- Schwaches SNR-Regime: Default minSnrDb = 2 dB deckt Detektion ab ~2 dB
  (CW/Chirp) bzw. ~3-4 dB (breitbandige PSK/QAM) ab. Noch darunter:
  minSnrDb weiter senken, gleichzeitig tCpPeak/tRsProm herabsetzen und
  die QAM-Dichte-Fallback-Schwelle qamDenseSnr an das reale SNR-Regime
  anpassen.
- Kürzere Pulse als 400 Samples: minPulseLen senken und die Testsuite
  um derartige Szenarien erweitern (die Kalibrierung gilt für die
  getesteten Pulslängen).
- Andere Oversamplingfaktoren: Der Symboltiming-Kandidatenbereich ist
  generisch (1.6..n/12 Samples/Symbol); die Kalibrierung prüfte 4/8/16.
