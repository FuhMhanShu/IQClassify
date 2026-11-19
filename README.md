# IQClassify

Modulationsklassifikation von Pulsen in IQ-Signalen (C++17, statische
Bibliothek).

IQClassify durchsucht einen komplexen Abtaststrom (I/Q, Basisband) nach
Pulsen, trennt gleichzeitig auftretende Pulse in Frequenz und Zeit und
klassifiziert die Modulationsart jedes Pulses in zwei Stufen:

- Stufe 1 (grob): Ton, Frequenzmodulation, Phasenmodulation, flaechiges
  Signal.
- Stufe 2 (fein): Unmodulated (CW), BPSK, QPSK, FSK2, Chirp, OFDM,
  QAM16, QAM64, QAM256.

## Merkmale

- Keine Abhaengigkeiten ausser PocketFFT (header-only, vendored unter
  `third_party/pocketfft/`).
- Interne Rechnung in Doppelpraezision; Eingabe als `float` oder `double`.
- Mehrere Pulse pro Puffer: zeitlich getrennt und gleichzeitig in
  verschiedenen Frequenzbaendern.
- Lokalisierung: jeder detektierte Puls wird mit samplegenauer Zeitlage
  (`startSampleFine`/`numSamplesFine`) und Sub-Bin-Mittenfrequenz
  (`centerFreqRefined`) ausgegeben; `extractPulse` liefert den Puls als
  Basisband-Slice direkt aus dem IQ-Puffer (Wiederfinden im Signal).
- Klassenkapselung: eine Fassade (`ModulationClassifier`), alle
  Schwellwerte konfigurierbar (`ClassifierConfig`).
- Deterministische Testsuite ohne externes Testframework: 106 Szenarien
  aus der Qualitaetsmatrix (SNR-Raster, Oversampling, Frequenzlagen,
  Chirp-Varianten, Pulslaengen, Mehrpuls, Seed-Streuung) plus
  Haertungsfaelle (Nullzeiger, leere und NaN-beladene Eingaben,
  Dropouts im Puls).
- Robust gegen Amplitudenschwankungen (Fading): Ein Puls mit schwankender
  Amplitude bleibt durch Fade-Merge ein Puls (solange der Einbruch nicht
  bis auf den Rauschboden faellt), inkl. korrekter Klassifikation und
  Lokalisierung.
- Haertung der Eingabepfade: NaN/Inf-Samples verhalten sich wie
  Schweigephasen und vergiften weder FFT-Fenster noch Statistiken
  (siehe docs/METHODOLOGY.md).

## Bauen mit CMake

    cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build
    ctest --test-dir build        # Testsuite (Option IQCLASSIFY_BUILD_TESTS)

## Bauen mit qmake

    qmake IQClassify.pro
    make                          # erzeugt libIQClassify.a

Testsuite:

    cd tests
    qmake IQClassifyTests.pro && make && ./iqclassify_tests

Einbindung in ein qmake-Konsumentenprojekt:

    include(IQClassify.pri)

## Verwendung

    #include "iqclassify/iqclassify.hpp"

    iqclassify::ClassifierConfig cfg;              // Optionen/Schwellwerte
    iqclassify::ModulationClassifier clf(cfg);
    std::vector<std::complex<float>> iq = /* ... */;
    for (const auto& pulse : clf.classifyF(iq.data(), iq.size())) {
        // pulse.modulation, pulse.confidence, pulse.centerFreqNorm, ...
        // Puls zeitlich/frequenztechnisch wiederfinden (fs = 1 normiert):
        //   t0 = pulse.startSampleFine, len = pulse.numSamplesFine,
        //   f0 = pulse.centerFreqRefined
        auto slice = clf.extractPulseF(iq.data(), iq.size(), pulse);
        // slice[k] == iq[pulse.startSampleFine + k], auf f0 verschoben
    }

Alle Frequenzen sind auf die Abtastrate normiert (fs = 1).

## Bilder

![Pipeline-Blockschaltbild](docs/images/pipeline.png)

![Detektion und Lokalisierung](docs/images/detection.png)

Zeit-Huellkurve und Detektionsspektrum eines Demo-Puffers (QPSK + CW + FSK2
bei 20 dB SNR) mit den lokalisierten Pulsen: feine Zeitfenster (gelb),
detektierte Baender (blau) und Sub-Bin-Mittenfrequenzen (farbige Linien).

![Spektrogramm mit lokalisierten Pulsen](docs/images/spectrogram.png)

Spektrogramm desselben Puffers mit den Puls-Rechtecken (Zeit x Frequenz) aus
den Klassifikationsergebnissen.

Die Bilder werden deterministisch aus einem Demo-Puffer und den echten
Klassifikationsergebnissen erzeugt (aus dem Repo-Root):

    cmake --build build --target iqclassify_readme_images
    ./build/iqclassify_readme_images docs/images

## Dokumentation

- `docs/REQUIREMENTS.md` Anforderungen und Leistungsumfang
- `docs/METHODOLOGY.md` Methodik, Schwellwerte, Validierung, Grenzen
- `docs/API.md`   Programmierschnittstelle

## Verzeichnisstruktur

    include/iqclassify/   oeffentliche Header (Schirmheader: iqclassify.hpp)
    src/                  Implementierung
    tests/                 Testsuite (Signalgeneratoren, Szenarien)
    tools/                 Bildgenerator fuer die README-Bilder (docs/images/)
    third_party/pocketfft/ vendorte FFT (BSD-3-Klausel, Max-Planck-Gesellschaft)
    docs/                  Doku (Anforderung, Methodik, API) und Bilder (images/)

## Lizenz

MIT, siehe `LICENSE`. PocketFFT unterliegt seiner eigenen Lizenz, siehe
`third_party/pocketfft/LICENSE.md`.
