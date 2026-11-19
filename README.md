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
- Klassenkapselung: eine Fassade (`ModulationClassifier`), alle
  Schwellwerte konfigurierbar (`ClassifierConfig`).
- Deterministische Testsuite ohne externes Testframework.

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
    }

Alle Frequenzen sind auf die Abtastrate normiert (fs = 1).

## Dokumentation

- `docs/REQUIREMENTS.md` Anforderungen und Leistungsumfang
- `docs/METHODOLOGY.md` Methodik, Schwellwerte, Validierung, Grenzen
- `docs/API.md`   Programmierschnittstelle

## Verzeichnisstruktur

    include/iqclassify/   oeffentliche Header (Schirmheader: iqclassify.hpp)
    src/                  Implementierung
    tests/                 Testsuite (Signalgeneratoren, Szenarien)
    third_party/pocketfft/ vendorte FFT (BSD-3-Klausel, Max-Planck-Gesellschaft)
    docs/                  Anforderungs-, Methodik- und API-Dokumentation

## Lizenz

MIT, siehe `LICENSE`. PocketFFT unterliegt seiner eigenen Lizenz, siehe
`third_party/pocketfft/LICENSE.md`.
