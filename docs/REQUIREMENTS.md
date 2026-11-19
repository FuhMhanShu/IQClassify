# Anforderungen und Leistungsumfang

## Aufgabe

Klassifikation der Modulationsart von Pulsen in komplexen IQ-Abtastströmen
(zweistufig):

- Stufe 1, grobe Klasse: Ton, Frequenzmodulation (Chirp/FSK),
  Phasenmodulation (PSK/QAM), flächiges rauschartiges Signal (OFDM).
- Stufe 2, Modulationsart: Unmodulated (unmodulierter Träger, CW), BPSK,
  QPSK, FSK2, Chirp, OFDM, QAM16, QAM64, QAM256. Nicht entscheidbare Pulse
  werden als Unknown mit niedriger Konfidenz ausgegeben.

## Mehrpuls-Szenarien

Ein Puffer kann mehrere Pulse enthalten:

1. zeitlich getrennte Pulse (sequenziell, mit Schweigen dazwischen),
2. mehrere Pulse gleichzeitig in verschiedenen Frequenzbändern
   (frequenzgetrennt, überlagert in der Zeit).

Sind zwei Pulse sowohl zeitlich als auch spektral überlappt (gleiche
Frequenzlage), sind sie physikalisch nicht trennbar und werden als ein
Puls klassifiziert; dies ist eine dokumentierte Grenze, kein Fehler.

## Lokalisierung

Jeder erkannte Puls muss im IQ-Signal wiederfindbar sein:

- samplegenaue Zeitlage (`startSampleFine`/`numSamplesFine` in Samples,
  `t = startSample / fs`),
- Sub-Bin-Mittenfrequenz (`centerFreqRefined`, fraktional zur Abtastrate,
  zyklisch in [0,1), `f = centerFreqRefined * fs`),
- `extractPulse`/`extractPulseF`: Basisband-Slice des Pulses aus dem
  Originalpuffer (Zeitfenster der feinen Lage, Bandpass um die feine
  Mittenfrequenz), sodass der Puls ohne erneute Suche weiterverarbeitet
  werden kann.

## Low-SNR-Detektion

Auch Signale unterhalb von 6 dB SNR müssen erkannt werden, Ziel ~2 dB:

- SNR-Gate (`minSnrDb`) Default 2 dB; Banddetektor brückt Lücken zwischen
  Bandläufen, die noch deutlich über dem Rauschboden liegen
  (`bandBridgeDb`), damit schwache Breitbandsignale nicht fragmentieren.
- Reines Rauschen darf auch bei niedrigem SNR keine Pulse erzeugen.
-Geeignete IQ-Testsignale (CW, Chirp, BPSK, QPSK, QAM16, FSK2, OFDM) bei
  0-4 dB SNR in der Testsuite verankert.

## Analyse-API je erkanntem Puls

- PSK (BPSK/QPSK): `synchronizePsk` liefert den phasenrichtig
  aufsynchronisierten Pulsanfang — Symboltakt, Trägerfrequenz und
  Trägerphase am Pulsanfang (auf das M-PSK-Raster aufgelöst) sowie die
  derotierten Symbole. Phasenkorrektur ohne Einschwingvorgang eines
  Clock-Synchronizers (Batch-Schätzung über den ganzen Puls).
- Chirp: `analyzeChirp` berechnet alle Parameter — Anfangs-/Endfrequenz,
  Mittenfrequenz, Hub, Steigung, Richtung, Bestimmtheitsmaß r², SNR.
- Symbol-Schnittstellen: `extractSymbols` liefert für PSK die Symbole und
  für FSK2 Frequenz/Ton-Entscheidung pro Symbol (inkl. Symboltakt), als
  Ausgangspunkt für die spätere Symbolgewinnung.

## Rahmenbedingungen

- Einbindung als statische Bibliothek, unterstützt via CMake
  (IQClassify STATIC, CTest-Integration) und qmake (CONFIG += staticlib,
  IQClassify.pri für Konsumentenprojekte).
- Keine Abhängigkeiten außer PocketFFT (header-only, unter
  third_party/pocketfft/ mitgeführt).
- C++17, keine RTTI-/Ausnahme-Anforderungen über die Standardbibliothek
  hinaus.
- Klassenkapselung: alle Verarbeitungsschritte in eigenen Klassen, eine
  öffentliche Fassade, alle Schwellwerte in einer Konfigurationsstruktur.
- Dokumentation von Anforderungen, Methodik und API in docs/.
- Deterministische Testsuite ohne externes Testframework.

## Nicht-Funktionsanforderungen

- Interne Rechnung in double; Eingabe zusätzlich als float komfortabel.
- Vorhersehbarer Speicherbedarf: Eingabepuffer werden nicht kopiert
  außer für fensterweise Spektralanalysen.
- Konfigurierbarkeit: Jeder Schwellwert der Referenzvalidierung bleibt
  einstellungsfähig (ClassifierConfig), siehe docs/METHODOLOGY.md.
