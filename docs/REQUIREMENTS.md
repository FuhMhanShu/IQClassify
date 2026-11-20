# Anforderungen und Leistungsumfang

## Aufgabe

Klassifikation der Modulationsart von Pulsen in komplexen IQ-Abtastströmen
(zweistufig):

- Stufe 1, grobe Klasse: Ton, Frequenzmodulation (Chirp/FSK),
  Phasenmodulation (PSK/QAM), flächiges rauschartiges Signal (OFDM).
- Stufe 2, Modulationsart: Unmodulated (unmodulierter Träger, CW), BPSK,
  QPSK, FSK2, Chirp, OFDM, QAM16, QAM64, QAM256, FrankCode (polyphase
  Frank-Code, Radar), PCode (P1/P2; P3/P4 erscheinen als Chirp) und
  BarkerCode (BPSK mit Barker-Folge 7/11/13). Nicht entscheidbare Pulse
  werden als Unknown mit
  niedriger Konfidenz ausgegeben.

## Mehrpuls-Szenarien

Ein Puffer kann mehrere Pulse enthalten:

1. zeitlich getrennte Pulse (sequenziell, mit Schweigen dazwischen),
2. mehrere Pulse gleichzeitig in verschiedenen Frequenzbändern
   (frequenzgetrennt, überlagert in der Zeit).

Sind zwei Pulse sowohl zeitlich als auch spektral überlappt (gleiche
Frequenzlage), sind sie physikalisch nicht trennbar und werden als ein
Puls klassifiziert; dies ist eine dokumentierte Grenze, kein Fehler.

## Spektrale Kanalisierung (erster Pipelineschritt)

Die Pipeline beginnt mit einer spektralen Kanalisierung: das Spektrum wird in
Kanäle zerlegt, die sich nicht gleichverteilt, sondern um die Bereiche im
Spektrum legen, in denen Energie detektiert wurde (adaptive Kanäle, je
Energiebereich ein Kanal um dessen Energieschwerpunkt). Jeder Kanal wird
einzeln untersucht: extrahiert (Bandpass), auf seine Mitte heruntergemischt
und optional resamplet (dezimiert mit Anti-Alias-Filter). Die Folgemodule
(Zeitsegmentierung, Banddetektion, Merkmale, Klassifikation) laufen pro
Kanal; alle Ergebnisse beziehen sich auf das Originalsignal (Zeit in
Original-Samples, Frequenzen zyklisch in [0,1) der Originalrate). Über
`Channelizer` sind die Kanal-Basisbänder auch direkt zugreifbar
(`channelize`/`extract`/`resample`).

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

## Frank-Codes (polyphase Radar-Codes)

- Frank-Code der Ordnung M (Laenge N = M^2 Chips, Chip (i,j) mit Phase
  2*pi*i*j/M) wird als `FrankCode` klassifiziert. Erkannt wird die
  Treppenstruktur der Momentanfrequenz: M Stufen, gleichmaessig beabstandet,
  jede lange am Stueck (Run-Laenge), Stufen-Mitten linear gegen den Index
  (R^2). Abgrenzung: Chirp (kontinuierlicher Sweep, kein langer Run pro
  Stufe), FSK2 (2 Stufen), PSK (kein Sweep).
- Analyse-API `analyzeFrank`: Ordnung M, Code-Laenge N = M^2, Chip-Rate,
  Mittenfrequenz, belegtes Band, SNR, R^2 der Stufen-Mitten und die
  geschaetzten Code-Phasen (N Werte, [0, 2pi)).
- Testgenerator `iqtest::genFrank(rng, M, sps)`; Testsuite Szenario 21.
- Dokumentierte Grenze: das Kamm-Spektrum hoeherer Ordnungen kann bei
  niedrigem SNR in Subbaender zerfallen (dann mehrere FrankCode-Pulse);
  direkt an einen starken CW grenzende kurze Pulse teilen sich ein
  Zeitsegment und werden als ein Puls behandelt (bekannte Segment-Grenze).

## P-Codes (P1..P4) und Barker-Codes

- P1/P2-Code der Ordnung M (Laenge N = M^2 Chips, Phasen-Gesetz
  `-(pi/M)(M-2i)(M-2j)` bzw. `-(pi/M)(M-2i)(M-2j-1)`) wird als `PCode`
  klassifiziert. Erkannt ueber die `f_inst`-Treppe (wie Frank) plus
  Matched-Filter-Vergleich der Phasen-Gesetze gegen Frank/P1/P2
  (`polyLawMf`, Chip-Rate fein um n/M^2 gesucht); ohne erkennbare Treppe
  (kleine Ordnung, starke Impulsformung) greift ein strengerer Fallback
  (Score > 0.10 und Ratio > 1.8), damit FSK2 nicht als PCode durchgeht.
- P3/P4 (Chirp-artig, Laenge M) werden als `Chirp` klassifiziert; P3
  erscheint durch Faltung am Chip-Rand als Saegezahn in `f_inst` und wird
  ueber den gefalteten Chirp-Test erkannt.
- Analyse-API `analyzePcode`: Subtyp (1 = P1, 2 = P2), Ordnung M,
  Code-Laenge N = M^2, Chip-Rate, Mittenfrequenz, SNR, R^2, lawErr und die
  N geschaetzten Code-Phasen.
- Barker-Code (BPSK-Puls mit Barker-Folge 7/11/13, wiederholt) wird als
  `BarkerCode` fein-klassifiziert: Barker-Match direkt auf dem Slice
  (Symboltakt fein gesucht, alle Rotationen und Vorzeichen, >= 95 %
  Uebereinstimmung). Repetitive Barker-Codes haben ein Linienspektrum und
  CP-aehnliche Periodizitaet — die OFDM-Tests sind deshalb gegen 2-PSK-
  Faltung abgesichert. Analyse-API `analyzeBarker`: Code, Laenge,
  Wiederholungen, sps, Traeger.
- Testgeneratoren `iqtest::genPcode(rng, p, M, sps)` (p = 1..4) und
  `iqtest::genBarker(rng, len, sps, repeats)`; Testsuite Szenarien 22 und 23
  (inkl. Negativtests: FSK2 bleibt FSK2, zufaelliger BPSK bleibt BPSK).
- Doku mit Bildern: docs/MODULATIONEN.md.

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
