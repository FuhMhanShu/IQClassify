# Methodik

Alle Algorithmen und Schwellwerte wurden in einer ausführbaren
Referenzimplementierung (Python) gegen simulierte Signale kalibriert und
anschließend elementweise nach C++ portiert. Die Testsuite (tests/) deckt
die Kalibrierungskonfigurationen ab.

Verallgemeinerung: alle Frequenzen sind fraktional zur Abtastrate (fs = 1).
Die Pipeline läuft in fünf Phasen.

## Phase 0: Pulssegmentierung (Energie)

- Blockenergie über 64-Sample-Fenster.
- Rauschboden per Minimum-Statistik: Median je 16-Block-Gruppe, Minimum
  über alle Gruppen. Gruppen mit exakter digitaler Stille
  (Median <= 1e-20) werden übersprungen; sie würden den Boden auf null
  drücken.
- Hysterese 8 dB (Ein) / 5 dB (Aus) über dem Boden, Sicherheitsrand 64
  Samples, Mindestpuls 400 Samples. Direkt anschließende Segmente werden
  verschmolzen (Flickerschutz).

## Phase 1: Banddetektion und Basisbandextraktion

Detektion: gleitende Hann-gefensterte Detektions-FFTs (Fensterlänge bis
4096, 50 % Überlappung), per-Bin-Maximum über alle Teilfenster. So
erscheint jeder Puls mit voller Stärke, unabhängig von seiner Position im
Segment. Bänder sind zusammenhängende Bereiche über dem Rauschboden in
zyklischen Koordinaten ab dem jeweils stillsten Bin.

Bandkriterien: Flanke 12 dB, Peak 16 dB, Mittel 9 dB über Boden;
Mindestbreite max(8, nft/512) Bins; Lücken <= max(4, nft/512) Bins
verschmelzen. Schmale starke Bänder (Peak > Boden + 31 dB) umgehen die
Breiten-Grenzen (CW-Ton).

Satellitenverschmelzung: Phasen-/Frequenz-Keying erzeugt spektrale
Satelliten etwa eine Symbolrate neben dem Hauptzipfel; tiefe Einbrüche
zerlegen sie zu eigenständigen Läufen. Verschmolzen wird bei Lücke <=
5 % der größeren Breite, oder bei kleiner Absolutlücke (max(4, nft/64))
an eine mindestens sechsmal breitere Band. Zwei schmale Nachbarn bleiben
getrennt (eng benachbarte Töne überleben).

SNR (Spektraldichte): Verhältnis der mittleren in-band-Leistungsdichte zur
Außen-Leistungsdichte. Säume anderer Bänder werden (mit Breiten-Wächter)
aus der Außen-Referenz ausgeschlossen; im Mehrpuls-Fall bricht die
SNR-Schätzung sonst zusammen. Flachheit:
Standardabweichung der stark geglätteten in-band-dB-Kurve.

Extraktion: ungefensterte FFT des gesamten Segments, Bandmaske mit Wächter
(breite/12, mindestens 3 Bins; eine exakt anliegende Maske würde den
Übergang abschneiden und Ringing-ISI injizieren), zyklische Verschiebung
zur DC, IFFT. Der resultierende Basisband-Slice enthält alle Pulse der
Band in der Zeit.

Zeitliche Crops: Ein-/Aus-Hysterese der Slice-Blockenergie mit relativen
Toren gegen das Slice-Maximum (12/20 dB). Hat der Slice kein eigenes
Schweigen (ein Puls füllt das Segment), wird der gesamte Slice als ein
Puls ausgegeben.

## Lokalisierung (Wiederfinden des Pulses im IQ-Signal)

Jedes `PulseResult` trägt neben dem groben Blockraster-Fenster eine feine
Lage, direkt aus den Zwischengroessen der Pipeline abgeleitet:

- Zeit (samplegenau): Hüllkurve des (bandgefilterten, sanitisierten) Crops;
  erstes und letztes Sample oberhalb 25 % des Crop-Maximums. Die Raender des
  Crops (Sicherheitsrand, Blockraster der Segmentierung) fallen damit weg.
  Umrechnung: `t = startSampleFine / fs`.
- Frequenz (Sub-Bin): leistungsgewichteter Schwerpunkt des In-Band-PSD
  (Max-Hold-Detektionsspektrum, lineare Leistung) über den detektierten
  Baendern, zyklisch ab `band.lo`; bei schmalen Tönen (CW) liegt der
  Schwerpunkt auf der Tonbin, bei breiten Baendern in der Bandmitte.
  Umrechnung: `f = centerFreqRefined * fs` (zyklisch in [0,1)).
- `ModulationClassifier::extractPulse` extrahiert den Puls anschließend
  direkt aus dem Originalpuffer: Zeitfenster `[startSampleFine,
  startSampleFine + numSamplesFine)`, Bandmaske um `centerFreqRefined` mit
  `bandwidthNorm` (über `BandExtractor::extract`, inkl. Wächter), IFFT und
  Rückverschiebung auf die Bandmitte. `slice[k]` entspricht
  `iq[startSampleFine + k] * e^{-j*2pi*f0*k}`.

## Phase 2: Merkmale

Auf jedem (Band, Crop):
- Hüllkurve: CV, PAPR, fenstergemitteltes PAPR (Fenster 16).
- Trägerkohärenz carrCoh = |mean(z * e^{-j*2pi*f_res*t})| / mean(|z|):
  Ein reiner Ton hält den Wert nahe 1, jede Modulation mittelt den Zeiger
  weg. Entscheidend ist die vorherige Entfernung des Restträgers
  (Median der Momentanfrequenz): die Bandmitte liegt auf dem quantisierten
  FFT-Raster, ein sub-Bin-Restoffset würde den Zeiger über lange Crop-
  Längen wegmitten (Sinc-Auslöschung) und die CW-Erkennung pulslängen-
  und positionsabhängig machen. Ohne diese Entfernung versagte der
  Ton-Test beispielsweise bei einem 1000-Sample-CW-Puls und bei CW neben
  einem Chirp. Trennt weiterhin BPSK (Phasensprünge mitteln den Zeiger
  weg) vom unmodulierten Träger.
- Momentanfrequenz: kontinuierliche Phase, Median-Restträger-Entfernung.
- Chirp-Fit: lineare Regression auf f_inst, hüllkurven²-gewichtet,
  mittlere 80 % des Pulses; Bestimmtheitsmaß r² und Hub/Bandbreite.
- f_inst-Moden: Histogramm (24 Bins, nur Hüllkurve > 0.4 * Mittel;
  Randnormalisierung /3 wie Referenz), zusammenhängende Bins über 35 % des
  Maximums bilden eine Mode; Konzentration der Masse um die Top-2-Moden.

## Phase 3: Stufe 1 (grobe Klasse)

Reihenfolge (jede Regel schließt ab):
1. SNR-Gate: SNR < minSnrDb (6 dB) überspringt den Abschnitt.
2. Ton: Anteil der Samples mit |f_inst| < 0.004 über 70 % und
   Trägerkohärenz > 0.85. Der Test ist unabhängig von der detektierten
   Bandbreite; der Splatter benachbarter Pulse darf einen CW-Puls nicht in
   die Phasenanalyse drängen.
3. Schmalband (<= 14 Bins): Frequenz-Pfad.
4. Chirp: r² > 0.70 und Hub > 0.50 der Bandbreite: Frequenz-Pfad.
5. FSK: mindestens 2 f_inst-Moden mit Konzentration > 0.55.
6. Flächig: fenstergemitteltes PAPR > 6 dB und CV > 0.40 (OFDM-Gruppe).
7. sonst: Phasen-Pfad.

## Phase 4: Stufe 2

Ton → Unmodulated (Konfidenz 0.9).

Frequenz-Pfad: Chirp bei r²/Hub-Kriterien, sonst FSK2 bei >= 2 Moden,
Unmodulated bei genau einer Mode, sonst Unknown.

Phasen-Pfad (PSK/QAM/OFDM):
1. OFDM über zyklisches Präfix: normalisierte zirkulare
   Autokorrelation der Leistung, Kandidaten-Nutzteil 16..min(n/4, 1024);
   Spitze > 0.09, Nutzteil >= 48, Kontrast > 3 gegen Nachbarschaft,
   Talreinheit >= 0.80 (Anteil der Lags im +-48-Fenster unterhalb 30 % der
   Spitze; echte CP-Spitzen ragen aus sauberen Tälern, Einzelträger-
   Zufallsspitzen sitzen auf abfallenden Rücken). Bei Reinheit >= 0.95
   akzeptiert, sonst dieselbe Lag +-4 in beiden Puls-Slices bestätigt
   (Halbseiten-Konsistenz).
2. OFDM-Fallback: fenstergemitteltes PAPR > 6 dB bei glattem Spektrum
   (Flachheit < 2.5 + 1.5 dB).
3. Symboltiming: Symbolratlinie im zyklischen Leistungsspektrum
   (Kandidatenbereich Bandbreite/1.5 bis /1.05), parabolische Interpolation
   und Sub-Bin-Verfeinerung über die zyklische Kohärenz (biasfrei; ein
   Ratenfehler von 1e-4/Sample verdriftet das Gitter über wenige tausend
   Symbole um ein volles Sample). Augenschwerpunkt aus dem überabgetasteten
   Leistungsprofil (Schritt 0.125 Sample, parabolisch verfeinert).
4. Empfangsfilter: gematchte RRC-Filterung (Rolloff 0.22) vor der
   Symbolabtastung; gewinnt die volle Symbol-SNR zurück (gegenüber direkter
   Augenabtastung bei sps=8 rund 8 dB). Symbole unter 30 % des Median
   (Rampe am Crop-Rand) werden verworfen.
5. Restträger aus der Symbolfolge per z⁴-Spektrum (Pow-2-Puffer), parabolisch
   + Kohärenz-Ternärsuche auf 1e-7/Symbol.
6. Kumulant k42 = -C42/C21² (positiv, SNR-korrigiert), Symbolamplituden-CV,
   Phasenfaltung f2/f4.
   - BPSK: f2 > 0.75 and k42 > 1.5.
   - QPSK: f4 > 0.75 and 0.75 < k42 < 1.35 (sowie CV-/Falt-Fallbacks).
   - CV < 0.18: konstanter Modulus ohne Faltstruktur: Unknown.
   - sonst QAM-Ordnung.
7. QAM-Ordnung: rotationsinvarianter Rasterfit. Die I-Projektion wird
   gegen die analytischen Levelgitter von QAM16/64/256 geprüft
   (Drehung -45 bis +45 Grad, Schrittweite 2); der Quantisierungsfehler
   ist auf die Gitterkonstante normiert, irrelevante Gitter erzielen ~1/3,
   score = 1/3 - Fehler. Akzeptiert nur bei absolutem Score > 0.08 und
   Vorsprung > 0.03 vor dem Zweitplatzierten. Ohne auflösbare
   Levelstruktur greift der Dichte-Fallback: QAM256 ab SNR 20 dB (die
   QAM256-Levelabstände brauchen jedes dB), sonst QAM64.

OFDM-Gruppe (flächiges Signal): OFDM mit konfidenzabhängiger Entscheidung.

## Konfiguration

Jeder genannte Wert steht in `ClassifierConfig` (include/Config.hpp) und
ist dort mit seiner Bedeutung kommentiert. Eine Änderung erfordert eine
 erneute Validierung gegen die Testsuite.

## Validierung (Referenzimplementierung, C++ portiert ohne Änderungen)

Die Referenzimplementierung führt die C++-Pipeline elementweise aus
(äquivalente Arithmetik) und wurde gegen eine erweiterte Qualitätsmatrix
aus 168 Szenarien gefahren. Die C++-Testsuite (tests/test_main.cpp,
100 Szenarien) deckt den verifizierten robusten Bereich ab; jedes Szenario
der Suite wurde in der Referenz über mindestens drei Seeds bestätigt.

- Basiskalibrierung, SNR 20 dB, sps 8: 9/9 korrekt (BPSK, QPSK, FSK2,
  Chirp, CW, OFDM, QAM16, QAM64, QAM256).
- SNR-Raster 30/20/15/12/10/8/6/3/0 dB je Typ:
  - 15 dB: alle 9 Typen korrekt.
  - 10 dB: 8/9 (QAM256 fällt auf QAM64, dokumentierte Grenze; FSK2
    erst ab ~12 dB über alle Seeds stabil).
  - 6 dB: BPSK/QPSK/QAM16/QAM64/Chirp/CW/OFDM korrekt; darunter
    zunehmende Fragmentierung des Banddetektors (unterhalb der
    Spezifikation von 6 dB Mindest-SNR).
- Oversampling: BPSK/QPSK/QAM16 bei sps 4/8/16/32 bestanden; FSK2 bei
  sps 4/8 bestanden.
- Frequenzlagen fc 0.00 bis 0.45 und -0.30 (QPSK, FSK2, Chirp, CW):
  bestanden.
- Chirp-Varianten: up/down, Bandbreiten 0.10/0.25/0.40, SNR 20 und
  10 dB: 12/12 bestanden.
- Pulslängen: QPSK 150 bis 1000 Symbole (1200 bis 8000 Samples); CW 500
  bis 8000 Samples: bestanden.
- Mehrpuls: Sequenz QPSK/FSK2/CW; zwei BPSK hintereinander mit
  vertauschten Frequenzen; zeitgleich überlagerte Pulse BPSK+QAM16,
  QPSK+FSK2, CW+QPSK, Chirp+CW, OFDM+QPSK (frequenzgetrennt): bestanden.
- Seed-Streuung: QPSK bei 10 dB und QAM16 bei 20 dB, je 5 Seeds ohne
  Ausreißer (Referenz); OFDM (nfft 256, cp 64, 40 Symbole, Last 0.3)
  bei 10/15/20 dB über 3 Seeds stabil.
- Härtungsfälle: Nullzeiger, leere und Miniatureingaben, reines Rauschen,
  Null- und DC-Puffer, vollständig NaN-/Inf-beladene Puffer sowie
  NaN-/Inf-/Null-Dropouts mitten in einem Puls: definierte Ergebnisse,
  keine Abstürze, Dropout verhält sich wie Stille (Puls bleibt
  klassifizierbar).

## Härtung

Eingabesicherung (defense in depth, mehrschichtig):
- ModulationClassifier::classify: Nullzeiger und Länge 0 liefern leer.
- Banddetektion/Extraktion: nicht-endliche IQ-Samples werden vor der
  FFT als Ausfall (Null) behandelt; ein einzelnes NaN vergiftet sonst
  über die Butterfly-Operation das gesamte Spektrum des Fensters und
  damit den ganzen Puls. Der Effekt ist identisch zu Null-Dropouts.
- Der Klassifikations-Crop wird zusätzlich sanitisiert (NaN/Inf -> 0).
- Median: nicht-endliche Werte werden gefiltert; std::sort mit NaN
  verletzt die strikte schwache Ordnung und wäre undefiniertes Verhalten
  bis hin zum Absturz.
- f_inst-Histogramm: isfinite-Schutz vor der Indexumrechnung.
- Frühabrutzzweig der Banddetektion (Fenster < 16) füllt die PSD mit
  -300 dB, damit die SNR-Schätzung nie auf ein leeres PSD stößt
  (Modulo-Division durch null).
- estimateSnrFlat: leere PSD oder degenerierte Bandbreite liefern
  SNR -99 dB, Flachheit 0.

## Grenzen

- Mindestpuls 400 Samples, Mindest-SNR 6 dB (Vollauslastung aller
  Typen ab ~12 dB; QAM256 ab ~15-20 dB; FSK2 ab ~12-15 dB stabil).
- Oversamplingfaktor: sps 4 bis 8 für alle Typen; BPSK/QPSK/QAM bis
  sps 32. sps 2 erzeugt pulsförmige Breitbandsignale und wird nicht
  unterstützt; FSK2 zerfällt ab sps 16 in Tonsatelliten.
- OFDM: zyklisches Präfix muss mindestens 48 Samples lang sein
  (Kandidaten-Nutzteil beginnt dort); cp 16 ist demnach nicht
  erkennbar. Traglast 0.3 benötigt mindestens ~40 Symbole (bei 20
  Symbolen liegt die CP-Autokorrelationsspitze zu nah an der Schwelle);
  Last 1.0 mit nfft >= 512 füllt das gesamte Nyquist-Band und fällt
  aus dem Bandmodell.
- Gleichzeitige Pulse in derselben Frequenzlage sind nicht trennbar;
  frequenzgetrennte Überlagerungen benötigen einen Mittenabstand von
  grob 0.10 oder mehr (halbierende Bandbreiten dürfen sich nicht
  überlappen).
- Die Mindeststatistik der Pulssegmentierung benötigt mindestens eine
  ruhige Teilphase im Puffer: ein durchgehendes Signal ohne jede
  Schweigephase liefert keinen Rauschboden und damit keinen Puls.
- Symbolraten außerhalb 1.6..n/12 Samples/Symbol werden nicht verfolgt.
