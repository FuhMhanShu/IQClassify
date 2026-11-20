// IQClassify - Konfiguration aller Schwellwerte und Parameter.
//
// Jede Groesse entspricht einem validierten Parameter des
// Referenzprototyps; die Standardwerte sind gegen simulierte Pulsen
// (BPSK, QPSK, FSK2, Chirp, CW, OFDM, QAM16/64/256) bei 5..20 dB SNR
// und Oversamplingfaktoren 4/8/16 kalibriert. Aenderungen erfordern eine
// erneute Validierung gegen die Testsuite (tests/).

#pragma once

#include <cstddef>

namespace iqclassify {

struct ClassifierConfig {
    // ---- Spektrale Kanalisierung (erster Pipelineschritt)
    bool channelize = true;   // Spektrum zuerst in Kanaele zerlegen, die sich
                              // um die detektierten Energiebereiche legen (je
                              // Band ein Kanal um dessen Energieschwerpunkt);
                              // erst danach laeuft je Kanal die weitere Pipeline
                              // (Zeitsegmentierung, Crops, Merkmale, Stufe 1/2).
                              // false = Legacy-Pfad (erst Zeitsegmentierung, dann
                              // Banddetektion je Segment)
    int upsample = 0;         // Upsampling des extrahierten Basisband-Slices fuer
                              // sehr kurze Pulse (z. B. 180-Sample-Puffer mit 13
                              // BPSK-Bit bei sps=8): 0 = automatisch (Slice < 1024
                              // Samples -> U = clamp(ceil(1024/len),1,8)), 1 = aus,
                              // >1 = fester Faktor U. Sample-basierte Schwellwerte
                              // werden mit U skaliert, Bandbreite/Bins der Merkmale
                              // umgerechnet, Ergebniszeiten auf das Original
                              // zurueckgerechnet (/U); Frequenzen bleiben in der
                              // Originalrate (Kanalisierung laeuft unveraendert)
    int threads = 0;          // Parallelverarbeitung der Kanaele/Segmente:
                              // 0 = automatisch (min(hardware_concurrency, 8)),
                              // 1 = seriell. Wirkt ab ~32k Samples Pufferlaenge;
                              // Ergebnisse sind threadzahl-unabhaengig identisch
    int channelDecimate = 1;  // Resample-Faktor je Kanal (1 = volle Abtastrate);
                              // wird pro Kanal auf <= floor(1/Kanalbreite)
                              // begrenzt; sample-basierte Schwellwerte der
                              // Folgestufen werden entsprechend skaliert und die
                              // Ergebnisse zurueckgerechnet

    // ---- Pulssegmentierung (Energie im Block)
    std::size_t segWin = 64;        // Blocklaenge / Energie-Fenster
    std::size_t segFloorGroup = 16; // Bloecke pro Rauschboden-Gruppe
    double segOnDb = 8.0;           // Hysterese: Ein-Schwelle ueber Boden
    double segOffDb = 5.0;          // Hysterese: Aus-Schwelle ueber Boden
    std::size_t segMargin = 64;      // beidseitig angelagerte Samples
    std::size_t minPulseLen = 400;   // minimale akzeptierte Pulslaenge
    std::size_t minAnalyzeLen = 256; // Mindestlaenge eines Crops fuer Symbolanalyse
                                     // und CP-Autokorrelation (Stufe 2)
    double cropOnRelDb = 12.0;       // Crop-Ein: innerhalb dieses Bereichs des Slice-Maximums
    double cropOffRelDb = 20.0;      // Crop-Aus: unterhalb dieses Bereichs des Slice-Maximums
    double cropFadeMergeDb = 5.0;    // Fade-Merge: zwei Crops verschmelzen, wenn die
                                     // Luecke nie unter (cropFadeMergeDb ueber dem
                                     // slice-eigenen Minimum) fiel, also kein
                                     // Schweigen, sondern ein Amplituden-Einbruch ist

    // ---- Banddetektion
    std::size_t detectWin = 4096;   // Teilfensterlaenge der Detektions-FFTs
    double bandGuardDb = 12.0;   // Schwelle ueber Boden fuer Bandstart
    double bandPeakDb = 16.0;    // erforderlicher Bandpeak ueber Boden
    double bandMeanDb = 9.0;     // erforderliches Bandmittel ueber Boden
    double bandBridgeDb = 6.0;   // Low-SNR-Bruecke: Luecke zwischen zwei Laeufen
                                 // verschmelzen, wenn ihr mittleres PSD noch
                                 // bandBridgeDb ueber dem Boden liegt (Signal
                                 // in der Luecke, kein Rauschen)
    double bandStrongDb = 31.0;  // schmale starke Toene umgehen die Breitengrenzen
    double bandMergeFrac = 0.05; // Satellitenmerge: Luecke <= Anteil * Breite
    double satMergeWidthRatio = 6.0;  // relatives Breitenverhaeltnis fuer Absolut-Merge
    long satMergeMin = 4;               // Mindestluecke fuer Satelliten-Anhaengung
    long satMergeDenom = 64;            // Absolutluecke = max(4, nft/64) Bins

    // ---- SNR / Filterung
    double minSnrDb = 2.0;       // SNR-Gate Stufe 1; Default erlaubt Detektion
                                 // bis herunter zu ~2 dB SNR (CW/Chirp), breit-
                                 // bandige Typen brauchen ~3-4 dB fuer korrekte
                                 // Klassifikation (siehe docs/METHODOLOGY.md)
    double tFlatDb = 2.5;        // in-band Spektralglaette (Std, dB)

    // ---- Stufe 1
    std::size_t tBwCwBins = 14;  // Bin-Anzahl fuer den Schmalband-Pfad
    double tCwFthr = 0.004;      // |f_inst|-Toleranz der Tonentscheidung
    double tCwFrac = 0.70;       // Anteil der Samples innerhalb der Toleranz
    double tCwCoh = 0.85;        // Notwendige Traegerkohaerenz fuer einen Ton
    double tChirpR2 = 0.70;      // Chirp-Fitqualitaet
    double tChirpSweep = 0.50;   // Chirp-Hubanteil der Bandbreite
    double tModeConc = 0.55;     // f_inst-Modenkonzentration (FSK)
    // ---- Frank-Code (polyphase Radar-Code): Treppe der Momentanfrequenz
    int tFrankMinModes = 3;      // Mindestanzahl Stufen (Ordnung M >= 3)
    int tFrankMaxModes = 24;     // maximale Stufenzahl (<= histBins)
    double tFrankR2 = 0.95;      // R^2 der Stufen-Mitten gegen Stufen-Index
    double tFrankMassEq = 6.0;   // erlaubtes max/min der Stufen-Massen
    double tFrankCover = 0.40;   // Anteil der gated Samples in Stufen
    double tFrankMono = 0.80;    // Anteil der monotonen (zyklischen) Stufenfolge
    int tBarkerMinLen = 7;       // kleinste Barker-Laenge, die als BarkerCode
                                 // gelabelt wird (kuerzere Folgen sind von
                                 // Zufalls-BPSK nicht unterscheidbar)
    double tPaprOfdmDb = 6.0;    // Fenstergemitteltes (16) PAPR der OFDM-Gruppe
    double tCvOfdm = 0.40;       // Hoellkurven-CV der OFDM-Gruppe

    // ---- OFDM: zyklisches Praefix
    long ofdmTauMin = 16;        // kleinste Kandidatenlaenge des Nutzteils
    long ofdmTauMax = 1024;      // groesste Kandidatenlaenge des Nutzteils
    std::size_t ofdmMinBandBins = 32;
    double tCpPeak = 0.09;       // CP-Autokorrelationsspitze
    long tCpMinTau = 48;         // minimale Nutzteillaenge (nfft >= 48)
    double tCpContrast = 3.0;    // Spitze dominiert ihre Nachbarschaft
    double tCpClean = 0.80;      // Talreinheit um die Spitze
    double tCpDeepFrac = 0.30;   // Anteil unterhalb dessen ein Lag als Tal gilt
    double tCpCleanSure = 0.95; // Sauberkeit, die ohne Halbtest akzeptiert
    long cpCleanWin = 48;        // +-Lagfenster des Sauberkeitstests
    long cpHalfTol = 4;          // Lagtoleranz zwischen den Haelften
    double cpHalfRhoRel = 0.5;   // Mindestspitzenverhaeltnis der Haelften
    long cpIsolation = 8;        // Lags nahe der Spitze zaehlen nicht als zweite

    // ---- Symboltiming
    std::size_t timingMinSyms = 24;
    double tRsProm = 5.0;        // Prominenz der zyklischen Leistungslinie
    double rsBwDivLo = 1.5;      // untere K-Grenze der Symbolratensuche:
                                 // K ab bw/rsBwDivLo (fr = K/nft); kurze Pulse
                                 // haben ein relativ breites Band -> dort weiter
                                 // nach unten oeffnen (siehe shortPulseProfile)
    double rsBwDivHi = 1.05;     // obere K-Grenze: K bis bw/rsBwDivHi
    double mfAlpha = 0.22;       // Empfangs-RRC-Rolloff fuer Symbolabtastung
    std::size_t mfSpan = 8;      // Empfangsfilterlaenge in Symbolen
    double mfMinRel = 0.30;      // Symbolamplituden-Gate relativ zum Median

    // ---- Stufe 2: Phasenanalyse
    double tFold2 = 0.75;        // BPSK-Faltkonzentration
    double tFold4 = 0.75;        // QPSK-Faltkonzentration
    double tCvPskHi = 0.35;      // CV-Decke fuer PSK-Fallbacks
    double tCvConstMod = 0.18;    // darunter: konstanter Modulus, keine Entscheidung
    double k42Qam64 = 0.616;
    double k42Qam256 = 0.60;

    // ---- QAM-Ordnung (Rasterfit)
    double qamScoreStrong = 0.08;  // Score zur Akzeptanz des Argmax
    double qamScoreMargin = 0.03;   // erforderlicher Vorsprung vor Rang 2
    double qamDenseSnr = 20.0;       // SNR, ab dem unaufgeloeste Konstellationen
                                     // als QAM256 bewertet werden

    // ---- f_inst-Histogramm / Moden
    int histBins = 24;
    double histEnvRel = 0.4;     // Hoellkurven-Gate relativ zum Mittel
    double modeBinThr = 0.35;    // Bin muss diesen Anteil des Maximums erreichen
    int modeSplit = 2;           // Binluecke, die zwei Moden trennt

    // ---- Dwell (informativ)
    int tDwellRun = 4;

    // Profil fuer sehr kurze Pulse (z. B. ~180-Sample-Puffer mit 13
    // BPSK-Bit bei sps=8): feinere Blocksegmentierung, kleinere
    // Mindestpulslaenge, niedrigere Symbolschwelle und kleinere
    // Analyse-Mindestlaenge. Default-Profil bleibt auf laengere Pulse
    // kalibriert (minPulseLen 400).
    static ClassifierConfig shortPulseProfile() {
        ClassifierConfig c;
        c.segWin = 16;
        c.segFloorGroup = 4;
        c.segMargin = 16;
        c.minPulseLen = 96;
        c.detectWin = 256;
        c.timingMinSyms = 8;
        c.minAnalyzeLen = 64;
        c.tRsProm = 3.0;
        c.rsBwDivLo = 2.0;  // kurze Pulse: Band relativ breit gegen Symbolrate
        // Das fenstergemittelte PAPR kurzer Pulse ist durch die Pulsflanken
        // natuerlich erhoeht; die OFDM-Schwelle wird daher angehoben, damit
        // kurze PSK-Pulse nicht als NoiseLike (OFDM) durchgereicht werden.
        c.tPaprOfdmDb = 9.0;
        return c;
    }
};

}  // namespace iqclassify
