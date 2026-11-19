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
    // ---- Pulssegmentierung (Energie im Block)
    std::size_t segWin = 64;        // Blocklaenge / Energie-Fenster
    std::size_t segFloorGroup = 16; // Bloecke pro Rauschboden-Gruppe
    double segOnDb = 8.0;           // Hysterese: Ein-Schwelle ueber Boden
    double segOffDb = 5.0;          // Hysterese: Aus-Schwelle ueber Boden
    std::size_t segMargin = 64;      // beidseitig angelagerte Samples
    std::size_t minPulseLen = 400;   // minimale akzeptierte Pulslaenge
    double cropOnRelDb = 12.0;       // Crop-Ein: innerhalb dieses Bereichs des Slice-Maximums
    double cropOffRelDb = 20.0;      // Crop-Aus: unterhalb dieses Bereichs des Slice-Maximums

    // ---- Banddetektion
    std::size_t detectWin = 4096;   // Teilfensterlaenge der Detektions-FFTs
    double bandGuardDb = 12.0;   // Schwelle ueber Boden fuer Bandstart
    double bandPeakDb = 16.0;    // erforderlicher Bandpeak ueber Boden
    double bandMeanDb = 9.0;     // erforderliches Bandmittel ueber Boden
    double bandStrongDb = 31.0;  // schmale starke Toene umgehen die Breitengrenzen
    double bandMergeFrac = 0.05; // Satellitenmerge: Luecke <= Anteil * Breite
    double satMergeWidthRatio = 6.0;  // relatives Breitenverhaeltnis fuer Absolut-Merge
    long satMergeMin = 4;               // Mindestluecke fuer Satelliten-Anhaengung
    long satMergeDenom = 64;            // Absolutluecke = max(4, nft/64) Bins

    // ---- SNR / Filterung
    double minSnrDb = 6.0;
    double tFlatDb = 2.5;        // in-band Spektralglaette (Std, dB)

    // ---- Stufe 1
    std::size_t tBwCwBins = 14;  // Bin-Anzahl fuer den Schmalband-Pfad
    double tCwFthr = 0.004;      // |f_inst|-Toleranz der Tonentscheidung
    double tCwFrac = 0.70;       // Anteil der Samples innerhalb der Toleranz
    double tCwCoh = 0.85;        // Notwendige Traegerkohaerenz fuer einen Ton
    double tChirpR2 = 0.70;      // Chirp-Fitqualitaet
    double tChirpSweep = 0.50;   // Chirp-Hubanteil der Bandbreite
    double tModeConc = 0.55;     // f_inst-Modenkonzentration (FSK)
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
};

}  // namespace iqclassify
