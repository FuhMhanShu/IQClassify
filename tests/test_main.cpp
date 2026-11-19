// IQClassify - Testsuite.
//
// Szenarien (jeweils gegen die Python-Referenzmatrix mit mehreren Seeds
// abgesichert; siehe docs/METHODOLOGY.md):
//   1. Basiskalibrierung: 9 Modulationen, SNR 20 dB, sps 8
//   2. SNR 15 dB: alle 9 Typen
//   3. SNR 10 dB: 8 Typen (ohne FSK2: ab ~12 dB erst zuverlaessig),
//      QAM256 als QAM64-Fallback dokumentiert
//   4. Sequenzen: drei Pulse hintereinander; 2x BPSK mit vertauschten fc
//   5. Gleichzeitige Pulse in getrennten Frequenzbaendern (5 Kombinationen)
//   6. Oversampling: BPSK/QPSK/QAM16 bei sps 4, 16, 32; FSK2 bei sps 4
//   7. Frequenzlagen: fc 0.0 bis 0.45 und -0.30
//   8. Chirp up/down, Bandbreiten 0.10/0.25/0.40, SNR 20/10 dB
//   9. Pulslaengen: QPSK 150/1000 Symbole, CW 500/8000 Samples
//  10. Seed-Streuung: QPSK @10 dB und QAM16 @20 dB, je 3 Seeds
//  11. Haertung: Nullzeiger, leere/miniature Eingaben, reines Rauschen,
//      Nullen-, DC-, NaN- und Inf-Eingaben, Dropouts,.float-Pfad,
//      kurzer CW-Puls
//  12. Lokalisierung: samplegenaue Zeitlage + Sub-Bin-Frequenz der
//      detektierten Pulse (QPSK/CW/Chirp/FSK2/OFDM, 2-Puls-Szene) und
//      extractPulse/extractPulseF (Wiederfinden im IQ-Puffer)
//  13. Amplitudenschwankende Pulse: Fade-Merge (ein Puls bleibt ein Puls),
//      korrekte Klassifikation und Lokalisierung trotz schwankender Amplitude
//  14. Low-SNR bis 2 dB: lange Pulse werden erkannt, lokalisiert und
//      klassifiziert (Default-Gate minSnrDb=2, Low-SNR-Bruecke im Banddetektor)
//  15. Analyse-API: PSK-Synchronisation (phasenrichtiger Pulsanfang ohne
//      Clocksync-Einschwingen), Chirp-Parameter, Symbol-Schnittstelle PSK/FSK
//  16. Spektrale Kanalisierung (erster Pipelineschritt): adaptive Kanaele um
//      die Energiebereiche, Extract/Mix/Resample je Kanal, Aequivalenz zum
//      Legacy-Pfad, Resampling in classify, breiter Chirp bleibt ein Puls
//  17. Skalierung/Performance: 2 Mio Samples in einem Durchlauf klassifizieren
//      (Laufzeit wird ausgegeben, nur gegen Pathologien abgesichert)
//  18. Haertung II: Extremamplituden (1e150/1e-300), Nyquist-Muster, Puls am
//      Pufferanfang/-ende, 60 dB SNR, Mini-Puffer (n=1..70), Puffer vollstaendig
//      vom Puls ausgefuellt (dokumentierte Segmenter-Grenze: 0 Pulse)
//  19. Parallelverarbeitung: threads=1 vs. threads=4 liefern identische
//      Ergebnisse (kanalisierte Pipeline + Channelizer)
//  20. Sehr kurze Pulse (~180 Samples, 13 BPSK-Bit bei sps=8): Detektion und
//      Lokalisierung via automatischem Slice-Upsampling (cfg.upsample) +
//      shortPulseProfile; Rauschboden bei wenigen Bloecken per Perzentil;
//      CW/FSK2/QPSK ebenso; Default-Profil erwartet laengere Pulse (0 Pulse)
//  21. Frank-Codes (polyphase Radar-Codes, Ordnung M, Laenge M^2): Detektion
//      als FrankCode via Treppenstruktur der Momentanfrequenz (Stufen gleich-
//      abstaendig, lange Runs, R^2 der Stufen-Mitten), M=8/16 bei 20/15 dB,
//      Mehrpuls-Szene Frank+CW, Analyse-API analyzeFrank (Ordnung, Chip-Rate,
//      Code-Phasen); Chirp/FSK2/QPSK duerfen nicht als FrankCode durchgehen
//
// Ein Szenario besteht genau dann den Test, wenn die sortierte Labelmenge
// der erkannten Pulse der erwarteten entspricht. Fuer
// Robustheitsfaelle gilt: definiertes Ergebnis ohne Absturz,bei
// Bedeutungslosigkeit der Eingabe exakt keine Pulse.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "SignalGenerator.hpp"
#include "iqclassify/iqclassify.hpp"

namespace {

int g_total = 0;
int g_ok = 0;

std::string joinLabels(std::vector<std::string> labels) {
    std::sort(labels.begin(), labels.end());
    std::string out;
    for (std::size_t i = 0; i < labels.size(); ++i) {
        if (i) out += "/";
        out += labels[i];
    }
    return out;
}

void printPulse(const iqclassify::PulseResult& r) {
    std::printf("        start=%zu n=%zu fc=%.5f bw=%.5f snr=%.1f "
                "cv=%.3f papr=%.1f flat=%.2f modes=%d conf=%.2f "
                "fine=[%zu,+%zu) fcRef=%.6f\n",
                r.startSample, r.numSamples, r.centerFreqNorm,
                r.bandwidthNorm, r.snrDb, r.envelopeCv, r.paprWindowedDb,
                r.spectralFlatnessDb, r.freqModes, r.confidence,
                r.startSampleFine, r.numSamplesFine, r.centerFreqRefined);
}

// Zyklischer Abstand zweier normierter Frequenzen in [0,1).
double cycDist(double a, double b) {
    double d = std::fmod(a - b, 1.0);
    if (d > 0.5) d -= 1.0;
    if (d < -0.5) d += 1.0;
    return std::abs(d);
}

// Winkeldistanz in [-pi, pi).
double wrapPi(double a) {
    a = std::fmod(a, 2.0 * M_PI);
    if (a > M_PI) a -= 2.0 * M_PI;
    if (a < -M_PI) a += 2.0 * M_PI;
    return a;
}

// Winkeldistanz modulo des M-PSK-Rasters 2*pi/M (Restambiguität).
double modGridDist(double a, double b, int M) {
    double d = wrapPi(a - b);
    const double g = 2.0 * M_PI / M;
    d = std::fmod(d, g);
    if (d > g / 2) d -= g;
    if (d < -g / 2) d += g;
    return std::abs(d);
}

bool scenario(const char* name, const std::vector<iqtest::cd>& buf,
              const std::string& expect) {
    ++g_total;
    iqclassify::ClassifierConfig cfg;
    iqclassify::ModulationClassifier clf(cfg);
    const std::vector<iqclassify::PulseResult> res =
        clf.classify(buf.data(), buf.size());
    std::vector<std::string> labels;
    labels.reserve(res.size());
    for (const auto& r : res)
        labels.push_back(iqclassify::toString(r.modulation));
    const std::string got = joinLabels(labels);
    const bool ok = got == expect && res.size() >= 1;
    if (ok) {
        ++g_ok;
        std::printf("[OK ] %-24s expect=%-24s n_pulses=%zu\n", name,
                    expect.c_str(), res.size());
    } else {
        std::printf("[FAIL] %-24s expect=%-24s got=%-24s n_pulses=%zu\n",
                    name, expect.c_str(), got.c_str(), res.size());
        for (const auto& r : res) printPulse(r);
    }
    return ok;
}

// Lokalisierungs-Szenario: ein einzelner Puls an bekannter Position.
// Prueft feine Zeitlage, Sub-Bin-Frequenz und extractPulse (Wiederfinden
// des Pulses im IQ-Puffer).
bool scenarioLocate(const char* name, const std::vector<iqtest::cd>& buf,
                    const char* expectLabel, double trueFc,
                    std::size_t trueAt, std::size_t trueLen) {
    ++g_total;
    iqclassify::ClassifierConfig cfg;
    iqclassify::ModulationClassifier clf(cfg);
    const std::vector<iqclassify::PulseResult> res =
        clf.classify(buf.data(), buf.size());
    bool ok = res.size() == 1 &&
              std::string(iqclassify::toString(res[0].modulation)) ==
                  expectLabel;
    if (ok) {
        const iqclassify::PulseResult& p = res[0];
        const double fd = cycDist(p.centerFreqRefined, trueFc);
        const std::ptrdiff_t sd = static_cast<std::ptrdiff_t>(
                                      p.startSampleFine) -
                                  static_cast<std::ptrdiff_t>(trueAt);
        const std::ptrdiff_t ld = static_cast<std::ptrdiff_t>(
                                      p.numSamplesFine) -
                                  static_cast<std::ptrdiff_t>(trueLen);
        const bool inRange = p.centerFreqRefined >= 0.0 &&
                             p.centerFreqRefined < 1.0;
        ok = fd <= 0.01 && inRange && std::abs(sd) <= 160 &&
             std::abs(ld) <= 320;
        // extractPulse: Slice deckt das feine Fenster, behaelt den Grossteil
        // der Energie und liegt spektral auf DC (Basisband).
        const std::vector<std::complex<double>> slice =
            clf.extractPulse(buf.data(), buf.size(), p);
        double eSlice = 0.0;
        for (const auto& v : slice) eSlice += std::norm(v);
        double eRaw = 0.0;
        for (std::size_t i = p.startSampleFine;
             i < p.startSampleFine + p.numSamplesFine && i < buf.size();
             ++i)
            eRaw += std::norm(buf[i]);
        bool exOk = slice.size() == p.numSamplesFine && eRaw > 0.0 &&
                    eSlice > 0.8 * eRaw;
        if (exOk && slice.size() >= 16) {
            std::vector<iqclassify::detail::cd> sp = slice;
            iqclassify::detail::fftInPlace(sp, false);
            std::size_t am = 0;
            double mv = -1.0;
            for (std::size_t k = 0; k < sp.size(); ++k) {
                const double v = std::norm(sp[k]);
                if (v > mv) {
                    mv = v;
                    am = k;
                }
            }
            double fpeak = static_cast<double>(am) /
                           static_cast<double>(sp.size());
            if (fpeak > 0.5) fpeak -= 1.0;
            exOk = std::abs(fpeak) <= 0.5 * p.bandwidthNorm + 0.02;
        }
        // float-Pfad und Nullzeiger-Haertung von extractPulse
        std::vector<std::complex<float>> fbuf(buf.size());
        for (std::size_t i = 0; i < buf.size(); ++i)
            fbuf[i] = std::complex<float>(static_cast<float>(buf[i].real()),
                                          static_cast<float>(buf[i].imag()));
        const std::vector<std::complex<float>> fslice =
            clf.extractPulseF(fbuf.data(), fbuf.size(), p);
        exOk = exOk && fslice.size() == p.numSamplesFine &&
               clf.extractPulse(nullptr, 0, p).empty() &&
               clf.extractPulseF(nullptr, 0, p).empty();
        ok = ok && exOk;
        if (ok) {
            ++g_ok;
            std::printf("[OK ] %-24s %s fine=[%zu,+%zu) fcRef=%.6f "
                        "(true at=%zu fc=%.4f)\n",
                        name, expectLabel, p.startSampleFine,
                        p.numSamplesFine, p.centerFreqRefined, trueAt,
                        trueFc);
        } else {
            std::printf("[FAIL] %-24s label=%s fine=[%zu,+%zu) fcRef=%.6f "
                        "(true at=%zu len=%zu fc=%.4f) extract=%s\n",
                        name,
                        res.empty() ? "-"
                                    : iqclassify::toString(res[0].modulation),
                        res.empty() ? 0 : res[0].startSampleFine,
                        res.empty() ? 0 : res[0].numSamplesFine,
                        res.empty() ? 0.0 : res[0].centerFreqRefined, trueAt,
                        trueLen, trueFc, exOk ? "ok" : "BAD");
            for (const auto& r : res) printPulse(r);
        }
    } else {
        std::printf("[FAIL] %-24s expect=%s got=%zu pulses\n", name,
                    expectLabel, res.size());
        for (const auto& r : res) printPulse(r);
    }
    return ok;
}

// Erwartet die Erkennung von exakt null Pulsen (Robustheitsfaelle).
bool scenarioEmpty(const char* name, const std::vector<iqtest::cd>& buf) {
    ++g_total;
    iqclassify::ClassifierConfig cfg;
    iqclassify::ModulationClassifier clf(cfg);
    const std::vector<iqclassify::PulseResult> res =
        clf.classify(buf.data(), buf.size());
    const bool ok = res.empty();
    if (ok) {
        ++g_ok;
        std::printf("[OK ] %-24s expect=0 pulses\n", name);
    } else {
        std::printf("[FAIL] %-24s expect=0 pulses        n_pulses=%zu\n",
                    name, res.size());
        for (const auto& r : res) printPulse(r);
    }
    return ok;
}

// Modulationsfamilie mit wahlweisem SNR und OFDM/Symbolparametern.
// Der OFDM-Referenzfall nutzt 40 Symbole: bei 20 Symbolen und Traglast
// 0.3 liegt der CP-Autokorrelationspeak zu knapp an der Schwelle
// (dokumentierte Grenze, siehe METHODOLOGY.md).
void runFamily(const char* family, std::mt19937& rng, double snrDb,
               bool includeQam256, bool includeFsk2, int ofdmNsym,
               double ofdmLoad, std::size_t qamNsym, int sps) {
    std::vector<iqtest::cd> bpsk = iqtest::genPsk(rng, 2, 400, sps);
    std::vector<iqtest::cd> qpsk = iqtest::genPsk(rng, 4, 400, sps);
    std::vector<iqtest::cd> fsk = iqtest::genFsk2(rng, 400, sps);
    std::vector<iqtest::cd> chirp = iqtest::genChirp(4000, 0.20);
    std::vector<iqtest::cd> cw = iqtest::genCw(rng, 3000);
    std::vector<iqtest::cd> ofdm =
        iqtest::genOfdm(rng, 256, 64, ofdmNsym, ofdmLoad);
    std::vector<iqtest::cd> qam16 = iqtest::genQam(rng, 16, qamNsym, sps);
    struct Item {
        const char* name;
        const std::vector<iqtest::cd>* sig;
        const char* label;
        double fc;
    };
    std::vector<Item> items = {
        {"BPSK", &bpsk, "BPSK", 0.05},
        {"QPSK", &qpsk, "QPSK", -0.05},
        {"CHIRP", &chirp, "Chirp", -0.10},
        {"CW", &cw, "Unmodulated", 0.15},
        {"OFDM", &ofdm, "OFDM", -0.15},
        {"QAM16", &qam16, "QAM16", 0.07},
    };
    if (includeFsk2) items.insert(items.begin() + 2,
                                  Item{"FSK2", &fsk, "FSK2", 0.05});
    for (const Item& it : items) {
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(*it.sig, it.fc, snrDb);
        rec.silence(2000);
        scenario((std::string(family) + "-" + it.name).c_str(),
                 rec.buffer(), it.label);
    }
    if (includeQam256) {
        std::vector<iqtest::cd> qam64 = iqtest::genQam(rng, 64, qamNsym, sps);
        std::vector<iqtest::cd> qam256 =
            iqtest::genQam(rng, 256, qamNsym, sps);
        iqtest::Recorder rec64;
        rec64.silence(2000);
        rec64.add(qam64, 0.07, snrDb);
        rec64.silence(2000);
        scenario((std::string(family) + "-QAM64").c_str(),
                 rec64.buffer(), "QAM64");
        iqtest::Recorder rec256;
        rec256.silence(2000);
        rec256.add(qam256, 0.07, snrDb);
        rec256.silence(2000);
        // QAM256 ist erst ab ~15-20 dB sicher von QAM64 trennbar; bei
        // 10 dB ist der QAM64-Fallback die dokumentierte Grenze.
        const char* expect256 = snrDb >= 15.0 ? "QAM256" : "QAM64";
        scenario((std::string(family) + "-QAM256").c_str(),
                 rec256.buffer(), expect256);
    }
}

}  // namespace

int main() {
    // 1. Basiskalibrierung (SNR 20 dB)
    {
        std::mt19937 rng(42);
        runFamily("base", rng, 20.0, true, true, 40, 0.3, 4000, 8);
    }
    // 2. SNR 15 dB: alle Typen
    {
        std::mt19937 rng(314);
        runFamily("snr15", rng, 15.0, true, true, 40, 0.3, 1000, 8);
    }
    // 3. SNR 10 dB: ohne FSK2 (erst ab ~12 dB zuverlaessig, siehe
    //    METHODOLOGY.md); QAM256-Fallback dokumentiert.
    {
        std::mt19937 rng(271);
        runFamily("snr10", rng, 10.0, true, false, 40, 0.3, 1000, 8);
    }
    // 4. Sequenzen
    {
        std::mt19937 rng(7);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::genPsk(rng, 4, 400, 8), 0.05, 20.0);
        rec.silence(2000);
        rec.add(iqtest::genFsk2(rng, 400, 8), -0.05, 20.0);
        rec.silence(2000);
        rec.add(iqtest::genCw(rng, 3000), 0.15, 20.0);
        rec.silence(2000);
        scenario("seq-3pulses", rec.buffer(),
                 "FSK2/QPSK/Unmodulated");
    }
    {
        // gleicher Typ, vertauschte Frequenzen
        std::mt19937 rng(6002);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::genPsk(rng, 2, 400, 8), 0.10, 20.0);
        rec.silence(2000);
        rec.add(iqtest::genPsk(rng, 2, 400, 8), -0.10, 20.0);
        rec.silence(2000);
        scenario("seq-2xbpsk", rec.buffer(), "BPSK/BPSK");
    }
    // 5. Gleichzeitige Pulse in getrennten Frequenzbaendern
    {
        std::mt19937 rng(11);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::genQam(rng, 16, 1000, 8), 0.15, 20.0, 2000);
        rec.add(iqtest::genPsk(rng, 2, 400, 8), -0.15, 20.0, 2000);
        rec.silence(6000);
        scenario("overlap-qam16-bpsk", rec.buffer(), "BPSK/QAM16");
    }
    {
        // QPSK und FSK2 gleichzeitig. Mittenabstand 0.30: FSK2 hat bei
        // sps 8 eine halbe Bandbreite von ~0.145, QPSK ~0.075; bei
        // Abstand 0.20 ueberlappen die halben Bandbreiten und die Baender
        // sind nach dokumentierter Grenze (METHODOLOGY.md) nicht trennbar.
        std::mt19937 rng(12);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::genPsk(rng, 4, 400, 8), 0.15, 20.0, 2000);
        rec.add(iqtest::genFsk2(rng, 400, 8), -0.15, 20.0, 2000);
        rec.silence(6000);
        scenario("overlap-qpsk-fsk2", rec.buffer(), "FSK2/QPSK");
    }
    {
        std::mt19937 rng(13);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::genCw(rng, 3200), 0.20, 20.0, 2000);
        rec.add(iqtest::genPsk(rng, 4, 400, 8), -0.20, 20.0, 2000);
        rec.silence(6000);
        scenario("overlap-cw-qpsk", rec.buffer(), "QPSK/Unmodulated");
    }
    {
        std::mt19937 rng(6100);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::genChirp(4000, 0.25), 0.0, 20.0, 2000);
        rec.add(iqtest::genCw(rng, 3200), 0.25, 20.0, 2000);
        rec.silence(6000);
        scenario("overlap-chirp-cw", rec.buffer(), "Chirp/Unmodulated");
    }
    {
        std::mt19937 rng(6104);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::genOfdm(rng, 256, 64, 40, 0.3), 0.18, 20.0, 2000);
        rec.add(iqtest::genPsk(rng, 4, 400, 8), -0.12, 20.0, 2000);
        rec.silence(10000);
        scenario("overlap-ofdm-qpsk", rec.buffer(), "OFDM/QPSK");
    }
    // 6. Oversampling (FSK2 nur sps 4: ab sps 16 zerfallen die
    //    Tonscheiben in Satelliten, dokumentierte Grenze)
    {
        std::mt19937 rng(2024);
        for (int sps : {4, 16, 32}) {
            for (auto& mc : std::vector<std::pair<int, const char*>>{
                     {2, "BPSK"}, {4, "QPSK"}}) {
                iqtest::Recorder rec;
                rec.silence(2000);
                rec.add(iqtest::genPsk(rng, mc.first, 400, sps), 0.05,
                        20.0);
                rec.silence(2000);
                scenario((std::string("sps") + std::to_string(sps) + "-" +
                          mc.second)
                             .c_str(),
                         rec.buffer(), mc.second);
            }
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genQam(rng, 16, 1000, sps), 0.05, 20.0);
            rec.silence(2000);
            scenario((std::string("sps") + std::to_string(sps) + "-QAM16")
                         .c_str(),
                     rec.buffer(), "QAM16");
        }
        {
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genFsk2(rng, 400, 4), 0.05, 20.0);
            rec.silence(2000);
            scenario("sps4-FSK2", rec.buffer(), "FSK2");
        }
    }
    // 7. Frequenzlagen (SNR 20 dB)
    {
        std::mt19937 rng(61);
        const double fcs[] = {0.0, 0.15, 0.30, 0.45, -0.30};
        for (double fc : fcs) {
            char name[64];
            std::snprintf(name, sizeof(name), "freq%+.2f-qpsk", fc);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genPsk(rng, 4, 400, 8), fc, 20.0);
            rec.silence(2000);
            scenario(name, rec.buffer(), "QPSK");

            std::snprintf(name, sizeof(name), "freq%+.2f-cw", fc);
            iqtest::Recorder recCw;
            recCw.silence(2000);
            recCw.add(iqtest::genCw(rng, 3000), fc, 20.0);
            recCw.silence(2000);
            scenario(name, recCw.buffer(), "Unmodulated");

            std::snprintf(name, sizeof(name), "freq%+.2f-chirp", fc);
            iqtest::Recorder recChirp;
            recChirp.silence(2000);
            recChirp.add(iqtest::genChirp(4000, 0.25), fc, 20.0);
            recChirp.silence(2000);
            scenario(name, recChirp.buffer(), "Chirp");
        }
    }
    // 8. Chirp-Varianten
    {
        std::mt19937 rng(62);
        for (int upI = 0; upI < 2; ++upI) {
            const bool up = upI == 0;
            for (double bw : {0.10, 0.25, 0.40}) {
                for (double snr : {20.0, 10.0}) {
                    char name[64];
                    std::snprintf(name, sizeof(name), "chirp-%s-bw%.2f-%ddB",
                                  up ? "up" : "dn", bw,
                                  static_cast<int>(snr));
                    iqtest::Recorder rec;
                    rec.silence(2000);
                    rec.add(iqtest::genChirp(4000, bw, up), -0.10, snr);
                    rec.silence(2000);
                    scenario(name, rec.buffer(), "Chirp");
                }
            }
        }
    }
    // 9. Pulslaengen (SNR 20 dB)
    {
        std::mt19937 rng(63);
        for (std::size_t nsym : {150, 1000}) {
            char name[64];
            std::snprintf(name, sizeof(name), "len-nsym%zu", nsym);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genPsk(rng, 4, nsym, 8), -0.05, 20.0);
            rec.silence(2000);
            scenario(name, rec.buffer(), "QPSK");
        }
        for (std::size_t nsamp : {500, 8000}) {
            char name[64];
            std::snprintf(name, sizeof(name), "len-cw%zu", nsamp);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genCw(rng, nsamp), 0.15, 20.0);
            rec.silence(2000);
            scenario(name, rec.buffer(), "Unmodulated");
        }
    }
    // 10. Seed-Streuung
    {
        for (int i = 0; i < 3; ++i) {
            std::mt19937 rng(7000 + i);
            char name[64];
            std::snprintf(name, sizeof(name), "seed%d-qpsk-10dB", i);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genPsk(rng, 4, 400, 8), -0.05, 10.0);
            rec.silence(2000);
            scenario(name, rec.buffer(), "QPSK");

            std::mt19937 rng2(7100 + i);
            std::snprintf(name, sizeof(name), "seed%d-qam16-20dB", i);
            iqtest::Recorder rec2;
            rec2.silence(2000);
            rec2.add(iqtest::genQam(rng2, 16, 1000, 8), 0.07, 20.0);
            rec2.silence(2000);
            scenario(name, rec2.buffer(), "QAM16");
        }
    }
    // 11. Haertungs- und Robustheitsfaelle
    {
        iqclassify::ClassifierConfig cfg;
        iqclassify::ModulationClassifier clf(cfg);
        // Nullzeiger und Laenge 0: definiert leer
        ++g_total;
        if (clf.classify(nullptr, 0).empty() &&
            clf.classify(nullptr, 1234).empty()) {
            ++g_ok;
            std::printf("[OK ] %-24s expect=0 pulses\n", "null-input");
        } else {
            std::printf("[FAIL] %-24s expect=0 pulses\n", "null-input");
        }
        // leere Eingabe, Miniatureingabe
        scenarioEmpty("empty-buffer", {});
        scenarioEmpty("tiny-buffer", {iqtest::cd(1.0, 2.0),
                                      iqtest::cd(3.0, 1.0),
                                      iqtest::cd(2.0, 2.0)});
        // reines Rauschen
        {
            std::mt19937 rng(404);
            std::normal_distribution<double> nd(0.0, 1.0);
            std::vector<iqtest::cd> buf(8192);
            for (auto& v : buf) v = iqtest::cd(nd(rng), nd(rng));
            scenarioEmpty("pure-noise", buf);
        }
        // alles Null / konstanter DC-Puffer
        scenarioEmpty("all-zero", std::vector<iqtest::cd>(8192));
        scenarioEmpty("dc-constant",
                      std::vector<iqtest::cd>(8192, iqtest::cd(1.0, 0.0)));
        // kompletter NaN/Inf-Puffer
        {
            const double qnan = std::numeric_limits<double>::quiet_NaN();
            const double inf = std::numeric_limits<double>::infinity();
            scenarioEmpty("all-nan",
                          std::vector<iqtest::cd>(4096,
                                                  iqtest::cd(qnan, qnan)));
            scenarioEmpty("all-inf",
                          std::vector<iqtest::cd>(4096,
                                                  iqtest::cd(inf, inf)));
        }
        // NaN/Inf/Dropout mitten in einem Puls: wirken wie Stille;
        // der Puls muss weiterhin erkannt werden.
        {
            const double qnan = std::numeric_limits<double>::quiet_NaN();
            const double inf = std::numeric_limits<double>::infinity();
            for (int k = 0; k < 3; ++k) {
                std::mt19937 rng(90 + k);
                iqtest::Recorder rec;
                rec.silence(2000);
                rec.add(iqtest::genPsk(rng, 4, 400, 8), -0.05, 20.0);
                rec.silence(2000);
                std::vector<iqtest::cd> buf = rec.buffer();
                const char* nm = "nan-dropout";
                for (std::size_t i = 2500; i < 2600 && i < buf.size(); ++i) {
                    if (k == 0) buf[i] = iqtest::cd(qnan, qnan);
                    if (k == 1) buf[i] = iqtest::cd(inf, 0.0);
                    if (k == 2) buf[i] = iqtest::cd(0.0, 0.0);
                }
                if (k == 1) nm = "inf-dropout";
                if (k == 2) nm = "zero-dropout";
                scenario(nm, buf, "QPSK");
            }
        }
        // NaN-Salat ueber den ganzen Puffer: darf nicht abstuerzen und
        // darf keine Geisterpulse erzeugen.
        {
            const double qnan = std::numeric_limits<double>::quiet_NaN();
            std::mt19937 rng(404);
            std::normal_distribution<double> nd(0.0, 1.0);
            std::vector<iqtest::cd> buf(8192);
            for (std::size_t i = 0; i < buf.size(); ++i) {
                if (i % 37 == 0)
                    buf[i] = iqtest::cd(qnan, qnan);
                else
                    buf[i] = iqtest::cd(nd(rng), nd(rng));
            }
            scenarioEmpty("nan-sprinkled-noise", buf);
        }
        // float-Pfad (classifyF) auf dem Basisszenario
        {
            std::mt19937 rng(4242);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genPsk(rng, 4, 400, 8), -0.05, 20.0);
            rec.silence(2000);
            const std::vector<iqtest::cd>& buf = rec.buffer();
            std::vector<std::complex<float>> fbuf(buf.size());
            for (std::size_t i = 0; i < buf.size(); ++i)
                fbuf[i] = std::complex<float>(
                    static_cast<float>(buf[i].real()),
                    static_cast<float>(buf[i].imag()));
            ++g_total;
            const std::vector<iqclassify::PulseResult> res =
                clf.classifyF(fbuf.data(), fbuf.size());
            std::vector<std::string> labels;
            for (const auto& r : res)
                labels.push_back(iqclassify::toString(r.modulation));
            const bool okf =
                joinLabels(labels) == "QPSK" && res.size() == 1;
            if (okf) {
                ++g_ok;
                std::printf("[OK ] %-24s expect=%-24s\n", "float-path",
                            "QPSK");
            } else {
                std::printf("[FAIL] %-24s expect=%-24s\n", "float-path",
                            "QPSK");
            }
        }
        // kurzer CW-Puls (1000 Samples): ohne Resttraeger-Entfernung
        // haette die Traegerkohaerenz auf dieser Laenge versagt.
        {
            std::mt19937 rng(909);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genCw(rng, 1000), 0.15, 20.0);
            rec.silence(2000);
            scenario("cw-short-1000", rec.buffer(), "Unmodulated");
        }
    }
    // 12. Lokalisierung: detektierte Pulse zeitlich (samplegenau) und
    //     frequenztechnisch (Sub-Bin) im IQ-Signal wiederfinden
    {
        std::mt19937 rng(500);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::genPsk(rng, 4, 400, 8), 0.15, 20.0);
        rec.silence(2000);
        scenarioLocate("locate-qpsk", rec.buffer(), "QPSK", 0.15, 2000,
                       3200);
    }
    {
        std::mt19937 rng(501);
        iqtest::Recorder rec;
        rec.silence(3000);
        rec.add(iqtest::genCw(rng, 2500), -0.20, 20.0);
        rec.silence(2000);
        scenarioLocate("locate-cw-neg-fc", rec.buffer(), "Unmodulated",
                       -0.20, 3000, 2500);
    }
    {
        std::mt19937 rng(503);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::genChirp(4000, 0.25), 0.05, 20.0);
        rec.silence(2000);
        scenarioLocate("locate-chirp", rec.buffer(), "Chirp", 0.05, 2000,
                       4000);
    }
    {
        std::mt19937 rng(504);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::genFsk2(rng, 400, 8), 0.20, 20.0);
        rec.silence(2000);
        scenarioLocate("locate-fsk2", rec.buffer(), "FSK2", 0.20, 2000,
                       3200);
    }
    {
        std::mt19937 rng(505);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::genOfdm(rng, 256, 64, 40, 0.3), -0.15, 20.0);
        rec.silence(2000);
        scenarioLocate("locate-ofdm", rec.buffer(), "OFDM", -0.15, 2000,
                       40 * 320);
    }
    {
        // zwei sequenzielle Pulse: beide muessen lokalisierbar sein,
        // die feinen Fenster duerfen sich nicht ueberlappen und muessen
        // den jeweiligen wahren Puls zu >= 80 % enthalten.
        std::mt19937 rng(502);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::genPsk(rng, 2, 300, 8), 0.10, 20.0);  // [2000,4400)
        rec.silence(2500);
        rec.add(iqtest::genCw(rng, 2000), -0.15, 20.0);       // [6900,8900)
        rec.silence(2000);
        ++g_total;
        iqclassify::ClassifierConfig cfg;
        iqclassify::ModulationClassifier clf(cfg);
        const std::vector<iqclassify::PulseResult> res =
            clf.classify(rec.buffer().data(), rec.buffer().size());
        bool ok = res.size() == 2;
        if (ok) {
            std::vector<iqclassify::PulseResult> sorted = res;
            std::sort(sorted.begin(), sorted.end(),
                      [](const iqclassify::PulseResult& a,
                         const iqclassify::PulseResult& b) {
                          return a.startSampleFine < b.startSampleFine;
                      });
            const auto& res2 = sorted;
            const std::pair<double, std::pair<std::size_t, std::size_t>>
                truth[2] = {{0.10, {2000, 4400}}, {-0.15, {6900, 8900}}};
            for (std::size_t i = 0; i < 2 && ok; ++i) {
                const iqclassify::PulseResult& p = res2[i];
                const double fd = cycDist(p.centerFreqRefined, truth[i].first);
                const std::size_t lo = std::max(p.startSampleFine,
                                                truth[i].second.first);
                const std::size_t hi = std::min(
                    p.startSampleFine + p.numSamplesFine,
                    truth[i].second.second);
                const double cover =
                    hi > lo ? static_cast<double>(hi - lo) /
                                  static_cast<double>(
                                      truth[i].second.second -
                                      truth[i].second.first)
                            : 0.0;
                ok = fd <= 0.01 && cover >= 0.8;
            }
            ok = ok && res2[0].startSampleFine + res2[0].numSamplesFine <=
                           res2[1].startSampleFine;
        }
        if (ok) {
            ++g_ok;
            std::printf("[OK ] %-24s 2 Pulse lokalisiert, Fenster "
                        "disjunkt\n",
                        "locate-2pulses");
        } else {
            std::printf("[FAIL] %-24s n_pulses=%zu\n", "locate-2pulses",
                        res.size());
            for (const auto& r : res) printPulse(r);
        }
    }
    // 13. Amplitudenschwankende Pulse (Fading): ein Puls mit schwankender
    //     Amplitude darf nicht in mehrere Pulse zerlegt werden (Fade-Merge)
    //     und muss korrekt klassifiziert/lokaliert bleiben.
    {
        std::mt19937 rng(600);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::applyEnv(iqtest::genPsk(rng, 4, 400, 8), 0, 0.7,
                                 rng),
                0.15, 20.0);
        rec.silence(2000);
        scenarioLocate("fade-qpsk-rand", rec.buffer(), "QPSK", 0.15, 2000,
                       3200);
    }
    {
        std::mt19937 rng(601);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::applyEnv(iqtest::genPsk(rng, 4, 400, 8), 2, 0.7,
                                 rng),
                -0.10, 20.0);
        rec.silence(2000);
        scenarioLocate("fade-qpsk-mid", rec.buffer(), "QPSK", -0.10, 2000,
                       3200);
    }
    {
        // tiefe Mitten-Fade (-20 dB) bei hohem SNR: Luecke bleibt ueber dem
        // Rauschboden -> Fade-Merge liefert genau einen Puls.
        std::mt19937 rng(602);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::applyEnv(iqtest::genCw(rng, 3000), 2, 0.9, rng),
                -0.20, 40.0);
        rec.silence(2000);
        scenarioLocate("fade-cw-deep-40dB", rec.buffer(), "Unmodulated",
                       -0.20, 2000, 3000);
    }
    {
        std::mt19937 rng(603);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::applyEnv(iqtest::genCw(rng, 3000), 1, 0.5, rng),
                0.20, 20.0);
        rec.silence(2000);
        scenarioLocate("fade-cw-sinus6dB", rec.buffer(), "Unmodulated",
                       0.20, 2000, 3000);
    }
    {
        // FSK2 mit langsamer 6-dB-Amplitudenschwankung (Sinus): die beiden
        // Tonscheiben bleiben als ein Band erkennbar.
        std::mt19937 rng(604);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::applyEnv(iqtest::genFsk2(rng, 400, 8), 1, 0.5, rng),
                0.05, 20.0);
        rec.silence(2000);
        scenario("fade-fsk2-sinus", rec.buffer(), "FSK2");
    }
    {
        // Rampe (leiser Anfang ist am Detektionslimit): Puls wird als
        // genau ein Puls erkannt und korrekt klassifiziert.
        std::mt19937 rng(605);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::applyEnv(iqtest::genPsk(rng, 2, 400, 8), 1, 0.8,
                                 rng),
                0.10, 20.0);
        rec.silence(2000);
        scenario("fade-bpsk-ramp", rec.buffer(), "BPSK");
    }
    // 14. Low-SNR (bis 2 dB): lange Pulse (mehr Symbole/Energie) werden
    //     erkannt, lokalisiert und klassifiziert; Default-Gate minSnrDb=2.
    {
        // Hilfspruefung: genau ein Puls, korrekte Frequenzlage und Zeitlage.
        auto lowSnrCheck = [](const char* name,
                              const std::vector<iqtest::cd>& buf,
                              const char* expectLabel, double trueFc,
                              std::size_t trueAt) {
            ++g_total;
            iqclassify::ClassifierConfig cfg;
            iqclassify::ModulationClassifier clf(cfg);
            const auto res = clf.classify(buf.data(), buf.size());
            bool ok = res.size() == 1 &&
                      std::string(iqclassify::toString(res[0].modulation)) ==
                          expectLabel &&
                      cycDist(res[0].centerFreqRefined, trueFc) <= 0.02 &&
                      std::abs(static_cast<std::ptrdiff_t>(
                                   res[0].startSampleFine) -
                               static_cast<std::ptrdiff_t>(trueAt)) <= 200;
            if (ok) {
                ++g_ok;
                std::printf("[OK ] %-24s %s f=%.4f @%zu (snr est %.1f dB)\n",
                            name, expectLabel, res[0].centerFreqRefined,
                            res[0].startSampleFine, res[0].snrDb);
            } else {
                std::printf("[FAIL] %-24s expect=%s@%.3f/%zu got n=%zu\n", name,
                            expectLabel, trueFc, trueAt, res.size());
                for (const auto& r : res) printPulse(r);
            }
            return ok;
        };
        {
            std::mt19937 rng(810);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genCw(rng, 8000), 0.20, 2.0);
            rec.silence(2500);
            lowSnrCheck("lowsnr-cw-2dB", rec.buffer(), "Unmodulated", 0.20,
                        2000);
        }
        {
            std::mt19937 rng(811);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genChirp(4000, 0.25), -0.10, 2.0);
            rec.silence(2500);
            lowSnrCheck("lowsnr-chirp-2dB", rec.buffer(), "Chirp", -0.10,
                        2000);
        }
        {
            std::mt19937 rng(813);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genPsk(rng, 4, 1000, 8), 0.15, 2.0);
            rec.silence(2500);
            lowSnrCheck("lowsnr-qpsk-2dB", rec.buffer(), "QPSK", 0.15, 2000);
        }
        {
            std::mt19937 rng(813);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genPsk(rng, 2, 1000, 8), 0.10, 2.0);
            rec.silence(2500);
            lowSnrCheck("lowsnr-bpsk-2dB", rec.buffer(), "BPSK", 0.10, 2000);
        }
        {
            std::mt19937 rng(817);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genQam(rng, 16, 2000, 8), 0.07, 4.0);
            rec.silence(2500);
            lowSnrCheck("lowsnr-qam16-4dB", rec.buffer(), "QAM16", 0.07,
                        2000);
        }
        {
            // zwei Pulse bei niedrigem SNR in getrennten Baendern
            std::mt19937 rng(815);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genCw(rng, 6000), 0.25, 3.0, 2000);
            rec.add(iqtest::genPsk(rng, 4, 800, 8), -0.20, 5.0, 2000);
            rec.silence(6000);
            scenario("lowsnr-2pulses", rec.buffer(), "QPSK/Unmodulated");
        }
    }
    // 15. Analyse-API: PSK-Synchronisation (phasenrichtiger Pulsanfang),
    //     Chirp-Parameter, Symbol-Schnittstelle (PSK/FSK)
    {
        // PSK-Signal mit bekanntem Traeger fc und Startphase phi0 aufbauen.
        auto makePsk = [](int m, int seed, double fc, double phi0,
                          std::size_t nsym, int sps) {
            std::mt19937 rng(seed);
            const std::vector<iqtest::cd> psk =
                iqtest::genPsk(rng, m, nsym, sps);
            std::vector<iqtest::cd> x(psk.size());
            for (std::size_t i = 0; i < psk.size(); ++i) {
                const double ph =
                    2.0 * M_PI * fc * static_cast<double>(2000 + i) + phi0;
                x[i] = psk[i] * iqtest::cd(std::cos(ph), std::sin(ph));
            }
            return x;
        };
        // Phasen-Faltung: Konzentration der M-ten Potenz der Symbole.
        // Liegt nahe 1, wenn alle Symbole phasenrichtig auf dem M-PSK-Raster
        // liegen (unabhaengig von etwa durch das MF-Gate verworfenen Symbolen).
        auto foldConc = [](const std::vector<std::complex<double>>& syms,
                           int m) {
            if (syms.empty())
                return 0.0;
            std::complex<double> acc(0.0, 0.0);
            for (const auto& v : syms) {
                std::complex<double> p = v;
                for (int i = 1; i < m; ++i)
                    p *= v;
                acc += p;
            }
            return std::abs(acc) / static_cast<double>(syms.size());
        };
        // PSK-Synchronisation: BPSK und QPSK mit bekanntem Traeger/Phase.
        for (int mc = 0; mc < 2; ++mc) {
            const int m = mc == 0 ? 2 : 4;
            const char* label = mc == 0 ? "BPSK" : "QPSK";
            const char* name = mc == 0 ? "sync-bpsk" : "sync-qpsk";
            const double fc = mc == 0 ? 0.12 : -0.08;
            const double phi0 = mc == 0 ? 0.7 : -1.1;
            const int seed = mc == 0 ? 4242 : 777;
            std::vector<iqtest::cd> x = makePsk(m, seed, fc, phi0, 400, 8);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(x, 0.0, 20.0, 2000);
            rec.silence(2000);
            ++g_total;
            iqclassify::ClassifierConfig cfg;
            iqclassify::ModulationClassifier clf(cfg);
            const auto res =
                clf.classify(rec.buffer().data(), rec.buffer().size());
            bool ok = res.size() == 1 &&
                      std::string(iqclassify::toString(res[0].modulation)) ==
                          label;
            if (ok) {
                const iqclassify::PskSync s =
                    clf.synchronizePsk(rec.buffer().data(),
                                       rec.buffer().size(), res[0]);
                const double expectedPhi =
                    phi0 + 2.0 * M_PI * fc *
                               static_cast<double>(res[0].startSampleFine);
                const double phiErr =
                    modGridDist(s.carrierPhaseAtStart, expectedPhi, m);
                const double fold = foldConc(s.symbols, m);
                ok = s.valid && std::abs(s.sps - 8.0) < 0.6 &&
                     s.timingOffset >= 0.0 && s.timingOffset < s.sps &&
                     cycDist(s.carrierFreqNorm, fc) < 0.01 && phiErr < 0.35 &&
                     s.symbols.size() >= 300 && fold > 0.9;
                if (ok) {
                    ++g_ok;
                    std::printf("[OK ] %-24s sps=%.2f fc=%.4f phiErr=%.3f "
                                "fold=%.2f nsym=%zu\n",
                                name, s.sps, s.carrierFreqNorm, phiErr, fold,
                                s.symbols.size());
                } else {
                    std::printf("[FAIL] %-24s valid=%d sps=%.3f off=%.2f "
                                "fc=%.5f phiErr=%.3f fold=%.2f nsym=%zu\n",
                                name, (int)s.valid, s.sps, s.timingOffset,
                                s.carrierFreqNorm, phiErr, fold,
                                s.symbols.size());
                }
            } else {
                std::printf("[FAIL] %-24s classify n=%zu\n", name, res.size());
                for (const auto& r : res) printPulse(r);
            }
        }
        // Chirp-Parameter (up und down)
        for (int up = 0; up < 2; ++up) {
            const char* name = up == 0 ? "chirp-params-up" : "chirp-params-dn";
            const double fc = up == 0 ? -0.10 : 0.05;
            const double bw = 0.25;
            std::mt19937 rng(62 + up);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genChirp(4000, bw, up == 0), fc, 20.0);
            rec.silence(2000);
            ++g_total;
            iqclassify::ClassifierConfig cfg;
            iqclassify::ModulationClassifier clf(cfg);
            const auto res =
                clf.classify(rec.buffer().data(), rec.buffer().size());
            bool ok = res.size() == 1 &&
                      res[0].modulation == iqclassify::ModulationType::Chirp;
            if (ok) {
                const iqclassify::ChirpParams c =
                    clf.analyzeChirp(rec.buffer().data(), rec.buffer().size(),
                                     res[0]);
                const double fStartTrue =
                    up == 0 ? fc - bw / 2 : fc + bw / 2;
                const double fEndTrue = up == 0 ? fc + bw / 2 : fc - bw / 2;
                const double slopeTrue = bw / 4000.0 * (up == 0 ? 1.0 : -1.0);
                ok = c.valid && c.up == (up == 0) && c.r2 > 0.95 &&
                     cycDist(c.fStartNorm, fStartTrue) < 0.01 &&
                     cycDist(c.fEndNorm, fEndTrue) < 0.01 &&
                     std::abs(c.bandwidthNorm - bw) < 0.03 &&
                     std::abs(c.slopeNorm / slopeTrue - 1.0) < 0.15 &&
                     std::abs(static_cast<std::ptrdiff_t>(c.startSample) -
                              2000) <= 200;
                if (ok) {
                    ++g_ok;
                    std::printf("[OK ] %-24s f0=%.4f f1=%.4f bw=%.4f "
                                "slope=%.3e r2=%.3f @%zu\n",
                                name, c.fStartNorm, c.fEndNorm, c.bandwidthNorm,
                                c.slopeNorm, c.r2, c.startSample);
                } else {
                    std::printf("[FAIL] %-24s valid=%d up=%d r2=%.3f f0=%.4f "
                                "f1=%.4f bw=%.4f slope=%.3e start=%zu\n",
                                name, (int)c.valid, (int)c.up, c.r2, c.fStartNorm,
                                c.fEndNorm, c.bandwidthNorm, c.slopeNorm,
                                c.startSample);
                }
            } else {
                std::printf("[FAIL] %-24s classify n=%zu\n", name, res.size());
                for (const auto& r : res) printPulse(r);
            }
        }
        // Symbol-Schnittstelle FSK2: Takt, Frequenzen, Ton-Entscheidungen
        {
            std::mt19937 rng(31);
            const std::vector<iqtest::cd> fsk = iqtest::genFsk2(rng, 400, 8);
            std::mt19937 rngTx(31);
            std::uniform_real_distribution<double> u01(0.0, 1.0);
            std::vector<int> bits(400);
            for (int k = 0; k < 400; ++k) bits[k] = u01(rngTx) < 0.5 ? 1 : 0;
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(fsk, 0.05, 20.0, 2000);
            rec.silence(2000);
            ++g_total;
            iqclassify::ClassifierConfig cfg;
            iqclassify::ModulationClassifier clf(cfg);
            const auto res =
                clf.classify(rec.buffer().data(), rec.buffer().size());
            bool ok = res.size() == 1 &&
                      res[0].modulation == iqclassify::ModulationType::FSK2;
            if (ok) {
                const iqclassify::SymbolStream ss =
                    clf.extractSymbols(rec.buffer().data(),
                                       rec.buffer().size(), res[0]);
                double m1 = 0, m0 = 0;
                int c1 = 0, c0 = 0;
                int best = 0;
                std::size_t tot = 0;
                for (std::size_t k = 0; k < ss.symbolFreqs.size(); ++k) {
                    if (ss.fskTone[k]) {
                        m1 += ss.symbolFreqs[k];
                        ++c1;
                    } else {
                        m0 += ss.symbolFreqs[k];
                        ++c0;
                    }
                }
                const double sep =
                    (c1 > 0 && c0 > 0) ? m1 / c1 - m0 / c0 : 0.0;
                for (int sh = -8; sh <= 8; ++sh) {
                    int match = 0;
                    std::size_t t = 0;
                    for (std::size_t k = 0; k < ss.fskTone.size(); ++k) {
                        const long b = static_cast<long>(k) + sh;
                        if (b >= 0 && b < 400) {
                            ++t;
                            if (ss.fskTone[k] == bits[static_cast<std::size_t>(b)])
                                ++match;
                        }
                    }
                    if (t > 0 && match > best) {
                        best = match;
                        tot = t;
                    }
                }
                const double bitMatch =
                    tot ? static_cast<double>(best) / static_cast<double>(tot)
                        : 0.0;
                ok = ss.valid && ss.type == iqclassify::ModulationType::FSK2 &&
                     std::abs(ss.sps - 8.0) < 0.6 && ss.symbolFreqs.size() >= 350 &&
                     std::abs(sep - 0.125) < 0.03 && bitMatch >= 0.9;
                if (ok) {
                    ++g_ok;
                    std::printf("[OK ] %-24s sps=%.2f nsym=%zu sep=%.4f "
                                "bitmatch=%.2f\n",
                                "symbols-fsk2", ss.sps, ss.symbolFreqs.size(),
                                sep, bitMatch);
                } else {
                    std::printf("[FAIL] %-24s valid=%d sps=%.3f nsym=%zu "
                                "sep=%.4f bitmatch=%.2f\n",
                                "symbols-fsk2", (int)ss.valid, ss.sps,
                                ss.symbolFreqs.size(), sep, bitMatch);
                }
            } else {
                std::printf("[FAIL] %-24s classify n=%zu\n", "symbols-fsk2",
                            res.size());
                for (const auto& r : res) printPulse(r);
            }
        }
        // Symbol-Schnittstelle PSK: Symbole direkt verfuegbar
        {
            std::vector<iqtest::cd> x = makePsk(4, 777, -0.08, -1.1, 400, 8);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(x, 0.0, 20.0, 2000);
            rec.silence(2000);
            ++g_total;
            iqclassify::ClassifierConfig cfg;
            iqclassify::ModulationClassifier clf(cfg);
            const auto res =
                clf.classify(rec.buffer().data(), rec.buffer().size());
            bool ok = res.size() == 1 &&
                      res[0].modulation == iqclassify::ModulationType::QPSK;
            if (ok) {
                const iqclassify::SymbolStream ss =
                    clf.extractSymbols(rec.buffer().data(),
                                       rec.buffer().size(), res[0]);
                const double fold = foldConc(ss.symbols, 4);
                ok = ss.valid &&
                     ss.type == iqclassify::ModulationType::QPSK &&
                     ss.symbols.size() >= 300 && fold > 0.9;
                if (ok) {
                    ++g_ok;
                    std::printf("[OK ] %-24s nsym=%zu fold=%.2f sps=%.2f\n",
                                "symbols-qpsk", ss.symbols.size(), fold,
                                ss.sps);
                } else {
                    std::printf("[FAIL] %-24s valid=%d nsym=%zu fold=%.2f\n",
                                "symbols-qpsk", (int)ss.valid, ss.symbols.size(),
                                fold);
                }
            } else {
                std::printf("[FAIL] %-24s classify n=%zu\n", "symbols-qpsk",
                            res.size());
            }
        }
        // Auf CW-Puls sind die Analyse-Methoden definiert invalid.
        {
            std::mt19937 rng(5);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genCw(rng, 3000), 0.20, 20.0);
            rec.silence(2000);
            ++g_total;
            iqclassify::ClassifierConfig cfg;
            iqclassify::ModulationClassifier clf(cfg);
            const auto res =
                clf.classify(rec.buffer().data(), rec.buffer().size());
            bool ok = res.size() == 1;
            if (ok) {
                ok = !clf.synchronizePsk(rec.buffer().data(),
                                         rec.buffer().size(), res[0])
                          .valid &&
                     !clf.analyzeChirp(rec.buffer().data(), rec.buffer().size(),
                                       res[0])
                          .valid &&
                     !clf.extractSymbols(rec.buffer().data(),
                                         rec.buffer().size(), res[0])
                          .valid;
            }
            if (ok) {
                ++g_ok;
                std::printf("[OK ] %-24s auf CW-Puls definiert invalid\n",
                            "api-invalid-on-cw");
            } else {
                std::printf("[FAIL] %-24s\n", "api-invalid-on-cw");
            }
        }
    }
    // 16. Spektrale Kanalisierung: adaptive Kanaele um die Energiebereiche.
    {
        // QPSK@0.15 und CW@-0.20 gleichzeitig: zwei Kanaele um die
        // Energiebereiche, Mitten nahe den wahren fc.
        std::mt19937 rng(11);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::genPsk(rng, 4, 400, 8), 0.15, 20.0);
        rec.add(iqtest::genCw(rng, 3000), -0.20, 20.0, 2000);
        rec.silence(2000);
        const std::vector<iqtest::cd>& buf = rec.buffer();
        ++g_total;
        {
            iqclassify::ClassifierConfig cfg;
            iqclassify::Channelizer chz(cfg);
            const auto chans =
                chz.channelize(buf.data(), buf.size(), 1);
            bool ok = chans.size() == 2;
            bool hasQ = false;
            bool hasC = false;
            for (const auto& ch : chans) {
                if (ch.y.size() != buf.size()) ok = false;
                if (cycDist(ch.centerNorm, 0.15) < 0.03) hasQ = true;
                if (cycDist(ch.centerNorm, 0.80) < 0.03) hasC = true;
            }
            ok = ok && hasQ && hasC;
            if (ok) {
                ++g_ok;
                std::printf("[OK ] %-24s nchan=%zu (%.4f, %.4f)\n",
                            "chan-centers", chans.size(), chans[0].centerNorm,
                            chans[1].centerNorm);
            } else {
                std::printf("[FAIL] %-24s nchan=%zu hasQ=%d hasC=%d\n",
                            "chan-centers", chans.size(), (int)hasQ, (int)hasC);
            }
        }
        // Extract/Mix: im CW-Kanal liegt der Ton nahe DC (Mittenfrequenz ~0).
        ++g_total;
        {
            iqclassify::ClassifierConfig cfg;
            iqclassify::Channelizer chz(cfg);
            const auto chans = chz.channelize(buf.data(), buf.size(), 1);
            bool ok = false;
            for (const auto& ch : chans) {
                if (cycDist(ch.centerNorm, 0.80) >= 0.03) continue;
                iqtest::cd acc(0.0, 0.0);
                for (std::size_t i = 2000; i + 1 < 5000 && i + 1 < ch.y.size();
                     ++i)
                    acc += ch.y[i + 1] * std::conj(ch.y[i]);
                const double f = std::arg(acc) / (2.0 * M_PI);
                if (std::fabs(f) < 0.02) {
                    ok = true;
                    ++g_ok;
                    std::printf("[OK ] %-24s CW im Kanal nahe DC (f_inst=%.5f)\n",
                                "chan-extract-mix", f);
                    break;
                }
            }
            if (!ok) {
                std::printf("[FAIL] %-24s Ton nicht nahe DC im Kanal\n",
                            "chan-extract-mix");
            }
        }
        // Resample: Kanal mit decimate=4 -> Laenge n/4, Puls weiterhin da.
        ++g_total;
        {
            iqclassify::ClassifierConfig cfg;
            iqclassify::Channelizer chz(cfg);
            const auto chans = chz.channelize(buf.data(), buf.size(), 1);
            bool ok = false;
            for (const auto& ch : chans) {
                if (cycDist(ch.centerNorm, 0.80) >= 0.03) continue;
                auto ch4 = chz.extract(buf.data(), buf.size(), ch.band,
                                      ch.nft, 4);
                double e = 0.0;
                for (std::size_t m = 500; m < 1250 && m < ch4.y.size(); ++m)
                    e += std::norm(ch4.y[m]);
                ok = ch4.decimate == 4 &&
                     ch4.y.size() == (buf.size() + 3) / 4 && e > 100.0;
                if (ok) break;
            }
            if (ok) {
                ++g_ok;
                std::printf("[OK ] %-24s decimate=4 -> Laenge n/4, Puls drin\n",
                            "chan-resample");
            } else {
                std::printf("[FAIL] %-24s Resample-Faktor/Laenge/Puls\n",
                            "chan-resample");
            }
        }
        // Aequivalenz: kanalisierter Default-Pfad == Legacy-Pfad (Labels).
        ++g_total;
        {
            std::mt19937 rng2(7);
            iqtest::Recorder rec2;
            rec2.silence(2000);
            rec2.add(iqtest::genPsk(rng2, 4, 400, 8), 0.05, 20.0);
            rec2.silence(2000);
            rec2.add(iqtest::genFsk2(rng2, 400, 8), -0.05, 20.0);
            rec2.silence(2000);
            rec2.add(iqtest::genCw(rng2, 3000), 0.15, 20.0);
            rec2.silence(2000);
            const std::vector<iqtest::cd>& b2 = rec2.buffer();
            iqclassify::ClassifierConfig cOn;   // Default: channelize = true
            iqclassify::ClassifierConfig cOff;
            cOff.channelize = false;
            auto lab = [&](const iqclassify::ClassifierConfig& c) {
                iqclassify::ModulationClassifier clf(c);
                const auto res = clf.classify(b2.data(), b2.size());
                std::vector<std::string> l;
                for (const auto& r : res)
                    l.push_back(iqclassify::toString(r.modulation));
                return joinLabels(l);
            };
            const std::string a = lab(cOn);
            const std::string b = lab(cOff);
            const bool ok = a == b && a == "FSK2/QPSK/Unmodulated";
            if (ok) {
                ++g_ok;
                std::printf("[OK ] %-24s channelized==legacy: %s\n",
                            "chan-eq-legacy", a.c_str());
            } else {
                std::printf("[FAIL] %-24s chan=%s legacy=%s\n",
                            "chan-eq-legacy", a.c_str(), b.c_str());
            }
        }
        // Resampling in classify (channelDecimate=2): gleiche Labels.
        ++g_total;
        {
            iqclassify::ClassifierConfig cfg;
            cfg.channelDecimate = 2;
            iqclassify::ModulationClassifier clf(cfg);
            const auto res = clf.classify(buf.data(), buf.size());
            std::vector<std::string> l;
            for (const auto& r : res)
                l.push_back(iqclassify::toString(r.modulation));
            const std::string got = joinLabels(l);
            const bool ok = got == "QPSK/Unmodulated" && res.size() == 2;
            if (ok) {
                ++g_ok;
                std::printf("[OK ] %-24s channelDecimate=2: %s\n",
                            "chan-decimate-classify", got.c_str());
            } else {
                std::printf("[FAIL] %-24s got=%s n=%zu\n",
                            "chan-decimate-classify", got.c_str(), res.size());
            }
        }
        // Breiter Chirp (Hub 0.40, wuerde gleichverteilte Kanaele zerteilen):
        // adaptiver Kanal deckt den ganzen Hub ab -> genau ein Puls.
        ++g_total;
        {
            std::mt19937 rng3(31);
            iqtest::Recorder rec3;
            rec3.silence(2000);
            rec3.add(iqtest::genChirp(4000, 0.40, true), 0.0, 20.0);
            rec3.silence(2000);
            iqclassify::ClassifierConfig cfg;
            iqclassify::ModulationClassifier clf(cfg);
            const auto res =
                clf.classify(rec3.buffer().data(), rec3.buffer().size());
            const bool ok = res.size() == 1 &&
                            res[0].modulation ==
                                iqclassify::ModulationType::Chirp;
            if (ok) {
                ++g_ok;
                std::printf("[OK ] %-24s breiter Chirp = 1 Puls (%s)\n",
                            "chan-chirp-wide", iqclassify::toString(res[0].modulation));
            } else {
                std::printf("[FAIL] %-24s n=%zu\n", "chan-chirp-wide",
                            res.size());
                for (const auto& r : res) printPulse(r);
            }
        }
    }
    // 17. Skalierung: 2 Mio Samples, 3 Pulse sequenziell -> ein classify-Lauf.
    {
        ++g_total;
        std::mt19937 rng(20240607);
        iqtest::Recorder rec;
        rec.silence(2000000);
        rec.add(iqtest::genPsk(rng, 4, 400, 8), 0.15, 20.0, 200000);
        rec.add(iqtest::genCw(rng, 3000), -0.20, 20.0, 900000);
        rec.add(iqtest::genFsk2(rng, 400, 8), 0.05, 20.0, 1600000);
        const std::vector<iqtest::cd>& buf = rec.buffer();
        iqclassify::ClassifierConfig cfg;
        iqclassify::ModulationClassifier clf(cfg);
        const auto t0 = std::chrono::steady_clock::now();
        const auto res = clf.classify(buf.data(), buf.size());
        const double ms =
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t0)
                .count();
        std::vector<std::string> l;
        for (const auto& r : res)
            l.push_back(iqclassify::toString(r.modulation));
        const std::string got = joinLabels(l);
        // Laufzeit nur gegen Pathologien absichern (CI-Rechner variieren).
        const bool ok = got == "FSK2/QPSK/Unmodulated" && res.size() == 3 &&
                        ms < 120000.0;
        if (ok) {
            ++g_ok;
            std::printf("[OK ] %-24s n=%zu -> %s in %.0f ms (%.2f MS/s)\n",
                        "scale-2M", buf.size(), got.c_str(), ms,
                        buf.size() / 1e6 / (ms / 1e3));
        } else {
            std::printf("[FAIL] %-24s got=%s n=%zu ms=%.0f\n", "scale-2M",
                        got.c_str(), res.size(), ms);
            for (const auto& r : res) printPulse(r);
        }
    }
    // 18. Haertung II: Extremfaelle mit definiertem Ergebnis.
    {
        auto check = [](const char* name, const std::vector<iqtest::cd>& buf,
                        const std::string& expect, std::size_t maxFine0 = 0) {
            ++g_total;
            iqclassify::ClassifierConfig cfg;
            iqclassify::ModulationClassifier clf(cfg);
            const auto res = clf.classify(buf.data(), buf.size());
            std::vector<std::string> l;
            for (const auto& r : res)
                l.push_back(iqclassify::toString(r.modulation));
            const std::string got = joinLabels(l);
            bool ok = got == expect;
            if (ok && maxFine0 > 0 && !res.empty())
                ok = res[0].startSampleFine < maxFine0;
            if (ok) {
                ++g_ok;
                std::printf("[OK ] %-24s n=%zu -> %s\n", name, buf.size(),
                            got.empty() ? "(leer)" : got.c_str());
            } else {
                std::printf("[FAIL] %-24s expect=%s got=%s n=%zu\\n", name,
                            expect.c_str(), got.c_str(), res.size());
            }
        };
        {   // Extremamplitude 1e150: skaliert korrekt (QPSK), kein inf/NaN.
            std::mt19937 rng(1);
            iqtest::Recorder rec;
            rec.silence(2000);
            auto p = iqtest::genPsk(rng, 4, 400, 8);
            for (auto& v : p) v *= 1e150;
            rec.add(p, 0.12, 20.0);
            rec.silence(2000);
            check("hard-amp-1e150", rec.buffer(), "QPSK");
        }
        {   // Extrem kleine Amplitude 1e-300: im Rauschen unsichtbar -> leer.
            std::mt19937 rng(2);
            iqtest::Recorder rec;
            rec.silence(2000);
            auto p = iqtest::genPsk(rng, 4, 400, 8);
            for (auto& v : p) v *= 1e-300;
            rec.add(p, 0.12, 20.0);
            rec.silence(2000);
            check("hard-amp-1e-300", rec.buffer(), "");
        }
        {   // Nyquist-Muster (+-1 alternierend): definiert, keine Pulse.
            std::vector<iqtest::cd> b(8000);
            for (std::size_t i = 0; i < b.size(); ++i)
                b[i] = iqtest::cd(i % 2 ? 1.0 : -1.0, 0.0);
            check("hard-nyquist-alt", b, "");
        }
        {   // Puls beginnt bei Sample 0.
            std::mt19937 rng(3);
            iqtest::Recorder rec;
            rec.silence(500);
            rec.add(iqtest::genPsk(rng, 4, 400, 8), 0.12, 20.0, 0);
            rec.silence(2000);
            check("hard-pulse-at-0", rec.buffer(), "QPSK", 200);
        }
        {   // Puls endet exakt am Pufferende.
            std::mt19937 rng(4);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genCw(rng, 3000), -0.20, 20.0, rec.buffer().size());
            check("hard-pulse-at-end", rec.buffer(), "Unmodulated");
        }
        {   // Sehr hohes SNR (60 dB).
            std::mt19937 rng(5);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genPsk(rng, 4, 400, 8), 0.12, 60.0);
            rec.silence(2000);
            check("hard-snr-60dB", rec.buffer(), "QPSK");
        }
        {   // Mini-Puffer n = 1..70 (Rauschen): definiert leer, kein Absturz.
            ++g_total;
            bool ok = true;
            for (std::size_t n : {1u, 2u, 3u, 15u, 16u, 17u, 63u, 64u, 65u, 70u}) {
                std::mt19937 rng(6);
                std::normal_distribution<double> g(0.0, 0.008);
                std::vector<iqtest::cd> b(n);
                for (auto& v : b) v = iqtest::cd(g(rng), g(rng));
                iqclassify::ClassifierConfig cfg;
                iqclassify::ModulationClassifier clf(cfg);
                if (!clf.classify(b.data(), b.size()).empty()) ok = false;
            }
            if (ok) {
                ++g_ok;
                std::printf("[OK ] %-24s n=1..70 Rauschen -> leer\\n",
                            "hard-mini-sweep");
            } else {
                std::printf("[FAIL] %-24s Mini-Puffer lieferte Pulse\\n",
                            "hard-mini-sweep");
            }
        }
        {   // Puffer vollstaendig vom Puls ausgefuellt (keine stille Referenz):
            // dokumentierte Segmenter-Grenze -> definiert 0 Pulse.
            std::mt19937 rng(7);
            auto p = iqtest::genPsk(rng, 4, 400, 8);
            std::vector<iqtest::cd> b(p.size());
            std::mt19937 r2(8);
            std::normal_distribution<double> g(0.0, 0.01);
            for (std::size_t i = 0; i < b.size(); ++i)
                b[i] = p[i] + iqtest::cd(g(r2), g(r2));
            check("hard-full-buffer", b, "");
        }
    }
    // 19. Parallelverarbeitung: identische Ergebnisse unabhaengig von der
    //     Threadanzahl (kanalisierte Pipeline und Channelizer).
    {
        // Puffer >= 32k mit mehreren Pulsen (aktiviert die Parallelpfade).
        std::mt19937 rng(99);
        iqtest::Recorder rec;
        rec.silence(60000);
        rec.add(iqtest::genPsk(rng, 4, 400, 8), 0.15, 20.0, 5000);
        rec.add(iqtest::genCw(rng, 3000), -0.20, 20.0, 20000);
        rec.add(iqtest::genFsk2(rng, 400, 8), 0.05, 20.0, 35000);
        rec.add(iqtest::genChirp(4000, 0.25, true), -0.05, 20.0, 50000);
        rec.silence(20000);
        const std::vector<iqtest::cd>& buf = rec.buffer();
        ++g_total;
        {
            iqclassify::ClassifierConfig c1;
            c1.threads = 1;
            iqclassify::ClassifierConfig c4;
            c4.threads = 4;
            iqclassify::ModulationClassifier f1(c1);
            iqclassify::ModulationClassifier f4(c4);
            const auto r1 = f1.classify(buf.data(), buf.size());
            const auto r4 = f4.classify(buf.data(), buf.size());
            bool ok = r1.size() == r4.size() && r1.size() == 4;
            for (std::size_t i = 0; ok && i < r1.size(); ++i) {
                ok = r1[i].modulation == r4[i].modulation &&
                     r1[i].startSampleFine == r4[i].startSampleFine &&
                     r1[i].numSamplesFine == r4[i].numSamplesFine &&
                     r1[i].centerFreqRefined == r4[i].centerFreqRefined;
            }
            if (ok) {
                ++g_ok;
                std::printf("[OK ] %-24s threads=1 == threads=4 (%zu Pulse)\n",
                            "par-identical", r1.size());
            } else {
                std::printf("[FAIL] %-24s n1=%zu n4=%zu\n", "par-identical",
                            r1.size(), r4.size());
            }
        }
        ++g_total;
        {
            iqclassify::ClassifierConfig c1;
            c1.threads = 1;
            iqclassify::ClassifierConfig c4;
            c4.threads = 4;
            iqclassify::Channelizer z1(c1);
            iqclassify::Channelizer z4(c4);
            const auto a = z1.channelize(buf.data(), buf.size(), 1);
            const auto b = z4.channelize(buf.data(), buf.size(), 1);
            bool ok = a.size() == b.size() && a.size() >= 1;
            for (std::size_t i = 0; ok && i < a.size(); ++i)
                ok = a[i].centerNorm == b[i].centerNorm &&
                     a[i].y.size() == b[i].y.size();
            if (ok) {
                ++g_ok;
                std::printf("[OK ] %-24s channelize seriell == parallel (%zu Kanaele)\n",
                            "par-channelize", a.size());
            } else {
                std::printf("[FAIL] %-24s n1=%zu n4=%zu\n", "par-channelize",
                            a.size(), b.size());
            }
        }
    }
    // 20. Sehr kurze Pulse (~180-Sample-Puffer, 13 BPSK-Bit bei sps=8):
    //     automatisches Upsampling des Basisband-Slices + shortPulseProfile.
    //     Default-Profil bleibt auf laengere Pulse kalibriert (0 Pulse).
    {
        // kind: 0 = PSK(m), 1 = FSK2, 2 = CW
        auto checkShort = [&](const char* name, int seed, int kind, int m,
                              const char* expectMod, double fcTrue) {  // expectMod may be nullptr
            ++g_total;
            std::mt19937 rng(seed);
            iqtest::Recorder rec;
            rec.silence(20);
            if (kind == 2) {
                rec.add(iqtest::genCw(rng, 120), fcTrue, 20.0);
            } else if (kind == 1) {
                rec.add(iqtest::genFsk2(rng, 13, 8), fcTrue, 20.0);
            } else {
                rec.add(iqtest::genPsk(rng, m, 13, 8), fcTrue, 20.0);
            }
            rec.silence(24);
            const auto& buf = rec.buffer();
            iqclassify::ClassifierConfig cfg =
                iqclassify::ClassifierConfig::shortPulseProfile();
            iqclassify::ModulationClassifier clf(cfg);
            const auto res = clf.classify(buf.data(), buf.size());
            // Default-Profil: minPulseLen=400 -> kurzer Puls wird nicht erwartet
            iqclassify::ModulationClassifier clfDef{iqclassify::ClassifierConfig()};
            const auto resDef =
                clfDef.classify(buf.data(), buf.size());
            const double fcExp = fcTrue;
            bool ok = res.size() == 1 && resDef.empty();
            if (ok) {
                const auto& r = res[0];
                // expectMod == nullptr: nur Detektion + Frequenzlage pruefen
                // (bei 13 Symbolen ist die Feinklasse PSK/QAM-Nachbarn).
                ok = (expectMod == nullptr ||
                      std::string(iqclassify::toString(r.modulation)) ==
                          expectMod) &&
                     cycDist(r.centerFreqRefined, fcExp) < 0.02 &&
                     r.startSampleFine < 156 &&
                     r.startSampleFine + r.numSamplesFine > 20;
            }
            if (ok) {
                ++g_ok;
                std::printf("[OK ] %-24s n=%zu -> %s fc=%.4f @%zu+%zu\n", name,
                            buf.size(),
                            iqclassify::toString(res[0].modulation),
                            res[0].centerFreqRefined, res[0].startSampleFine,
                            res[0].numSamplesFine);
            } else {
                std::printf("[FAIL] %-24s n=%zu nRes=%zu nDef=%zu", name,
                            buf.size(), res.size(), resDef.size());
                for (const auto& r : res)
                    std::printf(" %s(fc=%.4f @%zu+%zu)",
                                iqclassify::toString(r.modulation),
                                r.centerFreqRefined, r.startSampleFine,
                                r.numSamplesFine);
                std::printf("\n");
            }
        };
        // BPSK 13 Bit: Detektion + Frequenzlage; die Feinklassifikation bei
        // nur 13 Symbolen ist grenzwertig (PSK/QAM-Nachbarn), daher wird nur
        // die Detektion und Lokalisierung hart geprueft.
        ++g_total;
        {
            bool allDetected = true;
            double fcErr = 0.0;
            for (int seed : {123, 7, 42}) {
                std::mt19937 rng(seed);
                iqtest::Recorder rec;
                rec.silence(20);
                rec.add(iqtest::genPsk(rng, 2, 13, 8), 0.13, 20.0);
                rec.silence(24);
                iqclassify::ClassifierConfig cfg =
                    iqclassify::ClassifierConfig::shortPulseProfile();
                iqclassify::ModulationClassifier clf(cfg);
                const auto res = clf.classify(rec.buffer().data(),
                                             rec.buffer().size());
                if (res.size() != 1) {
                    allDetected = false;
                    break;
                }
                fcErr = std::max(fcErr, cycDist(res[0].centerFreqRefined, 0.13));
                if (res[0].startSampleFine >= 156 ||
                    res[0].startSampleFine + res[0].numSamplesFine <= 20)
                    allDetected = false;
            }
            if (allDetected && fcErr < 0.02) {
                ++g_ok;
                std::printf("[OK ] %-24s 3 Seeds -> 1 Puls, fc=0.13+-%.4f\n",
                            "short-bpsk13", fcErr);
            } else {
                std::printf("[FAIL] %-24s detected=%d fcErr=%.4f\n",
                            "short-bpsk13", (int)allDetected, fcErr);
            }
        }
        checkShort("short-qpsk13", 5, 0, 4, nullptr, -0.10);
        // FSK2-Generator erzeugt zwei feste Frequenzen um fc; die
        // Mittenfrequenz liegt bei fc (Hub symmetrisch um fc).
        checkShort("short-fsk2-13", 6, 1, 0, "FSK2", 0.20);
        checkShort("short-cw120", 7, 2, 0, "Unmodulated", 0.30);
        // Upsampling explizit aus: Slice bleibt kurz, Detektion klappt ueber
        // das Profil, aber die Symbolanalyse hat weniger Samples zur Verfuegung.
        ++g_total;
        {
            std::mt19937 rng(123);
            iqtest::Recorder rec;
            rec.silence(20);
            rec.add(iqtest::genPsk(rng, 2, 13, 8), 0.13, 20.0);
            rec.silence(24);
            iqclassify::ClassifierConfig cfg =
                iqclassify::ClassifierConfig::shortPulseProfile();
            cfg.upsample = 1;
            iqclassify::ModulationClassifier clf(cfg);
            const auto res = clf.classify(rec.buffer().data(),
                                         rec.buffer().size());
            const bool ok = res.size() == 1 &&
                            cycDist(res[0].centerFreqRefined, 0.13) < 0.02;
            if (ok) {
                ++g_ok;
                std::printf("[OK ] %-24s upsample=1 -> %s fc=%.4f @%zu+%zu\n",
                            "short-noupsample",
                            iqclassify::toString(res[0].modulation),
                            res[0].centerFreqRefined, res[0].startSampleFine,
                            res[0].numSamplesFine);
            } else {
                std::printf("[FAIL] %-24s nRes=%zu\n", "short-noupsample",
                            res.size());
            }
        }
    }
    // 21. Frank-Codes (polyphase Radar-Code): Treppe der Momentanfrequenz.
    {
        // M=8, 20 dB und 15 dB, je 3 Seeds: genau 1 Puls, FrankCode, fc-Lage.
        for (double snr : {20.0, 15.0}) {
            ++g_total;
            int okc = 0;
            double fcErr = 0.0;
            for (int seed = 1; seed <= 3; ++seed) {
                std::mt19937 rng(seed * 10 + 8);
                iqtest::Recorder rec;
                rec.silence(2000);
                rec.add(iqtest::genFrank(rng, 8, 8), 0.15, snr);
                rec.silence(2000);
                iqclassify::ClassifierConfig cfg;
                iqclassify::ModulationClassifier clf(cfg);
                const auto res =
                    clf.classify(rec.buffer().data(), rec.buffer().size());
                if (res.size() == 1 &&
                    res[0].modulation == iqclassify::ModulationType::FrankCode) {
                    ++okc;
                    fcErr = std::max(fcErr, cycDist(res[0].centerFreqRefined, 0.15));
                }
            }
            if (okc == 3 && fcErr < 0.03) {
                ++g_ok;
                std::printf("[OK ] %-24s M=8 snr=%.0f: 3/3 FrankCode fcErr=%.4f\n",
                            "frank-M8", snr, fcErr);
            } else {
                std::printf("[FAIL] %-24s M=8 snr=%.0f: %d/3 fcErr=%.4f\n",
                            "frank-M8", snr, okc, fcErr);
            }
        }
        // M=16, 20 dB: mindestens ein FrankCode-Puls (das Kamm-Spektrum kann
        // bei hoeherer Ordnung in Subbaender zerfallen).
        {
            ++g_total;
            std::mt19937 rng(77);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genFrank(rng, 16, 8), 0.15, 20.0);
            rec.silence(2000);
            iqclassify::ClassifierConfig cfg;
            iqclassify::ModulationClassifier clf(cfg);
            const auto res =
                clf.classify(rec.buffer().data(), rec.buffer().size());
            bool any = false;
            for (const auto& r : res)
                if (r.modulation == iqclassify::ModulationType::FrankCode)
                    any = true;
            if (any && !res.empty()) {
                ++g_ok;
                std::printf("[OK ] %-24s M=16 snr=20: %zu Pulse, mind. 1 FrankCode\n",
                            "frank-M16", res.size());
            } else {
                std::printf("[FAIL] %-24s M=16 snr=20: n=%zu\n", "frank-M16",
                            res.size());
                for (const auto& r : res) printPulse(r);
            }
        }
        // Mehrpuls: Frank + CW zeitlich getrennt in getrennten Baendern.
        {
            std::mt19937 rng(42);
            iqtest::Recorder rec;
            rec.silence(3000);
            rec.add(iqtest::genFrank(rng, 8, 8), 0.15, 20.0);
            rec.silence(800);
            rec.add(iqtest::genCw(rng, 4000), -0.35, 20.0);
            rec.silence(3000);
            scenario("frank-multi-cw", rec.buffer(), "FrankCode/Unmodulated");
        }
        // Analyse-API: Ordnung, Code-Laenge, Chip-Rate, R^2, Code-Phasen.
        {
            ++g_total;
            std::mt19937 rng(11);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genFrank(rng, 8, 8), 0.15, 20.0);
            rec.silence(2000);
            iqclassify::ClassifierConfig cfg;
            iqclassify::ModulationClassifier clf(cfg);
            const auto res =
                clf.classify(rec.buffer().data(), rec.buffer().size());
            bool ok = res.size() == 1 &&
                      res[0].modulation == iqclassify::ModulationType::FrankCode;
            iqclassify::FrankParams fp;
            if (ok) {
                fp = clf.analyzeFrank(rec.buffer().data(), rec.buffer().size(),
                                      res[0]);
                ok = fp.valid && fp.order >= 3 && fp.order <= 16 &&
                     fp.codeLength ==
                         static_cast<std::size_t>(fp.order) * fp.order &&
                     fp.r2 > 0.9 && fp.chipRateNorm > 0.01 &&
                     fp.chipRateNorm < 0.5 && fp.codePhases.size() == fp.codeLength;
            }
            if (ok) {
                ++g_ok;
                std::printf("[OK ] %-24s M=8: order=%d N=%zu chipRate=%.4f r2=%.3f\n",
                            "frank-analyze", fp.order, fp.codeLength,
                            fp.chipRateNorm, fp.r2);
            } else {
                std::printf("[FAIL] %-24s nRes=%zu valid=%d order=%d r2=%.3f\n",
                            "frank-analyze", res.size(), (int)fp.valid, fp.order,
                            fp.r2);
            }
        }
        // Negativ: Chirp und FSK2 duerfen nicht als FrankCode durchgehen
        // (wird zusaetzlich durch die gesamte uebrige Suite abgesichert).
        {
            std::mt19937 rng(5);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genChirp(8192, 0.25, true), 0.10, 20.0);
            rec.silence(2000);
            scenario("frank-neg-chirp", rec.buffer(), "Chirp");
        }
        {
            std::mt19937 rng(6);
            iqtest::Recorder rec;
            rec.silence(2000);
            rec.add(iqtest::genFsk2(rng, 1024, 8), 0.10, 20.0);
            rec.silence(2000);
            scenario("frank-neg-fsk2", rec.buffer(), "FSK2");
        }
    }
    std::printf("\nscore: %d/%d\n", g_ok, g_total);
    std::printf("%s\n", g_ok == g_total ? "ALL TESTS PASSED"
                                        : "TEST FAILURES PRESENT");
    return g_ok == g_total ? 0 : 1;
}
