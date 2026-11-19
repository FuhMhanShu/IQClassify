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
//
// Ein Szenario besteht genau dann den Test, wenn die sortierte Labelmenge
// der erkannten Pulse der erwarteten entspricht. Fuer
// Robustheitsfaelle gilt: definiertes Ergebnis ohne Absturz,bei
// Bedeutungslosigkeit der Eingabe exakt keine Pulse.

#include <algorithm>
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
                "cv=%.3f papr=%.1f flat=%.2f modes=%d conf=%.2f\n",
                r.startSample, r.numSamples, r.centerFreqNorm,
                r.bandwidthNorm, r.snrDb, r.envelopeCv, r.paprWindowedDb,
                r.spectralFlatnessDb, r.freqModes, r.confidence);
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
        std::mt19937 rng(12);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::genPsk(rng, 4, 400, 8), 0.10, 20.0, 2000);
        rec.add(iqtest::genFsk2(rng, 400, 8), -0.10, 20.0, 2000);
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
            rec.silence(2000);
            rec.add(iqtest::genCw(rng, 3000), fc, 20.0);
            rec.silence(2000);
            scenario(name, rec.buffer(), "Unmodulated");

            std::snprintf(name, sizeof(name), "freq%+.2f-chirp", fc);
            rec.silence(2000);
            rec.add(iqtest::genChirp(4000, 0.25), fc, 20.0);
            rec.silence(2000);
            scenario(name, rec.buffer(), "Chirp");
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
    std::printf("\nscore: %d/%d\n", g_ok, g_total);
    std::printf("%s\n", g_ok == g_total ? "ALL TESTS PASSED"
                                        : "TEST FAILURES PRESENT");
    return g_ok == g_total ? 0 : 1;
}
