// IQClassify - Testsuite.
//
// Szenarien:
//   1. Basisfall: 9 Modulationen, SNR 20 dB, sps 8 (Referenzkalibrierung)
//   2. Sequenz: drei Pulse hintereinander (QPSK, FSK2, CW)
//   3. Ueberlagerung: QAM16 und BPSK gleichzeitig in verschiedenen
//      Frequenzbaendern
//   4. SNR 10 dB: 8 Modulationen; QAM256 ist unterhalb der
//      Aufloesungsgrenze dokumentiert als QAM64-Fallback
//   5. Oversampling: BPSK/QPSK/FSK2 bei sps 4 und 16
//
// Ein Szenario besteht genau dann den Test, wenn die sortierte Labelmenge
// der erkannten Pulse der erwarteten entspricht.

#include <algorithm>
#include <cstdio>
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
        std::printf("[OK ] %-18s expect=%-28s n_pulses=%zu\n", name,
                     expect.c_str(), res.size());
    } else {
        std::printf("[FAIL] %-18s expect=%-28s got=%-28s n_pulses=%zu\n",
                    name, expect.c_str(), got.c_str(), res.size());
        for (const auto& r : res) printPulse(r);
    }
    return ok;
}

// Basisfall mit wahlweisem SNR. Liefert die Anzahl bestandener Faelle.
void runFamily(const char* family, std::mt19937& rng, double snrDb,
               bool includeQam256, std::size_t qamNsym, int sps) {
    // Signale in fester Reihenfolge ziehen (ein gemeinsamer Zufallsstrom
    // je Familie, wie in der Referenzkalibrierung)
    std::vector<iqtest::cd> bpsk = iqtest::genPsk(rng, 2, 400, sps);
    std::vector<iqtest::cd> qpsk = iqtest::genPsk(rng, 4, 400, sps);
    std::vector<iqtest::cd> fsk = iqtest::genFsk2(rng, 400, sps);
    std::vector<iqtest::cd> chirp = iqtest::genChirp(4000, 0.20);
    std::vector<iqtest::cd> cw = iqtest::genCw(rng, 3000);
    std::vector<iqtest::cd> ofdm = iqtest::genOfdm(rng, 256, 64, 20, 0.3);
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
        {"FSK2", &fsk, "FSK2", 0.05},
        {"CHIRP", &chirp, "Chirp", -0.10},
        {"CW", &cw, "Unmodulated", 0.15},
        {"OFDM", &ofdm, "OFDM", -0.15},
        {"QAM16", &qam16, "QAM16", 0.07},
    };
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
        const char* expect256 =
            snrDb >= 20.0 ? "QAM256" : "QAM64";  // dokumentierte Grenze
        scenario((std::string(family) + "-QAM256").c_str(),
                 rec256.buffer(), expect256);
    }
}

}  // namespace

int main() {
    // 1. Basiskalibrierung
    {
        std::mt19937 rng(42);
        runFamily("base", rng, 20.0, true, 4000, 8);
    }
    // 2. Sequenz dreier zeitlich getrennter Pulse
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
    // 3. Gleichzeitige Pulse in verschiedenen Frequenzbaendern
    {
        std::mt19937 rng(11);
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(iqtest::genQam(rng, 16, 4000, 8), 0.07, 20.0, 2000);
        rec.add(iqtest::genPsk(rng, 2, 4000, 8), -0.07, 20.0, 2000);
        rec.silence(6000);
        scenario("overlap-2pulses", rec.buffer(), "BPSK/QAM16");
    }
    // 4. SNR 10 dB (auesserhalb: QAM256-Fallback, dokumentierte Grenze)
    {
        std::mt19937 rng(123);
        runFamily("snr10", rng, 10.0, true, 4000, 8);
    }
    // 5. Oversampling-Varianten
    {
        std::mt19937 rng(2024);
        for (int sps : {4, 16}) {
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
            rec.add(iqtest::genFsk2(rng, 400, sps), 0.05, 20.0);
            rec.silence(2000);
            scenario((std::string("sps") + std::to_string(sps) + "-FSK2")
                         .c_str(),
                     rec.buffer(), "FSK2");
        }
    }
    std::printf("\nscore: %d/%d\n", g_ok, g_total);
    std::printf("%s\n", g_ok == g_total ? "ALL TESTS PASSED"
                                        : "TEST FAILURES PRESENT");
    return g_ok == g_total ? 0 : 1;
}
