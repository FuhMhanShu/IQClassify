// IQClassify - Performance-Benchmark.
//
// Misst die Laufzeit von classify() (kanalisierter Default-Pfad und
// Legacy-Pfad) sowie der Kanalisierung allein (Channelizer::channelize) auf
// deterministischen IQ-Puffern mit 10k .. 2M Samples und einer festen bzw.
// skalierenden Anzahl von Pulsen. Ausgabe: Tabelle auf stdout.
//
//   ./build/iqclassify_benchmark

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

#ifdef __unix__
#include <sys/resource.h>
#endif

#include "SignalGenerator.hpp"
#include "iqclassify/iqclassify.hpp"

namespace {

using iqtest::cd;
using clk = std::chrono::steady_clock;

double msSince(clk::time_point t0) {
    return std::chrono::duration<double, std::milli>(clk::now() - t0).count();
}

double peakRssMb() {
#ifdef __unix__
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) == 0)
        return static_cast<double>(ru.ru_maxrss) / 1024.0;  // KiB -> MiB
#endif
    return 0.0;
}

// Deterministischer Puffer: Rauschen ueberall, Pulse an festen Bruchteilen.
std::vector<cd> makeBuffer(std::size_t n, int npulse, unsigned seed) {
    std::mt19937 rng(seed);
    iqtest::Recorder rec;
    rec.silence(n);
    const int kinds = 5;
    for (int i = 0; i < npulse; ++i) {
        const std::size_t at =
            static_cast<std::size_t>((0.08 + 0.84 * (i + 0.5) / npulse) *
                                     static_cast<double>(n));
        const double fc = -0.30 + 0.6 * ((i * 0.37) - std::floor(i * 0.37));
        switch (i % kinds) {
            case 0:
                rec.add(iqtest::genPsk(rng, 4, 400, 8), fc, 20.0, at);
                break;
            case 1:
                rec.add(iqtest::genCw(rng, 3000), fc, 20.0, at);
                break;
            case 2:
                rec.add(iqtest::genFsk2(rng, 400, 8), fc, 20.0, at);
                break;
            case 3:
                rec.add(iqtest::genChirp(4000, 0.25, (i / 3) % 2 == 0), fc, 20.0,
                        at);
                break;
            default:
                rec.add(iqtest::genQam(rng, 16, 300, 8), fc, 20.0, at);
                break;
        }
    }
    return rec.buffer();
}

std::string labelsOf(const std::vector<iqclassify::PulseResult>& res) {
    std::vector<std::string> l;
    for (const auto& r : res)
        l.push_back(iqclassify::toString(r.modulation));
    std::sort(l.begin(), l.end());
    std::string s;
    for (const auto& x : l) {
        if (!s.empty()) s += "/";
        s += x;
    }
    return s;
}

struct Row {
    std::size_t n;
    int npulse;
    double msChan;
    double msLegacy;
    double msChannelize;
    std::size_t found;
    std::string labels;
};

template <typename F>
double timeMedian(F&& f, int reps) {
    std::vector<double> ts;
    for (int i = 0; i < reps; ++i) {
        const auto t0 = clk::now();
        f();
        ts.push_back(msSince(t0));
    }
    std::sort(ts.begin(), ts.end());
    return ts[ts.size() / 2];
}

}  // namespace

int main() {
    const std::size_t sizes[] = {10000, 50000, 100000, 250000,
                                 500000, 1000000, 2000000};
    std::printf("IQClassify Benchmark (fs = 1 normiert, 20 dB SNR, 3 Pulse bzw. "
                "n/8000 Pulse)\n");
    std::printf("%9s %7s | %12s %12s %12s %12s | %6s %8s | %s\n", "samples",
                "pulse", "auto ms", "1thr ms", "legacy ms", "kanal. ms", "found",
                "MS/s", "labels");
    std::printf("%9s %7s +%12s+%12s+%12s+%12s+%6s+%8s+%s\n", "-------", "-----",
                "-----------", "---------", "---------", "---------", "-----", "-----",
                "------");
    for (std::size_t n : sizes) {
        const int npulse = n < 100000 ? 3
                                     : static_cast<int>(std::max<std::size_t>(
                                           3, n / 8000));
        auto buf = makeBuffer(n, npulse, 20240607u);
        const int reps = n <= 100000 ? 5 : (n <= 500000 ? 3 : 1);
        iqclassify::ClassifierConfig cfg;  // Default: channelize = true, threads = auto
        iqclassify::ClassifierConfig cfg1;
        cfg1.threads = 1;  // seriell
        iqclassify::ClassifierConfig cfgL;
        cfgL.channelize = false;
        cfgL.threads = 1;
        iqclassify::ModulationClassifier clf(cfg);
        iqclassify::ModulationClassifier clf1(cfg1);
        iqclassify::ModulationClassifier clfL(cfgL);
        iqclassify::Channelizer chz(cfg);
        std::size_t found = 0;
        std::string labels;
        const double msC = timeMedian(
            [&] {
                const auto res = clf.classify(buf.data(), buf.size());
                if (found == 0) {
                    found = res.size();
                    labels = labelsOf(res);
                }
            },
            reps);
        const double ms1 = timeMedian(
            [&] { clf1.classify(buf.data(), buf.size()); }, reps);
        const double msL = timeMedian(
            [&] { clfL.classify(buf.data(), buf.size()); }, reps);
        const double msK = timeMedian(
            [&] { chz.channelize(buf.data(), buf.size(), 1); }, reps);
        std::printf("%9zu %7d | %12.2f %12.2f %12.2f %12.2f | %6zu %8.2f | %s\n", n,
                    npulse, msC, ms1, msL, msK, found, n / 1e6 / (msC / 1e3),
                    labels.c_str());
        std::fflush(stdout);
    }
    std::printf("\nPeak RSS: %.1f MiB\n", peakRssMb());
    std::printf("(auto = threads=auto (min(hw,8)), 1thr = threads=1, legacy = "
                "channelize=false + threads=1, kanal. = Channelizer::channelize; "
                "Median ueber Wiederholungen)\n");
    return 0;
}
