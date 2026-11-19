// IQClassify - interne DSP-Hilfsfunktionen (FFT, Statistik).
// Dieser Header ist ein Implementierungsdetail; Anwender binden
// iqclassify.hpp ein.

#pragma once

#include <atomic>
#include <complex>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

namespace iqclassify {
namespace detail {

using cd = std::complex<double>;

// Kleinstgroesse >= x, die eine Zweierpotenz ist.
std::size_t pow2ceil(std::size_t x);

// Fouriertransformierte in place, Laenge muss Zweierpotenz sein.
// forward: negativer Exponent; invert: positiver Exponent, 1/n-Skalierung.
// Ruft intern pocketfft auf (unterstuetzt beliebige Laengen, daher ohne
// Einschraenkung sicher); die Zweierpotenz-Anforderung haelt den
// Aufrufer einfach und die Plancache-Nutzung vorhersagbar).
void fftInPlace(std::vector<cd>& a, bool invert);

// Fouriertransformierte einer kopierten Quelle (Nullen angehaengt bis n).
std::vector<cd> fftPadded(const cd* x, std::size_t len, std::size_t n,
                          bool invert);

// 10*log10 mit Untergrenze gegen -300 dB.
double db10(double v);

// Median einer Liste (Mittel der mittleren beiden Elemente bei gerader Zahl).
double median(std::vector<double>& v);
double medianOf(const std::vector<double>& v);

// 2n auf Zweierpotenz gerundet (Pufferlaenge fuer zykische Korrelationen).
std::size_t protoPad(std::size_t n);

// Kontinuierliche Phase (Unwrapping +-pi) einer komplexen Folge; das erste
// Element ist phase(y[0]).
std::vector<double> unwrapPhase(const std::vector<cd>& y);

// Root-raised-cosine-Impulsantwort (Sender-/Empfangsfilter), summen-
// normiert auf 1. sps Samples pro Symbol, Rolloff alpha, Laenge
// span*sps+1 Taps. Analytische Stammfunktion an allen Singularitaeten.
std::vector<double> rrcTaps(long sps, double alpha, long span);

// Effektive Thread-Anzahl: requested <= 0 -> automatisch
// (min(hardware_concurrency, 8)); sonst min(requested, hardware_concurrency).
unsigned effectiveThreads(int requested);

// Fuehrt f(i) fuer i in [0, n) verteilt auf bis zu `threads` Threads aus
// (atomarer Work-Counter; threads <= 1 oder n == 0 -> seriell im aufrufenden
// Thread). pocketfft-Plancache und fftInPlace sind thread-sicher.
template <typename F>
void parallelFor(std::size_t n, unsigned threads, F&& f) {
    if (n == 0) return;
    if (threads <= 1) {
        for (std::size_t i = 0; i < n; ++i) f(i);
        return;
    }
    if (threads > n) threads = static_cast<unsigned>(n);
    std::atomic<std::size_t> next{0};
    const auto worker = [&]() {
        for (;;) {
            const std::size_t i = next.fetch_add(1);
            if (i >= n) return;
            f(i);
        }
    };
    std::vector<std::thread> pool;
    pool.reserve(threads - 1);
    for (unsigned t = 1; t < threads; ++t) pool.emplace_back(worker);
    worker();
    for (auto& th : pool) th.join();
}

}  // namespace detail
}  // namespace iqclassify
