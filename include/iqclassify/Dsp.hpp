// IQClassify - interne DSP-Hilfsfunktionen (FFT, Statistik).
// Dieser Header ist ein Implementierungsdetail; Anwender binden
// iqclassify.hpp ein.

#pragma once

#include <complex>
#include <cstddef>
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

}  // namespace detail
}  // namespace iqclassify
