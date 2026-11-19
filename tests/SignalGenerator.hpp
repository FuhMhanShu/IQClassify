// IQClassify - deterministische Testsignale fuer die Klassifikator-Testsuite.
//
// Die Generatoren entsprechen der Parameterwelt der Validierungsreferenz
// (RRC-Pulsformung alpha 0.25, FSK-Tonabstand 1/8, OFDM nfft 256, cp 64,
// Chirp-Hub 0.20). Alle Zufallsanteile laufen ueber std::mt19937 mit
// festen Seeds; die Suite ist dadurch voll reproduzierbar.

#pragma once

#include <complex>
#include <cstddef>
#include <random>
#include <vector>

namespace iqtest {

using cd = std::complex<double>;

// PSK mit M Phasenlagen (M = 2 -> BPSK, M = 4 -> QPSK).
std::vector<cd> genPsk(std::mt19937& rng, int m, std::size_t nsym, int sps,
                       double alpha = 0.25);

// Quadraturamplitude (16/64/256), gleichverteilig ueber das Levelgitter.
std::vector<cd> genQam(std::mt19937& rng, int m, std::size_t nsym, int sps,
                       double alpha = 0.25);

// 2-FSK mit Tonabstand sep (fraktional zur Abtastrate).
std::vector<cd> genFsk2(std::mt19937& rng, std::size_t nsym, int sps,
                        double sep = 1.0 / 8.0);

// Linearer Chirp ueber die Bandbreite bw (fraktional).
std::vector<cd> genChirp(std::size_t nsamp, double bw, bool up = true);

// Unmodulierter Traeger (CW) mit Zufallsphase.
std::vector<cd> genCw(std::mt19937& rng, std::size_t nsamp);

// Amplitudenschwankende Variante eines Pulses (Huellkurve multipliziert,
// anschliessend auf mittlere Leistung 1 normiert):
//   mode 0: Zufallsamplitude pro Sample in [1-depth, 1] (depth z. B. 0.7)
//   mode 1: langsamer Sinus-Fade zwischen 1-depth und 1
//   mode 2: Gauss'scher Einbruch in der Pulsmitte um den Faktor (1-depth)
std::vector<cd> applyEnv(const std::vector<cd>& sig, int mode, double depth,
                         std::mt19937& rng);

// OFDM: nfft Nutztraeger, cp Samples zyklisches Praefix, nsym Symbole,
// load Anteil belegter Traeger (symmetrisch um DC, QPSK-Beladung).
std::vector<cd> genOfdm(std::mt19937& rng, int nfft, int cp, int nsym,
                        double load = 1.0);

// Frequenzverschiebung um fc (fraktional, Vorzeichen wie die Referenz).
std::vector<cd> mix(const std::vector<cd>& sig, double fc);

// Additives weisses Gaussrauschen auf ein Signal mit 1 mittlerer Leistung.
std::vector<cd> noised(std::mt19937& rng, const std::vector<cd>& sig,
                       double snrDb);

// IQ-Strom-Recorder: platziert Pulse (sequenziell oder ueberlagert) und
// fuelle Zwischenraeume mit schwachem Schweberausch.
class Recorder {
public:
    static constexpr std::size_t kAtEnd =
        static_cast<std::size_t>(-1);

    explicit Recorder(double silenceSigma = 0.008)
        : sigma_(silenceSigma) {}

    void silence(std::size_t nsamp);

    // Fuegt einen Puls bei fc und SNR snrDb hinzu; at gibt die absolute
    // Startposition an (kAtEnd: anhaengen am Ende).
    void add(const std::vector<cd>& sig, double fc, double snrDb,
             std::size_t at = kAtEnd);

    const std::vector<cd>& buffer() const { return buf_; }

private:
    std::vector<cd> buf_;
    double sigma_;
    int silenceCount_ = 0;
    int pulseCount_ = 0;
};

}  // namespace iqtest
