// IQClassify - Generator fuer die README-Bilder (docs/images/*.png).
//
// Aufruf aus dem Repository-Root (CMake-Build):
//   ./build/iqclassify_readme_images [Ausgabeverzeichnis]   // Default: docs/images
//
// Erzeugt aus einem deterministischen Demo-Puffer (QPSK + CW + FSK2 bei
// 20 dB SNR) und den echten Klassifikationsergebnissen:
//   pipeline.png     - Blockschaltbild der Pipeline
//   detection.png    - Zeit-Huellkurve + Detektionsspektrum mit den
//                      lokalisierten Pulsen (feine Fenster, Baender, f_ref)
//   spectrogram.png  - Spektrogramm des Demo-Puffers mit ueberlagerten
//                      Puls-Rechtecken (Zeit x Frequenz)
//   codes.png        - f_inst-Treppen/Sweeps von Polyphase-Codes (Frank,
//                      P1, P2, P3, P4) und Barker13 mit echtem Label
//
// Reiner C++17-Code ohne externe Abhaengigkeiten: PNG (Palette, zlib
// Stored-Blocks) und 5x7-Bitmap-Font sind eingebettet.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "bitmap_font.h"
#include "iqclassify/iqclassify.hpp"
#include "SignalGenerator.hpp"

namespace {

// ---------------------------------------------------------------------------
// CRC32 (PNG) / Adler32 (zlib) / PNG-Schreiber (Farbpalette, stored blocks)
// ---------------------------------------------------------------------------
uint32_t crc32(const uint8_t* data, std::size_t len) {
    static uint32_t table[256];
    static bool init = false;
    if (!init) {
        for (uint32_t n = 0; n < 256; ++n) {
            uint32_t c = n;
            for (int k = 0; k < 8; ++k)
                c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
            table[n] = c;
        }
        init = true;
    }
    uint32_t c = 0xFFFFFFFFu;
    for (std::size_t i = 0; i < len; ++i)
        c = table[(c ^ data[i]) & 0xFF] ^ (c >> 8);
    return c ^ 0xFFFFFFFFu;
}

uint32_t adler32(const uint8_t* data, std::size_t len) {
    uint32_t a = 1, b = 0;
    for (std::size_t i = 0; i < len; ++i) {
        a = (a + data[i]) % 65521u;
        b = (b + a) % 65521u;
    }
    return (b << 16) | a;
}

void putU32(std::vector<uint8_t>& out, uint32_t v) {
    out.push_back(static_cast<uint8_t>((v >> 24) & 0xFF));
    out.push_back(static_cast<uint8_t>((v >> 16) & 0xFF));
    out.push_back(static_cast<uint8_t>((v >> 8) & 0xFF));
    out.push_back(static_cast<uint8_t>(v & 0xFF));
}

void pngChunk(std::vector<uint8_t>& out, const char* type,
              const std::vector<uint8_t>& data) {
    putU32(out, static_cast<uint32_t>(data.size()));
    const std::size_t tpos = out.size();
    for (int i = 0; i < 4; ++i) out.push_back(static_cast<uint8_t>(type[i]));
    out.insert(out.end(), data.begin(), data.end());
    putU32(out, crc32(out.data() + tpos, 4 + data.size()));
}

// ---------------------------------------------------------------------------
// Canvas: indiziertes Bild mit Palette + Zeichenprimitive
// ---------------------------------------------------------------------------
struct Canvas {
    int w = 0, h = 0;
    std::vector<uint8_t> px;               // Palette-Indizes
    std::vector<uint8_t> pal;              // RGB-Tripel
    std::map<uint32_t, uint8_t> palMap;

    Canvas(int width, int height, uint8_t bg)
        : w(width), h(height), px(static_cast<std::size_t>(width) * height, bg) {}

    uint8_t color(uint8_t r, uint8_t g, uint8_t b) {
        const uint32_t key = (uint32_t(r) << 16) | (uint32_t(g) << 8) | b;
        auto it = palMap.find(key);
        if (it != palMap.end()) return it->second;
        if (pal.size() / 3 >= 256) return 0;  // Palette voll -> Index 0
        const uint8_t idx = static_cast<uint8_t>(pal.size() / 3);
        pal.push_back(r);
        pal.push_back(g);
        pal.push_back(b);
        palMap[key] = idx;
        return idx;
    }

    void set(int x, int y, uint8_t c) {
        if (x >= 0 && x < w && y >= 0 && y < h)
            px[static_cast<std::size_t>(y) * w + x] = c;
    }
    void fillRect(int x0, int y0, int x1, int y1, uint8_t c) {
        for (int y = y0; y <= y1; ++y)
            for (int x = x0; x <= x1; ++x) set(x, y, c);
    }
    void rect(int x0, int y0, int x1, int y1, uint8_t c) {
        for (int x = x0; x <= x1; ++x) {
            set(x, y0, c);
            set(x, y1, c);
        }
        for (int y = y0; y <= y1; ++y) {
            set(x0, y, c);
            set(x1, y, c);
        }
    }
    void hline(int x0, int x1, int y, uint8_t c) {
        for (int x = x0; x <= x1; ++x) set(x, y, c);
    }
    void vline(int x, int y0, int y1, uint8_t c) {
        for (int y = y0; y <= y1; ++y) set(x, y, c);
    }
    void line(int x0, int y0, int x1, int y1, uint8_t c) {
        const int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
        const int dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
        int err = dx + dy;
        for (;;) {
            set(x0, y0, c);
            if (x0 == x1 && y0 == y1) break;
            const int e2 = 2 * err;
            if (e2 >= dy) {
                err += dy;
                x0 += sx;
            }
            if (e2 <= dx) {
                err += dx;
                y0 += sy;
            }
        }
    }
    int textWidth(const std::string& s, int scale) const {
        return static_cast<int>(s.size()) * 6 * scale;
    }
    void text(int x, int y, const std::string& s, int scale, uint8_t c) {
        int cx = x;
        for (char ch : s) {
            unsigned uc = static_cast<unsigned char>(ch);
            if (uc >= 97 && uc <= 122) uc -= 32;  // Lowercase -> Uppercase
            if (uc < 32 || uc > 126) uc = '?';
            const unsigned char* g = kFont[uc - 32];
            for (int r = 0; r < 7; ++r) {
                for (int col = 0; col < 5; ++col) {
                    if (g[r] & (1 << col)) {
                        for (int dy = 0; dy < scale; ++dy)
                            for (int dx = 0; dx < scale; ++dx)
                                set(cx + col * scale + dx, y + r * scale + dy, c);
                    }
                }
            }
            cx += 6 * scale;
        }
    }
    void textCentered(int cx, int y, const std::string& s, int scale,
                      uint8_t c) {
        text(cx - textWidth(s, scale) / 2, y, s, scale, c);
    }
    void arrow(int x0, int y0, int x1, int y1, uint8_t c) {
        line(x0, y0, x1, y1, c);
        const double ang = std::atan2(static_cast<double>(y1 - y0),
                                      static_cast<double>(x1 - x0));
        const int L = 8;
        line(x1, y1,
             x1 - static_cast<int>(L * std::cos(ang - 0.45)),
             y1 - static_cast<int>(L * std::sin(ang - 0.45)), c);
        line(x1, y1,
             x1 - static_cast<int>(L * std::cos(ang + 0.45)),
             y1 - static_cast<int>(L * std::sin(ang + 0.45)), c);
    }

    bool save(const std::string& path) const {
        std::vector<uint8_t> raw;
        raw.reserve(static_cast<std::size_t>(h) * (1 + w));
        for (int y = 0; y < h; ++y) {
            raw.push_back(0);  // Filter: none
            raw.insert(raw.end(), px.begin() + static_cast<std::size_t>(y) * w,
                       px.begin() + static_cast<std::size_t>(y + 1) * w);
        }
        // zlib-Strom mit gespeicherten (unkomprimierten) Deflate-Bloecken
        std::vector<uint8_t> z;
        z.push_back(0x78);
        z.push_back(0x01);
        std::size_t pos = 0;
        while (pos < raw.size() || (raw.empty() && pos == 0)) {
            const std::size_t chunk =
                std::min<std::size_t>(65535, raw.size() - pos);
            const bool final = (pos + chunk >= raw.size());
            z.push_back(final ? 0x01 : 0x00);
            const uint16_t len = static_cast<uint16_t>(chunk);
            z.push_back(static_cast<uint8_t>(len & 0xFF));
            z.push_back(static_cast<uint8_t>((len >> 8) & 0xFF));
            const uint16_t nlen = static_cast<uint16_t>(~len);
            z.push_back(static_cast<uint8_t>(nlen & 0xFF));
            z.push_back(static_cast<uint8_t>((nlen >> 8) & 0xFF));
            z.insert(z.end(), raw.begin() + pos, raw.begin() + pos + chunk);
            pos += chunk;
            if (raw.empty()) break;
        }
        putU32(z, adler32(raw.data(), raw.size()));

        std::vector<uint8_t> png;
        const uint8_t sig[8] = {0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
        png.insert(png.end(), sig, sig + 8);
        std::vector<uint8_t> ihdr;
        putU32(ihdr, static_cast<uint32_t>(w));
        putU32(ihdr, static_cast<uint32_t>(h));
        ihdr.push_back(8);  // Bitdepth
        ihdr.push_back(3);  // Colortype: Palette
        ihdr.push_back(0);
        ihdr.push_back(0);
        ihdr.push_back(0);
        pngChunk(png, "IHDR", ihdr);
        pngChunk(png, "PLTE", pal);
        pngChunk(png, "IDAT", z);
        pngChunk(png, "IEND", {});
        std::ofstream f(path, std::ios::binary);
        if (!f) return false;
        f.write(reinterpret_cast<const char*>(png.data()),
                static_cast<std::streamsize>(png.size()));
        return f.good();
    }
};

std::string fmt(double v, int prec) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", prec, v);
    return buf;
}

// ---------------------------------------------------------------------------
// Demo-Puffer + Klassifikation
// ---------------------------------------------------------------------------
struct Demo {
    std::vector<iqtest::cd> buf;
    std::vector<iqclassify::PulseResult> pulses;
    // (segment, band) je detektiertem Band, fuer die Spektrum-Ansicht
    struct BandView {
        long lo = 0, hi = 0;
        std::size_t nft = 0;
        double fcRef = 0.0;
    };
    std::vector<BandView> bands;
    std::vector<double> psdDb;  // globales Max-Hold-PSD (wie Banddetektion)
    std::size_t psdNft = 0;
};

Demo makeDemo() {
    Demo d;
    std::mt19937 rng(20240607);
    iqtest::Recorder rec;
    rec.silence(2000);
    rec.add(iqtest::genPsk(rng, 4, 400, 8), 0.15, 20.0);   // QPSK @ +0.15
    rec.silence(2000);
    rec.add(iqtest::genCw(rng, 2500), -0.20, 20.0);        // CW @ -0.20
    rec.silence(2000);
    rec.add(iqtest::genFsk2(rng, 400, 8), 0.05, 20.0);     // FSK2 @ +0.05
    rec.silence(2000);
    d.buf = rec.buffer();

    iqclassify::ClassifierConfig cfg;
    iqclassify::ModulationClassifier clf(cfg);
    d.pulses = clf.classify(d.buf.data(), d.buf.size());

    iqclassify::PulseSegmenter seg(cfg);
    iqclassify::BandExtractor bxs(cfg);
    for (const auto& ab : seg.detect(d.buf.data(), d.buf.size())) {
        auto det = bxs.bands(d.buf.data(), ab.first, ab.second);
        for (const auto& b : det.bands) {
            Demo::BandView bv;
            bv.lo = b.lo;
            bv.hi = b.hi;
            bv.nft = det.nft;
            d.bands.push_back(bv);
        }
    }
    // globales Max-Hold-PSD ueber den ganzen Puffer (Detektions-Charakteristik)
    {
        const std::size_t wlen = 4096;
        const std::size_t n = d.buf.size();
        const std::size_t nft = wlen;
        std::vector<double> psd(nft, 0.0);
        std::vector<double> win(wlen);
        for (std::size_t i = 0; i < wlen; ++i)
            win[i] = 0.5 - 0.5 * std::cos(2.0 * M_PI * i / wlen);
        for (std::size_t s = 0; s + wlen <= n; s += wlen / 2) {
            std::vector<iqclassify::detail::cd> b(nft, iqclassify::detail::cd(0, 0));
            for (std::size_t t = 0; t < wlen; ++t) {
                const auto v = d.buf[s + t];
                if (std::isfinite(v.real()) && std::isfinite(v.imag()))
                    b[t] = iqclassify::detail::cd(v.real() * win[t], v.imag() * win[t]);
            }
            iqclassify::detail::fftInPlace(b, false);
            for (std::size_t k = 0; k < nft; ++k) {
                const double p = std::norm(b[k]);
                if (p > psd[k]) psd[k] = p;
            }
        }
        d.psdNft = nft;
        d.psdDb.resize(nft);
        for (std::size_t k = 0; k < nft; ++k)
            d.psdDb[k] = iqclassify::detail::db10(psd[k]);
    }
    return d;
}

// ---------------------------------------------------------------------------
// Bild 1: Pipeline-Blockschaltbild
// ---------------------------------------------------------------------------
void renderPipeline(const std::string& outDir) {
    const int W = 1160, H = 400;
    Canvas c(W, H, 0);
    const uint8_t white = c.color(255, 255, 255);
    const uint8_t black = c.color(20, 20, 20);
    const uint8_t blue = c.color(41, 98, 180);
    const uint8_t lblue = c.color(219, 231, 248);
    const uint8_t green = c.color(46, 125, 50);
    const uint8_t lgreen = c.color(226, 240, 217);
    const uint8_t orange = c.color(198, 110, 20);
    const uint8_t lorange = c.color(252, 236, 214);
    const uint8_t gray = c.color(110, 110, 110);
    c.fillRect(0, 0, W - 1, H - 1, white);
    c.text(20, 14, "IQCLASSIFY - PIPELINE (ALLE FREQUENZEN AUF FS = 1 NORMIERT)", 2, black);
    c.hline(20, W - 20, 40, gray);

    struct Box {
        int x, y, w, h;
        const char* title;
        const char* sub;
        uint8_t border, fill;
    };
    const int bw = 236, bh = 74;
    const int y1 = 70, y2 = 230;
    const Box row1[] = {
        {20, y1, bw, bh, "IQ-PUFFER", "KOMPLEXER ABTASTSTROM I/Q", blue, lblue},
        {314, y1, bw, bh, "0 KANALISIERUNG", "KANAELE UM ENERGIEBEREICHE", blue, lblue},
        {608, y1, bw, bh, "1 SEGMENTIERUNG", "BLOCKENERGIE, HYSTERESE, BODEN", blue, lblue},
        {902, y1, bw, bh, "2 BAND+EXTRAKTION", "MAX-HOLD-PSD, MASKE, MIX, RESAMPLE", blue, lblue},
    };
    const Box row2[] = {
        {20, y2, bw, bh, "3 MERKMALE", "HUELLKURVE, F_INST, CHIRP-FIT", green, lgreen},
        {314, y2, bw, bh, "4 STUFE 1 (GROB)", "TON / FM / PM / FLACH (OFDM)", green, lgreen},
        {608, y2, bw, bh, "5 STUFE 2 (FEIN)", "CW BPSK QPSK FSK2 CHIRP OFDM QAM", green, lgreen},
        {902, y2, bw, bh, "PULSERESULT", "+ LOKALISIERUNG + ANALYSE-API", orange, lorange},
    };
    auto drawBox = [&](const Box& b) {
        c.fillRect(b.x, b.y, b.x + b.w - 1, b.y + b.h - 1, b.fill);
        c.rect(b.x, b.y, b.x + b.w - 1, b.y + b.h - 1, b.border);
        c.textCentered(b.x + b.w / 2, b.y + 12, b.title, 2, black);
        c.textCentered(b.x + b.w / 2, b.y + 40, b.sub, 1, gray);
    };
    for (const auto& b : row1) drawBox(b);
    for (const auto& b : row2) drawBox(b);
    for (int i = 0; i + 1 < 4; ++i)
        c.arrow(row1[i].x + row1[i].w + 4, y1 + bh / 2,
                row1[i + 1].x - 6, y1 + bh / 2, black);
    for (int i = 0; i + 1 < 4; ++i)
        c.arrow(row2[i].x + row2[i].w + 4, y2 + bh / 2,
                row2[i + 1].x - 6, y2 + bh / 2, black);
    // Uebergang Zeile 1 -> Zeile 2
    c.arrow(row1[3].x + bw / 2, y1 + bh + 4, row2[0].x + bw / 2, y2 - 6, black);
    c.text(20, 330, "JEDER PULS: MODULATION, KONFIDENZ, SNR, FEINES ZEITFENSTER [STARTSAMPLEFINE, +NUMSAMPLESFINE),", 1, gray);
    c.text(20, 346, "SUB-BIN-MITTENFREQUENZ CENTERFREQREFINED; EXTRACTPULSE LIEFERT DEN BASISBAND-SLICE AUS DEM IQ-PUFFER.", 1, gray);
    c.text(20, 372, "FADING-ROBUST: FADE-MERGE (AMPLITUDENSCHWANKUNG BLEIBT EIN PULS).", 1, gray);
    if (!c.save(outDir + "/pipeline.png"))
        std::fprintf(stderr, "Fehler beim Schreiben von pipeline.png\n");
}

// ---------------------------------------------------------------------------
// Bild 2: Detektion + Lokalisierung (Huellkurve + Spektrum)
// ---------------------------------------------------------------------------
void renderDetection(const Demo& d, const std::string& outDir) {
    const int W = 1160, H = 800;
    Canvas c(W, H, 0);
    const uint8_t white = c.color(255, 255, 255);
    const uint8_t black = c.color(20, 20, 20);
    const uint8_t gray = c.color(150, 150, 150);
    const uint8_t lgray = c.color(225, 225, 225);
    const uint8_t blue = c.color(41, 98, 180);
    const uint8_t red = c.color(200, 40, 40);
    const uint8_t green = c.color(46, 150, 60);
    const uint8_t lyellow = c.color(255, 246, 190);
    const uint8_t lblue = c.color(214, 232, 252);
    c.fillRect(0, 0, W - 1, H - 1, white);
    c.text(20, 12, "DETEKTION UND LOKALISIERUNG (DEMO-PUFFER, SNR 20 DB, FS = 1)", 2, black);
    c.hline(20, W - 20, 38, gray);

    const int ml = 90, mr = 30, mt = 60;
    const int pw = W - ml - mr;  // 1040

    // ---- oberes Panel: Huellkurve ----
    const int p1y0 = mt, p1y1 = 350;
    const int ph1 = p1y1 - p1y0;
    const std::size_t n = d.buf.size();
    // feine Fenster einfaerben
    for (const auto& p : d.pulses) {
        const int x0 = ml + static_cast<int>(p.startSampleFine * pw / n);
        const int x1 = ml + static_cast<int>((p.startSampleFine + p.numSamplesFine) * pw / n);
        c.fillRect(x0, p1y0, x1, p1y1 - 1, lyellow);
    }
    // Achsen
    c.rect(ml, p1y0, ml + pw - 1, p1y1 - 1, black);
    for (int i = 0; i <= 4; ++i) {
        const int x = ml + i * pw / 4;
        c.vline(x, p1y1 - 1, p1y1 + 4, black);
        c.textCentered(x, p1y1 + 8, std::to_string(i * n / 4), 1, black);
    }
    c.text(ml + pw / 2 - 60, p1y1 + 24, "SAMPLES", 1, black);
    const double dbLo = -50, dbHi = 0;
    for (int i = 0; i <= 5; ++i) {
        const double db = dbLo + (dbHi - dbLo) * i / 5;
        const int y = p1y1 - 1 - static_cast<int>((db - dbLo) / (dbHi - dbLo) * (ph1 - 1));
        c.hline(ml, ml + pw - 1, y, lgray);
        c.text(ml - 34, y - 3, fmt(db, 0) + "DB", 1, black);
    }
    // Huellkurve (pro Pixelspalte min/max als vertikaler Strich)
    double emax = 1e-30;
    for (const auto& v : d.buf) emax = std::max(emax, std::norm(v));
    const int colPerSample = std::max(1, static_cast<int>(n / pw));
    for (int x = 0; x < pw; ++x) {
        const std::size_t i0 = static_cast<std::size_t>(x) * n / pw;
        const std::size_t i1 = std::min(n, i0 + colPerSample + 1);
        double lo = 1e300, hi = 0.0;
        for (std::size_t i = i0; i < i1; ++i) {
            const double e = std::norm(d.buf[i]);
            lo = std::min(lo, e);
            hi = std::max(hi, e);
        }
        auto Y = [&](double e) {
            const double db = 10.0 * std::log10(e / emax);
            return p1y1 - 1 - static_cast<int>((db - dbLo) / (dbHi - dbLo) * (ph1 - 1));
        };
        const int y0 = std::min(ph1 - 1, std::max(0, Y(hi)));
        const int y1 = std::min(ph1 - 1, std::max(0, Y(lo)));
        c.vline(ml + x, p1y0 + y0, p1y0 + y1, blue);
    }
    // Labels der Pulse
    const uint8_t labCols[] = {red, green, blue};
    for (std::size_t i = 0; i < d.pulses.size(); ++i) {
        const auto& p = d.pulses[i];
        const int xc = ml + static_cast<int>((p.startSampleFine + p.numSamplesFine / 2) * pw / n);
        const uint8_t col = labCols[i % 3];
        std::string t1 = iqclassify::toString(p.modulation);
        std::string t2 = "F=" + fmt(p.centerFreqRefined, 3);
        c.textCentered(std::max(ml + 20, std::min(W - mr - 90, xc)), p1y0 + 6, t1, 2, col);
        c.textCentered(std::max(ml + 20, std::min(W - mr - 90, xc)), p1y0 + 24, t2, 1, col);
    }
    c.text(ml - 70, p1y0 - 14, "|IQ|", 1, black);

    // ---- unteres Panel: Spektrum ----
    const int p2y0 = 440, p2y1 = 720;
    const int ph2 = p2y1 - p2y0;
    // Baender einfaerben
    for (const auto& b : d.bands) {
        auto X = [&](long k) {
            double f = static_cast<double>(k) / b.nft;
            if (f >= 1.0) f -= 1.0;
            if (f < 0.0) f += 1.0;
            if (f > 0.5) f -= 1.0;
            return ml + static_cast<int>((f + 0.5) * pw);
        };
        c.fillRect(X(b.lo), p2y0, X(b.hi - 1), p2y1 - 1, lblue);
    }
    c.rect(ml, p2y0, ml + pw - 1, p2y1 - 1, black);
    double psdMin = 1e300, psdMax = -1e300;
    for (double v : d.psdDb) {
        psdMin = std::min(psdMin, v);
        psdMax = std::max(psdMax, v);
    }
    const double sLo = psdMin, sHi = psdMax;
    for (int i = 0; i <= 4; ++i) {
        const double f = -0.5 + i * 0.25;
        const int x = ml + static_cast<int>((f + 0.5) * pw);
        c.vline(x, p2y1 - 1, p2y1 + 4, black);
        c.textCentered(x, p2y1 + 8, fmt(f, 2), 1, black);
    }
    c.text(ml + pw / 2 - 30, p2y1 + 24, "F (FS = 1)", 1, black);
    for (int i = 0; i <= 4; ++i) {
        const double db = sLo + (sHi - sLo) * i / 4;
        const int y = p2y1 - 1 - static_cast<int>((db - sLo) / (sHi - sLo) * (ph2 - 1));
        c.hline(ml, ml + pw - 1, y, lgray);
        c.text(ml - 44, y - 3, fmt(db, 0), 1, black);
    }
    auto Y2 = [&](double db) {
        return p2y1 - 1 - static_cast<int>((db - sLo) / (sHi - sLo) * (ph2 - 1));
    };
    for (int x = 0; x < pw; ++x) {
        const std::size_t k = static_cast<std::size_t>(x) * d.psdNft / pw;
        c.set(ml + x, std::max(p2y0, std::min(p2y1 - 1, Y2(d.psdDb[k]))), black);
    }
    // Mittenfrequenzen markieren
    for (std::size_t i = 0; i < d.pulses.size(); ++i) {
        const auto& p = d.pulses[i];
        double f = p.centerFreqRefined;
        if (f > 0.5) f -= 1.0;
        const int x = ml + static_cast<int>((f + 0.5) * pw);
        const uint8_t col = labCols[i % 3];
        c.vline(x, p2y0, p2y1 - 1, col);
        c.textCentered(std::max(ml + 20, std::min(W - mr - 60, x)), p2y0 - 12,
                       fmt(p.centerFreqRefined, 3), 1, col);
    }
    c.text(ml - 40, p2y0 - 14, "PSD DB", 1, black);
    c.text(20, 760, "GELB: FEINES ZEITFENSTER (STARTSAMPLEFINE/NUMSAMPLESFINE)   BLAU: DETEKTIERTES BAND   FARBIGE LINIEN: CENTERFREQREFINED", 1, gray);
    if (!c.save(outDir + "/detection.png"))
        std::fprintf(stderr, "Fehler beim Schreiben von detection.png\n");
}

// ---------------------------------------------------------------------------
// Bild 3: Spektrogramm mit ueberlagerten Puls-Rechtecken
// ---------------------------------------------------------------------------
// Heatmap: t = 0 (laut) -> hell/gelb, t = 1 (leise) -> dunkelblau
uint8_t heat(Canvas& c, double t) {
    t = std::max(0.0, std::min(1.0, t));
    const double stops[5][3] = {
        {255.0, 252.0, 214.0},  // 0.00 laut: hellgelb
        {254.0, 196.0, 60.0},   // 0.25
        {240.0, 110.0, 30.0},   // 0.50 orange/rot
        {150.0, 40.0, 90.0},    // 0.75
        {12.0, 16.0, 70.0},     // 1.00 leise: dunkelblau
    };
    const double x = t * 4.0;
    const int i = std::min(3, static_cast<int>(x));
    const double f = x - i;
    const double r = stops[i][0] + f * (stops[i + 1][0] - stops[i][0]);
    const double g = stops[i][1] + f * (stops[i + 1][1] - stops[i][1]);
    const double b = stops[i][2] + f * (stops[i + 1][2] - stops[i][2]);
    return c.color(static_cast<uint8_t>(r), static_cast<uint8_t>(g),
                   static_cast<uint8_t>(b));
}

void renderSpectrogram(const Demo& d, const std::string& outDir) {
    const int W = 1160, H = 620;
    Canvas c(W, H, 0);
    const uint8_t white = c.color(255, 255, 255);
    const uint8_t black = c.color(20, 20, 20);
    const uint8_t gray = c.color(150, 150, 150);
    c.fillRect(0, 0, W - 1, H - 1, white);
    c.text(20, 12, "SPEKTROGRAMM DES DEMO-PUFFERS MIT LOKALISIERTEN PULSEN (ZEIT x FREQUENZ)", 2, black);
    c.hline(20, W - 20, 38, gray);

    const int ml = 90, mr = 30, mt = 60, mb = 46;
    const int pw = W - ml - mr, ph = H - mt - mb;

    // STFT: Fenster 256, Hop 64
    const std::size_t wlen = 256, hop = 64;
    const std::size_t n = d.buf.size();
    const std::size_t nframes = n >= wlen ? (n - wlen) / hop + 1 : 0;
    std::vector<double> win(wlen);
    for (std::size_t i = 0; i < wlen; ++i)
        win[i] = 0.5 - 0.5 * std::cos(2.0 * M_PI * i / wlen);
    std::vector<double> smag(nframes * (wlen / 2 + 1), 0.0);
    double gmax = 1e-30;
    for (std::size_t fr = 0; fr < nframes; ++fr) {
        std::vector<iqclassify::detail::cd> b(wlen, iqclassify::detail::cd(0, 0));
        for (std::size_t t = 0; t < wlen; ++t) {
            const auto v = d.buf[fr * hop + t];
            b[t] = iqclassify::detail::cd(v.real() * win[t], v.imag() * win[t]);
        }
        iqclassify::detail::fftInPlace(b, false);
        for (std::size_t k = 0; k <= wlen / 2; ++k) {
            const double p = std::norm(b[k]);
            smag[fr * (wlen / 2 + 1) + k] = p;
            gmax = std::max(gmax, p);
        }
    }
    const double dbRef = 10.0 * std::log10(gmax);
    const double dbSpan = 55.0;  // dB Dynamik
    // Bild: x = Frame (Zeit), y = Frequenz (oben +0.5, unten -0.5)
    for (int x = 0; x < pw; ++x) {
        const std::size_t fr = static_cast<std::size_t>(x) * nframes / pw;
        for (int y = 0; y < ph; ++y) {
            const double f = 0.5 - static_cast<double>(y) / ph;  // +0.5 .. -0.5
            std::size_t k = static_cast<std::size_t>(std::abs(f) * wlen);
            if (k > wlen / 2) k = wlen / 2;
            const double db = 10.0 * std::log10(smag[fr * (wlen / 2 + 1) + k] + 1e-30);
            double t = (dbRef - db) / dbSpan;  // 0 = laut, 1 = leise
            c.set(ml + x, mt + y, heat(c, t));
        }
    }
    // Achsen + Beschriftung
    c.rect(ml, mt, ml + pw - 1, mt + ph - 1, black);
    for (int i = 0; i <= 4; ++i) {
        const int x = ml + i * pw / 4;
        c.vline(x, mt + ph - 1, mt + ph + 4, black);
        c.textCentered(x, mt + ph + 8, std::to_string(i * n / 4), 1, black);
    }
    c.text(ml + pw / 2 - 30, mt + ph + 24, "SAMPLES", 1, black);
    for (int i = 0; i <= 4; ++i) {
        const double f = 0.5 - i * 0.25;
        const int y = mt + static_cast<int>((0.5 - f) * ph);
        c.hline(ml - 4, ml, y, black);
        c.text(ml - 40, y - 3, fmt(f, 2), 1, black);
    }
    c.text(20, mt + ph / 2 - 4, "F", 1, black);
    // Puls-Rechtecke: [startSampleFine, end) x [fc-bw/2, fc+bw/2]
    for (const auto& p : d.pulses) {
        const int x0 = ml + static_cast<int>(p.startSampleFine * pw / n);
        const int x1 = ml + static_cast<int>((p.startSampleFine + p.numSamplesFine) * pw / n);
        double f0 = p.centerFreqRefined;
        if (f0 > 0.5) f0 -= 1.0;
        const double bw = p.bandwidthNorm;
        auto Y = [&](double f) {
            return mt + static_cast<int>((0.5 - f) * ph);
        };
        const int y0 = Y(f0 + bw / 2), y1 = Y(f0 - bw / 2);
        c.rect(x0, y0, x1, y1, black);
        c.rect(x0 + 1, y0 + 1, x1 - 1, y1 - 1, white);
        const std::string label = std::string(iqclassify::toString(p.modulation)) +
                                  " F=" + fmt(p.centerFreqRefined, 3);
        const int tw = c.textWidth(label, 1);
        const int lx = std::max(ml + 4, std::min(ml + pw - tw - 6, x0 + 4));
        const int ly = std::max(mt + 4, std::min(mt + ph - 12, y0 - 12));
        c.fillRect(lx - 2, ly - 2, lx + tw + 2, ly + 9, white);
        c.text(lx, ly, label, 1, black);
    }
    if (!c.save(outDir + "/spectrogram.png"))
        std::fprintf(stderr, "Fehler beim Schreiben von spectrogram.png\n");
}

// ---------------------------------------------------------------------------
// Bild 4: Polyphase-Codes und Barker-Code mit echtem Klassifikations-Label
// ---------------------------------------------------------------------------
struct CodeDemo {
    std::string title;
    std::vector<iqtest::cd> sig;
    double fc;
};

void renderCodes(const std::string& outDir) {
    std::vector<CodeDemo> demos;
    {
        std::mt19937 rng(9);
        demos.push_back({"FRANK M=8 (SPS=8, F=0.15)", iqtest::genFrank(rng, 8, 8), 0.15});
    }
    for (int pp : {1, 2}) {
        std::mt19937 rng(50 + pp);
        char t[64];
        std::snprintf(t, sizeof(t), "P%d M=8 (SPS=8, F=0.15)", pp);
        demos.push_back({t, iqtest::genPcode(rng, pp, 8, 8), 0.15});
    }
    {
        std::mt19937 rng(70);
        demos.push_back({"P3 M=64 (SPS=8, F=0.15)", iqtest::genPcode(rng, 3, 64, 8), 0.15});
    }
    {
        std::mt19937 rng(71);
        demos.push_back({"P4 M=64 (SPS=8, F=0.15)", iqtest::genPcode(rng, 4, 64, 8), 0.15});
    }
    {
        std::mt19937 rng(100);
        demos.push_back({"BARKER13 x8 (SPS=8, F=0.15)", iqtest::genBarker(rng, 13, 8, 8), 0.15});
    }

    const int W = 1160, H = 920;
    Canvas c(W, H, 0);
    const uint8_t white = c.color(255, 255, 255);
    const uint8_t black = c.color(20, 20, 20);
    const uint8_t gray = c.color(150, 150, 150);
    const uint8_t lgray = c.color(228, 228, 228);
    const uint8_t blue = c.color(41, 98, 180);
    const uint8_t red = c.color(200, 40, 40);
    const uint8_t green = c.color(46, 150, 60);
    c.fillRect(0, 0, W - 1, H - 1, white);
    c.text(20, 12, "POLYPHASE-CODES UND BARKER-CODE: F_INST DES PULSES UND KLASSIFIKATIONS-LABEL (SNR 20 DB)", 2, black);
    c.hline(20, W - 20, 38, gray);

    iqclassify::ClassifierConfig cfg;
    iqclassify::ModulationClassifier clf(cfg);

    const int cols = 2, rows = 3;
    const int ml = 20, mt = 56, gap = 16;
    const int pw = (W - 2 * ml - (cols - 1) * gap) / cols;   // 554
    const int ph = (H - mt - (rows - 1) * gap - 30) / rows;  // ~272

    for (std::size_t i = 0; i < demos.size(); ++i) {
        const auto& dm = demos[i];
        const int x0 = ml + static_cast<int>(i % cols) * (pw + gap);
        const int y0 = mt + static_cast<int>(i / cols) * (ph + gap);
        // Klassifikation des echten Puffers
        iqtest::Recorder rec;
        rec.silence(2000);
        rec.add(dm.sig, dm.fc, 20.0);
        rec.silence(2000);
        const auto res = clf.classify(rec.buffer().data(), rec.buffer().size());
        std::string label = res.empty() ? "KEIN PULS" :
                            iqclassify::toString(res[0].modulation);
        std::string sub;
        if (!res.empty() && res[0].modulation == iqclassify::ModulationType::PCode) {
            const auto pp = clf.analyzePcode(rec.buffer().data(), rec.buffer().size(), res[0]);
            if (pp.valid) {
                char b[64];
                std::snprintf(b, sizeof(b), " SUBTYP P%d, M=%d, N=%zu", pp.subtype, pp.order, pp.codeLength);
                sub = b;
            }
        }
        if (!res.empty() && res[0].modulation == iqclassify::ModulationType::BarkerCode) {
            const auto bp = clf.analyzeBarker(rec.buffer().data(), rec.buffer().size(), res[0]);
            if (bp.valid) {
                char b[64];
                std::snprintf(b, sizeof(b), " LEN=%d, REPS=%d, SPS=%.1f", bp.length, bp.repeats, bp.sps);
                sub = b;
            }
        }
        // Panel
        c.rect(x0, y0, x0 + pw - 1, y0 + ph - 1, black);
        c.text(x0 + 8, y0 + 6, dm.title, 1, black);
        const std::string lab = "-> " + label + sub;
        c.text(x0 + 8, y0 + 20, lab, 1, red);
        // f_inst des Basisband-Pulses
        const auto& sig = dm.sig;
        const std::size_t n = sig.size();
        std::vector<double> up = iqclassify::detail::unwrapPhase(sig);
        std::vector<double> fi(n >= 1 ? n - 1 : 0);
        for (std::size_t k = 0; k + 1 < n; ++k)
            fi[k] = (up[k + 1] - up[k]) / (2.0 * M_PI);
        double fmin = 1e300, fmax = -1e300;
        for (double v : fi) {
            fmin = std::min(fmin, v);
            fmax = std::max(fmax, v);
        }
        if (fmax - fmin < 1e-9) fmax = fmin + 1e-9;
        const int gx0 = x0 + 10, gy0 = y0 + 44, gw = pw - 20, gh = ph - 60;
        c.rect(gx0, gy0, gx0 + gw - 1, gy0 + gh - 1, lgray);
        const double mid = 0.5 * (fmin + fmax);
        c.hline(gx0, gx0 + gw - 1, gy0 + gh / 2, gray);
        c.text(gx0 - 0, gy0 + gh + 4, "SAMPLES", 1, black);
        c.text(gx0, gy0 - 10, fmt(fmax, 4), 1, gray);
        c.text(gx0, gy0 + gh - 2, fmt(fmin, 4), 1, gray);
        c.text(gx0 + gw - 34, gy0 - 10, "F_INST", 1, gray);
        const int colPerSample = std::max(1, static_cast<int>(n / gw));
        for (int x = 0; x < gw; ++x) {
            const std::size_t k0 = static_cast<std::size_t>(x) * n / gw;
            const std::size_t k1 = std::min(fi.size(), k0 + colPerSample + 1);
            double lo = 1e300, hi = -1e300;
            for (std::size_t k = k0; k < k1; ++k) {
                lo = std::min(lo, fi[k]);
                hi = std::max(hi, fi[k]);
            }
            if (lo > hi) continue;
            auto Y = [&](double v) {
                return gy0 + gh - 1 - static_cast<int>((v - fmin) / (fmax - fmin) * (gh - 1));
            };
            c.vline(gx0 + x, std::max(gy0, Y(hi)), std::min(gy0 + gh - 1, Y(lo)), blue);
        }
        (void)green;
    }
    c.text(20, H - 18, "BLAU: MOMENTANFREQUENZ F_INST DES GENERIRTEN PULSES   ROT: LABEL DES KLASSIFIKATORS (ECHTE PIPELINE)", 1, gray);
    if (!c.save(outDir + "/codes.png"))
        std::fprintf(stderr, "Fehler beim Schreiben von codes.png\n");
}

}  // namespace

int main(int argc, char** argv) {
    const std::string outDir = argc > 1 ? argv[1] : "docs/images";
    Demo d = makeDemo();
    std::printf("Demo-Puffer: %zu Samples, %zu Pulse erkannt:", d.buf.size(),
                d.pulses.size());
    for (const auto& p : d.pulses)
        std::printf(" %s(f=%.3f @%zu)", iqclassify::toString(p.modulation),
                    p.centerFreqRefined, p.startSampleFine);
    std::printf("\n");
    renderPipeline(outDir);
    renderDetection(d, outDir);
    renderSpectrogram(d, outDir);
    renderCodes(outDir);
    std::printf("Bilder geschrieben nach %s\n", outDir.c_str());
    return 0;
}
