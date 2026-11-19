# Performance

Messung mit dem Benchmark-Tool des Repos (deterministische IQ-Puffer, 20 dB
SNR, Pulse vom Typ QPSK/CW/FSK2/Chirp/QAM16; bis 2 Mio Samples):

    cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
    ./build/iqclassify_benchmark

## Gemessene Laufzeiten nach der Optimierung (Release, Linux x86-64, Median)

Puffer mit `max(3, n/8000)` sequenziellen Pulsen; Vergleich direkt auf derselben
Maschine gegen den Stand vor der Optimierung (Commit 6d56fea):

| Samples | Pulse | classify vorher | classify nachher | legacy nachher | Speedup | MS/s nachher |
|--------:|------:|----------------:|-----------------:|---------------:|--------:|------------:|
|     10 000 |     3 |    7,7 ms |    5,7 ms |    5,1 ms |  1,35× |   1,7 |
|     50 000 |     3 |    6,0 ms |    5,2 ms |    3,7 ms |  1,15× |   9,6 |
|    100 000 |    12 |   23,7 ms |   19,5 ms |   16,2 ms |  1,22× |   5,1 |
|    250 000 |    31 |   66,2 ms |   49,5 ms |   41,7 ms |  1,34× |   5,1 |
|    500 000 |    62 |  121,9 ms |   97,2 ms |   81,5 ms |  1,25× |   5,1 |
|  1 000 000 |   125 |  256,8 ms |  200,0 ms |  162,2 ms |  1,28× |   5,0 |
|  2 000 000 |   250 |  496,7 ms |  389,0 ms |  328,1 ms |  1,28× |   5,1 |

Alle Pulse werden gefunden (250/250 bei 2M). Peak-RSS ~195 MiB (O(n)).

Spärliche Puffer profitieren stärker vom Energie-Gate der Detektion:
2 Mio Samples mit nur 3 Pulsen: **~20 ms (~100 MSamples/s**, vorher ~39 ms).

## Was optimiert wurde

- **Energie-Gate in der Banddetektion** (`BandExtractor::bands`): Fenster, deren
  Energie unter `1.5 × Rauschboden × Fensterlänge` liegt, werden nicht per FFT
  gerechnet (tragen zum Max-Hold nichts bei). Referenz ist der globale
  Puffer-Rauschboden (`PulseSegmenter::noiseFloor`). Spart bei spärlichen
  Puffern ~99 % der Detektionsfenster.
- **Schnelle zyklische Kohärenz** in der Sub-Bin-Symbolratensuche
  (`rsFromCyclic`/`rsLineFull`): inkrementelle Drehung (4 Mults/Sample statt
  sin/cos) — das war der dominierende Kostenblock der Stufe 2.
- **Kein Allokieren pro FFT-Aufruf** (`detail::fftInPlace`): shape/stride/axes
  thread-lokal gecacht; Detektionspuffer wird wiederverwendet.
- **SNR-Schätzung ohne `std::set`** (`estimateSnrFlat`): Byte-Markierung O(nft)
  statt O(Bänder × Breite × log) pro Aufruf.

## Einordnung

- **Skalierung ~linear (O(n))**: ab ~100k Samples konstant ~5 MSamples/s.
  Detektion arbeitet mit gleitenden Fenstern fester Länge (≤ 4096), Extraktion
  und Merkmale laufen pro Puls-Segment, nicht über den ganzen Puffer.
- **Kanalisierung kostet ~10-15 %** gegenüber dem Legacy-Pfad (zusätzliche
  Detektion auf dem Gesamtpuffer + Kanal-Routing; durch das Gate bei
  spärlichen Puffern nahezu kostenlos).
- **Speicher**: O(n), Eingabepuffer wird nicht kopiert.
- Die Analyse-API arbeitet nur auf dem extrahierten Puls-Slice (O(Pulslänge)).
- Resampling (`channelDecimate > 1`) reduziert die Rate pro Kanal und damit
  die Kosten der Folgestufen entsprechend.

## Parallelverarbeitung der Kanäle (threads)

Die Kanalverarbeitung ist parallelisiert (`ClassifierConfig::threads`,
Default 0 = automatisch min(hardware_concurrency, 8), 1 = seriell, wirkt ab
~32k Samples Pufferlänge): Segment-Detektion und Extraktion/Klassifikation
laufen pro (Segment, Kanalband) auf eigenen Threads, `Channelizer::channelize`
pro Kanal. pocketfft-Plancache und `fftInPlace` sind thread-sicher.
Die Ergebnisse sind **identisch unabhängig von der Threadanzahl**
(pro Thread eigene Buckets, danach deterministisch sortiert; Test Abschnitt 19).

Messung (5 Kerne, Median, 2M Samples / 250 Pulse):

| Modus | 100k | 500k | 1M | 2M (250 Pulse) | MS/s bei 2M |
|-------|-----:|-----:|---:|---------------:|------------:|
| threads=1 (seriell) | 19,3 ms | 101 ms | 199 ms | **393 ms** | 5,1 |
| threads=auto (5) | 8,7 ms | 39 ms | 76 ms | **~150 ms** | ~13 |

→ **~2,5-3× Speedup** durch die Parallelverarbeitung bei 2M Samples.

## Skalierungs-Test

Die Testsuite enthält Abschnitt 17 (`scale-2M`): 2 000 000 Samples mit 3
Pulsen werden in einem `classify`-Durchlauf korrekt klassifiziert.
