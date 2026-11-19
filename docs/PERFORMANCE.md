# Performance

Messung mit dem Benchmark-Tool des Repos (deterministische IQ-Puffer, 20 dB
SNR, Pulse vom Typ QPSK/CW/FSK2/Chirp/QAM16; bis 2 Mio Samples):

    cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
    ./build/iqclassify_benchmark

## Gemessene Laufzeiten (Release-Build, Linux x86-64, Median)

Puffer mit `max(3, n/8000)` sequenziellen Pulsen:

| Samples | Pulse | classify (kanalisiert) | legacy (`channelize=false`) | nur Kanalisierung | gefunden | MSamples/s |
|--------:|------:|----------------------:|----------------------------:|------------------:|---------:|-----------:|
|     10 000 |     3 |   7.4 ms |  7.0 ms |   2.7 ms |   3/3 |   1.3 |
|     50 000 |     3 |   5.9 ms |  4.8 ms |  12.0 ms |   3/3 |   8.5 |
|    100 000 |    12 |  23.7 ms | 22.1 ms |   9.9 ms |  12/12 |   4.2 |
|    250 000 |    31 |  61.9 ms | 56.7 ms |  21.6 ms |  31/31 |   4.0 |
|    500 000 |    62 | 119.9 ms | 112.1 ms |  44.1 ms |  62/62 |   4.2 |
|  1 000 000 |   125 | 255.3 ms | 235.4 ms |  99.6 ms | 125/125 |   3.9 |
|  2 000 000 |   250 | 494.5 ms | 453.9 ms | 230.3 ms | 250/250 |   4.0 |

Peak-RSS des Benchmark-Prozesses (inkl. 2M-Puffer, 32 MiB): ~193 MiB.

Ein Puffer mit nur wenigen Pulsen ist entsprechend schneller: 2 Mio Samples
mit 3 Pulsen klassifiziert die Testsuite in ~39 ms (~51 MSamples/s).

## Einordnung

- **Skalierung ~linear (O(n))**: ab ~100k Samples konstant ~4 MSamples/s.
  Ursache: die Detektion arbeitet mit gleitenden Fenstern fester Länge
  (≤ 4096, 50 % Überlappung), die Extraktion und Merkmalsextraktion laufen pro
  Puls-Segment (nicht über den ganzen Puffer).
- **Kanalisierung kostet ~5-10 %** gegenüber dem Legacy-Pfad (eine zusätzliche
  Detektion auf dem Gesamtpuffer plus Kanal-Routing); die reine Kanalisierung
  (`Channelizer::channelize`, volle Puffer-FFT je Kanal) liegt bei ~0.1 s für
  2M Samples.
- **Speicher**: O(n) für den Eingabepuffer (nicht kopiert) + O(n) Arbeits-
  puffer in der Kanalisierung/Extraktion; ~193 MiB Peak bei 2M Samples.
- Die Analyse-API (`synchronizePsk`, `analyzeChirp`, `extractSymbols`) arbeitet
  nur auf dem extrahierten Puls-Slice (O(Pulslänge)), nicht auf dem Puffer.
- Resampling (`channelDecimate > 1`) reduziert die Rate pro Kanal und damit
  die Kosten der Folgestufen entsprechend (bei gleicher Pulsanzahl).

## Skalierungs-Test

Die Testsuite enthält Abschnitt 17 (`scale-2M`): 2 000 000 Samples mit 3
Pulsen werden in einem `classify`-Durchlauf korrekt klassifiziert
(Laufzeit wird ausgegeben, gegen Pathologien abgesichert).
