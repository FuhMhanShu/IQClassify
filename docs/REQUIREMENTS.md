# Anforderungen und Leistungsumfang

## Aufgabe

Klassifikation der Modulationsart von Pulsen in komplexen IQ-Abtastströmen
(zweistufig):

- Stufe 1, grobe Klasse: Ton, Frequenzmodulation (Chirp/FSK),
  Phasenmodulation (PSK/QAM), flächiges rauschartiges Signal (OFDM).
- Stufe 2, Modulationsart: Unmodulated (unmodulierter Träger, CW), BPSK,
  QPSK, FSK2, Chirp, OFDM, QAM16, QAM64, QAM256. Nicht entscheidbare Pulse
  werden als Unknown mit niedriger Konfidenz ausgegeben.

## Mehrpuls-Szenarien

Ein Puffer kann mehrere Pulse enthalten:

1. zeitlich getrennte Pulse (sequenziell, mit Schweigen dazwischen),
2. mehrere Pulse gleichzeitig in verschiedenen Frequenzbändern
   (frequenzgetrennt, überlagert in der Zeit).

Sind zwei Pulse sowohl zeitlich als auch spektral überlappt (gleiche
Frequenzlage), sind sie physikalisch nicht trennbar und werden als ein
Puls klassifiziert; dies ist eine dokumentierte Grenze, kein Fehler.

## Rahmenbedingungen

- Einbindung als statische Bibliothek, unterstützt via CMake
  (IQClassify STATIC, CTest-Integration) und qmake (CONFIG += staticlib,
  IQClassify.pri für Konsumentenprojekte).
- Keine Abhängigkeiten außer PocketFFT (header-only, unter
  third_party/pocketfft/ mitgeführt).
- C++17, keine RTTI-/Ausnahme-Anforderungen über die Standardbibliothek
  hinaus.
- Klassenkapselung: alle Verarbeitungsschritte in eigenen Klassen, eine
  öffentliche Fassade, alle Schwellwerte in einer Konfigurationsstruktur.
- Dokumentation von Anforderungen, Methodik und API in docs/.
- Deterministische Testsuite ohne externes Testframework.

## Nicht-Funktionsanforderungen

- Interne Rechnung in double; Eingabe zusätzlich als float komfortabel.
- Vorhersehbarer Speicherbedarf: Eingabepuffer werden nicht kopiert
  außer für fensterweise Spektralanalysen.
- Konfigurierbarkeit: Jeder Schwellwert der Referenzvalidierung bleibt
  einstellungsfähig (ClassifierConfig), siehe docs/METHODOLOGY.md.
