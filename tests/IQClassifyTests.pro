# IQClassify - qmake-Projektdatei der Testsuite.
#
#   qmake IQClassifyTests.pro && make && ./iqclassify_tests
#
# Kompiliert die Bibliotheksquellen direkt hinein (kein Link-Schritt
# gegen eine zuvor gebaute Bibliothek noetig).

TARGET      = iqclassify_tests
TEMPLATE    = app
CONFIG     += console c++17
CONFIG     -= app_bundle
QT         -= core gui

HEADERS += \
    SignalGenerator.hpp \
    ../include/iqclassify/iqclassify.hpp \
    ../include/iqclassify/Types.hpp \
    ../include/iqclassify/Config.hpp \
    ../include/iqclassify/Dsp.hpp \
    ../include/iqclassify/PulseSegmenter.hpp \
    ../include/iqclassify/BandExtractor.hpp \
    ../include/iqclassify/FeatureExtractor.hpp \
    ../include/iqclassify/Classifier.hpp \
    ../third_party/pocketfft/pocketfft_hdronly.h

SOURCES += \
    SignalGenerator.cpp \
    test_main.cpp \
    ../src/Dsp.cpp \
    ../src/PulseSegmenter.cpp \
    ../src/BandExtractor.cpp \
    ../src/FeatureExtractor.cpp \
    ../src/Classifier.cpp

INCLUDEPATH += . ../include ../third_party/pocketfft
