# IQClassify - qmake-Projektdatei fuer die statische Bibliothek.
#
#   qmake IQClassify.pro && make          # Bibliothek bauen
#   (cd tests && qmake IQClassifyTests.pro && make && ./iqclassify_tests)
#
# Reine C++17-Bibliothek ohne Qt-Modulabhaengigkeit.

TARGET      = IQClassify
TEMPLATE    = lib
CONFIG     += staticlib
CONFIG     += c++17
QT         -= core gui

HEADERS += \
    include/iqclassify/iqclassify.hpp \
    include/iqclassify/Types.hpp \
    include/iqclassify/Config.hpp \
    include/iqclassify/Dsp.hpp \
    include/iqclassify/PulseSegmenter.hpp \
    include/iqclassify/BandExtractor.hpp \
    include/iqclassify/Channelizer.hpp \
    include/iqclassify/FeatureExtractor.hpp \
    include/iqclassify/Classifier.hpp \
    third_party/pocketfft/pocketfft_hdronly.h

SOURCES += \
    src/Dsp.cpp \
    src/PulseSegmenter.cpp \
    src/BandExtractor.cpp \
    src/Channelizer.cpp \
    src/FeatureExtractor.cpp \
    src/Classifier.cpp

INCLUDEPATH += include third_party/pocketfft
