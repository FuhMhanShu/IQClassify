// IQClassify - Schirmheader: bindet die komplette Bibliothek ein.
//
//   #include "iqclassify/iqclassify.hpp"
//   iqclassify::ClassifierConfig cfg;               // Schwellwerte
//   iqclassify::ModulationClassifier clf(cfg);
//   auto pulses = clf.classify(iq.data(), iq.size());

#pragma once

#include "Config.hpp"
#include "Types.hpp"
#include "PulseSegmenter.hpp"
#include "BandExtractor.hpp"
#include "FeatureExtractor.hpp"
#include "Classifier.hpp"
