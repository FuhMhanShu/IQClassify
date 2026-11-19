cmake_minimum_required(VERSION 3.10)
project(IQClassify
    VERSION 1.0.0
    DESCRIPTION "Modulationsklassifikation von Pulsen in IQ-Signalen"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(IQCLASSIFY_BUILD_TESTS "Testsuite von IQClassify bauen" ON)

# Statische Bibliothek. PocketFFT (header-only, third_party) ist die
# einzige Abhaengigkeit und wird privat in src/Dsp.cpp eingebunden.
add_library(IQClassify STATIC
    src/BandExtractor.cpp
    src/Classifier.cpp
    src/Dsp.cpp
    src/FeatureExtractor.cpp
    src/PulseSegmenter.cpp)

target_include_directories(IQClassify
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/pocketfft)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(IQClassify PRIVATE -Wall -Wextra)
endif()

if(IQCLASSIFY_BUILD_TESTS)
    enable_testing()
    add_executable(iqclassify_tests
        tests/SignalGenerator.cpp
        tests/test_main.cpp)
    target_include_directories(iqclassify_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    target_link_libraries(iqclassify_tests PRIVATE IQClassify)
    add_test(NAME iqclassify_tests COMMAND iqclassify_tests)
endif()

# Generator fuer die README-Bilder (docs/images/*.png), aus dem Repo-Root:
#   ./build/iqclassify_readme_images
option(IQCLASSIFY_BUILD_TOOLS "README-Bildgenerator bauen" ON)
if(IQCLASSIFY_BUILD_TOOLS)
    add_executable(iqclassify_readme_images
        tools/readme_images.cpp
        tests/SignalGenerator.cpp)
    target_include_directories(iqclassify_readme_images PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/tools
        ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    target_link_libraries(iqclassify_readme_images PRIVATE IQClassify)
endif()

# Installations-Support (optional)
include(GNUInstallDirs)
install(TARGETS IQClassify ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/iqclassify DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
