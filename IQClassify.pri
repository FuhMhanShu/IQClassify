# IQClassify - Einbindung in ein qmake-Konsumentenprojekt.
#
#   include(IQClassify.pri)
#
# Voraussetzung: die statische Bibliothek wurde gebaut (qmake
# IQClassify.pro && make) und liegt neben dieser Datei (.pri erwartet
# libIQClassify.a in $$PWD).

INCLUDEPATH += $$PWD/include
DEPENDPATH  += $$PWD/include

unix {
    LIBS += -L$$PWD -lIQClassify
}
win32 {
    LIBS += -L$$PWD -lIQClassify
}
