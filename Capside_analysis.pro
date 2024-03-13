TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += xdrfile-1.1.4/include/

SOURCES += main.cpp \
    xdrfile-1.1.4/src/xdrfile.c \
    xdrfile-1.1.4/src/xdrfile_c_test.c \
    xdrfile-1.1.4/src/xdrfile_trr.c \
    xdrfile-1.1.4/src/xdrfile_xtc.c

HEADERS += \
    diagram.h \
    dodecahedron.h \
    histogram.h \
    pentamermap.h \
    lammpscapsidparam.h \
    xdrfile-1.1.4/include/xdrfile.h \
    xdrfile-1.1.4/include/xdrfile_trr.h \
    xdrfile-1.1.4/include/xdrfile_xtc.h \
    xdrfile-1.1.4/config.h \
    xdrfile-1.1.4/include/xdrfile.h \
    diagram.h \
    dodecahedron.h \
    lammpscapsidparam.h \
    pentamermap.h \
    xtcanalysis.h \
    welford.h

