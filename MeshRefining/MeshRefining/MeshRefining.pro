QT += core
QT -= gui

# CONFIG += c++11

TARGET = MeshRefining
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    dataclass.cpp \
    vector.cpp \
    dataio.cpp \
    refinefunctions.cpp \
    boundary.cpp \
    miscellaneous.cpp \
    openfoamfile.cpp \
    boundarycondition.cpp \
    reflection.cpp

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

HEADERS += \
    globaldefine.h \
    dataclass.h \
    vector.h \
    dataio.h \
    refinefunctions.h \
    boundary.h \
    miscellaneous.h \
    openfoamfile.h \
    boundarycondition.h \
    lookup_table.h \
    reflection.h

INCLUDEPATH += /usr/local/include/
LIBS += -L/usr/local/lib -L/usr/lib64/mpich-3.2/lib/ -lmpi
LIBS += -L/usr/local/lib/ -lmetis
INCLUDEPATH +=/usr/lib/mpich/include #$$PWD/include
INCLUDEPATH +=/usr/include/mpich-3.2-x86_64
unix:!macx: LIBS += -lmpi -lcgns -lhdf5 -ldl -fPIC -pg -g
INCLUDEPATH +=$$PWD/include
# LIBS += -L$$PWD/lib/ -lreflection
# LIBS += -L$$PWD/lib/ -lrefletion
LIBS += -L$$PWD/lib/ -lgeometry
