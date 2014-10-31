TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -llapack -lblas -larmadillo
SOURCES += main.cpp \
    plot.cpp \
    schemes.cpp

HEADERS += \
    plot.h \
    schemes.h
