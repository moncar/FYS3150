TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -llapack -lblas -larmadillo
SOURCES += main.cpp \
    plot.cpp \
    schemes.cpp \
    schemes2d.cpp \
    markov_chain.cpp

HEADERS += \
    plot.h \
    schemes.h \
    schemes2d.h \
    markov_chain.h
