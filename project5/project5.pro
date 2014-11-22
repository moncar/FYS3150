TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += cpp11
LIBS += -llapack -lblas -larmadillo
SOURCES += main.cpp \
    plot.cpp \
    schemes.cpp \
    schemes2d.cpp \
    markov_chain.cpp \
    /home/filiphl/Desktop/FYS3150/project5/cppLibrary/lib.cpp \
    random.cpp

HEADERS += \
    plot.h \
    schemes.h \
    schemes2d.h \
    markov_chain.h \
    random.h
