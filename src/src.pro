TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
        particle_filter.cpp

HEADERS += \
    helper_functions.h \
    json.hpp \
    map.h \
    particle_filter.h
