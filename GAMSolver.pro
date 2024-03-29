#-------------------------------------------------
#
# Project created by QtCreator 2019-10-23T13:59:55
#
#-------------------------------------------------

QT       += core gui concurrent charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = GAMSolver
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11

SOURCES += \
        abstaractsolver.cpp \
        additionalsolver.cpp \
        additionalsolverforco2.cpp \
        additionalsolverforoxygen.cpp \
        argonsolver.cpp \
        chart.cpp \
        co22tsolver.cpp \
        co22tsolverK.cpp \
        co23tsolver.cpp \
        co2solver.cpp \
        global.cpp \
        main.cpp \
        nitrogensolver.cpp \
        oxygensolver.cpp \
        widget.cpp

HEADERS += \
        abstaractsolver.h \
        additionalSolverForOxygen.h \
        additionalsolver.h \
        additionalsolverforco2.h \
        argonsolver.h \
        chart.h \
        co22tsolver.h \
        co22tsolverK.h \
        co23tsolver.h \
        co2solver.h \
        global.h \
        nitrogensolver.h \
        oxygensolver.h \
        widget.h

FORMS += \
        widget.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
