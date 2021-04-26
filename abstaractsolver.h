#ifndef ABSTARACTSOLVER_H
#define ABSTARACTSOLVER_H

#include <QObject>
#include <QThread>
#include <QMutexLocker>
#include <QProgressDialog>
#include <QFutureWatcher>
#include <QtConcurrent>
#include "global.h"
#include "additionalsolver.h"
#include "additionalsolverforco2.h"
class AbstaractSolver : public QThread
{
    Q_OBJECT
public:
    AbstaractSolver(QObject *parent = nullptr);

    virtual void solve() = 0;
    virtual void prepareSolving()   = 0;
    void run() override {solve();}
    bool breaksolve = false;
    bool pauseSolve = false;
    macroParam leftParam;
    macroParam rightParam;
    solverParams solParam;
    int typeBC = 0;
    int typeTimeStepSolution = 0;
    int typeBulkVisc = 0;
    int typeShareVisc = 0;
    int typeEnergy = 0;
    AdditionalSolver additionalSolver;
    additionalSolverForCO2 additionalSolverCo2;

    QVector<double> R, P, Q_v, Q_t, T, Tv, Ent, Ent2;
    Matrix U1, U2, U3, U4, pres;
public slots:
    void pause();
    void breakSolver();
    void setTypePlot(int i);
protected:
    Matrix F1, F2, F3, F4;
    QVector <double> x;
    QVector<int> vectorForParallelSolving;
    double delta_h;
    Matrix timeSolvind;

    Matrix left_density;
    Matrix right_density;
    Matrix left_velocity;
    Matrix right_velocity;
    Matrix left_pressure;
    Matrix right_pressure;
    Matrix left_Tv;
    Matrix right_Tv;
    Matrix Tl,Tr;
    QVector<macroParam> rezultAfterPStart;
    QMutex mutex;
    QList<double> CvibrMass;
    double CVibrStartTemp;
    double CVibrStepTemp;
    QVector<double> EnergyVibr;
    double energyVibrStartTemp;
    double energyVibrStepTemp;

    QVector<double> Energy;
    double energyStartTemp;
    double energyStepTemp;


    void prepareVectors();

signals:
   void updateGraph(QVector<double> x, QVector<double> y, double lambda = 1);
   void updateAdditionalGraph(QVector<double> x, QVector<double> y, double lambda = 1);
   void updateAdditional2Graph(QVector<double> x, QVector<double> y, double lambda = 1);
   void updateTime(double time, double error = 0);
};

#endif // ABSTARACTSOLVER_H
