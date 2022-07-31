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
    virtual void prepareSolving();
    void run() override {solve();}
    virtual void calcRiemanPStar();
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

    Matrix R, P, Q_v, Q_t, T, Tv, Ent, Ent2, R_1, R_2, T12, T3, Q_v3, B_v, E_Z, PR;
    Matrix U1, U2, U3, U4,U5, pres;
public slots:
    void pause();
    void breakSolver();
    void setTypePlot(int i);
protected:
    Matrix F1, F2, F3, F4, F5;
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
    Matrix left_Tv, right_Tv;
    Matrix Tl,Tr;
    Matrix T12L,T12R, T3L, T3R;
    QVector<macroParam> rezultAfterPStart;
    QMutex mutex;
    QList<double> CvibrMass;
    double CVibrStartTemp;
    double CVibrStepTemp;
    QVector<double> EnergyVibr;
    double energyVibrStartTemp;
    double energyVibrStepTemp;
    QVector<double> EnergyVibr12;
    double energyVibrStartTemp12;
    double energyVibrStepTemp12;
    QVector<double> EnergyVibr3;
    double energyVibrStartTemp3;
    double energyVibrStepTemp3;

    QVector<double> Energy;
    double energyStartTemp;
    double energyStepTemp;


    void prepareVectors();

signals:
   void updateGraph(QVector<double> x, QVector<double> y, double lambda = 1);
   void updateAdditionalGraph(QVector<double> x, QVector<double> y, double lambda = 1);
   void updateAdditional2Graph(QVector<double> x, QVector<double> y, double lambda = 1);
   void updateTime(double time, double error = 0);
   void checkCalc();
};

#endif // ABSTARACTSOLVER_H
