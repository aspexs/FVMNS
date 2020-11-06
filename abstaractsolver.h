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
public slots:
    void pause();
    void breakSolver();
    virtual void setTypePlot(int i) = 0;
protected:
    Matrix U1, U2, U3, U4;
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
    QVector<macroParam> rezultAfterPStart;
    QMutex mutex;
signals:
   void updateGraph(QVector<double> x, QVector<double> y, double lambda = 1);
   void updateAdditionalGraph(QVector<double> x, QVector<double> y, double lambda = 1);
   void updateTime(double time);
};

#endif // ABSTARACTSOLVER_H