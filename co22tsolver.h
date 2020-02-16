#ifndef CO22TSOLVER_H
#define CO22TSOLVER_H
#include "global.h"
#include "abstaractsolver.h"
#include "additionalsolver.h"


class Co22TSolver: public AbstaractSolver
{
public:
    Co22TSolver(QObject *parent = nullptr);
    void prepareSolving() override;
    void solveFlux(Matrix U1L, Matrix U2L, Matrix pressureL,Matrix TvL, Matrix U1R, Matrix U2R, Matrix pressureR, Matrix TvR);
    void calcRiemanPStar();
    void calcFliux();
public slots:
    void solve() override;
    void setTypePlot(int i) override;

private:
    QList<double> CvibrMass;
    QList<double> EnergyVibr;
    QList<double> EnergyTr_Rot;
    double energyStartTemp;
    double energyStepTemp;

    double CVibrStartTemp;
    double CVibrStepTemp;

    double getEnergyVibrTemp(double energy);

    double getVibrTemp(double CVibr);
};

#endif // CO22TSOLVER_H
