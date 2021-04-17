#ifndef CO2Solver2T2_H
#define CO2Solver2T2_H

#include "global.h"
#include "abstaractsolver.h"
#include "additionalsolver.h"

class CO2Solver2T : public AbstaractSolver
{
     Q_OBJECT
public:
    CO2Solver2T(QObject *parent = nullptr);
    void prepareSolving() override;
    void solveFlux(Matrix U1L, Matrix U2L, Matrix pressureL, Matrix U1R, Matrix U2R, Matrix pressureR, Matrix TvR, Matrix TvL);
    void calcRiemanPStar();
    void calcFliux();

public slots:
    void solve() override;
    void setTypePlot(int i) override;
private:
     Matrix T, Tv;
       QVector<double> EnergyVibr, pres;
    QList<double> CvibrMass;
    //QList<double> Energy;
    double energyVibrStartTemp;
    double energyVibrStepTemp;

    double CVibrStartTemp;
    double CVibrStepTemp;

    //double getEnergyTemp(double energy);
    double getVibrTemp(double CVibr);
      double getEnergyVibrTemp(double energy);
};

#endif // CO2Solver2T2_H
