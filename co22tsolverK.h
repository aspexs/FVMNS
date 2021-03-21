#ifndef CO22TSOLVERК_H
#define CO22TSOLVERК_H
#include "global.h"
#include "abstaractsolver.h"
#include "additionalsolver.h"

class Co22TSolverK: public AbstaractSolver
{
public:
    Co22TSolverK(QObject *parent = nullptr);
    void prepareSolving() override;
    void solveFlux(const Matrix& U1, const Matrix& U2, const Matrix& U3, const Matrix& U4);
    void calcRiemanPStar();
    void calcFliux();
public slots:
    void solve() override;
    void setTypePlot(int i) override;

private:
    QList<double> CvibrMass;
    QVector<double> EnergyVibr;

    QList<double> EnergyTr_Rot;
    double energyStartTemp;
    double energyStepTemp;
    Matrix EnergyFullL, EnergyFullR;

    double CVibrStartTemp;
    double CVibrStepTemp;
    double dt;
     double leftEnergy;
    QVector<double> pres, R, T,Tv;
    Matrix F11, F22, F33,F44;
    void calcHLLE(const Matrix &U1, const Matrix &U2, const Matrix &U3, const Matrix &U4);
    void calcR(const Matrix &U1, const Matrix &U2, const Matrix &U3, const Matrix &U4);

    double getEnergyVibrTemp(double energy);
    double getVibrTemp(double CVibr);
};

#endif // CO22TSOLVERК_H
