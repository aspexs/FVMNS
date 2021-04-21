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
public slots:
    void solve() override;

private:
    Matrix EnergyFullL, EnergyFullR;
    double dt;
     double leftEnergy;
    void calcHLLE(const Matrix &U1, const Matrix &U2, const Matrix &U3, const Matrix &U4);
    void calcR(const Matrix &U1, const Matrix &U2, const Matrix &U3, const Matrix &U4);

    double getEnergyVibrTemp(double energy);
    double getVibrTemp(double CVibr);
};

#endif // CO22TSOLVERК_H
