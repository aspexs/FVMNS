#ifndef CO2SOLVER_H
#define CO2SOLVER_H

#include "global.h"
#include "abstaractsolver.h"
#include "additionalsolver.h"

class CO2Solver : public AbstaractSolver
{
     Q_OBJECT
public:
    CO2Solver(QObject *parent = nullptr);
    void prepareSolving() override;
    void solveFlux(Matrix U1L, Matrix U2L, Matrix pressureL, Matrix U1R, Matrix U2R, Matrix pressureR);
    void calcRiemanPStar();
    void calcFliux();

public slots:
    void solve() override;
    void setTypePlot(int i) override;
private:
    double getEnergyTemp(double energy);
    double getVibrTemp(double CVibr);
};

#endif // CO2SOLVER_H
