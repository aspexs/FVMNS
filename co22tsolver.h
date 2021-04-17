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
    void solveFlux(Matrix U1L, Matrix U2L, Matrix pressureL, Matrix TvL, Matrix Tl, Matrix U1R, Matrix U2R, Matrix pressureR, Matrix TvR, Matrix Tr, Matrix EnergyFull);
    void calcRiemanPStar();
    void calcFliux();
public slots:
    void solve() override;
    void setTypePlot(int i) override;

private:
    double dt, error;
    QVector<double> pres, R, P, Q_v, Q_t;
    Matrix T, Tv, Tl,Tr,E__;
    void calcR(const Matrix &U1, const Matrix &U2, const Matrix &U3, const Matrix &U4);

    double getEnergyVibrTemp(double energy);
    double getVibrTemp(double CVibr);
};

#endif // CO22TSOLVER_H
