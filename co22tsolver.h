#ifndef CO22TSOLVER_H
#define CO22TSOLVER_H
#include "global.h"
#include "abstaractsolver.h"
#include "additionalsolver.h"


class Co22TSolver: public AbstaractSolver
{
public:
    Co22TSolver(QObject *parent = nullptr);
    void solveFlux(Matrix U1L, Matrix U2L, Matrix pressureL, Matrix TvL, Matrix Tl, Matrix U1R, Matrix U2R, Matrix pressureR, Matrix TvR, Matrix Tr);
    void calcFliux();
public slots:
    void solve() override;
private:
    double dt, error;
    void calcR(const Matrix &U1, const Matrix &U2, const Matrix &U3, const Matrix &U4);

    double getEnergyVibrTemp(double energy);
    QVector<double> errors;
};

#endif // CO22TSOLVER_H
