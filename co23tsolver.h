#ifndef CO23TSOLVER_H
#define CO23TSOLVER_H

#include "abstaractsolver.h"

class Co23TSolver : public AbstaractSolver
{
public:
    explicit Co23TSolver(QObject *parent = nullptr);
     void prepareSolving() override;
     void solveFlux();
     void calcFliux();
public slots:
    void solve() override;
private:
    double dt, error;
    void calcR(const Matrix &U1, const Matrix &U2, const Matrix &U3, const Matrix &U4);
    double getEnergyVibr12Temp(double E);
    double getEnergyVibr3Temp(double E);
    void calcR(const Matrix &U1, const Matrix &U2, const Matrix &U3, const Matrix &U4, const Matrix &U5);
};

#endif // CO23TSOLVER_H
