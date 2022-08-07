#ifndef MIXTURECO2AR_H
#define MIXTURECO2AR_H

#include "abstaractsolver.h"

class MixtureCo2Ar : public AbstaractSolver
{
public:
    explicit MixtureCo2Ar(QObject *parent = nullptr);
    void prepareSolving() override;

    void solveFlux();
    void calcFliux();

public slots:
    void solve() override;
private:
    double dt, error;
    Matrix localTemp;
    void calcR(const Matrix &U1, const Matrix &U2, const Matrix &U3, const Matrix &U4,const Matrix &U5, const Matrix &U6 );

};

#endif // MIXTURECO2AR_H
