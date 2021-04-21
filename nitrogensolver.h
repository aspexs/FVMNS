#ifndef NITROGENSOLVER_H
#define NITROGENSOLVER_H
#include "global.h"
#include "abstaractsolver.h"
#include "additionalsolver.h"

class NitrogenSolver : public AbstaractSolver
{
    Q_OBJECT
public:
    NitrogenSolver(QObject *parent = nullptr);
    void prepareSolving() override;
    void solveFlux(Matrix U1L, Matrix U2L, Matrix U3L, Matrix U1R, Matrix U2R, Matrix U3R);
    void calcRiemanPStar();
    void calcFliux();

public slots:
    void solve() override;
};

#endif // NITROGENSOLVER_H
