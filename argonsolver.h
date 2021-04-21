#ifndef ARGONSOLVER_H
#define ARGONSOLVER_H
#include "global.h"
#include "abstaractsolver.h"
#include "additionalsolver.h"

class ArgonSolver : public AbstaractSolver
{
    Q_OBJECT
public:
    ArgonSolver(QObject *parent = nullptr);
    void prepareSolving() override;
    void solveFlux(Matrix U1L, Matrix U2L, Matrix U3L, Matrix U1R, Matrix U2R, Matrix U3R);
    void calcRiemanPStar();
    void calcFliux();
public slots:
    void solve() override;
    void setTypePlote(int i);

};

#endif // ARGONSOLVER_H
