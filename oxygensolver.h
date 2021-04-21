#ifndef OXYGENSOLVER_H
#define OXYGENSOLVER_H

#include <QObject>
#include "abstaractsolver.h"
#include "additionalSolverForOxygen.h"

class OxygenSolver : public AbstaractSolver
{
    Q_OBJECT
public:
    explicit OxygenSolver(QObject *parent = nullptr);

     void prepareSolving() override;
public slots:
    void solve() override;
private:
    additionalSolverForOxygen as;
    void solveFlux(Matrix U1L, Matrix U2L, Matrix pressureL,Matrix TvL,
                   Matrix U1R, Matrix U2R, Matrix pressureR,Matrix TvR);
    void calcRiemanPStar();
    void calcFliux();
    double dt;
signals:

public slots:

};

#endif // OXYGENSOLVER_H
