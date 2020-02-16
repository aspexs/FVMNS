#include "abstaractsolver.h"
AbstaractSolver::AbstaractSolver(QObject *parent) : QThread(parent)
{

}

void AbstaractSolver::pause()
{
    pauseSolve = !pauseSolve;
}

void AbstaractSolver::breakSolver()
{
    breaksolve = true;
}
