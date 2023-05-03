#include "global.h"

///////////////////////////////////////////////////////////////////////////////
/// class MixtureData
///////////////////////////////////////////////////////////////////////////////

double Mixture::mass(const int& i)
{
    switch (i)
    {
        case 0:
            return 31.9980 * ATOMIC_MASS_UNIT;
        break;
        case 1:
            return 15.9994 * ATOMIC_MASS_UNIT;
        break;
        default:
            return 0;
    }
}
double Mixture::sigma(const int& i)
{
    switch (i)
    {
        case 0:
            return 3.458e-10;
        break;
        case 1:
            return 2.750e-10;
        break;
        default:
            return 0;
    }
}
double Mixture::epsilon(const int& i)
{
    switch (i)
    {
        case 0:
            return 107.4;
        break;
        case 1:
            return 80.0;
        break;
        default:
            return 0;
    }
}
double Mixture::vEnergy(const int& i)
{
    return W_W * 1580.9 * i;
}
double Mixture::reducedMass(const int& i, const int& j)
{
    return mass(i) * mass(j) / (mass(i) + mass(j));
}

double Mixture::tauVTO2O2(const double& t, const double& rho_O2)
{
    double x = qPow(t, -1.0 / 3.0);
    double p = rho_O2 * K_BOLTZMANN * t / Mixture::mass(0) / 101325;
    return qExp(2.2036 * (61.03 * x - 10.19)) / p;
}
double Mixture::tauVTO2O(const double& t, const double& rho_O)
{
    double x = qPow(t, -1.0 / 3.0);
    double p = rho_O * K_BOLTZMANN * t / Mixture::mass(1) / 101325;
    return qExp(2.2036 * (-8.805 * x - 6.647)) / p;
}

///////////////////////////////////////////////////////////////////////////////
/// class MacroParam
///////////////////////////////////////////////////////////////////////////////

MacroParam::MacroParam()
{
    rho = {0.0, 0.0};
    p   = 0.0;
    v   = 0.0;
    t   = 0.0;
    tv  = 0.0;
}
MacroParam::MacroParam(const double& p, const double& v, const double& t,
                       const double& x_O2)
{
    this->rho = {0.0, 0.0};
    initialize(p, v, t, x_O2);
}
void MacroParam::computeRho(const double& x_O2)
{
    rho[0] = x_O2 * Mixture::mass(0) * p / (K_BOLTZMANN * t);
    rho[1] = (1.0 - x_O2) * Mixture::mass(1) * p / (K_BOLTZMANN * t);
}
void MacroParam::initialize(const double& p, const double& v, const double& t,
                            const double& x_O2)
{
    this->v   = v;
    this->t   = t;
    this->tv  = t;
    this->p   = p;
    this->computeRho(x_O2);
}

///////////////////////////////////////////////////////////////////////////////
/// class ProgressBar
///////////////////////////////////////////////////////////////////////////////

ProgressBar::ProgressBar()
{
    dx = 0.0;
    n0 = 0.0;
    n1 = 0.0;
}
void ProgressBar::initialize(const double& t)
{
    dx = 100.0 / t;
    for (int i = 0; i < 100; ++i)
    {
        std::cout << '/';
    }
    std::cout << " -- [100%]\n";
}
void ProgressBar::update(const double& dt)
{
    n1 += dt * dx;
    if (qFloor(n1) > qFloor(n0))
    {
        n0 = n1;
        std::cout << '/';
    }
}
