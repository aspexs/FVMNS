#include "global.h"

///////////////////////////////////////////////////////////////////////////////
/// class MixtureData
///////////////////////////////////////////////////////////////////////////////

double Mixture::mass(const int& i)
{
    switch (i)
    {
        case 0:
            return 44.01 * ATOMIC_MASS_UNIT;
        break;
        case 1:
            return 39.948 * ATOMIC_MASS_UNIT;
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
            return 3.763e-10;
        break;
        case 1:
            return 3.35e-10;
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
            return 244.0;
        break;
        case 1:
            return 141.5;
        break;
        default:
            return 0;
    }
}
double Mixture::vEnergy(const int& i, const int& j,
                                  const int& k)
{
    return W_W * (1345.04 * i + 667.25 * j + 2361.71 * k);
}
double Mixture::reducedMass(const int& i, const int& j)
{
    return mass(i) * mass(j) / (mass(i) + mass(j));
}
double Mixture::tauVTCO2(const double& t, const double& p)
{
    double a = -18.19;
    double b = 40.47;
    double c = 0.0;
    double d = 0.00423;
    return qExp(a + b * qPow(t, -1.0 / 3.0) + c * qPow(t, -2.0 / 3.0) +
                d * qPow(t, 1.0 / 3.0)) / p * 101325;
}
double Mixture::tauVVCO2(const double& t, const double& p)
{
    double a = -26.85;
    double b = 173.22;
    double c = -539.74;
    double d = 0.09645;
    double t13 = pow(t, -1.0 / 3.0);
    return qExp(a + b * t13 + c * qPow(t13, 2.0) + d / t13) / p * 101325;
}

///////////////////////////////////////////////////////////////////////////////
/// class MacroParam
///////////////////////////////////////////////////////////////////////////////

MacroParam::MacroParam()
{
    rho = {0.0, 0.0};
    p    = 0.0;
    v    = 0.0;
    t    = 0.0;
    t12  = 0.0;
    t3   = 0.0;
}

MacroParam::MacroParam(const double& p, const double& v, const double& t,
                       const double& x_CO2)
{
    rho = {0.0, 0.0};
    initialize(p, v, t, x_CO2);
}

void MacroParam::computeRho(const double& x_CO2)
{
    rho[0] = x_CO2 * Mixture::mass(0) * p / (K_BOLTZMANN * t);
    rho[1] = (1.0 - x_CO2) * Mixture::mass(1) * p / (K_BOLTZMANN * t);
}

void MacroParam::initialize(const double& p, const double& v, const double& t,
                            const double& x_CO2)
{
    this->v   = v;
    this->t   = t;
    this->t12 = t;
    this->t3  = t;
    this->p   = p;
    this->computeRho(x_CO2);
}

MacroParam MacroParam::proceed(const MacroParam& p0)
{
    MacroParam p2;
    p2.rho[0] = 2.0 * rho[0] - p0.rho[0];
    p2.rho[1] = 2.0 * rho[1] - p0.rho[1];
    p2.p = 2.0 * p - p0.p;
    p2.v = 2.0 * v - p0.v;
    p2.t = 2.0 * t - p0.t;
    p2.t12 = 2.0 * t12 - p0.t12;
    p2.t3 = 2.0 * t3 - p0.t3;
    return p2;
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

void ProgressBar::initialize(const double& n)
{
    dx = 80.0 / n;
    for (int i = 0; i < 80; ++i)
    {
        std::cout << '/';
    }
    std::cout << " -- [100%]\n";
}

void ProgressBar::update()
{
    n1 += dx;
    if (qFloor(n1) > qFloor(n0))
    {
        n0 = n1;
        std::cout << '/';
    }
}
