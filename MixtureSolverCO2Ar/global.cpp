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
double Mixture::vEnergy(const int& i, const int& j, const int& k)
{
    return W_W * (1345.04 * i + 667.25 * j + 2361.71 * k);
}
double Mixture::reducedMass(const int& i, const int& j)
{
    return mass(i) * mass(j) / (mass(i) + mass(j));
}

double Mixture::tauVTCO2CO2(const double& t, const double& rho_CO2)
{
    double x = qPow(t, -1.0 / 3.0);
    double p = rho_CO2 * K_BOLTZMANN * t / Mixture::mass(0) / 101325;
    return qExp(2.3025 * (-10.327 + 57.31 * x - 156.7 * qPow(x, 2.0))) / p;
}
double Mixture::tauVTCO2Ar(const double& t, const double& rho_Ar)
{
    double x = qPow(t, -1.0 / 3.0);
    double p = rho_Ar * K_BOLTZMANN * t / Mixture::mass(1) / 101325;
    return qExp(2.3025 * (-10.011 + 49.40 * x - 77.3 * qPow(x, 2.0))) / p;
}
double Mixture::tauVVCO2CO2(const double& t, const double& rho_CO2)
{
    double x = qPow(t, -1.0 / 3.0);
    double p = rho_CO2 * K_BOLTZMANN * t / Mixture::mass(0) / 101325;
    return qExp(2.3025 * (-12.662 + 88.87 * x - 272.5 * qPow(x, 2.0))) / p;
}
double Mixture::tauVVCO2Ar(const double& t, const double& rho_Ar)
{
    double x = qPow(t, -1.0 / 3.0);
    double p = rho_Ar * K_BOLTZMANN * t / Mixture::mass(1) / 101325;
    return qExp(2.3025 * (-11.511 + 77.67 * x - 209.7 * qPow(x, 2.0))) / p;
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
    t12 = 0.0;
    t3  = 0.0;
}

MacroParam::MacroParam(const double& p, const double& v, const double& t,
                       const double& x_CO2)
{
    this->rho = {0.0, 0.0};
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
