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
double Mixture::o2(const int& i)
{
    switch (i)
    {
        case 0:
            return 1345.04e2;
        break;
        case 1:
            return 667.25e2;
        break;
        case 2:
            return 2361.71e2;
        break;
        default:
            return 0;
    }
}
double Mixture::ox(const int& i)
{
    switch (i)
    {
        case 0:
            return -3.63e2;
        break;
        case 1:
            return 3.44e2;
        break;
        case 2:
            return -19.28e2;
        break;
        case 3:
            return -0.635e2;
        break;
        case 4:
            return -12.51e2;
        break;
        case 5:
            return -12.56e2;
        break;
        case 6:
            return 0.775e2;
        break;
        default:
            return 0;
    }
}
double Mixture::m(const int& i)
{
    switch (i)
    {
        case 0:
            return 1.46e-26;
        break;
        case 1:
            return 1.338e-26;
        break;
        case 2:
            return 1.46e-26;
        break;
        default:
            return 0;
    }
}
double Mixture::ni(const int& i)
{
    switch (i)
    {
        case 0:
            return o2(0) * C_LIGHT;
        break;
        case 1:
            return o2(1) * C_LIGHT;
        break;
        case 2:
            return o2(2) * C_LIGHT;
        break;
        default:
            return 0;
    }
}
double Mixture::a_SSH(const int& i)
{
    switch (i)
    {
        case 0:
            return 0.5;
        break;
        case 1:
            return 8.0 / 11.0;
        break;
        case 2:
            return 0.5;
        break;
        default:
            return 0;
    }
}
double Mixture::a1_SSH(const int& i)
{
    switch (i)
    {
        case 0:
            return 8.0 / 11.0;
        break;
        case 1:
            return 0.5;
        break;
        case 2:
            return 3.0 / 11.0;
        break;
        default:
            return 0;
    }
}

double Mixture::tauVTCO2CO2(const double& t, const double& rho_CO2)
{
    double a = -18.19;
    double b = 40.47;
    double c = 0.0;
    double d = 0.00423;
    double p = rho_CO2 * K_BOLTZMANN * t / Mixture::mass(0) / 101325;
    return qExp(a + b * qPow(t, -1.0 / 3.0) + c * qPow(t, -2.0 / 3.0) +
                d * qPow(t, 1.0 / 3.0)) / p;
}
double Mixture::tauVTCO2Ar(const double& t, const double& rho_Ar)
{
    double a = -18.19;
    double b = 40.47;
    double c = 0.0;
    double d = 0.00423;
    double p = rho_Ar * K_BOLTZMANN * t / Mixture::mass(1) / 101325;
    return qExp(a + b * qPow(t, -1.0 / 3.0) + c * qPow(t, -2.0 / 3.0) +
                d * qPow(t, 1.0 / 3.0)) / p;
}
double Mixture::tauVVCO2CO2(const double& t, const double& rho_CO2)
{
    double a = -26.85;
    double b = 173.22;
    double c = -539.74;
    double d = 0.09645;
    double t13 = pow(t, -1.0 / 3.0);
    double p = rho_CO2 * K_BOLTZMANN * t / Mixture::mass(0) / 101325;
    return qExp(a + b * t13 + c * qPow(t13, 2.0) + d / t13) / p;
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
