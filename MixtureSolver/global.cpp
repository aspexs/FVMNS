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

void MacroParam::computeRho(const double& x_CO2)
{
    rho[0] = x_CO2 * Mixture::mass(0) * p / (K_BOLTZMANN * t);
    rho[1] = (1.0 - x_CO2) * Mixture::mass(1) * p / (K_BOLTZMANN * t);
}
