#include "additionalSolverForOxygen.h"
#include <global.h>
#include <QtMath>
extern double sigma ;//=  3.458e-10;
extern double epsilonDevK;// = 107.4 ;
extern double molMass ;//= 32.0e-3;
extern double mass ;//= 32.0e-3/Nav;
static double h_nu_dev_k = 2276;

additionalSolverForOxygen::additionalSolverForOxygen()
{
   kb_dev_mass = kB/mass;
}

double additionalSolverForOxygen::getC_V(double T)
{
    const double temp = h_nu_dev_k/T;
    return pow(temp,2)* kb_dev_mass * exp(temp)/pow(exp(temp)-1,2);
}

double additionalSolverForOxygen::getVibrEnergy(double T)
{
    const double h_nu = h_nu_dev_k * kB;
    return 1/mass * h_nu /(exp(h_nu_dev_k/T)-1);
}

double additionalSolverForOxygen::getEtta(double T)
{
    return (5*kB*T) /(8*getOmega22(T));
}

double additionalSolverForOxygen::getZetta(double T, double etta)
{
    return 4.4/3 * M_PI * etta/ getF(T);
}

double additionalSolverForOxygen::getLambdaVibr(double T, double Tv)
{
    return (3.0*kB*T)/(8.0*getOmega11(T))*getC_V(Tv);
}

double additionalSolverForOxygen::getLambdaTr_Rot(double T)
{
    return (75.0*pow(kB,2)*T)/(32.0*mass*getOmega22(T)) + (3.0*kB*T)/(8.0*getOmega11(T))*kb_dev_mass ;
}

double additionalSolverForOxygen::getTauVibr(double T, double P)
{
     P = P/101325;
     float L = 1.11 - 407.1*pow(T, -1/3) + 6600.9*pow(T, -2/3) - 31307.9/T;
     return pow(10, L)/P;
}

double additionalSolverForOxygen::getTempFromEnergy(double E)
{
    double ii = h_nu_dev_k * kB /(E*mass);
    double i = ii + 1;
    double l =log(i);
   return h_nu_dev_k / l;
}

double additionalSolverForOxygen::getOmega22(double T)
{
    QVector<double> f = {-0.40811, -0.05086, 0.34010, 0.70375, -0.10699, 0.00763};
    double a22= 1.5;
    double x = (log(T/epsilonDevK))+a22;
    double omegaLD = pow(f[0] + f[1]/pow(x,2) + f[2]/x + f[3]*x + f[4]*pow(x,2)+f[5]*pow(x,3),-1);
    double omegaS = sqrt(kB*T/(M_PI*mass))*2*M_PI*pow(sigma,2);
    return omegaLD*omegaS;
}

double additionalSolverForOxygen::getF(double T)
{
    return 1 + pow(M_PI,3.0/2)/2 * pow(T/epsilonDevK,-0.5) +
              (pow(M_PI,2)/4 +2) * pow(T/epsilonDevK,-1) +
               pow(M_PI,3.0/2)*    pow(T/epsilonDevK,-1.5);
}

double additionalSolverForOxygen::getOmega11(double T)
{
    QVector<double> f = {-0.16845, -0.02258, 0.19779, 0.64373, -0.09267, 0.00711};
    double a11 = 1.4;

    double x = (log(T/epsilonDevK))+a11;
    double omegaLD = pow(f[0] + f[1]/pow(x,2) + f[2]/x + f[3]*x + f[4]*pow(x,2)+f[5]*pow(x,3),-1);
    double omegaS =sqrt(kB*T/(M_PI*mass))*M_PI*pow(sigma,2);
    return omegaLD*omegaS;
}
