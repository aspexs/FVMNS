#include "additionalsolverforco2.h"
#include <global.h>
#include <QtMath>
extern double sigma;
extern double epsilonDevK;
extern double molMass;
extern double mass;

additionalSolverForCO2::additionalSolverForCO2()
{

}

double additionalSolverForCO2::bulcViscosity(double T, double C_vibr, double C_tr, double C_rot, double pressure)
{
    double CV = C_rot + C_tr + C_vibr;

    pressure = 101325;
    double n = pressure / (kB*T);

    double E_av_full = (/*C_vibr / TauVibrNew(T, n, C_vibr) +*/ C_rot / tau_rot_CO2(T, pressure))*massaCO2 / (kB*n);
    return kB*T * pow((C_vibr + C_rot) / CV, 2) / E_av_full;
}

double additionalSolverForCO2::TauVibrNew(double T, double n, double cvibr)
{
    double divtau = 2. * kB*n*E_av_VT(T)/(massaCO2*cvibr);
    return 1./divtau;
}

double additionalSolverForCO2::E_av_VT(double T)
{
    double E_aver;
    double sigma;
    double S = 0;
    int i1, i2, i3;
    for (i1 = 0; i1 < 35; i1++)
    {
        for (i2 = 0; i2 < 65; i2++)
        {
            for (i3 = 0; i3 < 20; i3++)
            {
                S += S_VT2_CO2_d(i1, i2, i3, T);
            }
        }
    }
    sigma = M_PI*sigmaCO2CO2_mix2*sigmaCO2CO2_mix2;
    E_aver = pow(2. * M_PI*kB*T /(masRed_CO2_CO2), 1. / 2.) * sigma * S;
    return E_aver;
}

double additionalSolverForCO2::S_VT2_CO2_d(int v1, int v2, int v3, double T)
{
    double delta_E = (E_CO2(v1, v2, v3) - E_CO2(v1, v2 - 1, v3));// *1.602176e-22;
    double dkT =  delta_E / (kB*T);
    double ekT =  (E_CO2(v1, v2-1, v3) - E_CO2(0, 0, 0)) / (kB*T);
    return pow(v2 + 1, 1.)* pow(dkT, 2) * exp(-ekT) * P_VT2_CO2_d(v1, v2, v3, T) / pow(Zvibr_CO2(T), 1.);
}

double additionalSolverForCO2::E_CO2(int v1, int v2, int v3)
{
    double d[3] = { 1., 2., 1. };
    double ECO2 = hc * (o2[0] * (v1 + d[0] / 2.0) + o2[1] * (v2 + d[1] / 2.0) + o2[2] * (v3 + d[2] / 2.0) + ox[0] * (v1 + d[0] / 2.0) * (v1 + d[0] / 2.0) + ox[1] * (v1 + d[0] / 2.0) * (v2 + d[1] / 2.0) + ox[2] * (v1 + d[0] / 2.0) * (v3 + d[2] / 2.0) + ox[3] * (v2 + d[1] / 2.0) * (v2 + d[1] / 2.0) + ox[4] * (v2 + d[1] / 2.0) * (v3 + d[2] / 2.0) + ox[5] * (v3 + d[2] / 2.0) * (v3 + d[2] / 2.0));
    return ECO2;
}

double additionalSolverForCO2::P_VT2_CO2_d(int v1, int v2, int v3, double T)
{
    double delta_E = (E_CO2(v1, v2, v3) - E_CO2(v1, v2 - 1, v3));// *1.602176e-22;
    double v2_const = (pow(a_SSH[1] * alpha1_CO2, 2))*hPlank *(v2+1) / (8.0 * M_PI*M_PI*m[1] * ni[1]); //CO2(v1,v2+1,v3)+CO2<-->CO2(v1,v2,v3)+CO2
    double pVT2_CO2 = fsigm(delta_E, masRed_CO2_CO2, alpha1_CO2, T) * v2_const;
    if (pVT2_CO2 >= 1) pVT2_CO2 = 0.9;
    return pVT2_CO2;
}

double additionalSolverForCO2::fsigm(double delta_E, double masRed, double alpha, double T)
{
    double sigma1 = 0;
    if (delta_E > 0.0)
        sigma1 = 3. * pow((2. * pow(M_PI, 4) * pow(delta_E, 2)*masRed) / (pow(alpha, 2)*pow(hPlank, 2)*kB*T), 1. / 3.) + delta_E / (2. * kB*T);
    else
        sigma1 = 3. * pow((2. * pow(M_PI, 4) * pow(delta_E, 2)*masRed) / (pow(alpha, 2)*pow(hPlank, 2)*kB*T), 1. / 3.) - delta_E / (2. * kB*T);
    double sigma2 = pow(sigma1, 3. / 2.) * exp(-sigma1) / (1. - exp(-2.*sigma1 / 3.));
    double fsigma = sigma2 * 0.394 * pow((8. * pow(M_PI, 3) * delta_E * masRed) / (pow(alpha, 2)*pow(hPlank, 2)), 2);
    return fsigma;
}

double additionalSolverForCO2::Zvibr_CO2(double T)
{
    int i1, i2, i3;
    double S = 0;
    for (i1 = 0; i1 < 30; i1++)
    {
        for (i2 = 0; i2 < 65; i2++)
        {
            for (i3 = 0; i3 < 20; i3++)
            {
                S += (i2 + 1)*exp(-i1*e100 / kB / T)*exp(-i2*e010 / kB / T)*exp(-i3*e001 / kB / T);
            }
        }
    }
    return S;
}

double additionalSolverForCO2::tau_rot_CO2(double T, double P)
{
    double t_rot;
    double eta = 5. * kB*T / (8. * Omega_22_CO2_CO2(T) * 2. * Omega_rs_11_CO2_CO2(T));
    double p1 = P; //в паскалях или в атмосферах???
    t_rot = (ksi_rot_CO2(T)*eta*M_PI) / (4. * p1);
    return t_rot;
}

double additionalSolverForCO2::Omega_22_CO2_CO2(double T)
{
    double b = 7.75;
    double e0 = 19.911*1.602e-22;
    double a1 = (7.898524e-1) - (2.114115e-2)*b;
    double a2 = (-2.998325e-1) - (1.243977e-3)*b;
    double a3 = (7.077103e-1) + (3.583907e-2)*b;
    double a4 = (-8.946857e-1) - (2.473947e-2)*b;
    double a5 = (-2.958969) + (2.303358e-1)*b - (5.226562e-3)*b*b;
    double a6 = (4.348412) + (1.920321e-1)*b - (1.496557e-2)*b*b;
    double a7 = (2.205440) + (2.567027e-1)*b - (1.861359e-2)*b*b;
    double x = log(kB*T / e0);
    double lnOmega_22 = ((a1 + a2*x)*exp((x - a3) / a4)) / (exp((x - a3) / a4) + exp((a3 - x) / a4)) + (a5* exp((x - a6) / a7)) / (exp((x - a6) / a7) + exp((a6 - x) / a7));
    double Omega_22 = exp(lnOmega_22);
    return Omega_22;
}

double additionalSolverForCO2::Omega_11_CO2_CO(double T)
{
    double b = 7.87;
    double e0 = 15.418*1.602e-22;
    double a1 = (7.884756e-1) - (2.438494e-2)*b;
    double a2 = (-2.952759e-1) - (1.744149e-3)*b;
    double a3 = (5.020892e-1) + (4.316985e-2)*b;
    double a4 = (-9.042460e-1) - (4.017103e-2)*b;
    double a5 = (-3.373058) + (2.458538e-1)*b - (4.850047e-3)*b*b;
    double a6 = (4.161981) + (2.202737e-1)*b - (1.718010e-2)*b*b;
    double a7 = (2.462523) + (3.231308e-1)*b - (2.281072e-2)*b*b;
    double x = log(kB*T / e0);
    double lnOmega_11 = ((a1 + a2*x)*exp((x - a3) / a4)) / (exp((x - a3) / a4) + exp((a3 - x) / a4)) + (a5* exp((x - a6) / a7)) / (exp((x - a6) / a7) + exp((a6 - x) / a7));
    double Omega_11 = exp(lnOmega_11);

    return Omega_11;
}

double additionalSolverForCO2::Omega_rs_11_CO2_CO2(double T)
{
    double re = 4.119e-10;
    double xi_1 = 0.8002;
    double xi_2 = 0.049256;
    double b = 7.75;
    double x0 = xi_1* pow(b, xi_2);
    double s_sqr = x0*x0*re*re;
    double m_NN = (massaCO2*massaCO2) / (massaCO2 + massaCO2);
    double f = (kB*T) / (2 * 3.141592*m_NN);
    double Omega_rs = (pow(f, 0.5))*(3.141592*s_sqr);
    return Omega_rs;
}

double additionalSolverForCO2::ksi_rot_CO2(double T)
{
    double ksi_rot;
    double Fs;
    double eta = 5. * kB*T / (8. * Omega_22_CO2_CO2(T) * 2. * Omega_rs_11_CO2_CO2(T));
    double ksi_inf = 20.39;
    double kTe = T / 97.5; // e/k=97.5 for CO2
    Fs = 1. + (pow((M_PI), (1.5)) / 2.)*(pow(kTe, -0.5)) + (M_PI*M_PI / 4. + 2.)*(pow(kTe, -1)) + (pow((M_PI), (1.5)))*(pow(kTe, -1.5));
    ksi_rot = ksi_inf / Fs;
    return ksi_rot;
}
