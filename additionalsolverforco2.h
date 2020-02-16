#ifndef ADDITIONALSOLVERFORCO2_H
#define ADDITIONALSOLVERFORCO2_H


class additionalSolverForCO2
{
public:
    additionalSolverForCO2();

    static double bulcViscosity(double T, double C_vibr, double C_tr, double C_rot, double pressure);
    static double TauVibrNew(double T, double n, double cvibr);
    static double E_av_VT(double T);
    static double S_VT2_CO2_d(int v1, int v2, int v3, double T);
    static double E_CO2(int v1, int v2, int v3);
    static double P_VT2_CO2_d(int v1, int v2, int v3, double T);
    static double fsigm(double delta_E, double masRed, double alpha, double T);
    static double Zvibr_CO2(double T);

    static double  tau_rot_CO2(double T, double P);
    static double  Omega_22_CO2_CO2(double T);
    static double  Omega_11_CO2_CO(double T);
    static double  Omega_rs_11_CO2_CO2(double T); //Bruno et al
    static double ksi_rot_CO2(double T);
};

#endif // ADDITIONALSOLVERFORCO2_H
