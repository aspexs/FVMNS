#ifndef ADDITIONALSOLVERFOROXYGEN
#define ADDITIONALSOLVERFOROXYGEN


class additionalSolverForOxygen
{
public:
    additionalSolverForOxygen();
    double getC_V(double T);
    double getVibrEnergy(double T);
    double getEtta(double T);
    double getZetta(double T, double etta);
    double getLambdaVibr(double T, double Tv);
    double getLambdaTr_Rot(double T);
    double getTauVibr(double T, double P);
    double getTempFromEnergy(double E);
private:
    double getOmega11(double T);
    double getOmega22(double T);
    double getF(double T);
    double kb_dev_mass;
};

#endif // ADDITIONALSOLVERFOROXYGEN
