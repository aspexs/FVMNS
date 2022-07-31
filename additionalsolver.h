#ifndef ABSTRACTADDITIONALSOLVER_H
#define ABSTRACTADDITIONALSOLVER_H

#include <QObject>
#include <QVector>
#include <QMutex>
#include <global.h>

class AdditionalSolver : public QObject
{
    Q_OBJECT
public:
    enum Function
    {
        SHARE_VISC_SUPER_SIMPLE = 0,
        SHARE_VISC_SIMPLE,
        SHARE_VISC_OMEGA,
        BULC_VISC_SIMPLE,
        BULC_VISC_OLD,
        BULC_VISC_NEW,
        VIBR_ENERGY,
        ALL_ENERGY,
        C_VIBR,
        LAMBDA_TR,
        LAMBDA_VIBR,
        Z_VIBR,
        EVIBR12,
        EVIBR3,
        LAMBDA12,
        LAMBDA3,
        C_TR,
        C_ROT,
        LAMBDA_N2,
        FUNCTION_COUNT
    };
    enum BoundaryConditions
    {
        BC_RG,
        BC_RG2T,
        BC_PYTHON,
        BC_COUNT
    };
    enum TimeStepSolver
    {
        TS_SIMPLE,
        TS_FULL,
        TS_COUNT
    };
    enum ShareViscosity
    {
        SV_SS,
        SV_S,
        SV_O,
        SV_COUNT
    };

    enum BulkViscosity
    {
        BV_S,
        BV_OLD,
        BV_WITHOUT,
        BV_ONLY_RT_ROT,
        BV_COUNT
    };

    explicit AdditionalSolver(QObject *parent = nullptr);

    void solve();
    double startValue;
    double stopValue;
    double step;
    double pressure = 0;
    double density = 0;
    double gamma = 0;
    Function typeSolve;
    QVector <double> iterationVector;
    QVector <double> rezultVector;
    macroParam (*BoundaryCondition[BC_COUNT])(macroParam left, solverParams solParams) = {bondaryConditionRG2T,bondaryConditionPython };
    double (*TimeStepSolution[TS_COUNT])(Matrix velosity, Matrix density, Matrix pressure, double delta_h, solverParams solParams, Matrix energy) = {getTimeStep, getTimeStepFull};
    double (*shareViscosity[SV_COUNT])(double t1, double t2, double t3, double t4) = {shareViscositySuperSimple, shareViscositySimple, shareViscosityOmega};
    double (*bulkViscosity[BV_COUNT])(double t1, double t2, double t3, double t4) = {bulcViscositySimple, bulcViscosityOld2,bulcViscosityOnlyTRRot,bulcViscosityFalse, };

    macroParam ExacRiemanSolver(macroParam left, macroParam right, double Gamma);
    macroParam ExacRiemanSolverCorrect(macroParam left, macroParam right, double GammaL, double GammaR,double lambda = 0);

    QVector<QVector<double>> SolveEvolutionExplFirstOrder(Matrix F1,   Matrix F2,    Matrix F3, Matrix F4, Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old,
                                                          double dt, double delta_h);
    QVector<QVector<double>> SolveEvolutionExplFirstOrderForCO22(Matrix F1,   Matrix F2,    Matrix F3, Matrix F4, Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old,
                                                                 double dt, double delta_h, Matrix R = Matrix());
    QVector<QVector<double>> SolveEvolutionExplFirstOrderForO2(Matrix F1,   Matrix F2,    Matrix F3, Matrix F4, Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old,
                                                                 double dt, double delta_h);
    QVector<QVector<double>> SEEFOForCO2(Matrix F1,   Matrix F2,    Matrix F3, Matrix F4, Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old,
                                                                 double dt, double delta_h, Matrix F11, Matrix F22, Matrix F33, Matrix F44, Matrix R);
     QVector<Matrix> SolveMUSCL_UL_UR(Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old, double Gamma, double dt_dx,QVector<double>EnergyVibr, double energyStepTemp, double energyStartTemp,int LimType = 1, double omega = 0);

     QVector<QVector<double>> SEEFOForCO23T(Matrix F1,      Matrix F2,      Matrix F3,      Matrix F4,      Matrix F5,
                                            Matrix U1old,   Matrix U2old,   Matrix U3old,   Matrix U4old,   Matrix U5old,
                                            double dt, double delta_h, Matrix R_1, Matrix R_2);

signals:
    void completeSolution();
public:
    QMutex mutex;
    double (*func[FUNCTION_COUNT]) (double t1, double t2, double t3, double t4) = {shareViscositySuperSimple, shareViscositySimple, shareViscosityOmega, bulcViscositySimple,
                                                                                    bulcViscosityOld, bulcViscosityNew, vibrEnergy, fullEnergy, CVibrFunction, lambdaTr,
                                                                                    lambdaVibr, zVibr, EVibr12, EVibr3, Lambda12, Lambda3, Ctr, Crot,lambdaForNitrogen};
    static double shareViscositySuperSimple (double startT, double currentT, double density = 0, double pressure = 0);
    static double shareViscositySimple      (double startT, double currentT, double density = 0, double pressure = 0);
    static double shareViscosityOmega       (double startT, double currentT, double density = 0, double pressure = 0);
    static double bulcViscositySimple       (double startT, double currentT, double density = 0, double pressure = 0);
    static double bulcViscosityOld          (double startT, double currentT, double density = 0, double pressure = 0);

    static double bulcViscosityOld2          (double CVibr, double currentT, double density = 0, double pressure = 0);

    static double bulcViscosityNew          (double startT, double currentT, double density = 0, double pressure = 0);
    static double bulcViscosityOnlyTRRot    (double startT, double currentT, double density = 0, double pressure = 0);
    static double bulcViscosityFalse          (double startT, double currentT, double density = 0, double pressure = 0);
    static double vibrEnergy          (double startT, double currentT, double density = 0, double pressure = 0);
    static double fullEnergy          (double startT, double currentT, double density = 0, double pressure = 0);

    static double fullEnergy(double T, double Tv);
    static double CVibrFunction (double startT, double currentT, double density = 0, double pressure = 0);
    static double Crot(double startT = 0, double currentT = 0, double density = 0, double pressure = 0);
    static double Ctr(double startT = 0, double currentT = 0, double density = 0, double pressure = 0);
    static double CVibr(double T, double ZCO2Vibr);
    static double ZCO2Vibr(double T);
    static double DZDT(double T);
    static double Betta(double T, double CVibr);
    static double TauRot(double T,double P);
    static double TauVibr(double T, double P);

    static double lambda(double T, double CVibr);
    static double lambdaTr(double startT, double currentT, double density = 0, double pressure = 0);
    static double lambdaTr_Rot(double T);
    static double lambdaVibr2(double T, double Tv);
    static double lambdaVibr(double startT, double currentT, double density = 0, double pressure = 0);
    static double getOmega22(double T);
    static double getOmega11(double T);

    static double zVibr(double startT, double currentT, double density = 0, double pressure = 0);

    static double lambdaForNitrogen(double gamma,double T,double density = 0, double pressure = 0);

    static macroParam bondaryConditionRG(macroParam left, solverParams solParams);

    static macroParam bondaryConditionRG2T(macroParam left, solverParams solParams);
    static macroParam bondaryConditionPython(macroParam left, solverParams solParams);
    static double getTimeStep(Matrix velosity, Matrix density, Matrix pressure, double delta_h, solverParams solParams, Matrix gammas = Matrix());
    static double getTimeStepFull(Matrix velosity, Matrix density, Matrix pressure, double delta_h, solverParams solParams, Matrix energy = Matrix());
    static double getTimeStep2(Matrix velosity, Matrix density, Matrix pressure, double delta_h, solverParams solParams, double x, Matrix T);
    static QStringList runPython(macroParam left,int mlt);

    void MackCormack(Matrix U1, Matrix U2, Matrix U3, Matrix U4,
                     Matrix& F1, Matrix& F2, Matrix& F3, Matrix& F4,
                     Matrix& F11, Matrix& F22, Matrix& F33, Matrix& F44,Matrix& R,
                     QList<double>& EnergyVibr);
    double getEnergyVibrTemp(double energy, QList<double>& e);

    static double EVibr12(double ST, double T = 0, double r = 0, double t = 0);
    static double EVibr3(double sT, double T = 0, double r = 0, double t = 0);
    static double Lambda12(double ST, double T = 0, double r = 0, double t = 0);
    static double Lambda3(double sT, double T = 0, double r = 0, double t = 0);
    static double CVibr12(double T);
    static double CVibr3(double T);
    static double ZCO2Vibr12(double T);
    static double ZCO2Vibr3(double T);
    static double tauVibrVVLosev(double t, double P);
};

#endif // ABSTRACTADDITIONALSOLVER_H
