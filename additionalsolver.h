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
        FUNCTION_COUNT
    };
    enum BoundaryConditions
    {
        BC_RG,
        BC_RG2T,
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
        BV_COUNT
    };

    explicit AdditionalSolver(QObject *parent = nullptr);

    void solve();
    double startValue;
    double stopValue;
    double step;
    double pressure = 0;
    double density = 0;
    Function typeSolve;
    QVector <double> iterationVector;
    QVector <double> rezultVector;
    macroParam (*BoundaryCondition[BC_COUNT])(macroParam left, solverParams solParams) = {bondaryConditionRG, bondaryConditionRG2T};
    double (*TimeStepSolution[TS_COUNT])(Matrix velosity, Matrix density, Matrix pressure, double delta_h, solverParams solParams, Matrix energy) = {getTimeStep, getTimeStepFull};
    double (*shareViscosity[SV_COUNT])(double t1, double t2, double t3, double t4) = {shareViscositySuperSimple, shareViscositySimple, shareViscosityOmega};
    double (*bulkViscosity[BV_COUNT])(double t1, double t2, double t3, double t4) = {bulcViscositySimple, bulcViscosityOld2,bulcViscosityFalse};

    macroParam ExacRiemanSolver(macroParam left, macroParam right, double Gamma);
    macroParam ExacRiemanSolverCorrect(macroParam left, macroParam right, double GammaL, double GammaR);

    QVector<QVector<double>> SolveEvolutionExplFirstOrder(Matrix F1,   Matrix F2,    Matrix F3, Matrix U1old, Matrix U2old, Matrix U3old,
                                                          double dt, double delta_h);
    QVector<QVector<double>> SolveEvolutionExplFirstOrderForCO22(Matrix F1,   Matrix F2,    Matrix F3, Matrix F4, Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old,
                                                                 double dt, double delta_h);
signals:
    void completeSolution();
public:
    QMutex mutex;
    double (*func[FUNCTION_COUNT]) (double t1, double t2, double t3, double t4) = {shareViscositySuperSimple, shareViscositySimple, shareViscosityOmega,
                                                                                    bulcViscositySimple,bulcViscosityOld, bulcViscosityNew, vibrEnergy, fullEnergy, CVibrFunction, lambdaTr, lambdaVibr};
    static double shareViscositySuperSimple (double startT, double currentT, double density = 0, double pressure = 0);
    static double shareViscositySimple      (double startT, double currentT, double density = 0, double pressure = 0);
    static double shareViscosityOmega       (double startT, double currentT, double density = 0, double pressure = 0);
    static double bulcViscositySimple       (double startT, double currentT, double density = 0, double pressure = 0);
    static double bulcViscosityOld          (double startT, double currentT, double density = 0, double pressure = 0);

    static double bulcViscosityOld2          (double CVibr, double currentT, double density = 0, double pressure = 0);

    static double bulcViscosityNew          (double startT, double currentT, double density = 0, double pressure = 0);
    static double bulcViscosityFalse          (double startT, double currentT, double density = 0, double pressure = 0);
    static double vibrEnergy          (double startT, double currentT, double density = 0, double pressure = 0);
    static double fullEnergy          (double startT, double currentT, double density = 0, double pressure = 0);

    static double fullEnergy(double T, double Tv);
    static double CVibrFunction (double startT, double currentT, double density = 0, double pressure = 0);
    static double Crot();
    static double Ctr();
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



    static macroParam bondaryConditionRG(macroParam left, solverParams solParams);

    static macroParam bondaryConditionRG2T(macroParam left, solverParams solParams);
    static double getTimeStep(Matrix velosity, Matrix density, Matrix pressure, double delta_h, solverParams solParams, Matrix gammas = Matrix());
    static double getTimeStepFull(Matrix velosity, Matrix density, Matrix pressure, double delta_h, solverParams solParams, Matrix energy = Matrix());
};

#endif // ABSTRACTADDITIONALSOLVER_H
