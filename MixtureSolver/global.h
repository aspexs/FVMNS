#ifndef GLOBAL_H
#define GLOBAL_H

#include <QVector>
#include <QtMath>
#include <QDebug>
#include <QObject>
#include <QThread>
#include <QMutexLocker>
#include <QFutureWatcher>
#include <QtConcurrent>

// Постоянная Больцмана
#define K_BOLTZMANN 1.3805e-23

// Атомная единица массы
#define ATOMIC_MASS_UNIT 1.6605402e-27

// Постоянная Планка
#define H_PLANCK 6.6254E-34

// Число Авогадро
#define N_AVOGADRO 6.0221e+23

// Размерный коэффициент для расчета энергии
#define W_W 1.60219e-19 / 8065.47

// Кол-во сортов смеси
#define N_SORTS 2
#define N_POLYATOMIC_SORTS 1

// Характеристики внутренних степеней свободы CO2
#define N_VIBR_DEGREES 3
#define N_VIBR_L1 40
#define N_VIBR_L2 40
#define N_VIBR_L3 40

// Энергия диссоциации : CO2
#define DISS_ENERGY 64017

// Настройка метода Гаусса
#define N_MAX 30

// Число уравнений системы
#define SYSTEM_ORDER 6

/// Структура данных смеси СO2-Ar:
/// 1) Спектроскопические данные для молекулы СО2;
/// 2) Массы молекул сортов СO2, Ar;
/// 3) Газокинетические диаметры молекул сортов СO2, Ar;
/// 4) Глубина потенциальной ямы для молекул сортов : СO2, Ar.
class Mixture
{
public:

    // Методы расчета вспомогательных величин
    static double mass(const int& i);
    static double sigma(const int& i);
    static double epsilon(const int& i);
    static double vEnergy(const int& i, const int& j, const int& k);
    static double reducedMass(const int& i, const int& j);

    // Времена колебательной релаксации
    static double tauVTCO2(const double& t, const double& p);
    static double tauVVCO2(const double& t, const double& p);
};

// Основные макропараметры (rho[0] -> CO2, rho[1] -> Ar)
struct MacroParam
{
    // Основные макропараметры (не независимые)
    QVector<double> rho = {0.0, 0.0};
    double p    = 0.0;
    double v    = 0.0;
    double t    = 0.0;
    double t12  = 0.0;
    double t3   = 0.0;

    // Инициализация плотностей, зная {p, t, x_CO2 in [0, 1]}
    void computeRho(const double& x_CO2);
};

// Параметры решателя
struct SolverParams
{
    // Условия перед УВ и за УВ
    MacroParam lPoint;
    MacroParam rPoint;

    // Число внутренних ячеек и число итераций
    int nCell = 40;
    int nIter = 100;

    // Число Куранта и коэффициент сдвига давления
    double cfl = 0.95;
    double k   = 1e-3;
    double dx  = 0.003;

    // Настройки таблиц температур и энергий
    double t_min = 200.0;
    double t_max = 6000.0;
    int    t_n   = 5800;
};

#endif // GLOBAL_H
