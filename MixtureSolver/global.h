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
#include <iostream>

// Настройка используемых коэффицентов
#define USE_V_DIFF 1
#define USE_B_VISC 1

// Настройка температурного диапазона
#define T_MAX 5000.0
#define T_MIN 200.0
#define T_NUM 9600

// Малая величина
#define EPSILON 1e-6

// Число ячеек сетки решателя, число итераций
#define N_CELL 50
#define N_ITER 5000

// Число Куранта и длина одной ячейки
#define CFL 0.9
#define DX 0.0005

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
#define DISS_ENERGY 64017.0

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
class MacroParam
{
public:

    // Основные макропараметры (не независимые)
    QVector<double> rho = {0.0, 0.0};
    double p, v, t, t12, t3;

public:


    // Инициализация по умолчанию и для равновесного случая
    MacroParam();
    MacroParam(const double& p, const double& v, const double& t,
               const double& x_CO2);

    // Инициализация плотностей, зная {p, t, x_CO2 in [0, 1]}
    void computeRho(const double& x_CO2);

    // Инициализация для равновесного случая
    void initialize(const double& p, const double& v, const double& t,
                    const double& x_CO2);

    // Продолжение по гладкости
    MacroParam proceed(const MacroParam& p0);
};

/// ProgressBar - реализует шкалу прогресса в консоли
class ProgressBar
{
public:

    // Текущее и предыдущее заполнение, шаг
    double n0, n1, dx;

public:

    // Обнуление
    ProgressBar();

    // Инициализация для равновесного случая
    void initialize(const double& n);

    // Обновить шкалу
    void update();
};

#endif // GLOBAL_H
