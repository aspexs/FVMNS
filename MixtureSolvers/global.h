#ifndef GLOBAL_H
#define GLOBAL_H

#include <QVector>
#include <QtMath>
#include <QDebug>
#include <QObject>
#include <QThread>
#include <QMutexLocker>
#include <QProgressDialog>
#include <QFutureWatcher>
#include <QtConcurrent>

// Постоянная Больцмана
static const double K_BOLTZMANN = 1.3805e-23;

// Атомная единица массы
static const double ATOMIC_MASS_UNIT = 1.6605402e-27;

// Постоянная Планка
static const double H_PLANCK = 6.6254E-34;

// Число Авогадро
static const double N_AVOGADRO = 6.0221e+23;

// Размерный коэффициент для расчета энергии
static const double W_W = 1.60219e-19 / 8065.47;

// Кол-во сортов смеси
static const int N_SORTS = 2;
static const int N_POLYATOMIC_SORTS = 1;

// Характеристики внутренних степеней свободы CO2
static const int N_VIBR_DEGREES = 3;
static const int N_VIBR_L1 = 40;
static const int N_VIBR_L2 = 40;
static const int N_VIBR_L3 = 40;

// Энергия диссоциации : CO2
static const double DISS_ENERGY = 64017;

// Настройка метода Гаусса
const int N_MAX = 30;

// Число уравнений системы
const int SYSTEM_ORDER = 6;

/// Структура данных смеси СO2-Ar:
/// 1) Спектроскопические данные для молекулы СО2;
/// 2) Массы молекул сортов СO2, Ar;
/// 3) Газокинетические диаметры молекул сортов СO2, Ar;
/// 4) Глубина потенциальной ямы для молекул сортов : СO2, Ar.
class Mixture
{
public:

    // Конструктор и деструктор по умолчанию
    Mixture() = default;
    ~Mixture() = default;

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
    double cfl = 0.8;
    double k   = 1e-3;
    double dx  = 0.003;

    // Настройки таблиц температур и энергий
    double t_min = 200.0;
    double t_max = 5000.0;
    int    t_n   = 9600;
};

class Matrix
{
private:
    QVector<double> data;
public:
    Matrix (QVector<double> val)
    {
        data = val;
    }
    Matrix()
    {

    }
    void clear()
    {
        data.clear();
    }
    Matrix (int len, double val = 0)
    {
        data = QVector<double> (len, val);
    }
    QVector<double>::iterator begin()
    {
        return data.begin();
    }
    QVector<double>::iterator end()
    {
        return data.end();
    }
    int size()
    {
        return data.size();
    }
    operator QVector<double>() const
    {
        return data;
    }
    double& operator [](int i)
    {
        static double elseValue;
        if(i < data.size())
            return data[i];
        else
            return elseValue;
    }
    double first()
    {
        return data.first();
    }
    double last()
    {
        return data.last();
    }
    void push_back(double val)
    {
        data.push_back(val);
    }
    void push_front(double val)
    {
        data.push_front(val);
    }
    void removeLast()
    {
        data.removeLast();
    }
    void removeFirst()
    {
        data.removeFirst();
    }
    void resize(int i)
    {
        data.resize(i);
    }
    void fill(const double& t, int size = -1)
    {
        data.fill(t, size);
    }
    Matrix operator /(QVector<double> div)
    {
        QVector<double> output;
        if(data.size() != div.size())
            return data;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] / div[i]);
        return output;
    }

    Matrix operator *(QVector<double> div)
    {
        QVector<double> output;
        if(data.size() != div.size())
            return data;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] * div[i]);
        return output;
    }
    Matrix operator +(QVector<double> div)
    {
        QVector<double> output;
        if(data.size() != div.size())
            return data;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] + div[i]);
        return output;
    }
    Matrix operator -(QVector<double> div)
    {
        QVector<double> output;
        if(data.size() != div.size())
            return data;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] - div[i]);
        return output;
    }
    Matrix operator /(double div)
    {
        QVector<double> output;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] / div);
        return output;
    }
    Matrix operator *(double div)
    {
        QVector<double> output;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] * div);
        return output;
    }
    Matrix operator +(double div)
    {
        QVector<double> output;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] + div);
        return output;
    }
    Matrix operator -(double div)
    {
        QVector<double> output;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] - div);
        return output;
    }
    static Matrix POW(QVector<double> div, double param)
    {
        QVector<double> output;
        for(int i = 0; i < div.size(); i ++)
             output.push_back(pow(div[i],param));
        return output;
    }
    static Matrix SQRT(QVector<double> div)
    {
        QVector<double> output;
        for(int i = 0; i < div.size(); i ++)
             output.push_back(sqrt(div[i]));
        return output;
    }
    static Matrix REVERSE(QVector<double> div)
    {
        QVector<double> output;
        for(int i = 0; i < div.size(); i ++)
             output.push_back(1.0/div[i]);
        return output;
    }
};

#endif // GLOBAL_H
