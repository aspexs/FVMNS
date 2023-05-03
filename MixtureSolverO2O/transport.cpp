///////////////////////////////////////////////////////////////////////////////
/// Автор : Баталов Семен, 2022
/// Модуль предназаначен для расчета коэффициентов переноса для смеси газов
/// CO2-Ar, а именно: коэффициентов диффузии, термодиффузии, объемной и
/// сдвиговой вязкостей, коэффициентов теплопроводности.
///////////////////////////////////////////////////////////////////////////////

#include "transport.h"

///////////////////////////////////////////////////////////////////////////////
/// class OmegaIntegrals
///////////////////////////////////////////////////////////////////////////////

// Инициализация
void OmegaIntegrals::initialize()
{
    // Подгонка размеров столбцов матриц омега-интегралов
    omega11_vv.resize(N_SORTS);
    omega12_vv.resize(N_SORTS);
    omega13_vv.resize(N_SORTS);
    omega22_vv.resize(N_SORTS);
    aa_vv.resize(N_SORTS);
    bb_vv.resize(N_SORTS);
    cc_vv.resize(N_SORTS);

    // Подгонка размеров строк матриц омега-интегралов
    for (int i = 0; i < N_SORTS; ++i)
    {
        omega11_vv[i].fill(0.0, N_SORTS);
        omega12_vv[i].fill(0.0, N_SORTS);
        omega13_vv[i].fill(0.0, N_SORTS);
        omega22_vv[i].fill(0.0, N_SORTS);
        aa_vv[i].fill(0.0, N_SORTS);
        bb_vv[i].fill(0.0, N_SORTS);
        cc_vv[i].fill(0.0, N_SORTS);
    }
}

// Доступ к соответствующим полям
const QVector<QVector<double>>& OmegaIntegrals::omega11() const
{
    return omega11_vv;
}
const QVector<QVector<double>>& OmegaIntegrals::omega12() const
{
    return omega12_vv;
}
const QVector<QVector<double>>& OmegaIntegrals::omega13() const
{
    return omega13_vv;
}
const QVector<QVector<double>>& OmegaIntegrals::omega22() const
{
    return omega22_vv;
}
const QVector<QVector<double>>& OmegaIntegrals::aa() const
{
    return aa_vv;
}
const QVector<QVector<double>>& OmegaIntegrals::bb() const
{
    return bb_vv;
}
const QVector<QVector<double>>& OmegaIntegrals::cc() const
{
    return cc_vv;
}

///////////////////////////////////////////////////////////////////////////////
/// class OmegaIntegralsDc
///////////////////////////////////////////////////////////////////////////////

// Инициализация
void OmegaIntegralsDc::initialize()
{
    // Предварительная инициализация
    OmegaIntegrals::initialize();
}

///////////////////////////////////////////////////////////////////////////////
/// class OmegaIntegralsDcLjp
///////////////////////////////////////////////////////////////////////////////

// Инициализация
void OmegaIntegralsDcLjp::initialize()
{
    // Предварительная инициализация
    OmegaIntegralsDc::initialize();

    // Обнуление
    sig_ij = 0.0;
    e_ij = 0.0;
    s_ij = 0.0;
    tx = 0.0;
    x11 = 0.0;
}

// Расчет всех значений омега-интегралов
void OmegaIntegralsDcLjp::compute(const double& t)
{
    // Проход по возможным сочетаниям сортов
    for (int i = 0; i < N_SORTS; ++i)
    {
        for (int j = 0; j < N_SORTS; ++j)
        {
            // Инициализация вспомогательных переменных
            sig_ij = (Mixture::sigma(i) + Mixture::sigma(j)) / 2.0;
            e_ij = qSqrt(Mixture::epsilon(i) * Mixture::epsilon(j) *
                        qPow(Mixture::sigma(i) * 1e10, 6.0) *
                        qPow(Mixture::sigma(j) * 1e10, 6.0)) /
                        qPow(sig_ij * 1e10, 6.0);
            s_ij = M_PI * qPow(sig_ij, 2.0) *
                    qSqrt(K_BOLTZMANN * t / 2.0 / M_PI /
                          Mixture::reducedMass(i, j));
            tx = t / e_ij;

            // Расчет омега-интегралов в соответствии с
            // потенциалом Леннарда-Джонса
            x11 = qLn(tx) + 1.4;
            omega11_vv[i][j] = 1.0 / (-0.16845 - 2.25768e-2 / qPow(x11, 2.0) +
                                    0.19779 / x11 + 0.64373 * x11 -
                                    9.26718e-2 * qPow(x11, 2.0) +
                                    7.1131e-3 * qPow(x11, 3.0)) * s_ij;
            x11 = qLn(tx) + 1.5;
            omega22_vv[i][j] = 1.0 / (-0.40811 - 5.08552e-2 / qPow(x11, 2.0) +
                                    0.3401 / x11 + 0.70375 * x11 -
                                    0.10699 * qPow(x11, 2.0) +
                                    7.62686e-3 * qPow(x11, 3.0)) * s_ij * 2.0;

            x11 = qLn(tx) + 1.1;
            omega12_vv[i][j] = 1.0 / (0.40785 + 9.25303e-4 / qPow(x11, 2.0) +
                                    2.79680e-4 / x11 + 0.44739 * x11 -
                                    6.27242e-2 * qPow(x11, 2.0) +
                                    5.98567e-3 * qPow(x11, 3.0)) * s_ij * 3.0;
            x11 = qLn(tx) + 4.0;
            omega13_vv[i][j] = 1.0 / (25.04929 + 63.12444 / qPow(x11, 2.0) -
                                    65.87398 / x11 - 4.13758 * x11 +
                                    0.34999 * qPow(x11, 2.0) -
                                    1.0096e-2 * qPow(x11, 3.0)) * s_ij * 12.0;

            // Расчет комбинаций интегралов
            bb_vv[i][j] = (5.0 / 3.0 * omega12_vv[i][j] - 4.0 / 12.0 *
                           omega13_vv[i][j]) / omega11_vv[i][j];
            cc_vv[i][j] = 1.0 / omega11_vv[i][j] * omega12_vv[i][j] / 3.0;
            aa_vv[i][j] = 1.0 / omega11_vv[i][j] * omega22_vv[i][j] / 2.0;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
/// class SpecificHeat
///////////////////////////////////////////////////////////////////////////////

// Инициализация
void SpecificHeat::initialize()
{
    zv_s = 0.0;
    cv_s = 0.0;
}

// Доступ к соответствующим полям
const double& SpecificHeat::zv() const
{
    return zv_s;
}
const double& SpecificHeat::cv() const
{
    return cv_s;
}

///////////////////////////////////////////////////////////////////////////////
/// class SpecificHeatDc
///////////////////////////////////////////////////////////////////////////////

void SpecificHeatDc::computeZv(const double& tv)
{
    // Вспомогательные переменные
    double ee = 0.0;

    // Инициализация
    zv_s = 0.0;

    // Проходим по всем колебательным уровням O2
    for (int i = 0; i < N_VIBR_L; ++i)
    {
        ee = Mixture::vEnergy(i);
        zv_s += qExp(-ee / K_BOLTZMANN / tv);
    }
}
void SpecificHeatDc::computeCv(const double& tv)
{
    // Вспомогательные переменные
    double ee = 0.0;
    double s = 0.0;
    double ss = 0.0;
    double ppp = 0.0;

    // Проходим по всем колебательным уровням, мода 3
    for (int i = 0; i < N_VIBR_L; ++i)
    {
        ee = Mixture::vEnergy(i);
        ppp = qExp(-ee / K_BOLTZMANN / tv) / zv_s;
        s += ee / K_BOLTZMANN / tv * ppp;
        ss += qPow(ee / K_BOLTZMANN / tv, 2.0) * ppp;
    }

    // Вычисление удельных теплоемкостей
    cv_s = ss - qPow(s, 2.0);
}
void SpecificHeatDc::initialize()
{
    // Предварительная инициализация
    SpecificHeat::initialize();
}
void SpecificHeatDc::compute(const double& tv)
{
    // Порядок имеет значение
    computeZv(tv);
    computeCv(tv);
}

///////////////////////////////////////////////////////////////////////////////
/// class Energy
///////////////////////////////////////////////////////////////////////////////

void Energy::initialize()
{
    vE_s = 0.0;
    rE_s = 0.0;
    tE_s = 0.0;
    fullE_s = 0.0;
}
const double& Energy::vE() const
{
    return vE_s;
}
const double& Energy::rE() const
{
    return rE_s;
}
const double& Energy::tE() const
{
    return tE_s;
}
const double& Energy::fullE() const
{
    return fullE_s;
}

///////////////////////////////////////////////////////////////////////////////
/// class EnergyDc
///////////////////////////////////////////////////////////////////////////////

// Расчет соответствующих энергий
void EnergyDc::computeVe(const double& tv)
{
    double ee = 0.0;
    vE_s = 0.0;
    for (int i = 0; i < N_VIBR_L; i++)
    {
        ee = Mixture::vEnergy(i);
        vE_s += ee * qExp(-ee / K_BOLTZMANN / tv);
    }
    vE_s *= 1.0 / (heat_.zv() * Mixture::mass(0));
}
void EnergyDc::computeRe(const double& t)
{
    rE_s = K_BOLTZMANN * t / Mixture::mass(0);
}
void EnergyDc::computeTe(const MacroParam& param)
{
    tE_s = 1.5 * K_BOLTZMANN * param.t * (param.rho[0] / Mixture::mass(0) +
            param.rho[1] / Mixture::mass(1)) / (param.rho[0] + param.rho[1]);
}
void EnergyDc::computeFullE(const MacroParam& param)
{
    fullE_s = tE_s + param.rho[0] / (param.rho[0] + param.rho[1]) *
            (rE_s + vE_s);
}

// Инициализация
void EnergyDc::initialize()
{
    // Предварительная инициализация
    Energy::initialize();
    heat_.initialize();
}

// Расчет всех значений
void EnergyDc::compute(const MacroParam& param)
{
    // Порядок имеет значение
    heat_.compute(param.tv);
    computeVe(param.tv);
    computeRe(param.t);
    computeTe(param);
    computeFullE(param);
}
void EnergyDc::compute(const double& tv)
{
    // Порядок имеет значение
    heat_.compute(tv);
    computeVe(tv);
}
const SpecificHeat& EnergyDc::heat() const
{
    return heat_;
}

///////////////////////////////////////////////////////////////////////////////
/// class Temperature
///////////////////////////////////////////////////////////////////////////////

// Выделение памяти, обновление значений
void Temperature::initialize()
{
    vT_s = 0.0;
    tT_s = 0.0;
}

// Доступ к соответствующим полям
const double& Temperature::Tv() const
{
    return vT_s;
}
const double& Temperature::T() const
{
    return tT_s;
}

///////////////////////////////////////////////////////////////////////////////
/// class TemperatureNDc
///////////////////////////////////////////////////////////////////////////////

// Расчет всех значений температуры
void TemperatureNDc::calcTv(const double& ev)
{
    // Вспомогательные переменные
    int nL = 0;
    int nR = vE_v.size() - 1;
    int n = 0;

    // Двоичный поиск
    while (nR - nL > 1)
    {
        n = 0.5 * (nL + nR);
        if (vE_v[n] > ev)
        {
            nR = n;
        }
        else
        {
            nL = n;
        }
    }

    // Линейная интерполяция
    vT_s = vT_v[nL] + dT * (ev - vE_v[nL]) / (vE_v[nR] - vE_v[nL]);
}
void TemperatureNDc::calcT(const MacroParam& p, const double& e,
                           const double& ev)
{
    tT_s = ((p.rho[0] + p.rho[1]) * e - p.rho[0] * ev) /
            K_BOLTZMANN / (2.5 * p.rho[0] / Mixture::mass(0) +
            1.5 * p.rho[1] / Mixture::mass(1));
}

// Выделение памяти, обновление значений
void TemperatureNDc::initialize(const double& t0, const double& t1,
                                const int& n)
{
    Temperature::initialize();
    energy_.initialize();
    vT_v.fill(0.0, n);
    vE_v.fill(0.0, n);

    // Подготовка инструментов
    double tv = 0.0;
    dT = (t1 - t0) / n;

    // Заполняем таблицы значений
    for (int i = 0; i < n; ++i)
    {
        tv = t0 + dT * i;
        energy_.compute(tv);
        vT_v[i] = tv;
        vE_v[i] = energy_.vE();
    }
}

// Расчет всех значений
void TemperatureNDc::compute(const MacroParam& param, const double& ev,
                             const double& e)
{
    calcTv(ev);
    calcT(param, e, ev);
}

///////////////////////////////////////////////////////////////////////////////
/// class BracketIntegrals
//////////////////////////////////////////////////////////////////////////////

// Инициализация
void BracketIntegrals::initialize()
{
    // Подгонка размеров столбцов матриц
    lambda_vv.resize(N_SORTS);
    lambda00_vv.resize(N_SORTS);
    lambda01_vv.resize(N_SORTS);
    lambda11_vv.resize(N_SORTS);
    eta_vv.resize(N_SORTS);
    h00_vv.resize(N_SORTS);
    beta11_vv.resize(N_SORTS);
    beta01_vv.resize(N_SORTS);
    beta0011_v.resize(N_POLYATOMIC_SORTS);
    lambdaInt_v.fill(0.0, N_SORTS);

    // Подгонка размеров строк матриц омега-интегралов
    for (int i = 0; i < N_SORTS; ++i)
    {
        lambda_vv[i].fill(0.0, N_SORTS);
        lambda00_vv[i].fill(0.0, N_SORTS);
        lambda01_vv[i].fill(0.0, N_SORTS);
        lambda11_vv[i].fill(0.0, N_SORTS);
        eta_vv[i].fill(0.0, N_SORTS);
        h00_vv[i].fill(0.0, N_SORTS);
        beta11_vv[i].fill(0.0, N_SORTS);
        beta01_vv[i].fill(0.0, N_POLYATOMIC_SORTS);
    }
}

// Доступ к соответствующим полям
const QVector<QVector<double>>& BracketIntegrals::lambda() const
{
    return lambda_vv;
}
const QVector<QVector<double>>& BracketIntegrals::lambda00() const
{
    return lambda00_vv;
}
const QVector<QVector<double>>& BracketIntegrals::lambda01() const
{
    return lambda01_vv;
}
const QVector<QVector<double>>& BracketIntegrals::lambda11() const
{
    return lambda11_vv;
}
const QVector<QVector<double>>& BracketIntegrals::eta() const
{
    return eta_vv;
}
const QVector<QVector<double>>& BracketIntegrals::h00() const
{
    return h00_vv;
}
const QVector<QVector<double>>& BracketIntegrals::beta11() const
{
    return beta11_vv;
}
const QVector<QVector<double>>& BracketIntegrals::beta01() const
{
    return beta01_vv;
}
const QVector<double>& BracketIntegrals::beta0011() const
{
    return beta0011_v;
}
const QVector<double>& BracketIntegrals::lambdaInt() const
{
    return lambdaInt_v;
}

///////////////////////////////////////////////////////////////////////////////
/// class BracketIntegralsDc
///////////////////////////////////////////////////////////////////////////////

// Расчет phi
void BracketIntegralsDc::computePhi(const double& t)
{
    phi_v[0] = K_BOLTZMANN / 23.73 *
            (1.0 + qPow(M_PI, 1.5) / 2.0 * qSqrt(Mixture::epsilon(0) / t) +
            (qPow(M_PI, 2.0) / 4.0 + 2.0) * Mixture::epsilon(0) / t +
            qPow(M_PI * Mixture::epsilon(0) / t, 1.5));
}

// Расчет всех интегральных скобок
void BracketIntegralsDc::computeBrackets1(const double& t)
{
    for (int i = 0; i < N_SORTS; ++i)
    {
        for (int j = 0; j < N_SORTS; ++j)
        {
            lambda_vv[i][j] = 75.0 / 64.0 * t * K_BOLTZMANN /
                    Mixture::reducedMass(i, j) / omega_.omega22()[i][j] *
                    K_BOLTZMANN;
            eta_vv[i][j] = 5.0 / 8.0 * K_BOLTZMANN * t / omega_.omega22()[i][j];
        }
    }
}
void BracketIntegralsDc::computeBrackets2(const double& t,
                                          const QVector<double>& x)
{
    for (int i = 0; i < N_SORTS; ++i)
    {
        lambdaInt_v[i] = 0.0;
    }
    for (int i = 0; i < N_SORTS; ++i)
    {
        for (int j = 0; j < N_SORTS; ++j)
        {
            if (i == j)
            {
                lambda00_vv[i][j] = 0.0;
                lambda01_vv[i][j] = 0.0;
                lambda11_vv[i][j] = qPow(x[i], 2.0) / lambda_vv[i][i];
                h00_vv[i][j] = qPow(x[i], 2.0) / eta_vv[i][i];
                beta11_vv[i][j] = 4.0 * t / M_PI * qPow(x[i], 2.0) /
                        eta_vv[i][j] * phi_v[i];
                for (int k = 0; k < N_SORTS; ++k)
                {
                    if (k != i)
                    {
                        lambda00_vv[i][j] = lambda00_vv[i][j] + x[i] * x[k] /
                                lambda_vv[i][k] / 2.0 / omega_.aa()[i][k];
                        lambda01_vv[i][j] = lambda01_vv[i][j] - x[i] * x[k] /
                                lambda_vv[i][k] / 4.0 / omega_.aa()[i][k] *
                                Mixture::mass(k) /
                                (Mixture::mass(i) + Mixture::mass(k)) *
                                (6.0 * omega_.cc()[i][k] - 5.0);
                        lambda11_vv[i][j] = lambda11_vv[i][j] + x[i] * x[k] /
                                lambda_vv[i][k] / 2.0 / omega_.aa()[i][k] *
                                (15.0 / 2.0 * qPow(Mixture::mass(i), 2.0) +
                                 25.0 / 4.0 * qPow(Mixture::mass(k), 2.0) -
                                 3.0 * qPow(Mixture::mass(k), 2.0) *
                                 omega_.bb()[i][k] + 4.0 *
                                 Mixture::mass(i) * Mixture::mass(k) *
                                 omega_.aa()[i][k]) /
                                qPow(Mixture::mass(i) +
                                     Mixture::mass(k), 2.0);
                        h00_vv[i][j] = h00_vv[i][j] + 2.0 * x[i] * x[k] /
                                eta_vv[i][k] * Mixture::mass(i) *
                                Mixture::mass(k) /
                                qPow(Mixture::mass(i) +
                                     Mixture::mass(k), 2.0) *
                                (5.0 / 3.0 / omega_.aa()[i][k] +
                                 Mixture::mass(k) / Mixture::mass(i));
                        beta11_vv[i][j] = beta11_vv[i][j] + x[i] * x[k] /
                                eta_vv[i][k] * Mixture::mass(k) /
                                qPow(Mixture::mass(i) +
                                     Mixture::mass(k), 2.0) *
                                (5.0 * K_BOLTZMANN * t * Mixture::mass(i) /
                                 omega_.aa()[i][k] + 4.0 * t *
                                 Mixture::mass(k) / M_PI *
                                 (phi_v[i] + phi_v[k]));
                    }
                }
            }
            else
            {
                lambda00_vv[i][j] = -x[i] * x[j] / lambda_vv[i][j] / 2.0 /
                        omega_.aa()[i][j];
                lambda01_vv[i][j] = x[i] * x[j] / lambda_vv[i][j] / 4.0 /
                        omega_.aa()[i][j] * Mixture::mass(i) /
                        (Mixture::mass(i) + Mixture::mass(j)) *
                        (6.0 * omega_.cc()[i][j] - 5.0);
                lambda11_vv[i][j] = -x[i] * x[j] / lambda_vv[i][j] / 2.0 /
                        omega_.aa()[i][j] * Mixture::reducedMass(i, j) /
                        (Mixture::mass(i) + Mixture::mass(j)) *
                        (55.0 / 4.0 - 3.0 * omega_.bb()[i][j] - 4.0 *
                         omega_.aa()[i][j]);
                h00_vv[i][j] = -2.0 * x[i] * x[j] / eta_vv[i][j] *
                        Mixture::reducedMass(i, j) /
                        (Mixture::mass(i) + Mixture::mass(j)) *
                        (5.0 / 3.0 / omega_.aa()[i][j] - 1.0);
                beta11_vv[i][j] = x[i] * x[j] / eta_vv[i][j] *
                        Mixture::mass(i) * Mixture::mass(j) /
                        qPow(Mixture::mass(i) +
                             Mixture::mass(j), 2.0) *
                        (-5.0 * K_BOLTZMANN * t / omega_.aa()[i][j] +
                         4.0 * t / M_PI * (phi_v[i] + phi_v[j]));
            }
            lambdaInt_v[i] += x[j] * Mixture::reducedMass(i, j) *
                    omega_.omega11()[i][j];
        }
    }
}
void BracketIntegralsDc::computeBrackets3(const double& t,
                                          const QVector<double>& x)
{
    for (int i = 0; i < N_POLYATOMIC_SORTS; ++i)
    {
        beta0011_v[i] = 0.0;
    }
    for (int i = 0; i < N_SORTS; ++i)
    {
        for (int j = 0; j < N_POLYATOMIC_SORTS; ++j)
        {
            if (i == j)
            {
                beta01_vv[i][j] = -4.0 * t / M_PI * qPow(x[i], 2.0) /
                        eta_vv[i][i] * phi_v[i];
                for (int k = 0; k < N_SORTS; ++k)
                {
                    if (k != i)
                    {
                        beta01_vv[i][j] = beta01_vv[i][j] -
                                4.0 * t / M_PI * x[i] * x[k] / eta_vv[i][k] *
                                Mixture::mass(k) /
                                (Mixture::mass(i) + Mixture::mass(k)) *
                                phi_v[i];
                    }
                }
            }
            else
            {
                beta01_vv[i][j] = -4.0 * t / M_PI * x[i] * x[j] /
                        eta_vv[i][j] * Mixture::mass(j) /
                        (Mixture::mass(i) + Mixture::mass(j)) *
                        phi_v[j];
            }
            beta0011_v[j] += 4.0 * t / M_PI * x[i] * phi_v[j] * x[j] /
                    eta_vv[i][j];
        }
    }
}

// Инициализация
void BracketIntegralsDc::initialize()
{
    // Предварительная инициализация
    BracketIntegrals::initialize();
    omega_.initialize();

    // Подготовка массива Phi
    phi_v.fill(0.0, N_SORTS);
}

// Расчет всех значений
void BracketIntegralsDc::compute(const MacroParam& param)
{
    // Мольная доля
    QVector<double> n = {param.rho[0] / Mixture::mass(0),
                         param.rho[1] / Mixture::mass(1)};
    QVector<double> x = {n[0] / (n[0] + n[1]),
                         n[1] / (n[0] + n[1])};

    // Порядок имеет значение
    omega_.compute(param.t);
    computePhi(param.t);
    computeBrackets1(param.t);
    computeBrackets2(param.t, x);
    computeBrackets3(param.t, x);
}
const OmegaIntegrals& BracketIntegralsDc::omega() const
{
    return omega_;
}

///////////////////////////////////////////////////////////////////////////////
/// class TransportCoefficients
///////////////////////////////////////////////////////////////////////////////

// Выделение памяти, обновление значений
void TransportCoefficients::initialize()
{
    // Обнуление значений
    sViscosity_s = 0.0;
    bViscosity_s = 0.0;
    tLambda_s = 0.0;
    rLambda_s = 0.0;
    cLambda_s = 0.0;
    vLambda_s = 0.0;

    // Выделение памяти и обнуление значений
    diffusion_vv.resize(N_SORTS);
    tDiffusion_v.fill(0.0, N_SORTS);
    for (int i = 0; i < N_SORTS; ++i)
    {
        diffusion_vv[i].fill(0.0, N_SORTS);
    }
}

// Доступ к соответствующим полям
const double& TransportCoefficients::sViscosity() const
{
    return sViscosity_s;
}
const double& TransportCoefficients::bViscosity() const
{
    return bViscosity_s;
}
const double& TransportCoefficients::tLambda() const
{
    return tLambda_s;
}
const double& TransportCoefficients::rLambda() const
{
    return rLambda_s;
}
const double& TransportCoefficients::cLambda() const
{
    return cLambda_s;
}
const double& TransportCoefficients::vLambda() const
{
    return vLambda_s;
}
const QVector<QVector<double>>& TransportCoefficients::diffusion() const
{
    return diffusion_vv;
}
const QVector<double>& TransportCoefficients::tDiffusion() const
{
    return tDiffusion_v;
}

///////////////////////////////////////////////////////////////////////////////
/// class TransportCoefficientsDc
///////////////////////////////////////////////////////////////////////////////

// Заполнение матриц систем и матриц свободных членов
void TransportCoefficientsDc::fillBMLambda(const double& nTot,
                                           const double& rho,
                                           const QVector<double>& x)
{
    // Заполняем матрицу системы
    for (int i = 0; i < N_SORTS; ++i)
    {
        for (int j = 0; j < N_SORTS; ++j)
        {
            mLambda_vv[i][j] = bracket_.lambda00()[i][j];
        }
    }
    for (int i = 0; i < N_SORTS; ++i)
    {
        for (int j = N_SORTS; j < N_SORTS * 2; ++j)
        {
            mLambda_vv[i][j] = bracket_.lambda01()[i][j - N_SORTS];
        }
    }
    for (int i = N_SORTS; i < N_SORTS * 2; ++i)
    {
        for (int j = 0; j < N_SORTS; ++j)
        {
            mLambda_vv[i][j] = mLambda_vv[j][i];
        }
    }
    for (int i = N_SORTS; i < N_SORTS * 2; ++i)
    {
        for (int j = N_SORTS; j < N_SORTS * 2; ++j)
        {
            mLambda_vv[i][j] = bracket_.lambda11()[i - N_SORTS][j - N_SORTS];
        }
    }
    for (int j = 0; j < N_SORTS; ++j)
    {
        mLambda_vv[0][j] = x[j] * Mixture::mass(j) * nTot / rho;
    }
    for (int j = N_SORTS; j < N_SORTS * 2; ++j)
    {
        mLambda_vv[0][j] = 0.0;
    }

    // Заполняем матрицу свободных членов
    for (int i = 0; i < N_SORTS; ++i)
    {
        bLambda_vv[i][0] = 0.0;
    }
    for (int i = N_SORTS; i < N_SORTS * 2; ++i)
    {
        bLambda_vv[i][0] = 4.0 / 5.0 / K_BOLTZMANN * x[i - N_SORTS];
    }
}
void TransportCoefficientsDc::fillBMDiffusion(const double& nTot,
                                              const double& rho,
                                              const QVector<double>& x)
{
    // Вспомогательные переменные
    double delta = 0.0;

    // Заполняем матрицу системы
    for (int i = 0; i < N_SORTS; ++i)
    {
        for (int j = 0; j < N_SORTS; ++j)
        {
            mDiffusion_vv[i][j] = mLambda_vv[i][j];
        }
    }

    // Заполняем матрицу свободных членов
    for (int i = 0; i < N_SORTS; ++i)
    {
        for (int j = 0; j < N_SORTS; ++j)
        {
            if (i == j)
            {
                delta = 1.0;
            }
            else
            {
                delta = 0.0;
            }
            bDiffusion_vv[i][j] = 8.0 / 25.0 / K_BOLTZMANN *
                    (delta - Mixture::mass(i) * x[i] * nTot / rho);
        }
        bDiffusion_vv[0][i] = 0.0;
    }
}
void TransportCoefficientsDc::fillBMSViscosity(const double& t,
                                               const QVector<double>& x)
{
    // Заполняем матрицу системы
    for (int i = 0; i < N_SORTS; ++i)
    {
        for (int j = 0; j < N_SORTS; ++j)
        {
            mSViscosity_vv[i][j] = bracket_.h00()[i][j];
        }
    }

    // Заполняем матрицу свободных членов
    for (int i = 0; i < N_SORTS; ++i)
    {
        bSViscosity_vv[i][0] = 2.0 / K_BOLTZMANN / t * x[i];
    }
}
void TransportCoefficientsDc::fillBMBViscosity(const double& nTot,
                                               const double& rho,
                                               const QVector<double>& x)
{
    // Вспомогательные переменные
    double cu = 0.0;
    double cuT = 0.0;

    // Заполняем матрицу системы
    for (int i = 0; i < N_SORTS; ++i)
    {
        for (int j = 0; j < N_SORTS; ++j)
        {
            mBViscosity_vv[i][j] = bracket_.beta11()[i][j];
        }
    }
    for (int i = 0; i < N_SORTS; ++i)
    {
        for (int j = N_SORTS; j < N_SORTS + N_POLYATOMIC_SORTS; ++j)
        {
            mBViscosity_vv[i][j] = bracket_.beta01()[i][j - N_SORTS];
        }
    }
    for (int i = N_SORTS; i < N_SORTS + N_POLYATOMIC_SORTS; ++i)
    {
        for (int j = 0; j < N_SORTS; ++j)
        {
            mBViscosity_vv[i][j] = mBViscosity_vv[j][i];
        }
    }
    for (int i = N_SORTS; i < N_SORTS + N_POLYATOMIC_SORTS; ++i)
    {
        for (int j = N_SORTS; j < N_SORTS + N_POLYATOMIC_SORTS; ++j)
        {
            mBViscosity_vv[i][j] = 0.0;
        }
    }
    for (int i = N_SORTS; i < N_SORTS + N_POLYATOMIC_SORTS; ++i)
    {
        mBViscosity_vv[i][i] = bracket_.beta0011()[i - N_SORTS];
    }
    for (int j = 0; j < N_SORTS; ++j)
    {
        mBViscosity_vv[0][j] = x[j] * 3.0 / 2.0 * K_BOLTZMANN * nTot / rho;
    }
    for (int j = N_SORTS; j < N_SORTS + N_POLYATOMIC_SORTS; ++j)
    {
        mBViscosity_vv[0][j] = x[j - N_SORTS] * K_BOLTZMANN /
                Mixture::mass(j - N_SORTS);
    }

    // Заполняем матрицу свободных членов
    cu = K_BOLTZMANN * nTot / rho * (3.0 / 2.0 + x[0]);
    cuT = K_BOLTZMANN * nTot / rho * x[0];
    for (int i = 0; i < N_SORTS; ++i)
    {
        bBViscosity_vv[i][0] = -x[i] * cuT / cu;
    }
    for (int i = N_SORTS; i < N_SORTS + N_POLYATOMIC_SORTS; ++i)
    {
        bBViscosity_vv[i][0] = x[i - N_SORTS] * nTot / rho * K_BOLTZMANN / cu;
    }
    bBViscosity_vv[0][0] = 0.0;
}

// Решение СЛУ методом Гаусса
void TransportCoefficientsDc::gauss(QVector<QVector<double>>& m,
                                    QVector<QVector<double>>& b)
{
    // Вспомогательные переменные
    double big = 0.0;
    double temp = 0.0;
    double pivInv = 0.0;
    int iCol = 0;
    int iRow = 0;
    int mSize = m.size();
    int bSize = b[0].size();

    // Диагонализация матрицы системы
    for (int i = 0; i < mSize; ++i)
    {
        iPiv_v[i] = 0;
    }
    for (int i = 0; i < mSize; ++i)
    {
        big = 0.0;
        for (int j = 0; j < mSize; ++j)
        {
            if (iPiv_v[j] != 1)
            {
                for (int k = 0; k < mSize; ++k)
                {
                    if  (iPiv_v[k] == 0)
                    {
                        if (qAbs(m[j][k]) >= big)
                        {
                            big = qAbs(m[j][k]);
                            iRow = j;
                            iCol = k;
                        }
                    }
                    else if (iPiv_v[k] > 1)
                    {
                        throw std::invalid_argument("Singular matrix in "
                                                    "Gauss!");
                    }
                }
            }
        }
        iPiv_v[iCol] = iPiv_v[iCol] + 1;
        if (iRow != iCol)
        {
            for (int l = 0; l < mSize; ++l)
            {
                temp = m[iRow][l];
                m[iRow][l] = m[iCol][l];
                m[iCol][l] = temp;
            }
            for (int l = 0; l < bSize; ++l)
            {
                temp = b[iRow][l];
                b[iRow][l] = b[iCol][l];
                b[iCol][l] = temp;
            }
        }
        indexRow_v[i] = iRow;
        indexColumn_v[i] = iCol;
        if (m[iCol][iCol] == 0.0)
        {
            throw std::invalid_argument("Singular matrix in Gauss!");
        }
        pivInv = 1.0 / m[iCol][iCol];
        m[iCol][iCol] = 1.0;
        for (int l = 0; l < mSize; ++l)
        {
            m[iCol][l] = m[iCol][l] * pivInv;
        }
        for (int l = 0; l < bSize; ++l)
        {
            b[iCol][l] = b[iCol][l] * pivInv;
        }
        for (int ll = 0; ll < mSize; ++ll)
        {
            if (ll != iCol)
            {
                temp = m[ll][iCol];
                m[ll][iCol] = 0.0;
                for (int l = 0; l < mSize; ++l)
                {
                    m[ll][l] = m[ll][l] - m[iCol][l] * temp;
                }
                for (int l = 0; l < bSize; ++l)
                {
                    b[ll][l] = b[ll][l] - b[iCol][l] * temp;
                }
            }
        }
    }

    // Восстановление матрицы m
    for (int l = 1; l <= mSize; ++l)
    {
        if (indexRow_v[mSize - l] != indexColumn_v[mSize - l])
        {
            for (int k = 0; k < mSize; ++k)
            {
                temp = m[k][indexRow_v[mSize - l]];
                m[k][indexRow_v[mSize - l]] = m[k][indexColumn_v[mSize - l]];
                m[k][indexColumn_v[mSize - l]] = temp;
            }
        }
    }
}

// Расчет всех коэффициентов переноса, базируясь на решении СЛУ
void TransportCoefficientsDc::computeTransport(const double& t,
                                               const double& nTot,
                                               const QVector<double>& x)
{
    // Коэффициенты термодиффузии
    for (int i = 0; i < N_SORTS; ++i)
    {
        tDiffusion_v[i] = -1.0 / 2.0 / nTot * bLambda_vv[i][0];
    }

    // Коэффициент теплопроводности (tr, rot, tr + rot, vibr)
    tLambda_s = 0.0;
    for (int i = N_SORTS; i < N_SORTS * 2; ++i)
    {
        tLambda_s = tLambda_s + 5.0 / 4.0 * K_BOLTZMANN * x[i - N_SORTS] *
                bLambda_vv[i][0];
    }
    rLambda_s = 3.0 * K_BOLTZMANN * t * x[0] / 16.0 /
            bracket_.lambdaInt()[0] * K_BOLTZMANN;
    cLambda_s = tLambda_s + rLambda_s;
    vLambda_s = 3.0 * K_BOLTZMANN * t * x[0] / 16.0 /
            bracket_.lambdaInt()[0] * K_BOLTZMANN * heat_.cv();

    // Коэффициенты диффузии
    for (int i = 0; i < N_SORTS; ++i)
    {
        for (int j = 0; j < N_SORTS; ++j)
        {
            diffusion_vv[i][j] = 1.0 / 2.0 / nTot * bDiffusion_vv[i][j];
        }
    }

    // Сдвиговая и объемная вязкости
    sViscosity_s = 0.0;
    for (int i = 0; i < N_SORTS; ++i)
    {
        sViscosity_s = sViscosity_s + K_BOLTZMANN * t / 2.0 *
                bSViscosity_vv[i][0] * x[i];
    }
    bViscosity_s = 0.0;
    for (int i = 0; i < N_SORTS; ++i)
    {
        bViscosity_s = bViscosity_s - K_BOLTZMANN * t *
                bBViscosity_vv[i][0] * x[i];
    }
}

// Выделение памяти, обновление значений
void TransportCoefficientsDc::initialize()
{
    // Предварительная инициализация
    TransportCoefficients::initialize();
    heat_.initialize();
    bracket_.initialize();

    // Инициализация (строки)
    mLambda_vv.resize(N_SORTS * 2);
    mDiffusion_vv.resize(N_SORTS);
    mSViscosity_vv.resize(N_SORTS);
    mBViscosity_vv.resize(N_SORTS + N_POLYATOMIC_SORTS);
    bLambda_vv.resize(N_SORTS * 2);
    bDiffusion_vv.resize(N_SORTS);
    bSViscosity_vv.resize(N_SORTS);
    bBViscosity_vv.resize(N_SORTS + N_POLYATOMIC_SORTS);

    // Инициализация (столбцы)
    for (int i = 0; i < N_SORTS * 2; ++i)
    {
        mLambda_vv[i].fill(0.0, N_SORTS * 2);
        bLambda_vv[i].fill(0.0, 1);
    }
    for (int i = 0; i < N_SORTS; ++i)
    {
        mDiffusion_vv[i].fill(0.0, N_SORTS);
        mSViscosity_vv[i].fill(0.0, N_SORTS);
        bDiffusion_vv[i].fill(0.0, N_SORTS);
        bSViscosity_vv[i].fill(0.0, 1);
    }
    for (int i = 0; i < N_SORTS + N_POLYATOMIC_SORTS; ++i)
    {
        mBViscosity_vv[i].fill(0.0, N_SORTS + N_POLYATOMIC_SORTS);
        bBViscosity_vv[i].fill(0.0, 1);
    }

    // Инициализация вспомогательных полей
    indexColumn_v.fill(0.0, N_MAX);
    indexRow_v.fill(0.0, N_MAX);
    iPiv_v.fill(0.0, N_MAX);
}

// Расчет всех значений
void TransportCoefficientsDc::compute(const MacroParam& param)
{
    // Мольная доля
    QVector<double> n = {param.rho[0] / Mixture::mass(0),
                         param.rho[1] / Mixture::mass(1)};
    QVector<double> x = {n[0] / (n[0] + n[1]),
                         n[1] / (n[0] + n[1])};

    // Вспомогательные переменные
    double nTot = n[0] + n[1];
    double rho = param.rho[0] + param.rho[1];

    // Инициализация матриц систем, порядок имеет значение
    heat_.compute(param.tv);
    bracket_.compute(param);
    fillBMLambda(nTot, rho, x);
    fillBMDiffusion(nTot, rho, x);
    fillBMSViscosity(param.t, x);
    fillBMBViscosity(nTot, rho, x);

    // Решение систем методом Гаусса
    gauss(mLambda_vv, bLambda_vv);
    gauss(mDiffusion_vv, bDiffusion_vv);
    gauss(mSViscosity_vv, bSViscosity_vv);
    gauss(mBViscosity_vv, bBViscosity_vv);

    // Вычисление коэффициентов переноса
    computeTransport(param.t, nTot, x);
}
const SpecificHeat& TransportCoefficientsDc::heat() const
{
    return heat_;
}
const BracketIntegrals& TransportCoefficientsDc::bracket() const
{
    return bracket_;
}

///////////////////////////////////////////////////////////////////////////////
/// class FlowMembers
///////////////////////////////////////////////////////////////////////////////

// Выделение памяти, обновление значений
void FlowMembers::initialize()
{
    trQ_s = 0.0;
    vQ_s = 0.0;
    diffQ_s = 0.0;
    fullQ_s = 0.0;
    xxP_s = 0.0;
    h_v.fill(0.0, N_SORTS);
    diffV_v.fill(0.0, N_SORTS);
    d_v.fill(0.0, N_SORTS);
    flow_v.fill(0.0, SYSTEM_ORDER);
}

// Доступ к соответствующим полям
const double &FlowMembers::trQ() const
{
    return trQ_s;
}
const double &FlowMembers::vQ() const
{
    return vQ_s;
}
const double &FlowMembers::diffQ() const
{
    return diffQ_s;
}
const double &FlowMembers::fullQ() const
{
    return fullQ_s;
}
const double &FlowMembers::xxP() const
{
    return xxP_s;
}
const QVector<double> &FlowMembers::h() const
{
    return h_v;
}
const QVector<double> &FlowMembers::diffV() const
{
    return diffV_v;
}
const QVector<double> &FlowMembers::d() const
{
    return d_v;
}
const QVector<double> &FlowMembers::flow() const
{
    return flow_v;
}

///////////////////////////////////////////////////////////////////////////////
/// class FlowMembersDc
///////////////////////////////////////////////////////////////////////////////

// Расчет соответствующих потоков энергии, энтальпии, компоненты тензора
void FlowMembersDc::computeD(const MacroParam& param,
                             const QVector<double>& dx_dx,
                             const double& dlnp_dx)
{
    QVector<double> n = {param.rho[0] / Mixture::mass(0),
                         param.rho[1] / Mixture::mass(1)};
    for (int i = 0; i < N_SORTS; ++i)
    {
        d_v[i] = dx_dx[i] + (n[i] / (n[0] + n[1]) - param.rho[i] /
                (param.rho[0] + param.rho[1])) * dlnp_dx;
    }
}
void FlowMembersDc::computeDiffV(const double& dlnT_dx)
{
    for (int i = 0; i < N_SORTS; ++i)
    {
        diffV_v[i] -= transport_.tDiffusion()[i] * dlnT_dx;
        for (int j = 0; j < N_SORTS; ++j)
        {
            diffV_v[i] -= transport_.diffusion()[i][j] * d_v[j];
        }
    }
}
void FlowMembersDc::computeH(const MacroParam& param)
{
    h_v[0] = 2.5 * K_BOLTZMANN * param.t / Mixture::mass(0) + energy_.rE() +
            energy_.vE();
    h_v[1] = 2.5 * K_BOLTZMANN * param.t / Mixture::mass(1);
}
void FlowMembersDc::computeDiffQ(const MacroParam& param)
{
    diffQ_s = 0.0;
    for (int i = 0; i < N_SORTS; ++i)
    {
        diffQ_s += param.rho[i] * h_v[i] * diffV_v[i] -
                param.p * d_v[i] * transport_.tDiffusion()[i];
    }
}
void FlowMembersDc::computeFullQ(const double& dT_dx, const double& dTv_dx)
{
    trQ_s = -transport_.cLambda() * dT_dx;
    vQ_s = -transport_.vLambda() * dTv_dx;
    fullQ_s = trQ_s + vQ_s + diffQ_s;
}
void FlowMembersDc::computeXxP(const MacroParam& param, const double& dv_dx)
{
    xxP_s = param.p - (4.0 / 3.0 * transport_.sViscosity() +
                       USE_B_VISC * transport_.bViscosity()) * dv_dx;
}

// Выделение памяти, обновление значений
void FlowMembersDc::initialize()
{
    // Предварительная инициализация
    FlowMembers::initialize();
    energy_.initialize();
    transport_.initialize();
}

// Расчет всех значений
void FlowMembersDc::compute(const MacroParam& param,
                            const QVector<double>& dx_dx, const double& dlnp_dx,
                            const double& dT_dx, const double& dlnT_dx,
                            const double& dTv_dx, const double& dv_dx)
{
    // Подготовка данных
    energy_.compute(param);
    transport_.compute(param);
    computeD(param, dx_dx, dlnp_dx);
    computeDiffV(dlnT_dx);
    computeH(param);
    computeDiffQ(param);
    computeFullQ(dT_dx, dTv_dx);
    computeXxP(param, dv_dx);

    // Расчет вектора поточных членов
    flow_v[0] = param.rho[0] * param.v + USE_V_DIFF * param.rho[0] * diffV_v[0];
    flow_v[1] = param.rho[1] * param.v - USE_V_DIFF * param.rho[0] * diffV_v[0];
    flow_v[2] = (param.rho[0] + param.rho[1]) * qPow(param.v, 2.0) + xxP_s;
    flow_v[3] = param.v * (param.rho[0] + param.rho[1]) *
            (energy_.fullE() + qPow(param.v, 2.0) / 2.0) + fullQ_s +
            param.v * xxP_s;
    flow_v[4] = param.v * param.rho[0] * energy_.vE() + vQ_s;
}
const Energy& FlowMembersDc::energy() const
{
    return energy_;
}
const TransportCoefficients& FlowMembersDc::transport() const
{
    return transport_;
}
