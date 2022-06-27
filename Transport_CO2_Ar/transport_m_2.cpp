///////////////////////////////////////////////////////////////////////////////
/// Автор : Баталов Семен, 2022
/// Модуль предназаначен для расчета коэффициентов переноса для смеси газов
/// CO2-Ar, а именно: коэффициентов диффузии, термодиффузии, объемной и
/// сдвиговой вязкостей, коэффициентов теплопроводности.
///////////////////////////////////////////////////////////////////////////////

#include "transport_m_2.h"

///////////////////////////////////////////////////////////////////////////////
/// class MixtureData
///////////////////////////////////////////////////////////////////////////////

/// private ///////////////////////////////////////////////////////////////////

// Расчет колебательных энергий
void tc_2::MixtureData::computeVEnergy()
{
    for (size_t i = 0; i < N_VIBR_L1; ++i)
    {
        for (size_t j = 0; j < N_VIBR_L2; ++j)
        {
            for (size_t k = 0; k < N_VIBR_L3; ++k)
            {
                vEnergy_vvv[i][j][k] = w_v[0] * i + w_v[1] * j + w_v[2] * k;
            }
        }
    }
}

// Расчет приведенных масс
void tc_2::MixtureData::computeReducedMass()
{
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        for (size_t j = 0; j < N_SORTS; ++j)
        {
            reducedMass_vv[i][j] = mass_v[i] * mass_v[j] /
                    (mass_v[i] + mass_v[j]);
        }
    }
}

/// public ////////////////////////////////////////////////////////////////////

// Инициализация
void tc_2::MixtureData::initialize()
{
    // Подгонка размеров строк
    w_v.resize(N_VIBR_DEGREES);
    mass_v.resize(N_SORTS);
    sigma_v.resize(N_SORTS);
    epsilon_v.resize(N_SORTS);

    // Подгонка размеров матрицы энергий
    vEnergy_vvv.resize(N_VIBR_L1);
    for (size_t i = 0; i < N_VIBR_L1; ++i)
    {
        vEnergy_vvv[i].resize(N_VIBR_L2);
        for (size_t j = 0; j < N_VIBR_L2; ++j)
        {
            vEnergy_vvv[i][j].resize(N_VIBR_L3);
        }
    }

    // Подгонка размеров матрицы приведенных масс
    reducedMass_vv.resize(N_SORTS);
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        reducedMass_vv[i].resize(N_SORTS);
    }

    // Спектроскопические данные
    w_v[0] = 1345.04 * W_W;
    w_v[1] = 667.25 * W_W;
    w_v[2] = 2361.71 * W_W;
    dEnergy_s = 64017;

    // Массы сортов
    mass_v[0] = 44.01 * ATOMIC_MASS_UNIT;
    mass_v[1] = 39.948 * ATOMIC_MASS_UNIT;

    // Газокинетические диаметры
    sigma_v[0] = 3.763e-10;
    sigma_v[1] = 3.35e-10;

    // Глубина потенциальной ямы
    epsilon_v[0] = 244.0;
    epsilon_v[1] = 141.5;

    // Колебательная энергия и приведенные массы
    computeVEnergy();
    computeReducedMass();
}

// Методы доступа к полям
const double& tc_2::MixtureData::dEnergy() const
{
    return dEnergy_s;
}
const QVector<double>& tc_2::MixtureData::mass() const
{
    return mass_v;
}
const QVector<double>& tc_2::MixtureData::sigma() const
{
    return sigma_v;
}
const QVector<double>& tc_2::MixtureData::epsilon() const
{
    return epsilon_v;
}
const QVector<QVector<QVector<double>>>& tc_2::MixtureData::vEnergy() const
{
    return vEnergy_vvv;
}
const QVector<QVector<double>>& tc_2::MixtureData::reducedMass() const
{
    return reducedMass_vv;
}

///////////////////////////////////////////////////////////////////////////////
/// class OmegaIntegrals
///////////////////////////////////////////////////////////////////////////////

/// public ////////////////////////////////////////////////////////////////////

// Инициализация
void tc_2::OmegaIntegrals::initialize()
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
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        omega11_vv[i].resize(N_SORTS);
        omega12_vv[i].resize(N_SORTS);
        omega13_vv[i].resize(N_SORTS);
        omega22_vv[i].resize(N_SORTS);
        aa_vv[i].resize(N_SORTS);
        bb_vv[i].resize(N_SORTS);
        cc_vv[i].resize(N_SORTS);

        // Обнуление значений
        for (size_t j = 0; j < N_SORTS; ++j)
        {
            omega11_vv[i][j] = 0.0;
            omega12_vv[i][j] = 0.0;
            omega13_vv[i][j] = 0.0;
            omega22_vv[i][j] = 0.0;
            aa_vv[i][j] = 0.0;
            bb_vv[i][j] = 0.0;
            cc_vv[i][j] = 0.0;
        }
    }
}

// Доступ к соответствующим полям
const QVector<QVector<double>>& tc_2::OmegaIntegrals::omega11() const
{
    return omega11_vv;
}
const QVector<QVector<double>>& tc_2::OmegaIntegrals::omega12() const
{
    return omega12_vv;
}
const QVector<QVector<double>>& tc_2::OmegaIntegrals::omega13() const
{
    return omega13_vv;
}
const QVector<QVector<double>>& tc_2::OmegaIntegrals::omega22() const
{
    return omega22_vv;
}
const QVector<QVector<double>>& tc_2::OmegaIntegrals::aa() const
{
    return aa_vv;
}
const QVector<QVector<double>>& tc_2::OmegaIntegrals::bb() const
{
    return bb_vv;
}
const QVector<QVector<double>>& tc_2::OmegaIntegrals::cc() const
{
    return cc_vv;
}

///////////////////////////////////////////////////////////////////////////////
/// class OmegaIntegralsDc
///////////////////////////////////////////////////////////////////////////////

/// public ////////////////////////////////////////////////////////////////////

// Инициализация
void tc_2::OmegaIntegralsDc::initialize()
{
    // Предварительная инициализация
    tc_2::OmegaIntegrals::initialize();

    // nullptr
    mixture_p = nullptr;
}

// Привязка к внешним классам
void tc_2::OmegaIntegralsDc::link(const MixtureData& mixture)
{
    mixture_p = &mixture;
}

///////////////////////////////////////////////////////////////////////////////
/// class OmegaIntegralsDcLjp
///////////////////////////////////////////////////////////////////////////////

/// public ////////////////////////////////////////////////////////////////////

// Инициализация
void tc_2::OmegaIntegralsDcLjp::initialize()
{
    // Предварительная инициализация
    tc_2::OmegaIntegralsDc::initialize();

    // Обнуление
    sig_ij = 0.0;
    e_ij = 0.0;
    s_ij = 0.0;
    tx = 0.0;
    x11 = 0.0;
}

// Расчет всех значений омега-интегралов
void tc_2::OmegaIntegralsDcLjp::compute(const double& t)
{
    // Проход по возможным сочетаниям сортов
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        for (size_t j = 0; j < N_SORTS; ++j)
        {
            // Инициализация вспомогательных переменных
            sig_ij = (mixture_p->sigma()[i] + mixture_p->sigma()[j]) / 2.0;
            e_ij = qSqrt(mixture_p->epsilon()[i] * mixture_p->epsilon()[j] *
                        qPow(mixture_p->sigma()[i] * 1e10, 6.0) *
                        qPow(mixture_p->sigma()[j] * 1e10, 6.0)) /
                        qPow(sig_ij * 1e10, 6.0);
            s_ij = M_PI * qPow(sig_ij, 2.0) *
                    qSqrt(K_BOLTZMANN * t / 2.0 / M_PI /
                          mixture_p->reducedMass()[i][j]);
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

/// public ////////////////////////////////////////////////////////////////////

// Инициализация
void tc_2::SpecificHeat::initialize()
{
    zvT12_s = 0.0;
    zvT3_s = 0.0;
    cvT12_s = 0.0;
    cvT3_s = 0.0;
}

// Доступ к соответствующим полям
const double& tc_2::SpecificHeat::zvT12() const
{
    return zvT12_s;
}
const double& tc_2::SpecificHeat::zvT3() const
{
    return zvT3_s;
}
const double& tc_2::SpecificHeat::cvT12() const
{
    return cvT12_s;
}
const double& tc_2::SpecificHeat::cvT3() const
{
    return cvT3_s;
}

///////////////////////////////////////////////////////////////////////////////
/// class SpecificHeatDc
///////////////////////////////////////////////////////////////////////////////

/// protected /////////////////////////////////////////////////////////////////

// Расчет статистических сумм и удельных теплоемкостей
void tc_2::SpecificHeatDc::computeZv(const double& t12, const double& t3)
{
    // Вспомогательные переменные
    size_t g = 0;
    double ee = 0.0;
    double ee1 = mixture_p->vEnergy()[1][0][0];
    double ee2 = mixture_p->vEnergy()[0][1][0];
    double ee3 = mixture_p->vEnergy()[0][0][1];

    // Инициализация
    zvT12_s = 0;
    zvT3_s = 0;

    // Проходим по всем колебательным уровням, мода 12
    for (size_t i = 0; i < N_VIBR_L1; ++i)
    {
        for (size_t j = 0; j < N_VIBR_L2; ++j)
        {
            g = j + 1;
            ee = mixture_p->vEnergy()[i][j][0];
            if (ee < mixture_p->dEnergy() * K_BOLTZMANN)
            {
                zvT12_s = zvT12_s + g * qExp(-((i * ee1 + j * ee2) /
                                             K_BOLTZMANN / t12));
            }
        }
    }

    // Проходим по всем колебательным уровням, мода 3
    for (size_t k = 0; k < N_VIBR_L3; ++k)
    {
        ee = mixture_p->vEnergy()[0][0][k];
        if (ee < mixture_p->dEnergy() * K_BOLTZMANN)
        {
            zvT3_s = zvT3_s + qExp(-(k * ee3 / K_BOLTZMANN / t3));
        }
    }
}
void tc_2::SpecificHeatDc::computeCv(const double& t12, const double& t3)
{
    // Вспомогательные переменные
    size_t g = 0;
    double ee = 0.0;
    double ee1 = mixture_p->vEnergy()[1][0][0];
    double ee2 = mixture_p->vEnergy()[0][1][0];
    double ee3 = mixture_p->vEnergy()[0][0][1];
    double s12 = 0.0;
    double s3 = 0.0;
    double ss12 = 0.0;
    double ss3 = 0.0;
    double ppp = 0.0;

    // Проходим по всем колебательным уровням, мода 12
    for (size_t i = 0; i < N_VIBR_L1; ++i)
    {
        for (size_t j = 0; j < N_VIBR_L2; ++j)
        {
            g = j + 1;
            ee = mixture_p->vEnergy()[i][j][0];
            if (ee < mixture_p->dEnergy() * K_BOLTZMANN)
            {
                ppp = g * qExp(-(i * ee1 + j * ee2) / K_BOLTZMANN / t12) /
                        zvT12_s;
                s12 = s12 + (i * ee1 + j * ee2) / K_BOLTZMANN / t12 * ppp;
                ss12 = ss12 + qPow((i * ee1 + j * ee2) / K_BOLTZMANN /
                                   t12, 2.0) * ppp;
            }
        }
    }

    // Проходим по всем колебательным уровням, мода 3
    for (size_t k = 0; k < N_VIBR_L3; ++k)
    {
        ee = mixture_p->vEnergy()[0][0][k];
        if (ee < mixture_p->dEnergy() * K_BOLTZMANN)
        {
            ppp = qExp(-(k * ee3 / K_BOLTZMANN / t3)) / zvT3_s;
            s3 = s3 + k * ee3 / K_BOLTZMANN / t3 * ppp;
            ss3 = ss3 + qPow(k * ee3 / K_BOLTZMANN / t3, 2.0) * ppp;
        }
    }

    // Вычисление удельных теплоемкостей
    cvT12_s = ss12 - qPow(s12, 2.0);
    cvT3_s = ss3 - qPow(s3, 2.0);
}

/// public ////////////////////////////////////////////////////////////////////

// Инициализация
void tc_2::SpecificHeatDc::initialize()
{
    // Предварительная инициализация
    tc_2::SpecificHeat::initialize();

    // nullptr
    mixture_p = nullptr;
}

// Привязка к внешним классам
void tc_2::SpecificHeatDc::link(const MixtureData& mixture)
{
    mixture_p = &mixture;
}

// Расчет всех значений
void tc_2::SpecificHeatDc::compute(const double& t12, const double& t3)
{
    // Порядок имеет значение
    computeZv(t12, t3);
    computeCv(t12, t3);
}

///////////////////////////////////////////////////////////////////////////////
/// class BracketIntegrals
///////////////////////////////////////////////////////////////////////////////

/// public ////////////////////////////////////////////////////////////////////

// Инициализация
void tc_2::BracketIntegrals::initialize()
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
    lambdaInt_v.resize(N_SORTS);

    // Подгонка размеров строк матриц омега-интегралов
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        lambda_vv[i].resize(N_SORTS);
        lambda00_vv[i].resize(N_SORTS);
        lambda01_vv[i].resize(N_SORTS);
        lambda11_vv[i].resize(N_SORTS);
        eta_vv[i].resize(N_SORTS);
        h00_vv[i].resize(N_SORTS);
        beta11_vv[i].resize(N_SORTS);
        beta01_vv[i].resize(N_POLYATOMIC_SORTS);

        // Обнуление значений
        for (size_t j = 0; j < N_SORTS; ++j)
        {
            lambda_vv[i][j] = 0.0;
            lambda00_vv[i][j] = 0.0;
            lambda01_vv[i][j] = 0.0;
            lambda11_vv[i][j] = 0.0;
            eta_vv[i][j] = 0.0;
            h00_vv[i][j] = 0.0;
            beta11_vv[i][j] = 0.0;
        }
        for (size_t j = 0; j < N_POLYATOMIC_SORTS; ++j)
        {
            beta01_vv[i][j] = 0.0;
            beta0011_v[j] = 0.0;
        }
        lambdaInt_v[i] = 0.0;
    }
}

// Доступ к соответствующим полям
const QVector<QVector<double>>& tc_2::BracketIntegrals::lambda() const
{
    return lambda_vv;
}
const QVector<QVector<double>>& tc_2::BracketIntegrals::lambda00() const
{
    return lambda00_vv;
}
const QVector<QVector<double>>& tc_2::BracketIntegrals::lambda01() const
{
    return lambda01_vv;
}
const QVector<QVector<double>>& tc_2::BracketIntegrals::lambda11() const
{
    return lambda11_vv;
}
const QVector<QVector<double>>& tc_2::BracketIntegrals::eta() const
{
    return eta_vv;
}
const QVector<QVector<double>>& tc_2::BracketIntegrals::h00() const
{
    return h00_vv;
}
const QVector<QVector<double>>& tc_2::BracketIntegrals::beta11() const
{
    return beta11_vv;
}
const QVector<QVector<double>>& tc_2::BracketIntegrals::beta01() const
{
    return beta01_vv;
}
const QVector<double>& tc_2::BracketIntegrals::beta0011() const
{
    return beta0011_v;
}
const QVector<double>& tc_2::BracketIntegrals::lambdaInt() const
{
    return lambdaInt_v;
}

///////////////////////////////////////////////////////////////////////////////
/// class BracketIntegralsDc
///////////////////////////////////////////////////////////////////////////////

/// protected /////////////////////////////////////////////////////////////////

// Расчет phi
void tc_2::BracketIntegralsDc::computePhi(const double& t)
{
    phi_v[0] = K_BOLTZMANN / 23.73 *
            (1.0 + qPow(M_PI, 1.5) / 2.0 * qSqrt(mixture_p->epsilon()[0] / t) +
            (qPow(M_PI, 2.0) / 4.0 + 2.0) * mixture_p->epsilon()[0] / t +
            qPow(M_PI * mixture_p->epsilon()[0] / t, 1.5));
}

// Расчет всех интегральных скобок
void tc_2::BracketIntegralsDc::computeBrackets1(const double& t)
{
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        for (size_t j = 0; j < N_SORTS; ++j)
        {
            lambda_vv[i][j] = 75.0 / 64.0 * t * K_BOLTZMANN /
                    mixture_p->reducedMass()[i][j] / omega_p->omega22()[i][j] *
                    K_BOLTZMANN;
            eta_vv[i][j] = 5.0 / 8.0 * K_BOLTZMANN * t /
                    omega_p->omega22()[i][j];
        }
    }
}
void tc_2::BracketIntegralsDc::computeBrackets2(const double& t,
                                                const QVector<double>& x)
{
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        lambdaInt_v[i] = 0.0;
    }
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        for (size_t j = 0; j < N_SORTS; ++j)
        {
            if (i == j)
            {
                lambda00_vv[i][j] = 0.0;
                lambda01_vv[i][j] = 0.0;
                lambda11_vv[i][j] = qPow(x[i], 2.0) / lambda_vv[i][i];
                h00_vv[i][j] = qPow(x[i], 2.0) / eta_vv[i][i];
                beta11_vv[i][j] = 4.0 * t / M_PI * qPow(x[i], 2.0) /
                        eta_vv[i][j] * phi_v[i];
                for (size_t k = 0; k < N_SORTS; ++k)
                {
                    if (k != i)
                    {
                        lambda00_vv[i][j] = lambda00_vv[i][j] + x[i] * x[k] /
                                lambda_vv[i][k] / 2.0 / omega_p->aa()[i][k];
                        lambda01_vv[i][j] = lambda01_vv[i][j] - x[i] * x[k] /
                                lambda_vv[i][k] / 4.0 / omega_p->aa()[i][k] *
                                mixture_p->mass()[k] /
                                (mixture_p->mass()[i] + mixture_p->mass()[k]) *
                                (6.0 * omega_p->cc()[i][k] - 5.0);
                        lambda11_vv[i][j] = lambda11_vv[i][j] + x[i] * x[k] /
                                lambda_vv[i][k] / 2.0 / omega_p->aa()[i][k] *
                                (15.0 / 2.0 * qPow(mixture_p->mass()[i], 2.0) +
                                 25.0 / 4.0 * qPow(mixture_p->mass()[k], 2.0) -
                                 3.0 * qPow(mixture_p->mass()[k], 2.0) *
                                 omega_p->bb()[i][k] + 4.0 *
                                 mixture_p->mass()[i] * mixture_p->mass()[k] *
                                 omega_p->aa()[i][k]) /
                                qPow(mixture_p->mass()[i] +
                                     mixture_p->mass()[k], 2.0);
                        h00_vv[i][j] = h00_vv[i][j] + 2.0 * x[i] * x[k] /
                                eta_vv[i][k] * mixture_p->mass()[i] *
                                mixture_p->mass()[k] /
                                qPow(mixture_p->mass()[i] +
                                     mixture_p->mass()[k], 2.0) *
                                (5.0 / 3.0 / omega_p->aa()[i][k] +
                                 mixture_p->mass()[k] / mixture_p->mass()[i]);
                        beta11_vv[i][j] = beta11_vv[i][j] + x[i] * x[k] /
                                eta_vv[i][k] * mixture_p->mass()[k] /
                                qPow(mixture_p->mass()[i] +
                                     mixture_p->mass()[k], 2.0) *
                                (5.0 * K_BOLTZMANN * t * mixture_p->mass()[i] /
                                 omega_p->aa()[i][k] + 4.0 * t *
                                 mixture_p->mass()[k] / M_PI *
                                 (phi_v[i] + phi_v[k]));
                    }
                }
            }
            else
            {
                lambda00_vv[i][j] = -x[i] * x[j] / lambda_vv[i][j] / 2.0 /
                        omega_p->aa()[i][j];
                lambda01_vv[i][j] = x[i] * x[j] / lambda_vv[i][j] / 4.0 /
                        omega_p->aa()[i][j] * mixture_p->mass()[i] /
                        (mixture_p->mass()[i] + mixture_p->mass()[j]) *
                        (6.0 * omega_p->cc()[i][j] - 5.0);
                lambda11_vv[i][j] = -x[i] * x[j] / lambda_vv[i][j] / 2.0 /
                        omega_p->aa()[i][j] * mixture_p->reducedMass()[i][j] /
                        (mixture_p->mass()[i] + mixture_p->mass()[j]) *
                        (55.0 / 4.0 - 3.0 * omega_p->bb()[i][j] - 4.0 *
                         omega_p->aa()[i][j]);
                h00_vv[i][j] = -2.0 * x[i] * x[j] / eta_vv[i][j] *
                        mixture_p->reducedMass()[i][j] /
                        (mixture_p->mass()[i] + mixture_p->mass()[j]) *
                        (5.0 / 3.0 / omega_p->aa()[i][j] - 1.0);
                beta11_vv[i][j] = x[i] * x[j] / eta_vv[i][j] *
                        mixture_p->mass()[i] * mixture_p->mass()[j] /
                        qPow(mixture_p->mass()[i] +
                             mixture_p->mass()[j], 2.0) *
                        (-5.0 * K_BOLTZMANN * t / omega_p->aa()[i][j] +
                         4.0 * t / M_PI * (phi_v[i] + phi_v[j]));
            }
            lambdaInt_v[i] += x[j] * mixture_p->reducedMass()[i][j] *
                    omega_p->omega11()[i][j];
        }
    }
}
void tc_2::BracketIntegralsDc::computeBrackets3(const double& t,
                                                const QVector<double>& x)
{
    for (size_t i = 0; i < N_POLYATOMIC_SORTS; ++i)
    {
        beta0011_v[i] = 0.0;
    }
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        for (size_t j = 0; j < N_POLYATOMIC_SORTS; ++j)
        {
            if (i == j)
            {
                beta01_vv[i][j] = -4.0 * t / M_PI * qPow(x[i], 2.0) /
                        eta_vv[i][i] * phi_v[i];
                for (size_t k = 0; k < N_SORTS; ++k)
                {
                    if (k != i)
                    {
                        beta01_vv[i][j] = beta01_vv[i][j] -
                                4.0 * t / M_PI * x[i] * x[k] / eta_vv[i][k] *
                                mixture_p->mass()[k] /
                                (mixture_p->mass()[i] + mixture_p->mass()[k]) *
                                phi_v[i];
                    }
                }
            }
            else
            {
                beta01_vv[i][j] = -4.0 * t / M_PI * x[i] * x[j] /
                        eta_vv[i][j] * mixture_p->mass()[j] /
                        (mixture_p->mass()[i] + mixture_p->mass()[j]) *
                        phi_v[j];
            }
            beta0011_v[j] += 4.0 * t / M_PI * x[i] * phi_v[j] * x[j] /
                    eta_vv[i][j];
        }
    }
}

/// public ////////////////////////////////////////////////////////////////////

// Инициализация
void tc_2::BracketIntegralsDc::initialize()
{
    // Предварительная инициализация
    tc_2::BracketIntegrals::initialize();

    // nullptr
    mixture_p = nullptr;
    omega_p = nullptr;

    // Подготовка массива Phi
    phi_v.resize(N_SORTS);
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        phi_v[i] = 0;
    }
}

// Привязка к внешним классам
void tc_2::BracketIntegralsDc::link(const MixtureData& mixture)
{
    mixture_p = &mixture;
}
void tc_2::BracketIntegralsDc::link(const OmegaIntegrals& omega)
{
    omega_p = &omega;
}

// Расчет всех значений
void tc_2::BracketIntegralsDc::compute(const double& t,
                                       const QVector<double>& x)
{
    // Порядок имеет значение
    computePhi(t);
    computeBrackets1(t);
    computeBrackets2(t, x);
    computeBrackets3(t, x);
}

///////////////////////////////////////////////////////////////////////////////
/// class TransportCoefficients
///////////////////////////////////////////////////////////////////////////////

/// protected /////////////////////////////////////////////////////////////////

/// public ////////////////////////////////////////////////////////////////////

// Выделение памяти, обновление значений
void tc_2::TransportCoefficients::initialize()
{
    // Обнуление значений
    sViscosity_s = 0.0;
    bViscosity_s = 0.0;
    tLambda_s = 0.0;
    rLambda_s = 0.0;
    cLambda_s = 0.0;
    vLambdaT12_s = 0.0;
    vLambdaT3_s = 0.0;

    // Выделение памяти и обнуление значений
    diffusion_vv.resize(N_SORTS);
    tDiffusion_v.resize(N_SORTS);
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        diffusion_vv[i].resize(N_SORTS);
        for (size_t j = 0; j < N_SORTS; ++j)
        {
            diffusion_vv[i][j] = 0.0;
        }
    }
}

// Доступ к соответствующим полям
const double& tc_2::TransportCoefficients::sViscosity() const
{
    return sViscosity_s;
}
const double& tc_2::TransportCoefficients::bViscosity() const
{
    return bViscosity_s;
}
const double& tc_2::TransportCoefficients::tLambda() const
{
    return tLambda_s;
}
const double& tc_2::TransportCoefficients::rLambda() const
{
    return rLambda_s;
}
const double& tc_2::TransportCoefficients::cLambda() const
{
    return cLambda_s;
}
const double& tc_2::TransportCoefficients::vLambdaT12() const
{
    return vLambdaT12_s;
}
const double& tc_2::TransportCoefficients::vLambdaT3() const
{
    return vLambdaT3_s;
}
const QVector<QVector<double>>& tc_2::TransportCoefficients::diffusion() const
{
    return diffusion_vv;
}
const QVector<double>& tc_2::TransportCoefficients::tDiffusion() const
{
    return tDiffusion_v;
}

///////////////////////////////////////////////////////////////////////////////
/// class TransportCoefficientsDc
///////////////////////////////////////////////////////////////////////////////

/// protected /////////////////////////////////////////////////////////////////

// Заполнение матриц систем и матриц свободных членов
void tc_2::TransportCoefficientsDc::fillBMLambda(const double& nTot,
                                                 const double& rho,
                                                 const QVector<double>& x)
{
    // Заполняем матрицу системы
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        for (size_t j = 0; j < N_SORTS; ++j)
        {
            mLambda_vv[i][j] = bracket_p->lambda00()[i][j];
        }
    }
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        for (size_t j = N_SORTS; j < N_SORTS * 2; ++j)
        {
            mLambda_vv[i][j] = bracket_p->lambda01()[i][j - N_SORTS];
        }
    }
    for (size_t i = N_SORTS; i < N_SORTS * 2; ++i)
    {
        for (size_t j = 0; j < N_SORTS; ++j)
        {
            mLambda_vv[i][j] = mLambda_vv[j][i];
        }
    }
    for (size_t i = N_SORTS; i < N_SORTS * 2; ++i)
    {
        for (size_t j = N_SORTS; j < N_SORTS * 2; ++j)
        {
            mLambda_vv[i][j] = bracket_p->lambda11()[i - N_SORTS][j - N_SORTS];
        }
    }
    for (size_t j = 0; j < N_SORTS; ++j)
    {
        mLambda_vv[0][j] = x[j] * mixture_p->mass()[j] * nTot / rho;
    }
    for (size_t j = N_SORTS; j < N_SORTS * 2; ++j)
    {
        mLambda_vv[0][j] = 0.0;
    }

    // Заполняем матрицу свободных членов
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        bLambda_vv[i][0] = 0.0;
    }
    for (size_t i = N_SORTS; i < N_SORTS * 2; ++i)
    {
        bLambda_vv[i][0] = 4.0 / 5.0 / K_BOLTZMANN * x[i - N_SORTS];
    }
}
void tc_2::TransportCoefficientsDc::fillBMDiffusion(const double& nTot,
                                                    const double& rho,
                                                    const QVector<double>& x)
{
    // Вспомогательные переменные
    double delta = 0.0;

    // Заполняем матрицу системы
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        for (size_t j = 0; j < N_SORTS; ++j)
        {
            mDiffusion_vv[i][j] = mLambda_vv[i][j];
        }
    }

    // Заполняем матрицу свободных членов
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        for (size_t j = 0; j < N_SORTS; ++j)
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
                    (delta - mixture_p->mass()[i] * x[i] * nTot / rho);
        }
        bDiffusion_vv[0][i] = 0.0;
    }
}
void tc_2::TransportCoefficientsDc::fillBMSViscosity(const double& t,
                                                     const QVector<double>& x)
{
    // Заполняем матрицу системы
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        for (size_t j = 0; j < N_SORTS; ++j)
        {
            mSViscosity_vv[i][j] = bracket_p->h00()[i][j];
        }
    }

    // Заполняем матрицу свободных членов
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        bSViscosity_vv[i][0] = 2.0 / K_BOLTZMANN / t * x[i];
    }
}
void tc_2::TransportCoefficientsDc::fillBMBViscosity(const double& nTot,
                                                     const double& rho,
                                                     const QVector<double>& x)
{
    // Вспомогательные переменные
    double cu = 0.0;
    double cuT = 0.0;

    // Заполняем матрицу системы
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        for (size_t j = 0; j < N_SORTS; ++j)
        {
            mBViscosity_vv[i][j] = bracket_p->beta11()[i][j];
        }
    }
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        for (size_t j = N_SORTS; j < N_SORTS + N_POLYATOMIC_SORTS; ++j)
        {
            mBViscosity_vv[i][j] = bracket_p->beta01()[i][j - N_SORTS];
        }
    }
    for (size_t i = N_SORTS; i < N_SORTS + N_POLYATOMIC_SORTS; ++i)
    {
        for (size_t j = 0; j < N_SORTS; ++j)
        {
            mBViscosity_vv[i][j] = mBViscosity_vv[j][i];
        }
    }
    for (size_t i = N_SORTS; i < N_SORTS + N_POLYATOMIC_SORTS; ++i)
    {
        for (size_t j = N_SORTS; j < N_SORTS + N_POLYATOMIC_SORTS; ++j)
        {
            mBViscosity_vv[i][j] = 0.0;
        }
    }
    for (size_t i = N_SORTS; i < N_SORTS + N_POLYATOMIC_SORTS; ++i)
    {
        mBViscosity_vv[i][i] = bracket_p->beta0011()[i - N_SORTS];
    }
    for (size_t j = 0; j < N_SORTS; ++j)
    {
        mBViscosity_vv[0][j] = x[j] * 3.0 / 2.0 * K_BOLTZMANN * nTot / rho;
    }
    for (size_t j = N_SORTS; j < N_SORTS + N_POLYATOMIC_SORTS; ++j)
    {
        mBViscosity_vv[0][j] = x[j - N_SORTS] * K_BOLTZMANN /
                mixture_p->mass()[j - N_SORTS];
    }

    // Заполняем матрицу свободных членов
    cu = K_BOLTZMANN * nTot / rho * (3.0 / 2.0 + x[0]);
    cuT = K_BOLTZMANN * nTot / rho * x[0];
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        bBViscosity_vv[i][0] = -x[i] * cuT / cu;
    }
    for (size_t i = N_SORTS; i < N_SORTS + N_POLYATOMIC_SORTS; ++i)
    {
        bBViscosity_vv[i][0] = x[i - N_SORTS] * nTot / rho * K_BOLTZMANN / cu;
    }
    bBViscosity_vv[0][0] = 0.0;
}

// Решение СЛУ методом Гаусса
void tc_2::TransportCoefficientsDc::gauss(QVector<QVector<double>>& m,
                                          QVector<QVector<double>>& b)
{
    // Вспомогательные переменные
    double big = 0.0;
    double temp = 0.0;
    double pivInv = 0.0;
    size_t iCol = 0;
    size_t iRow = 0;
    size_t mSize = m.size();
    size_t bSize = b[0].size();

    // Диагонализация матрицы системы
    for (size_t i = 0; i < mSize; ++i)
    {
        iPiv_v[i] = 0;
    }
    for (size_t i = 0; i < mSize; ++i)
    {
        big = 0.0;
        for (size_t j = 0; j < mSize; ++j)
        {
            if (iPiv_v[j] != 1)
            {
                for (size_t k = 0; k < mSize; ++k)
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
            for (size_t l = 0; l < mSize; ++l)
            {
                temp = m[iRow][l];
                m[iRow][l] = m[iCol][l];
                m[iCol][l] = temp;
            }
            for (size_t l = 0; l < bSize; ++l)
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
        for (size_t l = 0; l < mSize; ++l)
        {
            m[iCol][l] = m[iCol][l] * pivInv;
        }
        for (size_t l = 0; l < bSize; ++l)
        {
            b[iCol][l] = b[iCol][l] * pivInv;
        }
        for (size_t ll = 0; ll < mSize; ++ll)
        {
            if (ll != iCol)
            {
                temp = m[ll][iCol];
                m[ll][iCol] = 0.0;
                for (size_t l = 0; l < mSize; ++l)
                {
                    m[ll][l] = m[ll][l] - m[iCol][l] * temp;
                }
                for (size_t l = 0; l < bSize; ++l)
                {
                    b[ll][l] = b[ll][l] - b[iCol][l] * temp;
                }
            }
        }
    }

    // Восстановление матрицы m
    for (size_t l = 1; l <= mSize; ++l)
    {
        if (indexRow_v[mSize - l] != indexColumn_v[mSize - l])
        {
            for (size_t k = 0; k < mSize; --k)
            {
                temp = m[k][indexRow_v[mSize - l]];
                m[k][indexRow_v[mSize - l]] = m[k][indexColumn_v[mSize - l]];
                m[k][indexColumn_v[mSize - l]] = temp;
            }
        }
    }
}

// Расчет всех коэффициентов переноса, базируясь на решении СЛУ
void tc_2::TransportCoefficientsDc::computeTransport(const double& t,
                                                     const double& nTot,
                                                     const QVector<double>& x)
{
    // Коэффициенты термодиффузии
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        tDiffusion_v[i] = -1.0 / 2.0 / nTot * bLambda_vv[i][0];
    }

    // Коэффициент теплопроводности (tr, rot, tr + rot, vibr)
    tLambda_s = 0.0;
    for (size_t i = N_SORTS; i < N_SORTS * 2; ++i)
    {
        tLambda_s = tLambda_s + 5.0 / 4.0 * K_BOLTZMANN * x[i - N_SORTS] *
                bLambda_vv[i][0];
    }
    rLambda_s = 3.0 * K_BOLTZMANN * t * x[0] / 16.0 /
            bracket_p->lambdaInt()[0] * K_BOLTZMANN;
    cLambda_s = tLambda_s + rLambda_s;
    vLambdaT12_s = 3.0 * K_BOLTZMANN * t * x[0] / 16.0 /
            bracket_p->lambdaInt()[0] * K_BOLTZMANN * heat_p->cvT12();
    vLambdaT3_s = 3.0 * K_BOLTZMANN * t * x[0] / 16.0 /
            bracket_p->lambdaInt()[0] * K_BOLTZMANN * heat_p->cvT3();

    // Коэффициенты диффузии
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        for (size_t j = 0; j < N_SORTS; ++j)
        {
            diffusion_vv[i][j] = 1.0 / 2.0 / nTot * bDiffusion_vv[i][j];
        }
    }

    // Сдвиговая и объемная вязкости
    sViscosity_s = 0.0;
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        sViscosity_s = sViscosity_s + K_BOLTZMANN * t / 2.0 *
                bSViscosity_vv[i][0] * x[i];
    }
    bViscosity_s = 0.0;
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        bViscosity_s = bViscosity_s - K_BOLTZMANN * t *
                bBViscosity_vv[i][0] * x[i];
    }
}

/// public ////////////////////////////////////////////////////////////////////

// Выделение памяти, обновление значений
void tc_2::TransportCoefficientsDc::initialize()
{
    // Предварительная инициализация
    tc_2::TransportCoefficients::initialize();

    // nullptr
    mixture_p = nullptr;
    bracket_p = nullptr;
    heat_p = nullptr;

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
    for (size_t i = 0; i < N_SORTS * 2; ++i)
    {
        mLambda_vv[i].resize(N_SORTS * 2);
        bLambda_vv[i].resize(1);
        for (size_t j = 0; j < N_SORTS * 2; ++j)
        {
            mLambda_vv[i][j] = 0.0;
        }
        bLambda_vv[i][0] = 0.0;
    }
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        mDiffusion_vv[i].resize(N_SORTS);
        mSViscosity_vv[i].resize(N_SORTS);
        bDiffusion_vv[i].resize(N_SORTS);
        bSViscosity_vv[i].resize(1);
        for (size_t j = 0; j < N_SORTS; ++j)
        {
            mDiffusion_vv[i][j] = 0.0;
            mSViscosity_vv[i][j] = 0.0;
            bDiffusion_vv[i][j] = 0.0;
        }
        bSViscosity_vv[i][0] = 0.0;
    }
    for (size_t i = 0; i < N_SORTS + N_POLYATOMIC_SORTS; ++i)
    {
        mBViscosity_vv[i].resize(N_SORTS + N_POLYATOMIC_SORTS);
        bBViscosity_vv[i].resize(1);
        for (size_t j = 0; j < N_SORTS + N_POLYATOMIC_SORTS; ++j)
        {
            mBViscosity_vv[i][j] = 0.0;
        }
        bBViscosity_vv[i][0] = 0.0;
    }

    // Инициализация вспомогательных полей
    indexColumn_v.resize(N_MAX);
    indexRow_v.resize(N_MAX);
    iPiv_v.resize(N_MAX);
    for (size_t i = 0; i < N_MAX; ++i)
    {
        indexColumn_v[i] = 0;
        indexRow_v[i] = 0;
        iPiv_v[i] = 0;
    }
}

// Привязка к внешним классам
void tc_2::TransportCoefficientsDc::link(const MixtureData& mixture)
{
    mixture_p = &mixture;
}
void tc_2::TransportCoefficientsDc::link(const BracketIntegrals& bracket)
{
    bracket_p = &bracket;
}
void tc_2::TransportCoefficientsDc::link(const SpecificHeat& heat)
{
    heat_p = &heat;
}

// Расчет всех значений
void tc_2::TransportCoefficientsDc::compute(const double& t,
                                            const QVector<double>& x,
                                            const double& p)
{
    // Вспомогательные переменные
    double nTot = 0.0;
    double rho = 0.0;

    // Подготовка вспомогательных величин
    nTot = p / K_BOLTZMANN / t;
    for (size_t i = 0; i < N_SORTS; ++i)
    {
        rho = rho + x[i] * mixture_p->mass()[i] * nTot;
    }

    // Инициализация матриц систем, порядок имеет значение
    fillBMLambda(nTot, rho, x);
    fillBMDiffusion(nTot, rho, x);
    fillBMSViscosity(t, x);
    fillBMBViscosity(nTot, rho, x);

    // Решение систем методом Гаусса
    gauss(mLambda_vv, bLambda_vv);
    gauss(mDiffusion_vv, bDiffusion_vv);
    gauss(mSViscosity_vv, bSViscosity_vv);
    gauss(mBViscosity_vv, bBViscosity_vv);

    // Вычисление коэффициентов переноса
    computeTransport(t, nTot, x);
}

///////////////////////////////////////////////////////////////////////////////
/// class DataDc
///////////////////////////////////////////////////////////////////////////////

/// public ////////////////////////////////////////////////////////////////////

// Выделение памяти, обновление значений
void tc_2::DataDc::initialize()
{
    mixture_.initialize();
}

// Доступ к данным смеси
const tc_2::MixtureData& tc_2::DataDc::mixture() const
{
    return mixture_;
}

///////////////////////////////////////////////////////////////////////////////
/// class Computer
///////////////////////////////////////////////////////////////////////////////

/// public ////////////////////////////////////////////////////////////////////

// Выделение памяти, обновление значений
void tc_2::ComputerDc::initialize()
{
    omega_.initialize();
    heat_.initialize();
    bracket_.initialize();
    transport_.initialize();
}

// Привязка к классу данных
void tc_2::ComputerDc::link(const DataDc& data)
{
    omega_.link(data.mixture());
    heat_.link(data.mixture());
    bracket_.link(data.mixture());
    bracket_.link(omega_);
    transport_.link(data.mixture());
    transport_.link(heat_);
    transport_.link(bracket_);
}

// Расчет всех коэффициентов
void tc_2::ComputerDc::compute(double t, double t12, double t3,
                               QVector<double> x, double p)
{
    // Порядок имеет значение
    omega_.compute(t);
    heat_.compute(t12, t3);
    bracket_.compute(t, x);
    transport_.compute(t, x, p);
}

// Доступ к коэффициентам переноса
const tc_2::OmegaIntegrals& tc_2::ComputerDc::omega() const
{
    return omega_;
}
const tc_2::SpecificHeat& tc_2::ComputerDc::heat() const
{
    return heat_;
}
const tc_2::BracketIntegrals& tc_2::ComputerDc::bracket() const
{
    return bracket_;
}
const tc_2::TransportCoefficients& tc_2::ComputerDc::transport() const
{
    return transport_;
}
