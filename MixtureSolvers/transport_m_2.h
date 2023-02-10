#ifndef TRANSPORT_M_2_H
#define TRANSPORT_M_2_H

///////////////////////////////////////////////////////////////////////////////
/// Автор : Баталов Семен, 2022
/// Модуль предназаначен для расчета коэффициентов переноса для смеси газов
/// CO2-Ar, а именно: коэффициентов диффузии, термодиффузии, объемной и
/// сдвиговой вязкостей, коэффициентов теплопроводности.
///////////////////////////////////////////////////////////////////////////////

#include <QtGlobal>
#include <QtMath>
#include <QVector>
#include <QDebug>
#include <exception>
#include "global.h"

// Transport Coefficients, mixture 2 -> CO2, Ar
namespace tc_2
{

// Постоянная Больцмана
const double K_BOLTZMANN = 1.3805e-23;

// Атомная единица массы
const double ATOMIC_MASS_UNIT = 1.6605402e-27;

// Постоянная Планка
const double H_PLANCK = 6.6254E-34;

// Число Авогадро
const double N_AVOGADRO = 6.0221e+23;

// Размерный коэффициент для расчета энергии
const double W_W = 1.60219e-19 / 8065.47;

// Кол-во сортов смеси
const size_t N_SORTS = 2;
const size_t N_POLYATOMIC_SORTS = 1;

// Характеристики внутренних степеней свободы CO2
const size_t N_VIBR_DEGREES = 3;
const size_t N_VIBR_L1 = 40;
const size_t N_VIBR_L2 = 40;
const size_t N_VIBR_L3 = 40;

// Настройка метода Гаусса
const size_t N_MAX = 30;

// Число уравнений системы
const size_t SYSTEM_ORDER = 6;

/// Структура данных смеси СO2-Ar:
/// 1) Спектроскопические данные для молекулы СО2;
/// 2) Массы молекул сортов СO2, Ar;
/// 3) Газокинетические диаметры молекул сортов СO2, Ar;
/// 4) Глубина потенциальной ямы для молекул сортов : СO2, Ar.
class MixtureData
{
private:

    // Энергия диссоциации : CO2
    double dEnergy_s;

    // Спектроскопические данные : CO2 -> (моды 1, 2, 3)
    QVector<double> w_v;

    // Массы молекул : СO2, Ar -> (1, 2)
    QVector<double> mass_v;

    // Газокинетические диаметры : СO2, Ar -> (1, 2)
    QVector<double> sigma_v;

    // Глубина потенциальной ямы : СO2, Ar -> (1, 2)
    QVector<double> epsilon_v;

    // Колебательные энергии : СО2 -> (1:l1, 1:l2, 1:l3)
    QVector<QVector<QVector<double>>> vEnergy_vvv;

    // Приведенные массы
    QVector<QVector<double>> reducedMass_vv;

    // Расчет колебательных энергий
    void computeVEnergy();

    // Расчет приведенных масс
    void computeReducedMass();

public:

    // Конструктор и деструктор по умолчанию
    MixtureData() = default;
    ~MixtureData() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Методы доступа к полям
    const double& dEnergy() const;
    const QVector<double>& mass() const;
    const QVector<double>& sigma() const;
    const QVector<double>& epsilon() const;
    const QVector<QVector<QVector<double>>>& vEnergy() const;
    const QVector<QVector<double>>& reducedMass() const;
};

/// OmegaIntegrals - основной родительский класс, содержит результаты
/// расчета омега-интегралов и их комбинаций, способ расчета определяется в
/// дочерних классах.
class OmegaIntegrals
{
protected:

    // Значения соответствующих омега-интегралов (11, 12; 21, 22)
    QVector<QVector<double>> omega11_vv;
    QVector<QVector<double>> omega12_vv;
    QVector<QVector<double>> omega13_vv;
    QVector<QVector<double>> omega22_vv;

    // Значения комбинаций омега-интегралов (11, 12; 21, 22)
    QVector<QVector<double>> aa_vv;
    QVector<QVector<double>> bb_vv;
    QVector<QVector<double>> cc_vv;

public:

    // Конструктор и деструктор по умолчанию
    OmegaIntegrals() = default;
    virtual ~OmegaIntegrals() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Расчет всех значений омега-интегралов
    virtual void compute(const double& t) = 0;

    // Доступ к соответствующим полям
    const QVector<QVector<double>>& omega11() const;
    const QVector<QVector<double>>& omega12() const;
    const QVector<QVector<double>>& omega13() const;
    const QVector<QVector<double>>& omega22() const;
    const QVector<QVector<double>>& aa() const;
    const QVector<QVector<double>>& bb() const;
    const QVector<QVector<double>>& cc() const;
};

/// OmegaIntegralsDc (Direct Computation) - дочерний класс, производит расчет
/// омега-интегралов и их комбинаций по определению, для этого дополнительно
/// использует класс MixtureData. Тип потенциала взаимодействия определяется
/// в дочерних классах.
class OmegaIntegralsDc : public OmegaIntegrals
{
protected:

    // Указатель на данные смеси
    const MixtureData* mixture_p;

public:

    // Конструктор и деструктор по умолчанию
    OmegaIntegralsDc() = default;
    virtual ~OmegaIntegralsDc() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Привязка к внешним классам
    void link(const MixtureData& mixture);

    // Расчет всех значений омега-интегралов
    virtual void compute(const double& t) = 0;
};

/// OmegaIntegralsDcLjp (Lennard-Jones Potential) - дочерний класс, производит
/// расчет омега-интегралов и их комбинаций в соответствии с потенциалом
/// Леннарда-Джонса.
class OmegaIntegralsDcLjp : public OmegaIntegralsDc
{
protected:

    // Вспомогательные переменные
    double sig_ij;
    double e_ij;
    double s_ij;
    double tx;
    double x11;

public:

    // Конструктор и деструктор по умолчанию
    OmegaIntegralsDcLjp() = default;
    virtual ~OmegaIntegralsDcLjp() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Расчет всех значений омега-интегралов
    virtual void compute(const double& t);
};

/// SpecificHeat - основной родительский класс, содержит результаты
/// расчета удельных теплоемкостей колебательных степеней свободы молекулы
/// CO2, способ расчета определяется в дочерних классах.
class SpecificHeat
{
protected:

    // Статистические суммы для колебательных мод : 12, 3
    double zvT12_s;
    double zvT3_s;

    // Удельные теплоемкости для колебательных мод : 12, 3
    double cvT12_s;
    double cvT3_s;

public:

    // Конструктор и деструктор по умолчанию
    SpecificHeat() = default;
    virtual ~SpecificHeat() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Расчет всех значений
    virtual void compute(const double& t12, const double& t3) = 0;

    // Доступ к соответствующим полям
    const double& zvT12() const;
    const double& zvT3() const;
    const double& cvT12() const;
    const double& cvT3() const;
};

/// SpecificHeatDc (Direct Computation) - дочерний класс, производит расчет
/// удельных теплоемкостей колебательных степеней свободы молекулы CO2 по
/// определению, для этого дополнительно использует класс MixtureData.
class SpecificHeatDc : public SpecificHeat
{
protected:

    // Указатель на данные смеси
    const MixtureData* mixture_p;

    // Расчет статистических сумм и удельных теплоемкостей
    void computeZv(const double& t12, const double& t3);
    void computeCv(const double& t12, const double& t3);

public:

    // Конструктор и деструктор по умолчанию
    SpecificHeatDc() = default;
    virtual ~SpecificHeatDc() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Привязка к внешним классам
    void link(const MixtureData& mixture);

    // Расчет всех значений
    virtual void compute(const double& t12, const double& t3);
};

/// Energy - основной родительский класс, содержит результаты
/// расчета полной и колебательных энергий смеси, способ расчета
/// определяется в дочерних классах.
class Energy
{
protected:

    // Колебательные, вращательная, поступательная и полная энергии
    double vE12_s;
    double vE3_s;
    double rE_s;
    double tE_s;
    double fullE_s;

public:

    // Конструктор и деструктор по умолчанию
    Energy() = default;
    virtual ~Energy() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Расчет всех значений
    virtual void compute(const macroParam& param) = 0;

    // Доступ к соответствующим полям
    const double& vE12() const;
    const double& vE3() const;
    const double& rE() const;
    const double& tE() const;
    const double& fullE() const;
};

/// EnergyDc (Direct Computation) - дочерний класс, производит расчет
/// колебательных, вращательной, поступательной энергий смеси по
/// определению, использует классы MixtureData и SpecificHeat.
class EnergyDc : public Energy
{
protected:

    // Указатель на данные смеси
    const MixtureData* mixture_p;

    // Указатель на инструмент расчета удельных теплопроводностей
    const SpecificHeat* heat_p;

    // Расчет соответствующих энергий
    void computeVe12(const macroParam& param);
    void computeVe3(const macroParam& param);
    void computeRe(const macroParam& param);
    void computeTe(const macroParam& param);
    void computeFullE(const macroParam& param);

public:

    // Конструктор и деструктор по умолчанию
    EnergyDc() = default;
    virtual ~EnergyDc() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Привязка к внешним классам
    void link(const MixtureData& mixture);
    void link(const SpecificHeat& heat);

    // Расчет всех значений
    virtual void compute(const macroParam& param);
};

/// BracketIntegrals - основной родительский класс, содержит результаты
/// расчета интегральных скобок, способ расчета определяется в дочерних
/// классах.
class BracketIntegrals
{
protected:

    // Интегральные скобки (N_SORTS x N_SORTS)
    QVector<QVector<double>> lambda_vv;
    QVector<QVector<double>> lambda00_vv;
    QVector<QVector<double>> lambda01_vv;
    QVector<QVector<double>> lambda11_vv;
    QVector<QVector<double>> eta_vv;
    QVector<QVector<double>> h00_vv;
    QVector<QVector<double>> beta11_vv;

    // Интегральные скобки (N_SORTS x N_POLYATOMIC_SORTS)
    QVector<QVector<double>> beta01_vv;

    // Интегральные скобки (N_POLYATOMIC_SORTS)
    QVector<double> beta0011_v;

    // Внутренние коэффициенты теплопроводности (N_SORTS)
    QVector<double> lambdaInt_v;

public:

    // Конструктор и деструктор по умолчанию
    BracketIntegrals() = default;
    virtual ~BracketIntegrals() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Расчет всех значений
    virtual void compute(const double& t, const QVector<double>& x) = 0;

    // Доступ к соответствующим полям
    const QVector<QVector<double>>& lambda() const;
    const QVector<QVector<double>>& lambda00() const;
    const QVector<QVector<double>>& lambda01() const;
    const QVector<QVector<double>>& lambda11() const;
    const QVector<QVector<double>>& eta() const;
    const QVector<QVector<double>>& h00() const;
    const QVector<QVector<double>>& beta11() const;
    const QVector<QVector<double>>& beta01() const;
    const QVector<double>& beta0011() const;
    const QVector<double>& lambdaInt() const;
};

/// OmegaIntegralsDc (Direct Computation) - дочерний класс, производит расчет
/// интегральных скобок по определению, для этого дополнительно использует
/// классы MixtureData и OmegaIntegrals.
class BracketIntegralsDc : public BracketIntegrals
{
protected:

    // Указатель на данные смеси
    const MixtureData* mixture_p;

    // Указатель на инструмент расчета омега-интегралов
    const OmegaIntegrals* omega_p;

    // Вспомогательная интегральная скобка (N_SORTS)
    QVector<double> phi_v;

    // Расчет phi
    void computePhi(const double& t);

    // Расчет всех интегральных скобок по группам (1, 2, 3)
    void computeBrackets1(const double& t);
    void computeBrackets2(const double& t, const QVector<double>& x);
    void computeBrackets3(const double& t, const QVector<double>& x);

public:

    // Конструктор и деструктор по умолчанию
    BracketIntegralsDc() = default;
    virtual ~BracketIntegralsDc() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Привязка к внешним классам
    void link(const MixtureData& mixture);
    void link(const OmegaIntegrals& omega);

    // Расчет всех значений
    virtual void compute(const double& t, const QVector<double>& x);
};

/// TransportCoefficients - основной родительский класс, содержит результаты
/// расчета коэффициентов переноса, способ расчета определяется в дочерних
/// классах.
class TransportCoefficients
{
protected:

    // Коэффициенты сдвиговой (shear) и объемной (bulk) вязкости
    double sViscosity_s;
    double bViscosity_s;

    // Коэффициенты теплопроводности (tr, rot, tr + rot, vibr)
    double tLambda_s;
    double rLambda_s;
    double cLambda_s;
    double vLambdaT12_s;
    double vLambdaT3_s;

    // Коэффициенты диффузии (11, 12; 21, 22)
    QVector<QVector<double>> diffusion_vv;

    // Коэффициенты термодиффузии (1, 2)
    QVector<double> tDiffusion_v;

public:

    // Конструктор и деструктор по умолчанию
    TransportCoefficients() = default;
    virtual ~TransportCoefficients() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Расчет всех значений
    virtual void compute(const double& t, const QVector<double>& x,
                         const double& p) = 0;

    // Доступ к соответствующим полям
    const double& sViscosity() const;
    const double& bViscosity() const;
    const double& tLambda() const;
    const double& rLambda() const;
    const double& cLambda() const;
    const double& vLambdaT12() const;
    const double& vLambdaT3() const;
    const QVector<QVector<double>>& diffusion() const;
    const QVector<double>& tDiffusion() const;
};

/// TransportCoefficientsDc (Direct Computation) - дочерний класс, производит
/// расчет всех коэффициентов переноса по определению, для этого дополнительно
/// использует классы MixtureData и BracketIntegrals.
class TransportCoefficientsDc : public TransportCoefficients
{
protected:

    // Указатель на данные смеси
    const MixtureData* mixture_p;

    // Указатель на инструмент расчета интегральных скобок
    const BracketIntegrals* bracket_p;

    // Указатель на инструмент расчета удельных теплопроводностей
    const SpecificHeat* heat_p;

    // Матрицы систем
    QVector<QVector<double>> mLambda_vv;
    QVector<QVector<double>> mDiffusion_vv;
    QVector<QVector<double>> mSViscosity_vv;
    QVector<QVector<double>> mBViscosity_vv;

    // Матрицы свободных членов
    QVector<QVector<double>> bLambda_vv;
    QVector<QVector<double>> bDiffusion_vv;
    QVector<QVector<double>> bSViscosity_vv;
    QVector<QVector<double>> bBViscosity_vv;

    // Вспомогательные поля
    QVector<size_t> indexColumn_v;
    QVector<size_t> indexRow_v;
    QVector<size_t> iPiv_v;

    // Заполнение матриц систем и матриц свободных членов
    void fillBMLambda(const double& nTot, const double& rho,
                      const QVector<double>& x);
    void fillBMDiffusion(const double& nTot, const double& rho,
                         const QVector<double>& x);
    void fillBMSViscosity(const double& t, const QVector<double>& x);
    void fillBMBViscosity(const double& nTot, const double& rho,
                          const QVector<double>& x);

    // Решение СЛУ методом Гаусса
    void gauss(QVector<QVector<double>>& m, QVector<QVector<double>>& b);

    // Расчет всех коэффициентов переноса, базируясь на решении СЛУ
    void computeTransport(const double& t, const double& nTot,
                          const QVector<double>& x);

public:

    // Конструктор и деструктор по умолчанию
    TransportCoefficientsDc() = default;
    virtual ~TransportCoefficientsDc() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Привязка к внешним классам
    void link(const MixtureData& mixture);
    void link(const BracketIntegrals& bracket);
    void link(const SpecificHeat& heat);

    // Расчет всех значений
    virtual void compute(const double& t, const QVector<double>& x,
                         const double& p);
};

/// FlowMembers - основной родительский класс, содержит результаты
/// расчета вектора поточных членов, способ расчета
/// определяется в дочерних классах.
class FlowMembers
{
protected:

    // Потоки энергии и компонента тензора напряжений
    double trQ_s;
    double vQ12_s;
    double vQ3_s;
    double diffQ_s;
    double fullQ_s;
    double xxP_s;

    // Энтальпия, скорость диффузии сортов, дифф-ая термодинамическая сила
    QVector<double> h_v;
    QVector<double> diffV_v;
    QVector<double> d_v;

    // Вектор поточных членов
    QVector<double> flow_v;

public:

    // Конструктор и деструктор по умолчанию
    FlowMembers() = default;
    virtual ~FlowMembers() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Расчет всех значений
    virtual void compute(const macroParam& param, const QVector<double>& dx_dx,
                         const double& dp_dx, const double& dT_dx,
                         const double& dT12_dx, const double& dT3_dx,
                         const double& dv_dx) = 0;

    // Доступ к соответствующим полям
    const double& trQ() const;
    const double& vQ12() const;
    const double& vQ3() const;
    const double& diffQ() const;
    const double& fullQ() const;
    const double& xxP() const;
    const QVector<double>& h() const;
    const QVector<double>& diffV() const;
    const QVector<double>& d() const;
    const QVector<double>& flow() const;
};

/// FlowMembersDc (Direct Computation) - дочерний класс, производит расчет
/// вектора поточных членов смеси по определению, использует классы
/// MixtureData, Energy, TransportCoefficients.
class FlowMembersDc : public FlowMembers
{
protected:

    // Указатели на необходимые инструменты
    const MixtureData* mixture_p;
    const Energy* energy_p;
    const TransportCoefficientsDc* transport_p;

    // Расчет соответствующих потоков энергии, энтальпии, компоненты тензора
    void computeD(const macroParam& param, const QVector<double>& dx_dx,
                  const double& dp_dx);
    void computeDiffV(const macroParam& param, const double& dT_dx);
    void computeH(const macroParam& param);
    void computeDiffQ(const macroParam& param);
    void computeFullQ(const double& dT_dx, const double& dT12_dx,
                      const double& dT3_dx);
    void computeXxP(const macroParam& param, const double& dv_dx);

public:

    // Конструктор и деструктор по умолчанию
    FlowMembersDc() = default;
    virtual ~FlowMembersDc() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Привязка к внешним классам
    void link(const MixtureData& mixture);
    void link(const Energy& energy);
    void link(const TransportCoefficientsDc& transport);

    // Расчет всех значений
    virtual void compute(const macroParam& param, const QVector<double>& dx_dx,
                         const double& dp_dx, const double& dT_dx,
                         const double& dT12_dx, const double& dT3_dx,
                         const double& dv_dx);
};

/// Computer - расчитывает коэффициенты переноса и поточные члены,
/// предоставляет к ним доступ, дополнительно использует MixtureData.
class Computer
{
private:

    // Инструменты для расчета вектора потоков
    OmegaIntegralsDcLjp omegaComputer_;
    SpecificHeatDc heatComputer_;
    EnergyDc energyComputer_;
    BracketIntegralsDc bracketComputer_;
    TransportCoefficientsDc transportComputer_;
    FlowMembersDc flowComputer_;

    // Указатели на необходимые инструменты
    const MixtureData* mixture_p;

public:

    // Выделение памяти, обновление значений
    void initialize();

    // Привязка к классу данных
    void link(const MixtureData& mixture);

    // Расчет всех коэффициентов
    void compute(const macroParam& param, const QVector<double>& dx_dx,
                 const double& dp_dx, const double& dT_dx,
                 const double& dT12_dx, const double& dT3_dx,
                 const double& dv_dx);

    // Доступ к соответствующим коэффициентам
    const OmegaIntegrals& omegaComputer() const;
    const SpecificHeat& heatComputer() const;
    const Energy& energyComputer() const;
    const BracketIntegrals& bracketComputer() const;
    const TransportCoefficients& transportComputer() const;
    const FlowMembers& flowComputer() const;
};

} // namespace tc_2

#endif // TRANSPORT_M_2_H
