#ifndef TRANSPORT_M_2_H
#define TRANSPORT_M_2_H

/// Автор : Баталов Семен, 2022
/// Модуль предназаначен для расчета коэффициентов переноса для смеси газов
/// O2-O, а именно: коэффициентов диффузии, термодиффузии, объемной и
/// сдвиговой вязкостей, коэффициентов теплопроводности.

#include "global.h"

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
public:

    // Конструктор и деструктор по умолчанию
    OmegaIntegralsDc() = default;
    virtual ~OmegaIntegralsDc() = default;

    // Выделение памяти, обновление значений
    void initialize();

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
/// O2, способ расчета определяется в дочерних классах.
class SpecificHeat
{
protected:

    // Статистическая сумма для колебательной моды O2
    double zv_s;

    // Удельная колебательная теплоемкость O2
    double cv_s;

public:

    // Конструктор и деструктор по умолчанию
    SpecificHeat() = default;
    virtual ~SpecificHeat() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Расчет всех значений
    virtual void compute(const double& tv) = 0;

    // Доступ к соответствующим полям
    const double& zv() const;
    const double& cv() const;
};

/// SpecificHeatDc (Direct Computation) - дочерний класс, производит расчет
/// удельных теплоемкостей колебательных степеней свободы молекулы O2 по
/// определению, для этого дополнительно использует класс MixtureData.
class SpecificHeatDc : public SpecificHeat
{
protected:

    // Расчет статистических сумм и удельных теплоемкостей
    void computeZv(const double& tv);
    void computeCv(const double& tv);

public:

    // Конструктор и деструктор по умолчанию
    SpecificHeatDc() = default;
    virtual ~SpecificHeatDc() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Расчет всех значений
    virtual void compute(const double& tv);
};

/// Energy - основной родительский класс, содержит результаты
/// расчета полной и колебательных энергий смеси, способ расчета
/// определяется в дочерних классах.
class Energy
{
protected:

    // Колебательные, вращательная, поступательная и полная энергии
    double vE_s;
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
    virtual void compute(const MacroParam& param) = 0;

    // Доступ к соответствующим полям
    const double& vE() const;
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

    // Инструмент расчета удельных теплопроводностей
    SpecificHeatDc heat_;

    // Расчет соответствующих энергий
    void computeVe(const double& tv);
    void computeRe(const double& t);
    void computeTe(const MacroParam& param);
    void computeFullE(const MacroParam& param);

public:

    // Конструктор и деструктор по умолчанию
    EnergyDc() = default;
    virtual ~EnergyDc() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Расчет всех значений
    virtual void compute(const MacroParam& param);
    void compute(const double& tv);
    const SpecificHeat& heat() const;
};

/// Temperature - основной родительский класс, содержит результаты
/// расчета температур смеси по ее энергиям, способ расчета
/// определяется в дочерних классах.
class Temperature
{
protected:

    // Колебательные и поступательная температуры
    double vT_s;
    double tT_s;

public:

    // Конструктор и деструктор по умолчанию
    Temperature() = default;
    virtual ~Temperature() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Расчет всех значений
    virtual void compute(const MacroParam& param, const double& ev,
                         const double& e) = 0;

    // Доступ к соответствующим полям
    const double& Tv() const;
    const double& T() const;
};

/// TemperatureNDc (Not Direct Computation) - дочерний класс, производит
/// расчет колебательных и поступательной температур смеси с помощью линейной
/// интерполяции, использует класс EnergyDc.
class TemperatureNDc : public Temperature
{
protected:

    // Таблицы для интерполяции
    QVector<double> vT_v;
    QVector<double> vE_v;
    EnergyDc energy_;
    double dT;

    // Расчет всех значений температуры
    void calcTv(const double& ev);
    void calcT(const MacroParam& p, const double& e, const double& ev);

public:

    // Конструктор и деструктор по умолчанию
    TemperatureNDc() = default;
    virtual ~TemperatureNDc() = default;

    // Выделение памяти, обновление значений
    void initialize(const double& t0, const double& t1, const int& n);

    // Расчет всех значений
    virtual void compute(const MacroParam& param, const double& ev,
                         const double& e);
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
    virtual void compute(const MacroParam& param) = 0;

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

    // Инструмент расчета омега-интегралов
    OmegaIntegralsDcLjp omega_;

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

    // Расчет всех значений
    virtual void compute(const MacroParam& param);
    const OmegaIntegrals& omega() const;
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
    double vLambda_s;

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
    virtual void compute(const MacroParam& param) = 0;

    // Доступ к соответствующим полям
    const double& sViscosity() const;
    const double& bViscosity() const;
    const double& tLambda() const;
    const double& rLambda() const;
    const double& cLambda() const;
    const double& vLambda() const;
    const QVector<QVector<double>>& diffusion() const;
    const QVector<double>& tDiffusion() const;
};

/// TransportCoefficientsDc (Direct Computation) - дочерний класс, производит
/// расчет всех коэффициентов переноса по определению, для этого дополнительно
/// использует классы MixtureData и BracketIntegrals.
class TransportCoefficientsDc : public TransportCoefficients
{
protected:

    // Инструмент расчета интегральных скобок
    BracketIntegralsDc bracket_;

    // Инструмент расчета удельных теплопроводностей
    SpecificHeatDc heat_;

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
    QVector<int> indexColumn_v;
    QVector<int> indexRow_v;
    QVector<int> iPiv_v;

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

    // Расчет всех значений
    virtual void compute(const MacroParam& param);
    const SpecificHeat& heat() const;
    const BracketIntegrals& bracket() const;
};

/// FlowMembers - основной родительский класс, содержит результаты
/// расчета вектора поточных членов, способ расчета
/// определяется в дочерних классах.
class FlowMembers
{
protected:

    // Потоки энергии и компонента тензора напряжений
    double trQ_s;
    double vQ_s;
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
    virtual void compute(const MacroParam& param, const QVector<double>& dx_dx,
                         const double& dlnp_dx, const double& dT_dx,
                         const double& dlnT_dx, const double& dTv_dx,
                         const double& dv_dx) = 0;

    // Доступ к соответствующим полям
    const double& trQ() const;
    const double& vQ() const;
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

    // Необходимые инструменты
    EnergyDc energy_;
    TransportCoefficientsDc transport_;

    // Расчет соответствующих потоков энергии, энтальпии, компоненты тензора
    void computeD(const MacroParam& param, const QVector<double>& dx_dx,
                  const double& dlnp_dx);
    void computeDiffV(const double& dlnT_dx);
    void computeH(const MacroParam& param);
    void computeDiffQ(const MacroParam& param);
    void computeFullQ(const double& dT_dx, const double& dTv_dx);
    void computeXxP(const MacroParam& param, const double& dv_dx);

public:

    // Конструктор и деструктор по умолчанию
    FlowMembersDc() = default;
    virtual ~FlowMembersDc() = default;

    // Выделение памяти, обновление значений
    void initialize();

    // Расчет всех значений
    virtual void compute(const MacroParam& param, const QVector<double>& dx_dx,
                         const double& dlnp_dx, const double& dT_dx,
                         const double& dlnT_dx, const double& dTv_dx,
                         const double& dv_dx);
    const Energy& energy() const;
    const TransportCoefficients& transport() const;
};

#endif // TRANSPORT_M_2_H
