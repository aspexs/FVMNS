#ifndef MIXTURECO2AR_H
#define MIXTURECO2AR_H

#include "transport_m_2.h"

/// BorderCondition - инструмент расчета макропараметров за ударной волной.
/// Требует информации о характеристиках потока перед УВ, использует для расчета
/// законы сохранения.
class BorderCondition
{
public:

    // Находит равновесные макропараметры за УВ
    void compute(const MacroParam& lP);

    // Возвращает равновесные макропараметры за УВ
    const MacroParam& rP() const;

private:

    // Макропараметры перед УВ и за УВ соответственно
    MacroParam lP_, rP_;

    // Иструмент расчета энергий
    EnergyDc energy_;

    // Общая энергия, массовые доли CO2 и Ar перед УВ соответственно
    double lE_, y0_, y1_, lRho_;

    // Вспомогательные величины
    double alpha_;

private:

    // Расчет вспомогательных величин
    void initialize(const MacroParam& lP);

    // Возвращает дискриминант
    double computeD(const double& t);

    // Возвращает плотность смеси по температуре
    double computeRho(const double& t);

    // Функция для нахождения равновесной температуры за УВ: F(T_eq) = 0
    double computeF(const double& t);

    // Находит равновесную температуру за УВ методом бисекции
    void computeT();

    // Находит равновесное давление за УВ
    void computeP();

    // Находит равновесные значения плотностей сортов за УВ
    void computeRho();

    // Находит равновесную скорость потока за УВ
    void computeV();
};

/// MixtureCo2Ar - предоставляет инструменты для решения задачи о
/// релаксационных процессах за ударной волной в смеси CO2-Ar
class MixtureCo2Ar
{
public:

    // Шаг по времени, суммарное абсолютное изменение вектора конс. переменных
    double dt, error;

    // Хранит номер последней итерации
    int currIter;

public:

    MixtureCo2Ar();

    void initialize(const MacroParam& lP);
    void solve();

    QVector<QVector<double>> saveMacroParams();
    QVector<QVector<double>> saveU();
    QVector<QVector<double>> saveF();
    QVector<QVector<double>> saveHlleF();
    QVector<QVector<double>> saveR();

private:

    // Вспомогательные структуры
    QMutex mutex;
    QVector<int> parAll_v;
    QVector<int> parIn_v;
    TemperatureNDc computeT;

    // Флаг выхода из итерационного процесса
    bool notFinished;

    // Все макропараметры течения во всех точках
    QVector<MacroParam> points;

    // Скорость звука и показатель адиабаты
    QVector<double> a, k;

    // Консервативные переменные, поточные члены, релаксационные члены
    QVector<QVector<double>> U, F, R;

    // Значения потока на границах ячеек по методу HLLE
    QVector<QVector<double>> hlleF;

private:

    // Расчет вектора потоков во всех ячейках
    void computeF();

    // Расчет релаксационных членов
    void computeR();

    // Расчет потоков на стыках ячеек методом HLLE
    void computeHlleF();

    // Производит один шаг итерационного процесса
    void step();

    // Возврат к макропараметрам Ui -> points
    void updateMacroParam();

    // Расчет показателей адиабаты 'k' и скоростей звука 'a'
    void updateAK();

    // Рассчитывает временной шаг по критерию Куранта-Фридрихса-Леви
    void updateDt();
};

#endif // MIXTURECO2AR_H
