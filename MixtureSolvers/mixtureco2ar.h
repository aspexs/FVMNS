#ifndef MIXTURECO2AR_H
#define MIXTURECO2AR_H

#include "global.h"
#include "transport_m_2.h"

/// MixtureCo2Ar - предоставляет инструменты для решения задачи о
/// релаксационных процессах за ударной волной в смеси CO2-Ar
class MixtureCo2Ar
{
public:

    void initialize(const SolverParams& init);
    void solve();

    MixtureCo2Ar();

private:

    // Вспомогательные структуры
    QMutex mutex;
    QVector<int> parAll_v;
    QVector<int> parIn_v;
    TemperatureNDc computeT;

    // Настройки решателя
    SolverParams solParam;

    // Шаг по времени, абсолютное суммарное изменение потоков
    double dt, error;

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

    // Обновляет условие на границе за ударной волной
    void updateBoundCond();
};

#endif // MIXTURECO2AR_H
