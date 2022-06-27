#include <QCoreApplication>
#include <thread>
#include <chrono>
#include "transport_m_2.h"

struct Data
{
public:

    // Температура, давление, концентрация
    double t_, t12_, t3_, p_;
    QVector<double> x_;

    // Инициализация
    void initialize(double t, double t12, double t3, double p,
                    QVector<double> x)
    {
        t_ = t;
        t12_ = t12;
        t3_ = t3;
        p_ = p;
        x_ = x;
    }
};

// Расчет всех коэффициентов переноса в диапазоне значений
void testCompute(size_t beg, size_t end, const QVector<Data>& data,
                 const tc_2::DataDc& dataDc)
{
    // Создаем объект - вычислитель
    tc_2::ComputerDc computer;

    // Инициализируем и привязываем к data
    computer.initialize();
    computer.link(dataDc);

    // Расчет значений в определенном диапазоне точек
    for (size_t i = beg; i < end; ++i)
    {
        computer.compute(data[i].t_, data[i].t12_, data[i].t3_, data[i].x_,
                         data[i].p_);
    }
}

// Пример работы класса вычислителя
void testComputeSeq()
{
    // Вспомогательные переменные
    double t = 2000;
    double t12 = 1000;
    double t3 = 2000;
    double p = 101300;
    QVector<double> x = {0.5, 0.5};

    // Данные смеси
    tc_2::DataDc dataDc;
    dataDc.initialize();

    // Создаем объект - вычислитель
    tc_2::ComputerDc computer;
    computer.initialize();
    computer.link(dataDc);

    // Расчет значений в определенной точке
    computer.compute(t, t12, t3, x, p);

    // Вывод данных на экран
    qDebug() << computer.transport().tLambda();
    qDebug() << computer.transport().rLambda();
    qDebug() << computer.transport().cLambda();
    qDebug() << computer.transport().vLambdaT12();
    qDebug() << computer.transport().vLambdaT3();
    qDebug() << computer.transport().sViscosity();
    qDebug() << computer.transport().bViscosity();
    qDebug() << computer.transport().diffusion();
    qDebug() << computer.transport().tDiffusion();
}

// Пример работы класса вычислителя (распараллеливание)
void testComputeParallel(size_t subSize, size_t num)
{
    // Вспомогательные переменные
    double stepT            = 10;
    double begT             = 200;
    double begT12           = 2000;
    double begT3            = 2000;
    double begP             = 101300;
    QVector<double> begX    = {0.5, 0.5};
    size_t mainSize         = subSize * num;

    // Данные смеси
    tc_2::DataDc dataDc;
    dataDc.initialize();

    // Инициализация тестового массива точек
    QVector<Data> data(mainSize);
    for (size_t i = 0; i < mainSize; ++i)
    {
        data[i].initialize(begT + i * stepT, begT12, begT3, begP, begX);
    }

    // Расчет коэффициентов во всех точках
    std::thread* threads = new std::thread[num];
    for (size_t i = 0; i < num; ++i)
    {
        threads[i] = std::thread(testCompute, i * subSize, (i + 1) * subSize,
                                 std::cref(data), std::cref(dataDc));
    }
    for (size_t i = 0; i < num; ++i)
    {
        threads[i].join();
    }
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    // Начало отсчета
    auto begin = std::chrono::steady_clock::now();

    // Пример параллельного расчета (5 потоков, по 200 точек на поток,
    // всего 1000 точек)
    testComputeParallel(200, 5);

    // Пример расчета коэффициентов переноса
    testComputeSeq();

    // Конец отсчета
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::
            milliseconds>(end - begin);

    // Вывод продолжительности работы
    qDebug() << "The time: " << elapsed_ms.count() << " ms\n";

    return 0;
}
