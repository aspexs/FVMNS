#include <QCoreApplication>
#include <QFile>
#include <QTextStream>
#include <QString>
#include <QRandomGenerator>
#include "mixture.h"

// Запись данных в файл
void writeToFile(const QString& name, const QVector<QVector<double>> table)
{
    // Открываем файл для записи
    QFile file1(name);
    file1.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out1(&file1);

    // Настройка параметров вывода
    out1.setFieldWidth(16);
    out1.setRealNumberPrecision(8);
    out1.setRealNumberNotation(QTextStream::ScientificNotation);

    // Запись таблицы в файл
    for (int i = 0; i < table[0].size(); ++i)
    {
        for (int j = 0; j < table.size(); ++j)
        {
            out1 << table[j][i];
        }
        out1 << '\n';
    }

    // Закрытие файла
    file1.close();
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    MixtureO2O solver;

    // Прямая УВ
    solver.initialize(MacroParam(66.6, 1370, 300, 0.30));
    solver.solve();
    writeToFile("data_O2O.txt", solver.getMacroParams());

    std::cout << "\n > Time        [ms]   : " << 1000 * solver.time;
    std::cout << "\n > dt          [ms]   : " << 1000 * solver.dt;
    std::cout << "\n > Iterations  [-]   : " << solver.currIter << " / "
              << MAX_N_ITER << "\n\n";

//    MacroParam point(66.6, 1370, 300, 0.30);
//    double n0 = point.rho[0] / Mixture::mass(0);
//    double n1 = point.rho[1] / Mixture::mass(1);
//    double x0 = n0 / (n0 + n1);
//    double d = x0 * Mixture::sigma(0) + (1.0 - x0) * Mixture::sigma(1);
//    double L = K_BOLTZMANN * point.t / (qSqrt(2) * M_PI * qPow(d, 2.0) * point.p);
//    qDebug() << L;

    return 0;
}
