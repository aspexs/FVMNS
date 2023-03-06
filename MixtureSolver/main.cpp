#include <QCoreApplication>
#include <QFile>
#include <QTextStream>
#include <QString>
#include "mixtureco2ar.h"

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

SolverParams createParam(double p, double v, double t, double x_CO2,
                         int nIter, double k)
{
    SolverParams param;

    // Left point
    param.lPoint.v   = v;
    param.lPoint.t   = t;
    param.lPoint.t12 = t;
    param.lPoint.t3  = t;
    param.lPoint.p   = p;
    param.lPoint.computeRho(x_CO2);

    // Right point
    param.rPoint.t   = t;
    param.rPoint.t12 = t;
    param.rPoint.t3  = t;
    param.rPoint.p   = 0.5 * (param.lPoint.rho[0] + param.lPoint.rho[1]) *
            qPow(param.lPoint.v, 2.0) + param.lPoint.p;
    param.rPoint.computeRho(x_CO2);
    param.rPoint.v   = (param.lPoint.rho[0] + param.lPoint.rho[1]) *
            param.lPoint.v / (param.rPoint.rho[0] + param.rPoint.rho[1]);

    // Grid settings
    param.nIter = nIter;
    param.k = k;
    return param;
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    MixtureCo2Ar solver;

    solver.initialize(createParam(1e3, 600, 300, 0.5, 1000, 2e-2));
    solver.solve();
    writeToFile("macro_params.txt", solver.saveMacroParams());
    writeToFile("u.txt", solver.saveU());
    writeToFile("f.txt", solver.saveF());
    writeToFile("hlle_f.txt", solver.saveHlleF());
    writeToFile("r.txt", solver.saveR());

    return 0;
}
