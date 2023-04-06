#include <QCoreApplication>
#include <QFile>
#include <QTextStream>
#include <QString>
#include <QRandomGenerator>
#include "mixtureco2ar.h"

//const double P2 = 1.195 * 101325;
//const double T2 = 955;
//const double P1 = 0.417 * 101325; //0.0851 * 101325;
const double P0 = 0.0851 * 101325;
const double T0 = 293;
const double P1 = 0.417 * 101325;
const double P2 = 1.195 * 101325;

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

// Функция ошибки (отклонения давления перед отраженной УВ)
double error(const double& v, const double& p0, const double& t0,
             const double& p1)
{
    BorderCondition bc;
    bc.compute(MacroParam(p0, v, t0, 0.02));
    return qAbs(bc.rP().p - p1);
}

// Находит минимум функции ошибки
MacroParam optimizer(const double& down, const double& up, const double& p0,
                     const double& t0, const double& p1)
{
    double v_new = down;
    double v_old = down;

    double newErr = error(v_new, p0, t0, p1);
    double oldErr = newErr;
    double dv = (up - down) / 500;

    for (int i = 0; i < 500; ++i)
    {
        v_new = down + dv * i;
        newErr = error(v_new, p0, t0, p1);

        if (newErr < oldErr)
        {
            v_old = v_new;
            oldErr = newErr;
        }
    }
    return MacroParam(p0, v_old, t0, 0.02);
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    MixtureCo2Ar solver;
    BorderCondition bc0, bc1;
    MacroParam data0 = optimizer(600, 700, P0, T0, P1);
    bc0.compute(data0);
    MacroParam data1 = optimizer(690, 790, P1, bc0.rP().t, P2);
    bc1.compute(data1);

    double v = data1.v + bc0.rP().v - data0.v;
    qDebug() << v << bc1.rP().v << bc1.rP().t << bc1.rP().p;

    // solver.initialize(MacroParam(66.6, 1370, 300, 0.999));
    solver.initialize(data1);
    solver.solve();
    writeToFile("macro_params.txt", solver.saveMacroParams());
    writeToFile("F.txt", solver.saveF());
    writeToFile("hlleF.txt", solver.saveHlleF());
    writeToFile("U.txt", solver.saveU());
    writeToFile("R.txt", solver.saveR());

    std::cout << "\n > Time           : " << solver.time;
    std::cout << "\n > dt             : " << solver.dt;
    std::cout << "\n > Shock position : " << solver.curShockPos;
    std::cout << "\n > Iterations     : " << solver.currIter << " / "
              << MAX_N_ITER << '\n';

    return 0;
}
