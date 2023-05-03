#include <QCoreApplication>
#include <QFile>
#include <QTextStream>
#include <QString>
#include <QRandomGenerator>
#include "mixture.h"

//const double P0 = 0.0835 * 101325;
//const double T0 = 297; //293;
//const double P1 = 0.415 * 101325;
//const double P2 = 1.195 * 101325; //1.195 * 101325;

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
//    BorderCondition bc0, bc1;
//    MacroParam data0 = optimizer(600, 700, P0, T0, P1);
//    bc0.compute(data0);
//    MacroParam data1 = optimizer(700, 800, P1, bc0.rP().t, P2);
//    bc1.compute(data1);

    // Прямая УВ
    solver.initialize(MacroParam(66.6, 1370, 300, 0.3), "energy_data.dat", 3);
    solver.solve();
    writeToFile("data_CO2AR_3T.txt", solver.saveMacroParams());

    std::cout << "\n > Time        [s]   : " << solver.time;
    std::cout << "\n > dt          [s]   : " << solver.dt;
    std::cout << "\n > Iterations  [-]   : " << solver.currIter << " / "
              << MAX_N_ITER << "\n\n";

//    TemperatureNDc temp;
//    MacroParam p(66.6, 1370, 4000, 0.0001);
//    EnergyDc energy;
//    p.t = 4000;
//    p.t12 = 4000;
//    energy.initialize();
//    energy.compute(p);
//    // temp.initialize(T_MIN, T_MAX, T_NUM, Y_NUM);
//    temp.readFromFile("energy_data.dat", T_NUM, Y_NUM);
//    // temp.writeToFile("energy_data.dat", T_NUM, Y_NUM);

//    temp.compute(p, energy.vE12(), energy.vE3(), energy.fullE());
//    qDebug() << temp.T() << temp.T12() << temp.T3();
//    temp.compute(p, energy.vE12() + energy.vE3(), energy.fullE());
//    qDebug() << temp.T() << temp.T12() << temp.T3();
//    temp.compute(p, energy.fullE());
//    qDebug() << temp.T() << temp.T12() << temp.T3();

    return 0;
}
