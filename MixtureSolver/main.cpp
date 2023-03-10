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

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    MixtureCo2Ar solver;

    solver.initialize(MacroParam(66.6, 1370, 300, 0.999));
    solver.solve();
    writeToFile("macro_params.txt", solver.saveMacroParams());

    std::cout << "\n > errMax : " << solver.errMax << '\n';

//    writeToFile("u.txt", solver.saveU());
//    writeToFile("f.txt", solver.saveF());
//    writeToFile("hlle_f.txt", solver.saveHlleF());
//    writeToFile("r.txt", solver.saveR());

    return 0;
}
