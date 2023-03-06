#include <QCoreApplication>
#include <QFile>
#include <QTextStream>
#include "transport_m_2.h"

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
    for (int i = 0; i < table.size(); ++i)
    {
        for (int j = 0; j < table[i].size(); ++j)
        {
            out1 << table[i][j];
        }
        out1 << '\n';
    }

    // Закрытие файла
    file1.close();
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    TemperatureNDc tComputer;
    tComputer.initialize(200, 1000, 800);
    QVector<QVector<double>> table;

    table.push_back(tComputer.vT12_v);
    table.push_back(tComputer.vE12_v);
    table.push_back(tComputer.vT3_v);
    table.push_back(tComputer.vE3_v);

    writeToFile("test_table.txt", table);

    return a.exec();
}
