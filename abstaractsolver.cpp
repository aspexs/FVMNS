#include "abstaractsolver.h"
AbstaractSolver::AbstaractSolver(QObject *parent) : QThread(parent)
{
    QFile fileCVibr(QDir::currentPath() + "\\CVibr.csv");
    QFile fileVibrEnergy(QDir::currentPath() + "\\vibrEnergy.csv");
    QFile fileEnergy(QDir::currentPath() + "\\Energy.csv");
    if( fileCVibr.open(QFile::ReadOnly) )
    {
        QTextStream outCVibr(&fileCVibr);
        QStringList line = outCVibr.readLine().split(";");
        CVibrStartTemp = line[0].toDouble();
        CvibrMass.push_back(line[1].toDouble());
        line = outCVibr.readLine().split(";");
        CVibrStepTemp = line[0].toDouble() - CVibrStartTemp;
        CvibrMass.push_back(line[1].toDouble());
        while (!outCVibr.atEnd())
        {
            QStringList line = outCVibr.readLine().split(";");
            if(line.size() != 2)
                break;
            CvibrMass.push_back(line[1].toDouble());
        }
    }
    else
    {
        qDebug() << "Нет файла с CVIBR";
    }
    if( fileVibrEnergy.open(QFile::ReadOnly) )
    {
        QTextStream outEnergy(&fileVibrEnergy);
        QStringList line = outEnergy.readLine().split(";");
        energyVibrStartTemp = line[0].toDouble();
        EnergyVibr.push_back(line[1].toDouble());
        line = outEnergy.readLine().split(";");
        energyVibrStepTemp = line[0].toDouble() - energyVibrStartTemp;
        EnergyVibr.push_back(line[1].toDouble());

        while (!outEnergy.atEnd())
        {
           QStringList line = outEnergy.readLine().split(";");
           if(line.size() != 2)
               break;

           EnergyVibr.push_back(line[1].toDouble());
        }
    }
    else
    {
        qDebug() << "Нет файла с VibrEnergy";
    }

    if( fileEnergy.open(QFile::ReadOnly) )
    {
        QTextStream outEnergy(&fileEnergy);
        QStringList line = outEnergy.readLine().split(";");
        energyStartTemp = line[0].toDouble();
        Energy.push_back(line[1].toDouble());
        line = outEnergy.readLine().split(";");
        energyStepTemp = line[0].toDouble() - energyStartTemp;
        Energy.push_back(line[1].toDouble());

        while (!outEnergy.atEnd())
        {
           QStringList line = outEnergy.readLine().split(";");
           if(line.size() != 2)
               break;

           EnergyVibr.push_back(line[1].toDouble());
        }
    }
    else
    {
        qDebug() << "Нет файла с Full Energy";
    }

    fileVibrEnergy.close();
    fileCVibr.close();
}

void AbstaractSolver::pause()
{
    pauseSolve = !pauseSolve;
}

void AbstaractSolver::breakSolver()
{
    breaksolve = true;
}
