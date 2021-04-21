#include "abstaractsolver.h"
AbstaractSolver::AbstaractSolver(QObject *parent) : QThread(parent)
{
    QFile fileCVibr(QDir::currentPath() + "\\CVibr.csv");
    QFile fileVibrEnergy(QDir::currentPath() + "\\vibrEnergy.csv");
    QFile fileEnergy(QDir::currentPath() + "\\allEnergy.csv");
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
           Energy.push_back(line[1].toDouble());
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
void AbstaractSolver::prepareVectors()
{
    double x_right =solParam.lambda*solParam.lambdaSol; //% правая граница
    delta_h = (x_right) / solParam.NumCell;
    x.clear();
    x.push_back(0+0.5*delta_h);
    for(auto i = 1; i < solParam.NumCell; i++)
        x.push_back(x[i-1] + delta_h);
    x.push_back(x_right);
    x.push_front(0);
    U1[0]=U1[1];
    U2[0]=solParam.typeLeftBorder*U2[1];
    U3[0]=U3[1];
    U4[0]=U4[1];
    U1[solParam.NumCell+1]=U1[solParam.NumCell];
    U2[solParam.NumCell+1]=solParam.typeRightBorder*U2[solParam.NumCell];
    U3[solParam.NumCell+1]=U3[solParam.NumCell];
    U4[solParam.NumCell+1]=U4[solParam.NumCell];

    timeSolvind.push_back(0);
    for(int i = 0 ; i<  solParam.NumCell+1; i++)
        vectorForParallelSolving.push_back(i);
    F1.resize(solParam.NumCell+1);
    F2.resize(solParam.NumCell+1);
    F3.resize(solParam.NumCell+1);
    F4.resize(solParam.NumCell+1);
    P.resize(solParam.NumCell+2);
    Q_v.resize(solParam.NumCell+2);
    Q_t.resize(solParam.NumCell+2);
    R.resize(solParam.NumCell+1);
    rezultAfterPStart.resize(solParam.NumCell+1);
}
