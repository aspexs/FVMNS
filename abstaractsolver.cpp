#include "abstaractsolver.h"
#include <QMessageBox>
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
        QMessageBox msgBox;
        msgBox.setWindowTitle("Ошибка");
        msgBox.setText("Нет файла с CVIBR");
        msgBox.exec();
        breaksolve = true;
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
        QMessageBox msgBox;
        msgBox.setWindowTitle("Ошибка");
        msgBox.setText("Нет файла с VibrEnergy");
        msgBox.exec();
        breaksolve = true;
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
        QMessageBox msgBox;
        msgBox.setWindowTitle("Ошибка");
        msgBox.setText("Нет файла с Full Energy");
        msgBox.exec();
        breaksolve = true;
    }
    fileEnergy.close();
    fileVibrEnergy.close();
    fileCVibr.close();
}

void AbstaractSolver::prepareSolving()
{
    U1.resize(solParam.NumCell+2);
    U2.resize(solParam.NumCell+2);
    U3.resize(solParam.NumCell+2);
    U4.resize(solParam.NumCell+2);
    R.resize(solParam.NumCell+2);

    leftParam.density = leftParam.pressure /(UniversalGasConstant/molMass * leftParam.temp);
    double Z = additionalSolver.ZCO2Vibr(leftParam.temp);
    double Cv = 5.0/2 * kB/mass+ additionalSolver.CVibr(leftParam.temp,Z);
    solParam.Gamma = (UniversalGasConstant/molMass + Cv)/Cv;
    leftParam.soundSpeed = sqrt(solParam.Gamma *UniversalGasConstant/molMass * leftParam.temp);
    leftParam.velocity = solParam.Ma*leftParam.soundSpeed;
    leftParam.tempIntr = leftParam.temp;
    QFile pythonFile(QDir::currentPath() + "\\Fun.py");
    if( !pythonFile.open(QFile::ReadOnly) )
    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Ошибка");
        msgBox.setText("Нет файла с питоновским расчетом граничных значений");
        msgBox.exec();
        breaksolve = true;
        return;
    }
    rightParam = additionalSolver.bondaryConditionPython(leftParam, solParam);
    double leftEvibr = additionalSolver.vibrEnergy(0,leftParam.temp);
    double leftFullEnergy = 5.0/2*kB*leftParam.temp/mass + leftEvibr;
    double rightEVibr = additionalSolver.vibrEnergy(0,rightParam.temp);
    double rightFullEnergy = 5.0/2*kB*rightParam.temp/mass + rightEVibr;

    for(auto i  = 1; i < solParam.NumCell+1; i++)
    {
        if(i < solParam.NumCell/3 +1)
        {
            U1[i] = leftParam.density;
            U2[i] = leftParam.density*leftParam.velocity;
            U3[i] = leftParam.density*(leftFullEnergy + pow(leftParam.velocity,2)/2);
            U4[i] = leftParam.density*leftEvibr;
        }

        else
        {
            U1[i] = rightParam.density;
            U2[i] = rightParam.density*rightParam.velocity;
            U3[i] = rightParam.density*(rightFullEnergy + pow(rightParam.velocity,2)/2);
            U4[i] = rightParam.density*rightEVibr;
        }
    }
    prepareVectors();
}

void AbstaractSolver::calcRiemanPStar()
{
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calcPStar = [this](int& i)
    {
        macroParam left;
        macroParam right;
        mutex.lock();
        left.density = left_density[i];
        left.velocity = left_velocity[i];
        left.pressure = left_pressure[i];
        right.density = right_density[i];
        right.velocity = right_velocity[i];
        right.pressure = right_pressure[i];
        mutex.unlock();
        rezultAfterPStart[i] = additionalSolver.ExacRiemanSolverCorrect(left,right, solParam.Gamma, solParam.Gamma);
        return;
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calcPStar));
    futureWatcher.waitForFinished();
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
    U5[0]=U5[1];
    U1[solParam.NumCell+1]=U1[solParam.NumCell];
    U2[solParam.NumCell+1]=solParam.typeRightBorder*U2[solParam.NumCell];
    U3[solParam.NumCell+1]=U3[solParam.NumCell];
    U4[solParam.NumCell+1]=U4[solParam.NumCell];
    U5[solParam.NumCell+1]=U5[solParam.NumCell];

    timeSolvind.push_back(0);
    for(int i = 0 ; i<  solParam.NumCell+1; i++)
        vectorForParallelSolving.push_back(i);
    F1.resize(solParam.NumCell+1);
    F2.resize(solParam.NumCell+1);
    F3.resize(solParam.NumCell+1);
    F4.resize(solParam.NumCell+1);
    F5.resize(solParam.NumCell+1);
    P.resize(solParam.NumCell+2);
    Q_v.resize(solParam.NumCell+2);
    Q_v3.resize(solParam.NumCell+2);
    Q_t.resize(solParam.NumCell+2);
    R.resize(solParam.NumCell+1);
    R_1.resize(solParam.NumCell+1);
    R_2.resize(solParam.NumCell+1);
    T.resize(solParam.NumCell +2);
    Tv.resize(solParam.NumCell +2);
    T12.resize(solParam.NumCell +2);
    T3.resize(solParam.NumCell +2);
    Ent.resize(solParam.NumCell +2);
    Ent2.resize(solParam.NumCell +2);
    B_v.resize(solParam.NumCell +2);
    rezultAfterPStart.resize(solParam.NumCell+1);
}
void AbstaractSolver::setTypePlot(int i)
{
    solParam.typePlot = i;
    QVector<double> values, additionalValues;
    //additionalValues.resize(x.size());
    switch (solParam.typePlot)
    {
    case 0: values = pres;break;
    case 1: values = U1; break;
    case 2: values = U2/U1;break;
    case 3: values = T; additionalValues = Tv; break;
    case 4: values = P;break;
    case 5:values = Q_t;break;
    case 6:values = Q_v; additionalValues = Q_v3; break;
    case 7: values= Ent; additionalValues = Ent2; break;
    case 8: values = T12; additionalValues = T3; break;
    case 9: values = B_v; break;
    default: values = U1; break;
    }
    emit updateGraph(x, values, solParam.lambda);
    emit updateAdditionalGraph(x, additionalValues, solParam.lambda);
}
