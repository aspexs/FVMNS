#include "co23tsolver.h"
#include <QMessageBox>
Co23TSolver::Co23TSolver(QObject *parent) : AbstaractSolver(parent)
{
    QFile fileVibrEnergy12(QDir::currentPath() + "\\vibrEnergy12.csv");
    QFile fileVibrEnergy3(QDir::currentPath() + "\\vibrEnergy3.csv");
    if( fileVibrEnergy12.open(QFile::ReadOnly) )
    {
        QTextStream outEnergy(&fileVibrEnergy12);
        QStringList line = outEnergy.readLine().split(";");
        energyVibrStartTemp12 = line[0].toDouble();
        EnergyVibr12.push_back(line[1].toDouble());
        line = outEnergy.readLine().split(";");
        energyVibrStepTemp12 = line[0].toDouble() - energyVibrStartTemp12;
        EnergyVibr12.push_back(line[1].toDouble());

        while (!outEnergy.atEnd())
        {
           QStringList line = outEnergy.readLine().split(";");
           if(line.size() != 2)
               break;

           EnergyVibr12.push_back(line[1].toDouble());
        }
    }
    else
    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Ошибка");
        msgBox.setText("Нет файла с VibrEnergy12");
        msgBox.exec();
        breaksolve = true;
    }

    if( fileVibrEnergy3.open(QFile::ReadOnly) )
    {
        QTextStream outEnergy(&fileVibrEnergy3);
        QStringList line = outEnergy.readLine().split(";");
        energyVibrStartTemp3 = line[0].toDouble();
        EnergyVibr3.push_back(line[1].toDouble());
        line = outEnergy.readLine().split(";");
        energyVibrStepTemp3 = line[0].toDouble() - energyVibrStartTemp3;
        EnergyVibr3.push_back(line[1].toDouble());

        while (!outEnergy.atEnd())
        {
           QStringList line = outEnergy.readLine().split(";");
           if(line.size() != 2)
               break;

           EnergyVibr3.push_back(line[1].toDouble());
        }
    }
    else
    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Ошибка");
        msgBox.setText("Нет файла с VibrEnergy3");
        msgBox.exec();
        breaksolve = true;
    }
    fileVibrEnergy12.close();
    fileVibrEnergy3.close();
}

void Co23TSolver::prepareSolving()
{
    U1.resize(solParam.NumCell+2);
    U2.resize(solParam.NumCell+2);
    U3.resize(solParam.NumCell+2);
    U4.resize(solParam.NumCell+2);
    U5.resize(solParam.NumCell+2);
    R_1.resize(solParam.NumCell+2);
    R_2.resize(solParam.NumCell+2);

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
    double leftEvibr12 = additionalSolver.EVibr12(0,leftParam.temp);
    double leftEvibr3 = additionalSolver.EVibr3(0,leftParam.temp);
    double leftFullEnergy = 5.0/2*kB*leftParam.temp/mass + leftEvibr12 + leftEvibr3;
    double rightEVibr12 = additionalSolver.EVibr12(0,rightParam.temp);
    double rightEVibr3 = additionalSolver.EVibr3(0,rightParam.temp);
    double rightFullEnergy = 5.0/2*kB*rightParam.temp/mass + rightEVibr12 + rightEVibr3;

    for(auto i  = 1; i < solParam.NumCell+1; i++)
    {
        if(i < solParam.NumCell/3 +1)
        {
            U1[i] = leftParam.density;
            U2[i] = leftParam.density*leftParam.velocity;
            U3[i] = leftParam.density*(leftFullEnergy + pow(leftParam.velocity,2)/2);
            U4[i] = leftParam.density*leftEvibr12;
            U5[i] = leftParam.density*leftEvibr3;
        }

        else
        {
            U1[i] = rightParam.density;
            U2[i] = rightParam.density*rightParam.velocity;
            U3[i] = rightParam.density*(rightFullEnergy + pow(rightParam.velocity,2)/2);
            U4[i] = rightParam.density*rightEVibr12;
            U5[i] = rightParam.density*rightEVibr3;
        }
    }
    prepareVectors();
}

void Co23TSolver::solveFlux()
{
    calcRiemanPStar();
    calcFliux();
    calcR(U1,U2,U3,U4,U5);
}

void Co23TSolver::calcFliux()
{
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calcFlux = [this](int& i)
    {
        mutex.lock();
        auto point = rezultAfterPStart[i];
        double tempL = Tl[i];
        double tempR = Tr[i];
        auto du_dx = (right_velocity[i] - left_velocity[i])/delta_h;
        auto tempLT12 = T12L[i];
        auto tempRT12 = T12R[i];
        auto tempLT3 = T3L[i];
        auto tempRT3 = T3R[i];
        double Tx =  point.pressure/(point.density*UniversalGasConstant/molMass);
        double T12 = tempLT12;
        double T3 = tempLT3;
        mutex.unlock();
        double energyVibr12 =  additionalSolver.EVibr12(0,T12);
        double energyVibr3 =  additionalSolver.EVibr3(0,T3);

        double etta = additionalSolver.shareViscosityOmega(0,Tx);
        double zetta =additionalSolver.bulcViscosityOnlyTRRot(0,Tx);

        double P = (4.0/3*etta + zetta)*du_dx;
        double dt_dx = (tempR - tempL)/delta_h;
        double dt12_dx = (tempRT12 - tempLT12)/delta_h;
        double dt3_dx = (tempRT3 - tempLT3)/delta_h;
        double qVibr12 = -additionalSolver.Lambda12(0,T12)* dt12_dx;
        double qVibr3 = -additionalSolver.Lambda3(0,T3)* dt3_dx;
        double qTr = -additionalSolver.lambdaTr_Rot(Tx)*dt_dx;
        double Etr_rot = 5.0/2*kB*Tx/mass;
        double entalpi = Etr_rot + energyVibr12 + energyVibr3 + point.pressure/point.density + pow(point.velocity,2)/2;

        F1[i] = point.density * point.velocity;
        F2[i] = point.density * point.velocity*point.velocity + point.pressure - P;
        F3[i] = point.density * point.velocity*entalpi - P*point.velocity + qVibr12 + qVibr3 + qTr;
        F4[i] = point.density * point.velocity*energyVibr12 + qVibr12;
        F5[i] = point.density * point.velocity*energyVibr3 + qVibr3;
        Q_v[i] = qVibr12;
        Q_t[i] = qTr;
        Q_v3[i] = qVibr3;
        this->P[i] = P;
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calcFlux));
    futureWatcher.waitForFinished();
}

void Co23TSolver::solve()
{
     prepareSolving();
     for(auto i  = 0; i < solParam.MaxIter; i++)
     {
         if(breaksolve)
             break;
         while(pauseSolve)
         {
             if(breaksolve)
                 break;
         }
         Matrix velosity = U2/U1;

         auto EVibr12 = U4/U1;
         auto EVibr3 = U5/U1;
         auto energyFull = U3/U1 - Matrix::POW(velosity,2)/2;
         auto E_tr_rot = energyFull - EVibr12 - EVibr3;
         T.clear();
         T12.clear();
         T3.clear();
         Tv.clear();
         pres.clear();
         Matrix gammas;
         double Cv = 5.0/2 * kB/mass;
         auto gamma = (UniversalGasConstant/molMass + Cv)/Cv;
         for (auto energy: EVibr12)
             T12.push_back(getEnergyVibr12Temp(energy));
         for (auto energy: EVibr3)
             T3.push_back(getEnergyVibr3Temp(energy));
         for (auto energy: E_tr_rot)
         {
             T.push_back(energy*2*mass/(5*kB));
             gammas.push_back(gamma);
         }
         Tv = (T12 + T3)/2;

         pres = U1*T*UniversalGasConstant/molMass;

         dt = additionalSolver.getTimeStepFull(velosity, U1, pres, delta_h,solParam,gammas);
         timeSolvind.push_back(dt+timeSolvind.last());

         auto U1L = U1; U1L.removeLast();
         auto U2L = U2; U2L.removeLast();
         auto U3L = U3; U3L.removeLast();
         auto pressureL = pres;
         auto T12L = T12;
         auto T3L = T3;
         auto Templ = T;
         Templ.removeLast();
         pressureL.removeLast();
         left_pressure = pressureL;
         T12L.removeLast();
         T3L.removeLast();
         left_density   = U1L;
         left_velocity  = U2L / left_density;
         this->T12L = T12L;
         this->T3L = T3L;
         this->Tl = Templ;

         auto U1R = U1; U1R.removeFirst();
         auto U2R = U2; U2R.removeFirst();
         auto U3R = U3; U3R.removeFirst();
         auto pressureR = pres;
         auto T12R = T12;
         auto T3R = T3;
         auto TempR = T;
         TempR.removeFirst();
         pressureR.removeFirst();
         T12R.removeFirst();
         T3R.removeFirst();
         right_pressure = pressureR;
         right_density  = U1R;
         right_velocity = U2R / right_density;
         this->T12R = T12R;
         this->T3R = T3R;
         this->Tr = TempR;

         solveFlux();

         auto res = additionalSolver.SEEFOForCO23T(F1, F2, F3, F4, F5, U1, U2, U3, U4, U5, dt, delta_h, R_1, R_2);
         auto copyU1 = U1;
         auto copyU2 = U2;
         auto copyU3 = U3;
         auto copyU4 = U4;
         auto copyU5 = U5;
         U1 = res[0];
         U2 = res[1];
         U3 = res[2];
         U4 = res[3];
         U5 = res[4];
         error = res[5].first();

         U1[0]=U1[1];                         U1[U1.size() - 1] = U1[U1.size() - 2];
         U2[0]=solParam.typeLeftBorder*U2[1]; U2[U2.size()-1]= solParam.typeRightBorder*U2[U2.size()-2];
         U4[0]=U4[1];                         U4[U4.size() - 1] = U4[U4.size() - 2];
         U5[0]=U5[1];                         U5[U5.size() - 1] = U5[U5.size() - 2];
         U3[0]=U3[1];
         auto T = rightParam.pressure* molMass/(UniversalGasConstant*U1.last());
         auto T12 = getEnergyVibr12Temp(U4.last()/U1.last());
         auto T3 = getEnergyVibr3Temp(U5.last()/U1.last());
         double Evibr12 = additionalSolver.EVibr12(0,T12);
         double Evibr3 = additionalSolver.EVibr3(0,T3);
         double rightFullEnergy = 5.0/2*kB*T/mass + Evibr12 + Evibr3;
         U3[U3.size() - 1] = U1.last()*(rightFullEnergy + pow(U2.last()/U1.last(),2)/2);

         Ent = U3/U1 + pres/U1;

         if(i % solParam.PlotIter == 0)
         {
             setTypePlot(solParam.typePlot);
             emit updateTime(timeSolvind.last(), error);

             QThread::msleep(200);
         }
         if(timeSolvind.last() > solParam.t_fin)
         {
             setTypePlot(solParam.typePlot);
             emit updateTime(timeSolvind.last(), error);
             QThread::msleep(200);
             break;
         }
     }
}

double Co23TSolver::getEnergyVibr12Temp(double E)
{
    for(auto i = 0 ; i < EnergyVibr12.size(); i++)
        if( E < EnergyVibr12[i])
            return (i-1)  * energyVibrStepTemp12 + energyVibrStartTemp12;
    return (EnergyVibr12.size()-1) * energyVibrStepTemp12 + energyVibrStartTemp12;
}

double Co23TSolver::getEnergyVibr3Temp(double E)
{
    for(auto i = 0 ; i < EnergyVibr3.size(); i++)
        if( E < EnergyVibr3[i])
            return (i-1)  * energyVibrStepTemp3 + energyVibrStartTemp3;
    return (EnergyVibr3.size()-1) * energyVibrStepTemp3 + energyVibrStartTemp3;
}

void Co23TSolver::calcR(const Matrix &U1, const Matrix &U2, const Matrix &U3, const Matrix &U4, const Matrix &U5)
{
    Matrix T122, T3, T,pressure;
    auto tempU1 = U1;
    auto tempU2 = U2;
    auto velosity = tempU2/U1;
    auto tempU3 = U3;
    auto tempU4 = U4;
    auto tempU5 = U5;
    auto EVibr12 = tempU4/U1;
    auto EVibr3 = tempU5/U1;
    auto energyFull = tempU3/U1 - Matrix::POW(velosity,2)/2;
    auto E_tr_rot = energyFull - EVibr12 - EVibr3;
    for (auto energy: EVibr12)
        T122.push_back(getEnergyVibr12Temp(energy));
    for (auto energy: EVibr3)
        T3.push_back(getEnergyVibr3Temp(energy));
    for (auto energy: E_tr_rot)
        T.push_back(energy*2*mass/(5*kB));

    //R12 =  rho*(EVibr12(T) - Evibr12(T12)) / tayVt2
    //+  rho*(EVibr12(T) - Evibr12(T12))/tayVv23 + rho*(EVibr12(T) - Evibr12(T12))/tayVv123
    //R3 =  rho*(EVibr3(T) - Evibr3(T3))/tayVv23 + rho*(EVibr3(T) - Evibr3(T3))/tayVv123
    pressure = tempU1*T*UniversalGasConstant/molMass;
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calc = [&](int& i)
    {
        double tauVibr = additionalSolver.TauVibr(T[i],pressure[i]);
        double tauVibr3 = additionalSolver.tauVibrVVLosev(T[i],pressure[i]);
        auto deltaE = (additionalSolver.EVibr12(0,T[i]) - additionalSolver.EVibr12(0,T122[i]));
        double r1 =tempU1[i] * (deltaE/tauVibr + deltaE/ tauVibr3);// + deltaE/ tauVibr3 );
        R_1[i]=r1;
        auto deltaE2 = (additionalSolver.EVibr3(0,T[i]) - additionalSolver.EVibr3(0,T3[i]));
        double r2 =tempU1[i]* (deltaE2/tauVibr3);// + deltaE2/tauVibr3);
        R_2[i]=r2;
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calc));
    futureWatcher.waitForFinished();
}
