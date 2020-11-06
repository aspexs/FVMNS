#include "co22tsolver.h"

Co22TSolver::Co22TSolver(QObject *parent): AbstaractSolver(parent)
{
    QFile fileEnergy(QDir::currentPath() + "\\vibrEnergy.csv");
    QFile fileCVibr(QDir::currentPath() + "\\CVibr.csv");
    QFile fileAllEnergy(QDir::currentPath() + "\\allEnergy.csv");
    if(fileEnergy.open(QFile::ReadOnly) && fileCVibr.open(QFile::ReadOnly) && fileAllEnergy.open(QFile::ReadOnly))
    {
        QTextStream outEnergy(&fileEnergy);
        QStringList line = outEnergy.readLine().split(";");
        energyStartTemp = line[0].toDouble();
        EnergyVibr.push_back(line[1].toDouble());
        line = outEnergy.readLine().split(";");
        energyStepTemp = line[0].toDouble() - energyStartTemp;
        EnergyVibr.push_back(line[1].toDouble());

        while (!outEnergy.atEnd())
        {
           QStringList line = outEnergy.readLine().split(";");
           if(line.size() != 2)
               break;

           EnergyVibr.push_back(line[1].toDouble());
        }
        QTextStream outCVibr(&fileCVibr);

        line = outCVibr.readLine().split(";");
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

        QTextStream outAllEnergy(&fileAllEnergy);
        line = outAllEnergy.readLine().split(";");
        energyStartTemp = line[0].toDouble();
        Energy.push_back(line[1].toDouble());
        line = outAllEnergy.readLine().split(";");
        energyStepTemp = line[0].toDouble() - energyStartTemp;
        Energy.push_back(line[1].toDouble());
        while (!outAllEnergy.atEnd())
        {
           QStringList line = outAllEnergy.readLine().split(";");
           if(line.size() != 2)
               break;

           Energy.push_back(line[1].toDouble());
        }
    }
    fileEnergy.close();
    fileCVibr.close();
}

void Co22TSolver::prepareSolving()
{
    U1.resize(solParam.NumCell+2);
    U2.resize(solParam.NumCell+2);
    U3.resize(solParam.NumCell+2);
    U4.resize(solParam.NumCell+2);

    leftParam.density = leftParam.pressure /(UniversalGasConstant/molMass * leftParam.temp);
    double Cv = 5.0/2 * kB/mass + CvibrMass[(leftParam.temp - CVibrStartTemp)/CVibrStepTemp];
    double gamma = (UniversalGasConstant/molMass + Cv)/Cv;
    leftParam.soundSpeed = sqrt(gamma*leftParam.pressure/leftParam.density);
    leftParam.velocity = solParam.Ma*leftParam.soundSpeed;
    leftParam.tempIntr = leftParam.temp;
    rightParam = additionalSolver.BoundaryCondition[typeBC+1](leftParam, solParam);
    rightParam.velocity = leftParam.density*leftParam.velocity/rightParam.density;

    double leftEvibr = additionalSolver.vibrEnergy(0,leftParam.tempIntr);
    double leftFullEnergy = 5.0/2*kB*leftParam.temp/mass + leftEvibr;

    double rightEVibr = additionalSolver.vibrEnergy(0,rightParam.tempIntr);
    double rightFullEnergy = 5.0/2*kB*rightParam.temp/mass + rightEVibr;

    for(auto i  = 1; i < solParam.NumCell+1; i++)
    {
        if(i < solParam.NumCell/2+1)
        {
            U1[i] = leftParam.density;
            U2[i] = leftParam.density*leftParam.velocity;
            U3[i] = leftParam.density*(leftFullEnergy + pow(leftParam.velocity,2)/2);
            U4[i] = leftParam.density*leftEvibr;
        }

        else
        {
            U1[i] = rightParam.density;
            U2[i] =  rightParam.density*rightParam.velocity;
            U3[i] = rightParam.density*(rightFullEnergy + pow(rightParam.velocity,2)/2);
            U4[i] = rightParam.density*rightEVibr;
        }
    }

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
    rezultAfterPStart.resize(solParam.NumCell+1);
}

void Co22TSolver::solveFlux(Matrix U1L, Matrix U2L, Matrix pressureL, Matrix TvL, Matrix U1R, Matrix U2R, Matrix pressureR, Matrix TvR)
{
    left_density   = U1L;
    right_density  = U1R;
    left_velocity  = U2L / left_density;
    right_velocity = U2R / right_density;
    left_pressure  = pressureL;
    right_pressure = pressureR;
    left_Tv = TvL;
    right_Tv = TvR;
    Txx.clear();
    Tvv.clear();
    pres.clear();
    Txx.resize(U1.size());
    Tvv.resize(U1.size());
    pres.resize(U1.size());
    calcRiemanPStar();
    calcFliux();
}

void Co22TSolver::calcRiemanPStar()
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
        auto tempL = left.pressure /(UniversalGasConstant/molMass * left.density);
        auto tempR = right.pressure /(UniversalGasConstant/molMass * right.density);
        mutex.lock();
        auto CvL = 5.0/2 * kB/mass + CvibrMass[(tempL - CVibrStartTemp)/CVibrStepTemp];
        auto CvR = 5.0/2 * kB/mass + CvibrMass[(tempR - CVibrStartTemp)/CVibrStepTemp];
        mutex.unlock();
        auto GammaL = (UniversalGasConstant/molMass + CvL)/CvL;
        auto GammaR = (UniversalGasConstant/molMass + CvR)/CvR;
        rezultAfterPStart[i] = additionalSolver.ExacRiemanSolverCorrect(left,right, GammaL, GammaR);
        return;
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calcPStar));
    futureWatcher.waitForFinished();
}

void Co22TSolver::calcFliux()
{
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calcFlux = [this](int& i)
    {
        mutex.lock();
        auto point = rezultAfterPStart[i];
        double tempL = left_pressure[i]/(left_density[i]*UniversalGasConstant/molMass);
        double tempR = right_pressure[i]/(right_density[i]*UniversalGasConstant/molMass);
        auto du_dx = (right_velocity[i] - left_velocity[i])/delta_h;
        auto tempLtv = left_Tv[i];//EnergyVibr[(left_Tv[i] - energyStartTemp)/energyStepTemp];
        auto tempRtv = right_Tv[i];//EnergyVibr[(right_Tv[i] - energyStartTemp)/energyStepTemp];
        double Tx = point.pressure/(point.density*UniversalGasConstant/molMass);
        double Tv = (tempRtv + tempLtv)/2;
        pres[i] = point.pressure;
        Txx[i] = Tx;
        Tvv[i] = Tv;
        mutex.unlock();


        double etta = additionalSolver.shareViscosity[typeShareVisc](leftParam.temp, Tx,0,0);
        double Cvibr = 0;
        double Evibr = 0;
        double Evibr2 = 0;
        double Etr_rot = 5.0/2*kB*Tx/mass;

        if((Tx - CVibrStartTemp)/CVibrStepTemp < CvibrMass.size())
        {
            Cvibr = CvibrMass[(Tx - CVibrStartTemp)/CVibrStepTemp];
            Evibr2 =  EnergyVibr[(Tv - energyStartTemp)/energyStepTemp];
            Evibr =  EnergyVibr[(Tx - energyStartTemp)/energyStepTemp];
        }
        else
        {
            Cvibr = CvibrMass.last();
            Evibr = EnergyVibr.last();
            Evibr2 =EnergyVibr.last();
        }
        double zetta =additionalSolver.bulkViscosity[typeBulkVisc](Cvibr,Tx,point.density, point.pressure);
        double P =  (4.0/3*etta + zetta)*du_dx;
        double dt_dx = (tempR - tempL)/delta_h;
        double dtv_dx = (tempRtv - tempLtv)/delta_h;
        double qVibr = -additionalSolver.lambdaVibr2(Tx,Tv) * dtv_dx;
        double qTr = -additionalSolver.lambdaTr_Rot(Tx)*dt_dx;
        double tauVibr = additionalSolver.TauVibr(Tv, point.pressure);

        double entalpi = Etr_rot + Evibr + point.pressure/point.density + pow(point.velocity,2)/2;
        F1[i] = (point.density * point.velocity);
        F2[i] = (F1[i]*point.velocity + point.pressure -P);
        F3[i] = (F1[i]*entalpi - P*point.velocity + qVibr+qTr);
        F4[i] =  qVibr - point.density/tauVibr*(Evibr2 - Evibr);
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calcFlux));
    futureWatcher.waitForFinished();
}

void Co22TSolver::solve()
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
        auto EVibr = U4/U1;
        auto energyFull = U3/U1 - Matrix::POW(velosity,2)/2;
        auto E_tr_rot = energyFull - EVibr;
        Matrix T, Tv;
        Matrix gammas;
        Matrix pressure;
        for (auto energy: EVibr)
            Tv.push_back(getEnergyVibrTemp(energy));

        for (auto energy: E_tr_rot)
        {
            T.push_back(energy*2*mass/(5*kB));
            double Cv = 5.0/2 * kB/mass + CvibrMass[(T.last() - CVibrStartTemp)/CVibrStepTemp];
            gammas.push_back((UniversalGasConstant/molMass + Cv)/Cv);
        }
        pressure = U1*T*UniversalGasConstant/molMass;
        double dt = additionalSolver.TimeStepSolution[AdditionalSolver::TimeStepSolver::TS_FULL](velosity, U1, pressure, delta_h,solParam,gammas);
        timeSolvind.push_back(dt);

        auto U1L = U1; U1L.removeLast();
        auto U2L = U2; U2L.removeLast();
        auto U3L = U3; U3L.removeLast();
        auto pressureL = pressure;
        auto TvL = Tv;
        pressureL.removeLast();
        TvL.removeLast();

        auto U1R = U1; U1R.removeFirst();
        auto U2R = U2; U2R.removeFirst();
        auto U3R = U3; U3R.removeFirst();
        auto pressureR = pressure;
        auto TvR = Tv;
        pressureR.removeFirst();
        TvR.removeFirst();

        solveFlux(U1L, U2L, pressureL,TvL, U1R, U2R,pressureR, TvR);
        auto res = additionalSolver.SolveEvolutionExplFirstOrderForCO22(F1, F2,F3,F4, U1, U2,U3, U4,dt,delta_h);

        auto copyU1 = U1;
        auto copyU2 = U2;
        auto copyU3 = U3;
        auto copyU4 = U4;
        U1 = res[0];
        U2 = res[1];
        U3 = res[2];
        U4 = res[3];
        U1[0]=U1[1];   U1[U1.size() - 1] =U1[U1.size() - 2];
        U2[0]=solParam.typeLeftBorder*U2[1];
        U2[U2.size()-1]= solParam.typeRightBorder*U2[U2.size()-2];
        U3[0]=U3[1];   U3[U3.size() - 1] =U3[U3.size() - 2];
        U4[0]=U4[1];   U4[U4.size() - 1] =U4[U4.size() - 2];

        Txx.last() =  Txx[Txx.size() -2];
        Tvv.last() =  Tv[Tv.size() -2];
        pres.last() =  pres[pres.size() -2];

        if(i % solParam.PlotIter == 0)
        {
            setTypePlot(solParam.typePlot);
            emit updateTime(timeSolvind.last());
            QThread::msleep(200);
        }
        if(timeSolvind.last() >solParam.t_fin)
            break;
    }
}

void Co22TSolver::setTypePlot(int i)
{
    solParam.typePlot = i;
    QVector<double> values, additionalValues;
    switch (solParam.typePlot)
    {
    case 0: values = pres;break;
    case 1: values = U1; break;
    case 2: values = Tvv;break;
    case 3: values = Txx; additionalValues = Tvv; break;

    default: values = U1; break;
    }
    emit updateGraph(x, values, solParam.lambda);
    emit updateAdditionalGraph(x, additionalValues, solParam.lambda);
}

double Co22TSolver::getEnergyVibrTemp(double energy)
{
    for(auto i = 0 ; i < EnergyVibr.size(); i++)
        if( energy < EnergyVibr[i])
            return i  * energyStepTemp + energyStartTemp;
    return (EnergyVibr.size()-1) * energyStepTemp + energyStartTemp;
}

double Co22TSolver::getEnergyTemp(double energy)
{
    for(auto i = 0 ; i < Energy.size(); i++)
        if( energy < Energy[i])
            return i  * energyStepTemp + energyStartTemp;
    return (Energy.size()-1) * energyStepTemp + energyStartTemp;
}


double Co22TSolver::getVibrTemp(double CVibr)
{
    for(auto i = 0 ; i < CvibrMass.size(); i++)
        if( CVibr < CvibrMass[i])
            return i  * CVibrStepTemp + CVibrStartTemp;
    return (CvibrMass.size()-1) * CVibrStepTemp + CVibrStartTemp;
}