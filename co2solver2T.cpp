#include "co2solver2T.h"

CO2Solver2T::CO2Solver2T(QObject *parent): AbstaractSolver(parent)
{
    QFile fileVibrEnergy(QDir::currentPath() + "\\vibrEnergy.csv");
     //QFile fileFullEnergy(QDir::currentPath() + "\\allEnergy.csv");
     QFile fileCVibr(QDir::currentPath() + "\\CVibr.csv");
     if(fileVibrEnergy.open(QFile::ReadOnly) && fileCVibr.open(QFile::ReadOnly) )
     {
         QTextStream outEnergy(&fileVibrEnergy);
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

     }
     fileVibrEnergy.close();
     fileCVibr.close();
}

void CO2Solver2T::prepareSolving()
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
    rightParam = additionalSolver.BoundaryCondition[typeBC](leftParam, solParam);
    rightParam.velocity = leftParam.density*leftParam.velocity/rightParam.density;
    double leftEvibr = additionalSolver.vibrEnergy(0,leftParam.temp);
    double rightEVibr = additionalSolver.vibrEnergy(0,rightParam.temp);
    double leftFullEnergy = 5.0/2*kB*leftParam.temp/mass + leftEvibr;
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
            U2[i] = rightParam.density*rightParam.velocity;
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

void CO2Solver2T::solveFlux(Matrix U1L, Matrix U2L, Matrix pressureL, Matrix U1R, Matrix U2R, Matrix pressureR,  Matrix TvR, Matrix TvL)
{
    left_density   = U1L;
    right_density  = U1R;
    left_velocity  = U2L / left_density;

    right_velocity = U2R / right_density;
    left_pressure  = pressureL;
    right_pressure = pressureR;
    left_Tv = TvL;
    right_Tv = TvR;
    calcRiemanPStar();
    calcFliux();
}

void CO2Solver2T::calcRiemanPStar()
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
        auto CvL = 5.0/2 * kB/mass;// + CvibrMass[(tempL - CVibrStartTemp)/CVibrStepTemp];
        auto CvR = 5.0/2 * kB/mass;// + CvibrMass[(tempR - CVibrStartTemp)/CVibrStepTemp];
        mutex.unlock();
        auto GammaL = (UniversalGasConstant/molMass + CvL)/CvL;
        auto GammaR = (UniversalGasConstant/molMass + CvR)/CvR;
        rezultAfterPStart[i] = additionalSolver.ExacRiemanSolverCorrect(left,right, GammaL, GammaR,0);
        return;
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calcPStar));
    futureWatcher.waitForFinished();
}

void CO2Solver2T::calcFliux()
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
        mutex.unlock();
        //pres[i] = point.pressure;
        double Tx = point.pressure/(point.density*UniversalGasConstant/molMass);
        double Tv = tempLtv;
        double Evibr =  AdditionalSolver::vibrEnergy(0,Tv);
        double etta = 0;//additionalSolver.shareViscosity[typeShareVisc](leftParam.temp, Tx,0,0);
        double Cvibr = 0;//additionalSolver.CVibr(Tx, additionalSolver.ZCO2Vibr(Tx));
        //double E = 0;
        //if((Tx - CVibrStartTemp)/CVibrStepTemp < CvibrMass.size())
        //{
        //    Cvibr = CvibrMass[(Tx - CVibrStartTemp)/CVibrStepTemp];
        //    E =  Energy[(Tx - CVibrStartTemp)/CVibrStepTemp];
        //}
        //else
        //{
        //    Cvibr = CvibrMass.last();
        //    E = Energy.last();
        //}
         double energyVibr1 = additionalSolver.vibrEnergy(0,Tx);
         double Etr_rot = 5.0/2*kB*Tx/mass;
        double zetta =0;//additionalSolver.bulkViscosity[typeBulkVisc](Cvibr,Tx,point.density, point.pressure);
        double P =  0;//(4.0/3*etta + zetta)*du_dx;
        double lambda = 0;//additionalSolver.lambda(Tx, Cvibr);
        double dt_dx = (tempR - tempL)/delta_h;
        double q =   -lambda*dt_dx;
        double entalpi = Etr_rot + Evibr + point.pressure/point.density + pow(point.velocity,2)/2;

        double a=-18.19;
        double b=40.47;
        //double c=0;
        double d=0.00423;
        double tauVibr = exp(a+b*pow(Tv,-1.0/3)+d/(pow(Tv,-1.0/3)))/(point.pressure/101325);
        auto deltaE = (energyVibr1 - Evibr);
        auto r = point.density/tauVibr*deltaE*delta_h;
        F1[i] = point.density * point.velocity;
        F2[i] = F1[i]*point.velocity + point.pressure -P;
        F3[i] = F1[i]*entalpi - P*point.velocity + q;
        F4[i] = F1[i]*Evibr - r;
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calcFlux));
    futureWatcher.waitForFinished();
}


void CO2Solver2T::solve()
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
        T.clear();
        Tv.clear();
        Matrix velosity = U2/U1;
        auto EVibr = U4/U1;
        auto energyReal = U3/U1 - Matrix::POW(velosity,2)/2; // real full energy
        auto E_tr_rot = energyReal - EVibr;
        Matrix gammas;
        Matrix temperature;
        pres.clear();
        for (auto energy: E_tr_rot)
        {
            //T.push_back(getEnergyTemp(energy));
            T.push_back(energy*2*mass/(5*kB));
            double Cv = 5.0/2 * kB/mass;// + CvibrMass[(temperature.last() - CVibrStartTemp)/CVibrStepTemp];
            gammas.push_back((UniversalGasConstant/molMass + Cv)/Cv);
        }
        for (auto energy: EVibr)
        {
            auto TempTv = getEnergyVibrTemp(energy);
            Tv.push_back(TempTv);
        }
        //pres.last() =  pres[pres.size() -2];
        pres = U1*T*UniversalGasConstant/molMass;
        double dt = additionalSolver.TimeStepSolution[AdditionalSolver::TimeStepSolver::TS_FULL](velosity, U1, pres, delta_h,solParam,gammas);
        timeSolvind.push_back(dt+timeSolvind.last());
        auto U1L = U1; U1L.removeLast();
        auto U2L = U2; U2L.removeLast();
        auto U3L = U3; U3L.removeLast();
        auto pressureL = pres;
        auto TvL = Tv;
        pressureL.removeLast();
        TvL.removeLast();

        auto U1R = U1; U1R.removeFirst();
        auto U2R = U2; U2R.removeFirst();
        auto U3R = U3; U3R.removeFirst();
        auto pressureR = pres;
        auto TvR = Tv;
        TvR.removeFirst();
        pressureR.removeFirst();

        solveFlux(U1L, U2L, pressureL, U1R, U2R,pressureR,TvR,TvL);

        auto res = additionalSolver.SolveEvolutionExplFirstOrder(F1, F2,F3, F4,U1, U2,U3,U4,dt,delta_h);
        auto copyU1 = U1;
        auto copyU2 = U2;
        auto copyU3 = U3;
        U1 = res[0];
        U2 = res[1];
        U3 = res[2];
        U4 = res[3];
        U1[0]=U1[1];   U1[U1.size() - 1] =U1[U1.size() - 2];
        U2[0]=solParam.typeLeftBorder*U2[1];
        U2[U2.size()-1]= solParam.typeRightBorder*U2[U2.size()-2];
        U3[0]=U3[1];   U3[U3.size() - 1] =U3[U3.size() - 2];
        U4[0]=U4[1];   U4[U4.size() - 1] =U4[U4.size() - 2];

        if(i % solParam.PlotIter == 0)
        {
            setTypePlot(solParam.typePlot);
            emit updateTime(timeSolvind.last());
            QThread::msleep(200);
        }
        if(timeSolvind.last() >solParam.t_fin)
            break;
    }
    breaksolve = false;
}

void CO2Solver2T::setTypePlot(int i)
{
    solParam.typePlot = i;
    QVector<double> values, additionalValues;
    switch (solParam.typePlot)
    {
    case 0: values = pres;break;
    case 1: values = U1; break;
    case 2: values = U2/U1;break;
    case 3: values = T; additionalValues = Tv; break;
    default: values = U1; break;
    }
    emit updateGraph(x, values, solParam.lambda);
     emit updateAdditionalGraph(x, additionalValues, solParam.lambda);
}

//double CO2Solver2T::getEnergyTemp(double energy)
//{
//    //for(auto i = 0 ; i < Energy.size(); i++)
//    //    if( energy < Energy[i])
//    //        return i  * energyStepTemp + energyStartTemp;
//    //return (Energy.size()-1) * energyStepTemp + energyStartTemp;
//}

double CO2Solver2T::getVibrTemp(double CVibr)
{
    for(auto i = 0 ; i < CvibrMass.size(); i++)
        if( CVibr < CvibrMass[i])
            return i  * CVibrStepTemp + CVibrStartTemp;
    return (CvibrMass.size()-1) * CVibrStepTemp + CVibrStartTemp;
}

double CO2Solver2T::getEnergyVibrTemp(double energy)
{
    for(auto i = 0 ; i < EnergyVibr.size(); i++)
        if( energy < EnergyVibr[i])
            return (i+1)  * energyStepTemp + energyStartTemp;
    return (EnergyVibr.size()-1) * energyStepTemp + energyStartTemp;
}
