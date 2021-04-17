#include "co2solver.h"

CO2Solver::CO2Solver(QObject *parent): AbstaractSolver(parent)
{
}

void CO2Solver::prepareSolving()
{
    AdditionalSolver().TauVibr(1200, 20);
    auto e = additionalSolver.ZCO2Vibr(1200);
    U1.resize(solParam.NumCell+2);
    U2.resize(solParam.NumCell+2);
    U3.resize(solParam.NumCell+2);
    leftParam.density = leftParam.pressure /(UniversalGasConstant/molMass * leftParam.temp);
    double Cv = 5.0/2 * kB/mass + CvibrMass[(leftParam.temp - CVibrStartTemp)/CVibrStepTemp];
    double gamma = (UniversalGasConstant/molMass + Cv)/Cv;
    leftParam.soundSpeed = sqrt(gamma*leftParam.pressure/leftParam.density);
    leftParam.velocity = solParam.Ma*leftParam.soundSpeed;

    rightParam = additionalSolver.bondaryConditionPython(leftParam, solParam);
    rightParam.velocity = leftParam.density*leftParam.velocity/rightParam.density;
    double leftFullEnergy = 5.0/2*kB*leftParam.temp/mass + additionalSolver.vibrEnergy(0,leftParam.temp);
    double rightFullEnergy = 5.0/2*kB*rightParam.temp/mass + additionalSolver.vibrEnergy(0,rightParam.temp);

    for(auto i  = 1; i < solParam.NumCell+1; i++)
    {
        if(i < solParam.NumCell/2+1)
        {
            U1[i] = leftParam.density;
            U2[i] = leftParam.density*leftParam.velocity;
            U3[i] = leftParam.density*(leftFullEnergy + pow(leftParam.velocity,2)/2);
        }

        else
        {
            U1[i] = rightParam.density;
            U2[i] = rightParam.density*rightParam.velocity;
            U3[i] = rightParam.density*(rightFullEnergy + pow(rightParam.velocity,2)/2);
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
    U1[solParam.NumCell+1]=U1[solParam.NumCell];
    U2[solParam.NumCell+1]=solParam.typeRightBorder*U2[solParam.NumCell];
    U3[solParam.NumCell+1]=U3[solParam.NumCell];
    timeSolvind.push_back(0);

    for(int i = 0 ; i<  solParam.NumCell+1; i++)
        vectorForParallelSolving.push_back(i);
    F1.resize(solParam.NumCell+1);
    F2.resize(solParam.NumCell+1);
    F3.resize(solParam.NumCell+1);
    rezultAfterPStart.resize(solParam.NumCell+1);

}

void CO2Solver::solveFlux(Matrix U1L, Matrix U2L, Matrix pressureL, Matrix U1R, Matrix U2R, Matrix pressureR)
{
    left_density   = U1L;
    right_density  = U1R;
    left_velocity  = U2L / left_density;
    right_velocity = U2R / right_density;
    left_pressure  = pressureL;
    right_pressure = pressureR;
    calcRiemanPStar();
    calcFliux();
}

void CO2Solver::calcRiemanPStar()
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
      // auto tempL = left.pressure /(UniversalGasConstant/molMass * left.density);
      // auto tempR = right.pressure /(UniversalGasConstant/molMass * right.density);
       mutex.lock();
       auto CvL = 5.0/2 * kB/mass;// + additionalSolver.CVibr(tempL, additionalSolver.ZCO2Vibr(tempL));
       auto CvR = 5.0/2 * kB/mass;// + additionalSolver.CVibr(tempR, additionalSolver.ZCO2Vibr(tempR));
       mutex.unlock();
       auto GammaL = (UniversalGasConstant/molMass + CvL)/CvL;
       auto GammaR = (UniversalGasConstant/molMass + CvR)/CvR;
        rezultAfterPStart[i] = additionalSolver.ExacRiemanSolverCorrect(left,right, GammaL,GammaR);
        return;
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calcPStar));
    futureWatcher.waitForFinished();
}

void CO2Solver::calcFliux()
{
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calcFlux = [this](int& i)
    {
        mutex.lock();
        auto point = rezultAfterPStart[i];
        double tempL = left_pressure[i]/(left_density[i]*UniversalGasConstant/molMass);
        double tempR = right_pressure[i]/(right_density[i]*UniversalGasConstant/molMass);
        auto du_dx = (right_velocity[i] - left_velocity[i])/delta_h*2;
        mutex.unlock();
        double Tx = point.pressure/(point.density*UniversalGasConstant/molMass);

        double etta = additionalSolver.shareViscosityOmega(leftParam.temp, Tx,0,0);
        double Cvibr = additionalSolver.CVibr(Tx, additionalSolver.ZCO2Vibr(Tx));
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
        double zetta =additionalSolver.bulcViscosityOld2(Cvibr,Tx,point.density, point.pressure);
        double P =  (4.0/3*etta + zetta)*du_dx;
        double lambda = additionalSolver.lambda(Tx, Cvibr);
        double dt_dx = (tempR - tempL)/delta_h*2;
        double q =   -lambda*dt_dx;
        double entalpi = Etr_rot + energyVibr1 + point.pressure/point.density + pow(point.velocity,2)/2;
        F1[i] = point.density * point.velocity;
        F2[i] = F1[i]*point.velocity + point.pressure -P;
        F3[i] = F1[i]*entalpi - P*point.velocity + q;
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calcFlux));
    futureWatcher.waitForFinished();
}


void CO2Solver::solve()
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
        auto energyReal = U3/U1 - Matrix::POW(velosity,2)/2; // real full energy
        Matrix gammas;
        Matrix pressure;
        Matrix temperature;
        for (auto energy: energyReal)
        {
            temperature.push_back(getEnergyTemp(energy));
            double Cv = 5.0/2 * kB/mass+ CvibrMass[(temperature.last() - CVibrStartTemp)/CVibrStepTemp];
            gammas.push_back((UniversalGasConstant/molMass + Cv)/Cv);
        }
        pressure = U1*temperature*UniversalGasConstant/molMass;
        double dt = additionalSolver.TimeStepSolution[AdditionalSolver::TimeStepSolver::TS_FULL](velosity, U1, pressure, delta_h,solParam,gammas);
            timeSolvind.push_back(dt+timeSolvind.last());
        auto U1L = U1; U1L.removeLast();
        auto U2L = U2; U2L.removeLast();
        auto U3L = U3; U3L.removeLast();
        auto pressureL = pressure;
        pressureL.removeLast();

        auto U1R = U1; U1R.removeFirst();
        auto U2R = U2; U2R.removeFirst();
        auto U3R = U3; U3R.removeFirst();
        auto pressureR = pressure;
        pressureR.removeFirst();
        solveFlux(U1L, U2L, pressureL, U1R, U2R,pressureR);
        auto res = additionalSolver.SolveEvolutionExplFirstOrderForO2(F1, F2,F3,F4,U1, U2,U3,U4,dt,delta_h);
        auto copyU1 = U1;
        auto copyU2 = U2;
        auto copyU3 = U3;
        U1 = res[0];
        U2 = res[1];
        U3 = res[2];
        U1[0]=U1[1];   U1[U1.size() - 1] =U1[U1.size() - 2];
        U2[0]=solParam.typeLeftBorder*U2[1];
        U2[U2.size()-1]= solParam.typeRightBorder*U2[U2.size()-2];

        U3[0]=U3[1];
        auto T = rightParam.pressure* molMass/(UniversalGasConstant*U1.last());
        double rightEVibr2 = additionalSolver.vibrEnergy(0,T);
        double rightFullEnergy2 = 5.0/2*kB*T/mass + rightEVibr2;
        U3[U3.size() - 1] = U1.last()*(rightFullEnergy2 + pow(U2.last()/U1.last(),2)/2);
        //U3[0]=U3[1];   U3[U3.size() - 1] =U3[U3.size() - 2];

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

void CO2Solver::setTypePlot(int i)
{
    solParam.typePlot = i;
    QVector<double> values;
    switch (solParam.typePlot)
    {
    case 0:
    {
        QVector<double> temp;
        Matrix temperature;
        auto energyReal = U3/U1 - Matrix::POW(U2/U1,2)/2;
        for (auto energy: energyReal)
            temperature.push_back(getEnergyTemp(energy));
        values = U1*temperature*UniversalGasConstant/molMass;break;
    }
    case 1: values = U1; break;
    case 2: values = U2/U1;break;
    case 3:
    {
        QVector<double> temp;
        Matrix temperature;
        auto energyReal = U3/U1 - Matrix::POW(U2/U1,2)/2;
        for (auto energy: energyReal)
            temperature.push_back(getEnergyTemp(energy));
        values = temperature;break;
    }
    default: values = U1; break;
    }
    emit updateGraph(x, values, solParam.lambda);
}

double CO2Solver::getEnergyTemp(double energy)
{
    for(auto i = 0 ; i < Energy.size(); i++)
        if( energy < Energy[i])
            return (i)  * energyStepTemp + energyStartTemp;
    return (Energy.size()-1) * energyStepTemp + energyStartTemp;
}

double CO2Solver::getVibrTemp(double CVibr)
{
    for(auto i = 0 ; i < CvibrMass.size(); i++)
        if( CVibr < CvibrMass[i])
            return i  * CVibrStepTemp + CVibrStartTemp;
    return (CvibrMass.size()-1) * CVibrStepTemp + CVibrStartTemp;
}
