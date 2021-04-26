#include "oxygensolver.h"

OxygenSolver::OxygenSolver(QObject *parent) : AbstaractSolver(parent)
{

}

void OxygenSolver::prepareSolving()
{
    U1.resize(solParam.NumCell+2);
    U2.resize(solParam.NumCell+2);
    U3.resize(solParam.NumCell+2);
    U4.resize(solParam.NumCell+2);
    leftParam.gas = "O2";
    leftParam.density = leftParam.pressure /(UniversalGasConstant/molMass * leftParam.temp);
    double Cv = 5.0/2 * kB/mass + as.getC_V(leftParam.temp);
    double gamma = (UniversalGasConstant/molMass + Cv)/Cv;
    leftParam.soundSpeed = sqrt(gamma*leftParam.pressure/leftParam.density);
    leftParam.velocity = solParam.Ma*leftParam.soundSpeed;
    leftParam.tempIntr = leftParam.temp;
    //rightParam = additionalSolver.BoundaryCondition[AdditionalSolver::BC_PYTHON](leftParam, solParam);
    rightParam = additionalSolver.BoundaryCondition[AdditionalSolver::BC_RG2T](leftParam, solParam);
    //rightParam.velocity = leftParam.density*leftParam.velocity/rightParam.density;

    double leftEvibr = as.getVibrEnergy(leftParam.temp);
    double leftFullEnergy = 5.0/2*kB*leftParam.temp/mass + leftEvibr;


    double rightEVibr = as.getVibrEnergy(rightParam.tempIntr);
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
    prepareVectors();
}

void OxygenSolver::solve()
{

    auto y =as.getVibrEnergy(300);
    auto e = as.getTempFromEnergy(y);
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
        T.clear();
        Tv.clear();
        pres.clear();
        Matrix gammas;
        for (auto energy: EVibr)
        {
            auto TempTv = as.getTempFromEnergy(energy);
            Tv.push_back(TempTv);
        }
        for (auto energy: E_tr_rot)
        {
            T.push_back(energy*2*mass/(5*kB));
            double Cv = 5.0/2 * kB/mass + as.getC_V(T.last());
            gammas.push_back((UniversalGasConstant/molMass + Cv)/Cv);
        }
        pres = U1*T*UniversalGasConstant/molMass;
        dt = additionalSolver.TimeStepSolution[AdditionalSolver::TimeStepSolver::TS_FULL](velosity, U1, pres, delta_h,solParam,gammas);
        timeSolvind.push_back(dt);

        auto U1L = U1; U1L.removeLast();
        auto U2L = U2; U2L.removeLast();
        auto U3L = U3; U3L.removeLast();
        auto pressureL = pres;
        auto TvL = Tv;
        auto Templ = T;
        Templ.removeLast();
        pressureL.removeLast();
        TvL.removeLast();

        auto U1R = U1; U1R.removeFirst();
        auto U2R = U2; U2R.removeFirst();
        auto U3R = U3; U3R.removeFirst();
        auto pressureR = pres;
        auto TvR = Tv;
        auto TempR = T;
        TempR.removeFirst();
        pressureR.removeFirst();
        TvR.removeFirst();

        solveFlux(U1L, U2L, pressureL,TvL, U1R, U2R,pressureR, TvR);
        auto res = additionalSolver.SolveEvolutionExplFirstOrderForO2(F1, F2,F3,F4, U1, U2,U3, U4,dt,delta_h);

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
        pres[pres.size() - 1] =  pres[pres.size() -2];
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

void OxygenSolver::solveFlux(Matrix U1L, Matrix U2L, Matrix pressureL,Matrix TvL,
                             Matrix U1R, Matrix U2R, Matrix pressureR,Matrix TvR)
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

void OxygenSolver::calcRiemanPStar()
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
        auto CvL = 5.0/2 * kB/mass + as.getVibrEnergy(tempL);
        auto CvR = 5.0/2 * kB/mass + as.getVibrEnergy(tempR);
        auto GammaL = (UniversalGasConstant/molMass + CvL)/CvL;
        auto GammaR = (UniversalGasConstant/molMass + CvR)/CvR;
        rezultAfterPStart[i] = additionalSolver.ExacRiemanSolverCorrect(left,right, GammaL, GammaR, dt/delta_h);
        return;
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calcPStar));
    futureWatcher.waitForFinished();
}

void OxygenSolver::calcFliux()
{
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calcFlux = [this](int& i)
    {
        mutex.lock();
        auto point = rezultAfterPStart[i];
        double tempL = left_pressure[i]/(left_density[i]*UniversalGasConstant/molMass);
        double tempR = right_pressure[i]/(right_density[i]*UniversalGasConstant/molMass);
        auto du_dx = (right_velocity[i] - left_velocity[i])/delta_h;
        auto tempLtv = left_Tv[i];
        auto tempRtv = right_Tv[i];
        double Tx =  point.pressure/(point.density*UniversalGasConstant/molMass);
        double Tv = (tempRtv+tempLtv)/2;
        pres[i] = point.pressure;
        mutex.unlock();

        double etta = as.getEtta(Tx);//additionalSolver.shareViscosity[typeShareVisc](leftParam.temp, Tx,0,0);
        double Evibr = as.getVibrEnergy(Tx);
        double Evibr2 = as.getVibrEnergy(Tv);
        double Etr_rot = 5.0/2*kB*Tx/mass;

        double zetta = as.getZetta(Tx,etta);//additionalSolver.bulkViscosity[AdditionalSolver::BulkViscosity::BV_ONLY_RT_ROT](Cvibr,Tx,point.density, point.pressure);
        double P =  (4.0/3*etta + zetta)*du_dx;
        double dt_dx = (tempR - tempL)/delta_h;
        double dtv_dx = (tempRtv - tempLtv)/delta_h;
        double qVibr = -as.getLambdaVibr(Tx, Tv) * dtv_dx;
        double qTr = -as.getLambdaTr_Rot(Tx)*dt_dx;

        double tauVibr = as.getTauVibr(Tv, point.pressure);
        double entalpi = Etr_rot + Evibr2 + point.pressure/point.density + pow(point.velocity,2)/2;
        F1[i] = (point.density * point.velocity);
        F2[i] = (F1[i]*point.velocity + point.pressure - P);
        F3[i] = (F1[i]*entalpi - P*point.velocity + qVibr + qTr);
        auto deltaE = Evibr - Evibr2;
        double n = point.pressure/(kB*Tx);
        auto r = point.density/tauVibr* deltaE;
        F4[i] =  F1[i]*Evibr2 + qVibr - r*dt;
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calcFlux));
    futureWatcher.waitForFinished();
}
