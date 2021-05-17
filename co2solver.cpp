#include "co2solver.h"

CO2Solver::CO2Solver(QObject *parent): AbstaractSolver(parent)
{
}

void CO2Solver::solveFlux(Matrix U1L, Matrix U2L, Matrix pressureL, Matrix U1R, Matrix U2R, Matrix pressureR, Matrix Tl, Matrix Tr)
{
    left_density   = U1L;
    right_density  = U1R;
    left_velocity  = U2L / left_density;
    right_velocity = U2R / right_density;
    left_pressure  = pressureL;
    right_pressure = pressureR;
    this->Tl = Tl;
    this->Tr = Tr;
    calcRiemanPStar();
    calcFliux();
}

void CO2Solver::calcFliux()
{
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calcFlux = [this](int& i)
    {
        mutex.lock();
        auto point = rezultAfterPStart[i];
        double tempL = Tl[i];
        double tempR = Tr[i];
        auto du_dx = (right_velocity[i] - left_velocity[i])/delta_h;
        mutex.unlock();
        double Tx = point.pressure/(point.density*UniversalGasConstant/molMass);
        double etta = additionalSolver.shareViscosityOmega(leftParam.temp, Tx,0,0);
        double Cvibr = additionalSolver.CVibr(Tx, additionalSolver.ZCO2Vibr(Tx));
        double energyVibr1 = additionalSolver.vibrEnergy(0,Tx);
        double Etr_rot = 5.0/2*kB*Tx/mass;
        double zetta = 0;//additionalSolver.bulcViscosityOld2(Cvibr,Tx,point.density, point.pressure);
        double P = (4.0/3*etta + zetta)*du_dx;
        double lambda = additionalSolver.lambda(Tx, Cvibr);
        double dt_dx = (tempR - tempL)/delta_h;
        double q =   -lambda*dt_dx;
        double entalpi = Etr_rot + energyVibr1 + point.pressure/point.density + pow(point.velocity,2)/2;
        F1[i] = point.density * point.velocity;
        F2[i] = F1[i]*point.velocity + point.pressure -P;
        F3[i] = F1[i]*entalpi - P*point.velocity + q;
         mutex.lock();
        Q_t[i] = q;
         B_v[i] = zetta;
          mutex.unlock();
        this->P[i] = P;
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
        double Cv = 5.0/2* kB/mass;
        T.clear();
        pres.clear();
        Ent.clear();
        for (auto energy: energyReal)
        {
            T.push_back(getEnergyTemp(energy));
            gammas.push_back((UniversalGasConstant/molMass + Cv)/Cv);
        }
        pres = U1*T*UniversalGasConstant/molMass;
        Ent = U3/U1 + pres/U1;
        double dt = additionalSolver.getTimeStepFull(velosity, U1, pres, delta_h,solParam,gammas);
            timeSolvind.push_back(dt+timeSolvind.last());
        auto U1L = U1; U1L.removeLast();
        auto U2L = U2; U2L.removeLast();
        auto U3L = U3; U3L.removeLast();
        auto pressureL = pres;
        pressureL.removeLast();
        auto Templ = T;
        Templ.removeLast();

        auto U1R = U1; U1R.removeFirst();
        auto U2R = U2; U2R.removeFirst();
        auto U3R = U3; U3R.removeFirst();
        auto pressureR = pres;
        pressureR.removeFirst();
        auto TempR = T;
        TempR.removeFirst();

        solveFlux(U1L, U2L, pressureL, U1R, U2R,pressureR,Templ,TempR);
        auto res = additionalSolver.SEEFOForCO2(F1, F2,F3,F4,U1, U2,U3,U4,dt,delta_h,{},{},{},{},R);
        auto copyU1 = U1;
        auto copyU2 = U2;
        auto copyU3 = U3;
        double error = res[5].first();
        U1 = res[0];
        U2 = res[1];
        U3 = res[2];
        double leftEvibr = additionalSolver.vibrEnergy(0,leftParam.temp);
        double leftFullEnergy = 5.0/2*kB*leftParam.temp/mass + leftEvibr;
        U1[0] = leftParam.density;
        U2[0] = leftParam.density*leftParam.velocity;
        U3[0] = leftParam.density*(leftFullEnergy + pow(leftParam.velocity,2)/2);
        U1[U1.size() - 1] =U1[U1.size() - 2];
        U2[U2.size()-1]= solParam.typeRightBorder*U2[U2.size()-2];
        auto T = rightParam.pressure* molMass/(UniversalGasConstant*U1.last());
        double rightEVibr2 = additionalSolver.vibrEnergy(0,T);
        double rightFullEnergy2 = 5.0/2*kB*T/mass + rightEVibr2;
        U3[U3.size() - 1] = U1.last()*(rightFullEnergy2 + pow(U2.last()/U1.last(),2)/2);

        Ent = U3/U1 + pres/U1;
        if(i % solParam.PlotIter == 0)
        {
            setTypePlot(solParam.typePlot);
            emit updateTime(timeSolvind.last(),error);
            QThread::msleep(200);
        }
        if(timeSolvind.last() >solParam.t_fin)
            break;
    }
    breaksolve = false;
}

double CO2Solver::getEnergyTemp(double energy)
{
    for(auto i = 0 ; i < Energy.size(); i++)
        if( energy < Energy[i])
            return (i-1)  * energyStepTemp + energyStartTemp;
    return (Energy.size()-1) * energyStepTemp + energyStartTemp;
}
