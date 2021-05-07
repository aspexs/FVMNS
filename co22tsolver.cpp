#include "co22tsolver.h"
Co22TSolver::Co22TSolver(QObject *parent): AbstaractSolver(parent)
{

}

void Co22TSolver::solveFlux(Matrix U1L, Matrix U2L, Matrix pressureL, Matrix TvL, Matrix Tl , Matrix U1R, Matrix U2R, Matrix pressureR, Matrix TvR, Matrix Tr)
{
    left_density   = U1L;
    right_density  = U1R;
    left_velocity  = U2L / left_density;
    right_velocity = U2R / right_density;
    left_pressure  = pressureL;
    right_pressure = pressureR;
    this->Tl = Tl;
    this->Tr = Tr;
    left_Tv = TvL;
    right_Tv = TvR;

    calcRiemanPStar();
    calcFliux();
    calcR(U1,U2,U3,U4);
}

void Co22TSolver::calcFliux()
{
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calcFlux = [this](int& i)
    {
        mutex.lock();
        auto point = rezultAfterPStart[i];
        double tempL = Tl[i];
        double tempR = Tr[i];
       auto du_dx = (right_velocity[i] - left_velocity[i])/delta_h;
        auto tempLtv = left_Tv[i];
        auto tempRtv = right_Tv[i];
        double Tx =  point.pressure/(point.density*UniversalGasConstant/molMass);
        double Tv = tempLtv;
        double energyVibr =  additionalSolver.vibrEnergy(0,Tv);
        mutex.unlock();

        double etta = additionalSolver.shareViscosityOmega(0,Tx);
        double zetta =additionalSolver.bulcViscosityOnlyTRRot(0,Tx);

        double P = (4.0/3*etta + zetta)*du_dx;
        double dt_dx = (tempR - tempL)/delta_h;
        double dtv_dx = (tempRtv - tempLtv)/delta_h;
        double qVibr = -additionalSolver.lambdaVibr2(Tx,Tv)* dtv_dx;
        double qTr = -additionalSolver.lambdaTr_Rot(Tx)*dt_dx;
        double Etr_rot = 5.0/2*kB*Tx/mass;
        double entalpi = Etr_rot + energyVibr + point.pressure/point.density + pow(point.velocity,2)/2;

        F1[i] = point.density * point.velocity;
        F2[i] = point.density * point.velocity*point.velocity + point.pressure - P;
        F3[i] = point.density * point.velocity*entalpi - P*point.velocity + qVibr + qTr;
        F4[i] = point.density * point.velocity*energyVibr + qVibr;
        Q_v[i] = qVibr;
        Q_t[i] = qTr;
        //Ent[i] = entalpi;
        this->P[i] = P;
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
        T.clear();
        Tv.clear();
        pres.clear();
        Matrix gammas;
        for (auto energy: EVibr)
        {
            auto TempTv = getEnergyVibrTemp(energy);
            Tv.push_back(TempTv);
        }
        double Cv = 5.0/2 * kB/mass;
        for (auto energy: E_tr_rot)
        {
            T.push_back(energy*2*mass/(5*kB));
            gammas.push_back((UniversalGasConstant/molMass + Cv)/Cv);
        }

        pres = U1*T*UniversalGasConstant/molMass;
        dt = additionalSolver.getTimeStepFull(velosity, U1, pres, delta_h,solParam,gammas);
        timeSolvind.push_back(dt+timeSolvind.last());

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

        solveFlux(U1L, U2L, pressureL,TvL,Templ, U1R, U2R,pressureR, TvR, TempR);
        auto res = additionalSolver.SEEFOForCO2(F1, F2,F3,F4, U1, U2,U3, U4,dt,delta_h,{},{},{},{},R);

        auto copyU1 = U1;
        auto copyU2 = U2;
        auto copyU3 = U3;
        auto copyU4 = U4;
        U1 = res[0];
        U2 = res[1];
        U3 = res[2];
        U4 = res[3];
        R = res[4];
        error = res[5].first();

        ////    Вариант расчета глупый, идет простой снос значений предпоследней ячеки
        //U1[0]=U1[1];   U1[U1.size() - 1] =U1[U1.size() - 2];
        //U2[0]=solParam.typeLeftBorder*U2[1];
        //U2[U2.size()-1]=solParam.typeRightBorder*U2[U2.size()-2];
        //U3[0]=U3[1];   U3[U3.size() - 1] =U3[U3.size() - 2];
        //U4[0]=U4[1];   U4[U4.size() - 1] =U4[U4.size() - 2];

        ////    Вариант расчета при фиксированном правом давлении
         U1[0]=U1[1];                         U1[U1.size() - 1] = U1[U1.size() - 2];
         U2[0]=solParam.typeLeftBorder*U2[1]; U2[U2.size()-1]= solParam.typeRightBorder*U2[U2.size()-2];
         U4[0]=U4[1];                         U4[U4.size() - 1] = U4[U4.size() - 2];

         U3[0]=U3[1];                         //U3[U3.size() - 1] = //U3[U3.size() - 2];

         auto T = rightParam.pressure* molMass/(UniversalGasConstant*U1.last());
         auto Tv = getEnergyVibrTemp(U4.last()/U1.last());

         double rightEVibr2 = additionalSolver.vibrEnergy(0,Tv);
         double rightFullEnergy2 = 5.0/2*kB*T/mass + rightEVibr2;
         U3[U3.size() - 1] = U1.last()*(rightFullEnergy2 + pow(U2.last()/U1.last(),2)/2);
         Ent = U3/U1 + pres/U1;

        ////    Вариант расчета при фиксированной правой части
        ///double rightEVibr = additionalSolver.vibrEnergy(0,rightParam.temp);
        ///double rightFullEnergy = 5.0/2*kB*rightParam.temp/mass + rightEVibr;
        ///U1[0]=U1[1];                         U1[U1.size() - 1] = rightParam.density;// U1[U1.size() - 2];
        ///U2[0]=solParam.typeLeftBorder*U2[1]; U2[U2.size() - 1]= rightParam.density* rightParam.velocity;// solParam.typeRightBorder*U2[U2.size()-2];
        ///U3[0]=U3[1];                         U3[U3.size() - 1] = rightParam.density*(rightFullEnergy + pow(rightParam.velocity,2)/2);
        ///U4[0]=U4[1];                         U4[U4.size() - 1] = rightParam.density*rightEVibr;
        ///R[0]=R[1];   R[R.size() - 1] =R[R.size() - 2];

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

void Co22TSolver::calcR(const Matrix &U1, const Matrix &U2, const Matrix &U3, const Matrix &U4)
{
    Matrix TV, T,pressure;
    auto tempU1 = U1;
    auto tempU2 = U2;
    auto velocity = tempU2/U1;
    auto tempU3 = U3;
    auto tempU4 = U4;
    auto EVibr = tempU4/U1;
    auto energyFull = tempU3/U1 - Matrix::POW(velocity,2)/2;
    auto E_tr_rot = energyFull - EVibr;
    for (auto i = 0; i < EVibr.size(); i++)
        TV.push_back(getEnergyVibrTemp(EVibr[i]));
    for (auto i = 0; i < E_tr_rot.size(); i++)
        T.push_back(E_tr_rot[i]*2*mass/(5*kB));

    pressure = tempU1*T*UniversalGasConstant/molMass;
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calc = [&](int& i)
    {
        double tauVibr = additionalSolver.TauVibr(T[i],pressure[i]);
        auto deltaE = (additionalSolver.vibrEnergy(0,T[i]) - EVibr[i]);
        R[i]=tempU1[i]/tauVibr*deltaE;
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calc));
    futureWatcher.waitForFinished();
}

double Co22TSolver::getEnergyVibrTemp(double energy)
{
    for(auto i = 0 ; i < EnergyVibr.size(); i++)
        if( energy < EnergyVibr[i])
            return (i-1)  * energyVibrStepTemp + energyVibrStartTemp;
    return (EnergyVibr.size()-1) * energyVibrStepTemp + energyVibrStartTemp;
}
