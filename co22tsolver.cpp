#include "co22tsolver.h"

Co22TSolver::Co22TSolver(QObject *parent): AbstaractSolver(parent)
{

}

void Co22TSolver::prepareSolving()
{
    U1.resize(solParam.NumCell+2);
    U2.resize(solParam.NumCell+2);
    U3.resize(solParam.NumCell+2);
    U4.resize(solParam.NumCell+2);
    R.resize(solParam.NumCell+2);

    leftParam.density = leftParam.pressure /(UniversalGasConstant/molMass * leftParam.temp);
    double Cv = 5.0/2 * kB/mass;
    double gamma = (UniversalGasConstant/molMass + Cv)/Cv;
    leftParam.soundSpeed = sqrt(gamma*leftParam.pressure/leftParam.density);
    leftParam.velocity = solParam.Ma*leftParam.soundSpeed;
    leftParam.tempIntr = leftParam.temp;

    rightParam = additionalSolver.bondaryConditionPython(leftParam, solParam);
    rightParam.velocity = leftParam.density*leftParam.velocity/rightParam.density;
//////////////////
// rightParam.temp = 500.5;
// rightParam.tempIntr = 500.5;
// rightParam.pressure = 30.259216085384180;
// rightParam.density = rightParam.pressure /(UniversalGasConstant/molMass * rightParam.temp);
// rightParam.velocity = leftParam.density*leftParam.velocity/rightParam.density;
    /////////////

    double leftEvibr = additionalSolver.vibrEnergy(0,leftParam.temp);
    double leftFullEnergy = 5.0/2*kB*leftParam.temp/mass + leftEvibr;
    double rightEVibr = additionalSolver.vibrEnergy(0,rightParam.temp);
    double rightFullEnergy = 5.0/2*kB*rightParam.temp/mass + rightEVibr;

    for(auto i  = 1; i < solParam.NumCell+1; i++)
    {
        if(i < solParam.NumCell/2 +1)
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
    P.resize(solParam.NumCell+2);
    Q_v.resize(solParam.NumCell+2);
    Q_t.resize(solParam.NumCell+2);

    E__.resize(U1.size());
    R.resize(solParam.NumCell+1);
    rezultAfterPStart.resize(solParam.NumCell+1);
}

void Co22TSolver::solveFlux(Matrix U1L, Matrix U2L, Matrix pressureL, Matrix TvL, Matrix Tl , Matrix U1R, Matrix U2R, Matrix pressureR, Matrix TvR, Matrix Tr, Matrix EnergyFull)
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
        right.density = right_density[i];
        right.velocity = right_velocity[i];
        left.pressure = left_pressure[i];
        right.pressure = right_pressure[i];
        left.tempIntr = left_Tv[i];
        right.tempIntr = right_Tv[i];
        mutex.unlock();
        rezultAfterPStart[i] = additionalSolver.ExacRiemanSolver(left,right,solParam.Gamma);
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
        this->P[i] = P;
        double dt_dx = (tempR - tempL)/delta_h;
        double dtv_dx = (tempRtv - tempLtv)/delta_h;

        double qVibr = -additionalSolver.lambdaVibr2(Tx,Tv)* dtv_dx;
        double qTr = -additionalSolver.lambdaTr_Rot(Tx)*dt_dx;

        Q_v[i] = qVibr;
        Q_t[i] = qTr;
        double Etr_rot = 5.0/2*kB*Tx/mass;
        double entalpi = Etr_rot + energyVibr + point.pressure/point.density + pow(point.velocity,2)/2;

        F1[i] = point.density * point.velocity;
        F2[i] = point.density * point.velocity*point.velocity + point.pressure - P;
        F3[i] = point.density * point.velocity*entalpi - P*point.velocity + qVibr + qTr;
        F4[i] = point.density * point.velocity*energyVibr + qVibr;
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
        Matrix gammas;
        Matrix pressure;
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
        pressure = U1*T*UniversalGasConstant/molMass;
        pres = pressure;
        dt = additionalSolver.getTimeStepFull(velosity, U1, pressure, delta_h,solParam,gammas);
        timeSolvind.push_back(dt+timeSolvind.last());
        //pres = pressure;
        auto U1L = U1; U1L.removeLast();
        auto U2L = U2; U2L.removeLast();
        auto U3L = U3; U3L.removeLast();
        auto pressureL = pressure;
        auto TvL = Tv;
        auto Templ = T;
        Templ.removeLast();
        pressureL.removeLast();
        TvL.removeLast();

        auto U1R = U1; U1R.removeFirst();
        auto U2R = U2; U2R.removeFirst();
        auto U3R = U3; U3R.removeFirst();
        auto pressureR = pressure;
        auto TvR = Tv;
        auto TempR = T;
        TempR.removeFirst();
        pressureR.removeFirst();
        TvR.removeFirst();

        solveFlux(U1L, U2L, pressureL,TvL,Templ, U1R, U2R,pressureR, TvR, TempR, energyFull);
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

void Co22TSolver::setTypePlot(int i)
{
    solParam.typePlot = i;
    QVector<double> values, additionalValues;
    additionalValues.resize(x.size());
    switch (solParam.typePlot)
    {
    case 0: values = pres;break;
    case 1: values = U1; break;
    case 2: values = U2/U1;break;
    case 3: values = T; additionalValues = Tv; break;
    case 4: values = P;break;
    case 5:values = Q_t;break;
    case 6:values = Q_v;break;

    default: values = U1; break;
    }
    emit updateGraph(x, values, solParam.lambda);
    emit updateAdditionalGraph(x, additionalValues, solParam.lambda);
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
            return (i)  * energyVibrStepTemp + energyVibrStartTemp;
    return (EnergyVibr.size()-1) * energyVibrStepTemp + energyVibrStartTemp;
}

double Co22TSolver::getVibrTemp(double CVibr)
{
    for(auto i = 0 ; i < CvibrMass.size(); i++)
        if( CVibr < CvibrMass[i])
            return i  * CVibrStepTemp + CVibrStartTemp;
    return (CvibrMass.size()-1) * CVibrStepTemp + CVibrStartTemp;
}
