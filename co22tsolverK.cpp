#include "co22tsolverK.h"

Co22TSolverK::Co22TSolverK(QObject *parent): AbstaractSolver(parent)
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
    }
    fileEnergy.close();
    fileCVibr.close();
}

void Co22TSolverK::prepareSolving()
{
    U1.resize(solParam.NumCell+2);
    U2.resize(solParam.NumCell+2);
    U3.resize(solParam.NumCell+2);
    U4.resize(solParam.NumCell+2);

    leftParam.density = leftParam.pressure /(UniversalGasConstant/molMass * leftParam.temp);
    double Cv = 5.0/2 * kB/mass;
    double gamma = (UniversalGasConstant/molMass + Cv)/Cv;
    solParam.Gamma = gamma;
    leftParam.soundSpeed = sqrt(gamma*leftParam.pressure/leftParam.density);
    leftParam.velocity = solParam.Ma*leftParam.soundSpeed;
    leftParam.tempIntr = leftParam.temp;

    rightParam = additionalSolver.BoundaryCondition[typeBC+1](leftParam, solParam);
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
    //leftEnergy = leftFullEnergy + pow(leftParam.velocity,2)/2  + leftParam.pressure/leftParam.density;
    //leftEnergy = rightFullEnergy + pow(rightParam.velocity,2)/2 + rightParam.pressure/rightParam.density;
    for(auto i  = 1; i < solParam.NumCell+1; i++)
    {
        if(i < 10)
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

    F11.resize(solParam.NumCell+1);
    F22.resize(solParam.NumCell+1);
    F33.resize(solParam.NumCell+1);
    F44.resize(solParam.NumCell+1);


    R.resize(solParam.NumCell+1);
    rezultAfterPStart.resize(solParam.NumCell+1);
}

void Co22TSolverK::solveFlux(const Matrix& U1, const Matrix& U2, const Matrix& U3, const Matrix& U4)
{
    calcHLLE(U1,U2,U3,U4);
    calcR(U1,U2,U3,U4);
}

void Co22TSolverK::calcRiemanPStar()
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
        auto Cv = 5.0/2 * kB/mass;
        auto Gamma = (UniversalGasConstant/molMass + Cv)/Cv;
        rezultAfterPStart[i] = additionalSolver.ExacRiemanSolver(left,right,Gamma);// GammaL, GammaR);
        return;
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calcPStar));
    futureWatcher.waitForFinished();
}

void Co22TSolverK::calcFliux()
{
//    QFutureWatcher<void> futureWatcher;
//    const std::function<void(int&)> calcFlux = [this](int& i)
//    {
//        mutex.lock();
//        auto point = rezultAfterPStart[i];
//        //point.temp = left_pressure[i]/(left_density[i]*UniversalGasConstant/molMass);
//        //point.pressure = left_pressure[i];
//        //point.density = left_density[i];
//        //point.density = leftParam.density*leftParam.velocity/point.velocity;

//        double tempL = Tl[i];// left_pressure[i]/(left_density[i]*UniversalGasConstant/molMass);
//        double tempR = Tr[i];//right_pressure[i]/(right_density[i]*UniversalGasConstant/molMass);
//        auto du_dx = (right_velocity[i] - left_velocity[i])/delta_h*2;
//        auto tempLtv = left_Tv[i];//EnergyVibr[(left_Tv[i] - energyStartTemp)/energyStepTemp];
//        auto tempRtv = right_Tv[i];//EnergyVibr[(right_Tv[i] - energyStartTemp)/energyStepTemp];
//        double Tx = point.pressure/(point.density*UniversalGasConstant/molMass);
//        //double Tx = left_pressure[i]/(left_density[i]*UniversalGasConstant/molMass);
//        double Tv = tempLtv;//tempRtv;//(tempRtv + tempLtv)/2;//getEnergyVibrTemp(energyVibr);
//        double energyVibr = AdditionalSolver::vibrEnergy(0,Tv);
//        Txx[i] = Tx;
//        pres[i] = point.pressure;
//        Tvv[i] = Tv;
//        T__[i] = point.velocity;
//        mutex.unlock();
//        //double omega11 = AdditionalSolver::getOmega11(Tx);
//        //double omega22 = AdditionalSolver::getOmega22(Tx);
//        //double zCO2Vibr = AdditionalSolver::ZCO2Vibr(Tv);
//        //double cVibr = AdditionalSolver::CVibr(Tv, zCO2Vibr);
//        //double etta = additionalSolver.shareViscosityOmega(0,Tx);//(5*kB*Tx) /(8*omega22);

//        double Evibr = AdditionalSolver::vibrEnergy(0,Tx);
//        double Etr_rot = 5.0/2*kB*Tx/mass;
//        //double zetta =additionalSolver.bulcViscosityOnlyTRRot(0,Tx);// (kB*Tx/betta)*0.16;
//        double P =0;//  (4.0/3*etta + zetta)*du_dx;
//       // double dt_dx = (tempR - tempL)/delta_h;
//        //double dtv_dx = (tempRtv - tempLtv)/delta_h;
//        //double qVibr = -additionalSolver.lambdaVibr2(Tx,Tv)* dtv_dx;// -(3.0*kB*Tx)/(8.0*omega11)*cVibr* dtv_dx;
//        //double qTr =-additionalSolver.lambdaTr_Rot(Tx)*dt_dx;//-((75.0*pow(kB,2)*Tx)/(32.0*mass*omega22) + (3.0*kB*Tx)/(8.0*omega11)*kB/mass)*dt_dx ;

//        double a=-18.19;
//        double b=40.47;
//        double c=0;
//        double d=0.00423;
//        //double tauVibr  = exp(a+b*pow(Tv,-1.0/3)+d/(pow(Tv,-1.0/3)))/(point.pressure/101325);
//        double tauVibr1 = exp(a+b*pow(Tx,-1.0/3)+d/(pow(Tx,-1.0/3)))/(point.pressure/101325);
//        double entalpi = Etr_rot + energyVibr + point.pressure/point.density + pow(point.velocity,2)/2;
//        auto deltaE = Evibr - energyVibr;
//        ////F1[i] = (point.density * point.velocity);
//        ////F2[i] = (F1[i]*point.velocity + point.pressure - P);
//        ////F3[i] = (F1[i]*entalpi - P*point.velocity + qVibr + qTr);

//        ////auto r = point.density/tauVibr1*deltaE*delta_h;
//        ////F4[i] =  F1[i]*energyVibr + qVibr - r;
//        F1[i]  = point.density * point.velocity;
//        F2[i]  = point.density * point.velocity*point.velocity + point.pressure;
//        F3[i]  = point.density * point.velocity*entalpi;
//        F4[i]  = point.density * point.velocity*energyVibr;
//        F11[i] = 0;
//        F22[i] = 0;//-P*point.velocity;
//        F33[i] = 0;//-additionalSolver.lambdaTr_Rot(Tx)*Tx - additionalSolver.lambdaVibr2(Tx,Tv)*Tv - P*pow(point.velocity,2);
//        F44[i] = 0;//-additionalSolver.lambdaVibr2(Tx,Tv)*Tv;
//          R[i] = point.density/tauVibr1*deltaE;

//    };
//    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calcFlux));
//    futureWatcher.waitForFinished();
}

void Co22TSolverK::solve()
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
        solveFlux(U1, U2, U3,U4);
        auto res = additionalSolver.SEEFOForCO2(F1, F2,F3,F4, U1, U2,U3, U4,dt,delta_h,F11,F22,F33,F44,R);

        U1 = res[0];
        U2 = res[1];
        U3 = res[2];
        U4 = res[3];


//        U1[0]=U1[1];                         U1[U1.size() - 1] = U1[U1.size() - 2];
//        U2[0]=solParam.typeLeftBorder*U2[1]; U2[U2.size()-1]= solParam.typeRightBorder*U2[U2.size()-2];
//        U3[0]=U3[1];                         U3[U3.size() - 1] = U3[U3.size() - 2];
//        U4[0]=U4[1];                         U4[U4.size() - 1] = U4[U4.size() - 2];

        //double rightEVibr = additionalSolver.vibrEnergy(0,rightParam.temp);
        //double rightFullEnergy = 5.0/2*kB*rightParam.temp/mass + rightEVibr;
          //U1[0]=U1[1];                         U1[U1.size() - 1] = rightParam.density;// U1[U1.size() - 2];
          //U2[0]=solParam.typeLeftBorder*U2[1]; U2[U2.size() - 1]= rightParam.density* rightParam.velocity;// solParam.typeRightBorder*U2[U2.size()-2];
          //U3[0]=U3[1];                         U3[U3.size() - 1] = rightParam.density*(rightFullEnergy + pow(rightParam.velocity,2)/2);
          //U4[0]=U4[1];                         U4[U4.size() - 1] = rightParam.density*rightEVibr;

          // сновис все, фиксируем давление. Пересчет энергии через фиксированное давление

          U1[0]=U1[1];                         U1[U1.size() - 1] = U1[U1.size() - 2];
          U2[0]=solParam.typeLeftBorder*U2[1]; U2[U2.size()-1]= solParam.typeRightBorder*U2[U2.size()-2];
          U4[0]=U4[1];                         U4[U4.size() - 1] = U4[U4.size() - 2];

          U3[0]=U3[1];                         //U3[U3.size() - 1] = //U3[U3.size() - 2];


                  auto T = rightParam.pressure* molMass/(UniversalGasConstant*U1.last());
                  auto Tv = getEnergyVibrTemp(U4.last()/U1.last());

                  double rightEVibr2 = additionalSolver.vibrEnergy(0,Tv);
                  double rightFullEnergy2 = 5.0/2*kB*T/mass + rightEVibr2;
                  U3[U3.size() - 1] = U1.last()*(rightFullEnergy2 + pow(U2.last()/U1.last(),2)/2);

        if(i % solParam.PlotIter == 0)
        {
            setTypePlot(solParam.typePlot);
            emit updateTime(timeSolvind.last());
            QThread::msleep(200);
        }
        if(timeSolvind.last() >solParam.t_fin)
        {
            setTypePlot(solParam.typePlot);
            emit updateTime(timeSolvind.last());
            QThread::msleep(200);
            break;
        }
    }
}

void Co22TSolverK::setTypePlot(int i)
{
    solParam.typePlot = i;
    QVector<double> values, additionalValues ;
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

void Co22TSolverK::calcHLLE(const Matrix& U1,const Matrix& U2,const Matrix& U3,const Matrix& U4)
{
    auto U1L = U1;              U1L.removeLast();
    auto U2L = U2;              U2L.removeLast();
    auto U3L = U3;              U3L.removeLast();
    auto U4L = U4;              U4L.removeLast();
    auto U1R = U1;              U1R.removeFirst();
    auto U2R = U2;              U2R.removeFirst();
    auto U3R = U3;              U3R.removeFirst();
    auto U4R = U4;              U4R.removeFirst();
    auto velocityL = U2L/U1L;
    auto velocityR = U2R/U1R;
    auto EVibrL = U4L/U1L;
    auto EVibrR = U4R/U1R;
    auto energyFullL = U3L/U1L - Matrix::POW(velocityL,2)/2;
    auto energyFullR = U3R/U1R - Matrix::POW(velocityR,2)/2;
    auto E_tr_rotL = energyFullL - EVibrL;
    auto E_tr_rotR = energyFullR - EVibrR;
    Matrix TVL, TVR, TL, TR,pressureL, pressureR;
    for (auto i = 0; i < EVibrL.size(); i++)
    {
        TVL.push_back(getEnergyVibrTemp(EVibrL[i]));
        TVR.push_back(getEnergyVibrTemp(EVibrR[i]));
    }
    for (auto i = 0; i < E_tr_rotL.size(); i++)
    {
        TL.push_back(E_tr_rotL[i]*2*mass/(5*kB));
        TR.push_back(E_tr_rotR[i]*2*mass/(5*kB));
    }
    pressureL = U1L*TL*UniversalGasConstant/molMass;
    pressureR = U1R*TR*UniversalGasConstant/molMass;
    auto coef = sqrt(solParam.Gamma*UniversalGasConstant/molMass);
    auto cL = Matrix::SQRT(TL)*coef;
    auto cR = Matrix::SQRT(TR)*coef; // скорости звука

    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calc = [&](int& i)
    {

         auto speed1 = {velocityL[i] + cL[i],velocityR[i] + cR[i], 0.0};
         auto speed2 = {velocityL[i] - cL[i],velocityR[i] - cR[i], 0.0};
         auto bR =*std::max_element(speed1.begin(), speed1.end());
         auto bL =*std::min_element(speed2.begin(), speed2.end());


         double du_dx = (velocityR[i] - velocityL[i])/delta_h;
         double dtv_dx = (TVR[i] - TVL[i])/delta_h;
         double dt_dx = (TR[i] - TL[i])/delta_h;

         double ettaL = additionalSolver.shareViscosityOmega(0,TL[i]);
         double zettaL =additionalSolver.bulcViscosityOnlyTRRot(0,TL[i]);
         float PL = (4.0/3*ettaL + zettaL)*du_dx;
         double Etr_rotL = 5.0/2*kB*TL[i]/mass;
         double entalpiL = Etr_rotL + EVibrL[i] + pressureL[i]/U1L[i] + pow(velocityL[i],2)/2;
         double qVibrL = -additionalSolver.lambdaVibr2(TL[i],TVL[i])* dtv_dx;
         double qTrL   =-additionalSolver.lambdaTr_Rot(TL[i])  *dt_dx;

         double ettaR = additionalSolver.shareViscosityOmega(0,TR[i]);
         double zettaR =additionalSolver.bulcViscosityOnlyTRRot(0,TR[i]);
         double PR = (4.0/3*ettaR + zettaR)*du_dx;
         double Etr_rotR = 5.0/2*kB*TR[i]/mass;
         double entalpiR = Etr_rotR + EVibrR[i] + pressureR[i]/U1R[i] + pow(velocityR[i],2)/2;
         double qVibrR = -additionalSolver.lambdaVibr2(TR[i],TVR[i])* dtv_dx;
         double qTrR   = -additionalSolver.lambdaTr_Rot(TR[i])  *dt_dx;


         auto F1L = U1L[i] * velocityL[i];
         auto F1R = U1R[i] * velocityR[i];
         auto F2L = F1L* velocityL[i] + pressureL[i] - PL;
         auto F2R = F1R* velocityR[i] + pressureR[i] - PR;
         auto F3L = F1L*entalpiL - PL*velocityL[i] + qVibrL + qTrL ;
         auto F3R = F1L*entalpiR - PR*velocityR[i] + qVibrR + qTrR ;
         auto F4L = F1L*EVibrL[i] + qVibrL;
         auto F4R = F1L*EVibrR[i] + qVibrR;


         auto res1 = (bR*F1L - bL*F1R + bR*bL*(U1R[i] - U1L[i]))/(bR-bL);
         auto res2 = (bR*F2L - bL*F2R + bR*bL*(U2R[i] - U2L[i]))/(bR-bL);
         auto res3 = (bR*F3L - bL*F3R + bR*bL*(U3R[i] - U3L[i]))/(bR-bL);
         auto res4 = (bR*F4L - bL*F4R + bR*bL*(U4R[i] - U4L[i]))/(bR-bL);
         F1[i] = res1;
         F2[i] = res2;
         F3[i] = res3;
         F4[i] = res4;
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calc));
    futureWatcher.waitForFinished();
}

void Co22TSolverK::calcR(const Matrix &U1, const Matrix &U2, const Matrix &U3, const Matrix &U4)
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

double Co22TSolverK::getEnergyVibrTemp(double energy)
{
    for(auto i = 0 ; i < EnergyVibr.size(); i++)
        if( energy < EnergyVibr[i])
            return i  * energyStepTemp + energyStartTemp;
    return (EnergyVibr.size()-1) * energyStepTemp + energyStartTemp;
}

double Co22TSolverK::getVibrTemp(double CVibr)
{
    for(auto i = 0 ; i < CvibrMass.size(); i++)
        if( CVibr < CvibrMass[i])
            return i  * CVibrStepTemp + CVibrStartTemp;
    return (CvibrMass.size()-1) * CVibrStepTemp + CVibrStartTemp;
}
