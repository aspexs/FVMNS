#include "mixtureco2ar.h"
#include <QMessageBox>
MixtureCo2Ar::MixtureCo2Ar(QObject *parent) : AbstaractSolver(parent)
{
    // тут надо подготовить и считать из файла все предварительно расчитанные данные,
    // если они есть
    localTemp.resize(solParam.NumCell+1);
}

void MixtureCo2Ar::prepareSolving()
{
    U1.resize(solParam.NumCell+2);
    U2.resize(solParam.NumCell+2);
    U3.resize(solParam.NumCell+2);
    U4.resize(solParam.NumCell+2);
    U5.resize(solParam.NumCell+2);
    U6.resize(solParam.NumCell+2);
    R_1.resize(solParam.NumCell+2);
    R_2.resize(solParam.NumCell+2);

    //структура leftParam содержит все что слева. ее надо корректно наполнить
    // возможно для смеси ее надо расширить.
    //leftParam.density ...
    double leftEvibr12 = additionalSolver.EVibr12(0,leftParam.temp);
    double leftEvibr3 = additionalSolver.EVibr3(0,leftParam.tempIntr);
    // особое вниамние тут, ибо полная энергия тут явно поменялась.
    double leftFullEnergy = 5.0/2*kB*leftParam.temp/mass + leftEvibr12 + leftEvibr3;

    // затем надо применить граничные условия для вычисления того, что справа.
    // этот шаг я возьму на себя в ближайшую неделю (Илья).

    // далее заполняем аналогичный параметр rightParam

    double rightEVibr12 = additionalSolver.EVibr12(0,rightParam.temp);
    double rightEVibr3 = additionalSolver.EVibr3(0,rightParam.temp);
    // тут тоже
    double rightFullEnergy = 5.0/2*kB*rightParam.temp/mass + rightEVibr12 + rightEVibr3;

    for(auto i  = 1; i < solParam.NumCell+1; i++)
    {
        if(i < solParam.NumCell/3 +1)
        {
            // подразумеваю что ты тут будешь заполнять через удельные концентрации
            //U1[i] = n_co2

            // Так же надо быть внимательным с плотностями, я так понимаю для U5 и U6 они другие
            // если жэжто так, поправь пожалуйста.
            U2[i] = leftParam.density;
            U3[i] = leftParam.density*leftParam.velocity;
            U4[i] = leftParam.density*(leftFullEnergy + pow(leftParam.velocity,2)/2);
            U5[i] = leftParam.density*leftEvibr12;
            U6[i] = leftParam.density*leftEvibr3;
        }

        else
        {
             //U1[i] = n_co2
            U2[i] = rightParam.density;
            U3[i] = rightParam.density*rightParam.velocity;
            U4[i] = rightParam.density*(rightFullEnergy + pow(rightParam.velocity,2)/2);
            U5[i] = rightParam.density*rightEVibr12;
            U6[i] = rightParam.density*rightEVibr3;
        }
    }
    prepareVectors();
}

void MixtureCo2Ar::solveFlux()
{
    calcRiemanPStar();
    calcFliux();
    calcR(U1,U2,U3,U4,U5,U6);
}

void MixtureCo2Ar::calcFliux()
{
    // наша главная функция и головная боль
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calcFlux = [this](int& i)
    {
        // забираем результаты из расчета разрыва потока.
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

        localTemp[i] = Tx;
        double dt_dx = (tempR - tempL)/delta_h;
        double dt12_dx = (tempRT12 - tempLT12)/delta_h;
        double dt3_dx = (tempRT3 - tempLT3)/delta_h;

        // тут блок который аналогичен тому, что было в трехтемпературном приближении
        double energyVibr12 =  additionalSolver.EVibr12(0,T12);
        double energyVibr3 =  additionalSolver.EVibr3(0,T3);

        // тут блок который изменился
        double etta     = 0; // FIXME additionalSolver.shareViscosityOmega(0,Tx);
        double zetta    = 0; // FIXME additionalSolver.bulcViscosityOld2(0,Tx);
        double lambdaTR = 0; // FIXME additionalSolver.lambdaTr_Rot(Tx);
        double qDiff    = 0; // FIXME дополнить
        double lambda12 = 0; // FIXME additionalSolver.Lambda12(Tx,T12);
        double lambda3  = 0; // FIXME additionalSolver.Lambda3(Tx,T3);
        double Etr_rot  = 0; // FIXME 5.0/2*kB*Tx/mass;
        // это уже конечные соотношения
        double P        = (4.0/3*etta + zetta)*du_dx;
        double qTr      = -lambdaTR*dt_dx;
        double qVibr12  = -lambda12* dt12_dx;
        double qVibr3   = -lambda3* dt3_dx;

        double E12_diff = 0; // Соответствует самой правой части уравнения 5
        double E3_diff  = 0; // Соответствует самой правой части уравнения 6
        double entalpi  = Etr_rot + energyVibr12 + energyVibr3 + point.pressure/point.density + pow(point.velocity,2)/2;
        F1[i] = 0; // дополнить
        F2[i] = point.density * point.velocity;
        F3[i] = point.density * point.velocity*point.velocity + point.pressure - P;
        F4[i] = point.density * point.velocity*entalpi - P*point.velocity + qVibr12 + qVibr3 + qTr + qDiff;
        F5[i] = point.density * point.velocity*energyVibr12 + qVibr12 - E12_diff;
        F6[i] = point.density * point.velocity*energyVibr3 + qVibr3 - E3_diff;
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calcFlux));
    futureWatcher.waitForFinished();
}

void MixtureCo2Ar::solve()
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

        auto EVibr12 = U5/U2;
        auto EVibr3 = U6/U2;
        auto energyFull = U4/U2 - Matrix::POW(velosity,2)/2;
        auto E_tr_rot = energyFull - EVibr12 - EVibr3;

        T.clear();
        T12.clear();
        T3.clear();
        Tv.clear();
        pres.clear();

        Matrix gammas;
        // FIXME Надо посмотреть как было в просто трехтемпературном приближении и выставить Т
        T = {0};
        gammas = {0};
        pres = U2*T*UniversalGasConstant/molMass;
        dt = additionalSolver.getTimeStepFull(velosity, U2, pres, delta_h,solParam, gammas);
        timeSolvind.push_back(dt+timeSolvind.last());

        // TODO проверить массивы, те ли значения берем
        auto U2L = U2; U2L.removeLast();
        auto U3L = U3; U3L.removeLast();
        auto U4L = U4; U4L.removeLast();
        auto pressureL = pres;
        auto T12L = T12;
        auto T3L = T3;
        auto Templ = T;
        Templ.removeLast();
        pressureL.removeLast();
        left_pressure = pressureL;
        T12L.removeLast();
        T3L.removeLast();
        left_density   = U2L;
        left_velocity  = U3L / left_density;
        this->T12L = T12L;
        this->T3L = T3L;
        this->Tl = Templ;

        auto U2R = U2; U2R.removeFirst();
        auto U3R = U3; U3R.removeFirst();
        auto U4R = U4; U4R.removeFirst();
        auto pressureR = pres;
        auto T12R = T12;
        auto T3R = T3;
        auto TempR = T;
        TempR.removeFirst();
        pressureR.removeFirst();
        T12R.removeFirst();
        T3R.removeFirst();
        right_pressure = pressureR;
        right_density  = U2R;
        right_velocity = U3R / right_density;
        this->T12R = T12R;
        this->T3R = T3R;
        this->Tr = TempR;

        solveFlux();
        auto res = additionalSolver.SEEFOForCO23T(F1, F2, F3, F4, F5, F6, U1, U2, U3, U4, U5, U5, dt, delta_h, R_1, R_2);

        auto copyU1 = U1;
        auto copyU2 = U2;
        auto copyU3 = U3;
        auto copyU4 = U4;
        auto copyU5 = U5;
        auto copyU6 = U6;
        U1 = res[0];
        U2 = res[1];
        U3 = res[2];
        U4 = res[3];
        U5 = res[4];
        U6 = res[5];
        error = res[6].first();

        U1[0]=U1[1];                         U1[U1.size() - 1] = U1[U1.size() - 2];
        U2[0]=solParam.typeLeftBorder*U2[1]; U2[U2.size()-1]= solParam.typeRightBorder*U2[U2.size()-2];
        U4[0]=U4[1];                         U4[U4.size() - 1] = U4[U4.size() - 2];
        U5[0]=U5[1];                         U5[U5.size() - 1] = U5[U5.size() - 2];
        U3[0]=U3[1];                         U6[U6.size() - 1] = U6[U6.size() - 2];

        auto T = rightParam.pressure* molMass/(UniversalGasConstant*U1.last());
        double Evibr12 = U5.last()/U1.last();
        double Evibr3 = U6.last()/U1.last();
        // FIXME исправить полную энергию.
        double rightFullEnergy = 5.0/2*kB*T/mass + Evibr12 + Evibr3;
        U4[U4.size() - 1] = U2.last()*(rightFullEnergy + pow(U3.last()/U2.last(),2)/2);

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

void MixtureCo2Ar::calcR(const Matrix &U1, const Matrix &U2, const Matrix &U3, const Matrix &U4, const Matrix &U5, const Matrix &U6)
{
    Matrix pressure;
    auto tempU2 = U2;
    auto tempU6 = U6;
    auto tempU5 = U5;
    auto EVibr12 = tempU5/U1;
    auto EVibr3 = tempU6/U1;
    // FIXME тут в идеале надо брать тесмпературу как в других приближениях, ожнако это затратно,
    // будем брать просто в ячейке
    pressure = tempU2*localTemp*UniversalGasConstant/molMass;
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calc = [&](int& i)
    {
        double tauVibr = additionalSolver.TauVibr(localTemp[i],pressure[i]);
        double tauVibr3 = additionalSolver.tauVibrVVLosev(localTemp[i],pressure[i]);
        auto deltaE = (additionalSolver.EVibr12(0,localTemp[i]) - EVibr12[i]);
        double r1 =tempU2[i] *(deltaE/tauVibr + 2*deltaE/ tauVibr3);

        auto deltaE2 = (additionalSolver.EVibr3(0,localTemp[i]) - EVibr3[i]);
        double r2 =tempU2[i]* (2*deltaE2/tauVibr3);
        mutex.lock();
        R_1[i]=r1;
        R_2[i]=r2;
        mutex.unlock();
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calc));
    futureWatcher.waitForFinished();
}
