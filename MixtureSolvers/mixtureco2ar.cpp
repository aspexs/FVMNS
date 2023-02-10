#include "mixtureco2ar.h"
#include <QMessageBox>

MixtureCo2Ar::MixtureCo2Ar(QObject *parent) : AbstaractSolver(parent)
{
    // тут надо подготовить и считать из файла все предварительно расчитанные данные,
    // если они есть
    localTemp.resize(solParam.NumCell+1);
    mixture.initialize();
}

void MixtureCo2Ar::prepareSolving()
{
    U1.resize(solParam.NumCell + 2);
    U2.resize(solParam.NumCell + 2);
    U3.resize(solParam.NumCell + 2);
    U4.resize(solParam.NumCell + 2);
    U5.resize(solParam.NumCell + 2);
    U6.resize(solParam.NumCell + 2);
    R_1.resize(solParam.NumCell + 2);
    R_2.resize(solParam.NumCell + 2);

    // Инициализация
    tc_2::SpecificHeatDc heat;
    tc_2::EnergyDc leftEnergy;
    tc_2::EnergyDc rightEnergy;
    heat.initialize();
    leftEnergy.initialize();
    rightEnergy.initialize();

    // Привязка объектов друг к другу
    heat.link(mixture);
    leftEnergy.link(mixture);
    leftEnergy.link(heat);
    rightEnergy.link(mixture);
    rightEnergy.link(heat);

    // Расчет энергий !Важно - на единицу объема, а не массы. + rho * v^2 / 2!
    heat.compute(leftParam.t12, leftParam.t3);
    leftEnergy.compute(leftParam);
    heat.compute(rightParam.t12, rightParam.t3);
    rightEnergy.compute(rightParam);

    // Исправил !Использую плотности каждого из сортов, а не их концентрации!
    for (auto i  = 1; i < solParam.NumCell+1; i++)
    {
        if (i < solParam.NumCell / 3 + 1)
        {
            U1[i] = leftParam.rho[0];
            U2[i] = leftParam.rho[1];
            U3[i] = (leftParam.rho[0] + leftParam.rho[1]) * leftParam.v;
            U4[i] = leftEnergy.fullE();
            U5[i] = leftEnergy.vE12();
            U6[i] = leftEnergy.vE3();
        }
        else
        {
            U1[i] = rightParam.rho[0];
            U2[i] = rightParam.rho[1];
            U3[i] = (rightParam.rho[0] + rightParam.rho[1]) * rightParam.v;
            U4[i] = rightEnergy.fullE();
            U5[i] = rightEnergy.vE12();
            U6[i] = rightEnergy.vE3();
        }
    }
    prepareVectors();
}

void MixtureCo2Ar::solveFlux()
{
    calcRiemanPStar();
    calcFliux();
    calcR(U1, U2, U3, U4, U5, U6);
}

void MixtureCo2Ar::calcFliux()
{
    // CHECK Внес изменения, не уверен, будет ли это все распараллеливаться...
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calcFlux = [this](int& i)
    {
        // забираем результаты из расчета разрыва потока.
//        mutex.lock();
//        auto point = rezultAfterPStart[i];
//        double lT = Tl[i];
//        double rT = Tr[i];
//        auto dv_dx = (right_velocity[i] - left_velocity[i]) / delta_h;
//        auto tempLT12 = T12L[i];
//        auto tempRT12 = T12R[i];
//        auto tempLT3 = T3L[i];
//        auto tempRT3 = T3R[i];
//        double Tx =  point.pressure/(point.density*UniversalGasConstant/molMass);
//        double T12 = tempLT12;
//        double T3 = tempLT3;
//        mutex.unlock();

        // TODO Нет информации по плотностямсортовслева и справа,
        // надо добавить, cм. dx_dx (формула (21) методички)
        macroParam point = rezultAfterPStart[i];
        double dv_dx = (right_velocity[i] - left_velocity[i]) / delta_h;
        double dT_dx = (Tr[i] - Tl[i]) / delta_h;
        double dT12_dx = (T12R[i] - T12L[i]) / delta_h;
        double dT3_dx = (T3R[i] - T3L[i]) / delta_h;
        double dp_dx = (right_pressure[i] - left_pressure[i]) / delta_h;
        localTemp[i] = point.p / (point.rho[0] / mixture.mass()[0] +
                point.rho[1]/ mixture.mass()[1]) / tc_2::K_BOLTZMANN;

        // Инициализация
        tc_2::Computer computer;
        computer.initialize();

        // Привязка и расчет вектора потоков
        computer.link(mixture);
        computer.compute(point, dx_dx, dp_dx, dT_dx, dT12_dx, dT3_dx, dv_dx);

        // Обновляем значения
        F1[i] = computer.flowComputer().flow()[0];
        F2[i] = computer.flowComputer().flow()[1];
        F3[i] = computer.flowComputer().flow()[2];
        F4[i] = computer.flowComputer().flow()[3];
        F5[i] = computer.flowComputer().flow()[4];
        F6[i] = computer.flowComputer().flow()[5];
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
        double rightFulleftEnergy = 5.0/2*kB*T/mass + Evibr12 + Evibr3;
        U4[U4.size() - 1] = U2.last()*(rightFulleftEnergy + pow(U3.last()/U2.last(),2)/2);

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

void MixtureCo2Ar::calcR(const Matrix &U1, const Matrix &U2, const Matrix &U3,
                         const Matrix &U4, const Matrix &U5, const Matrix &U6)
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
