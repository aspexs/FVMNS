#include "mixtureco2ar.h"

MixtureCo2Ar::MixtureCo2Ar()
{
    U.resize(SYSTEM_ORDER);
    F.resize(SYSTEM_ORDER);
    R.resize(SYSTEM_ORDER);
    hlleF.resize(SYSTEM_ORDER);
}

void MixtureCo2Ar::initialize(const SolverParams& init)
{
    // Обновляем настройки
    solParam = init;

    // Изменение размеров массивов
    for (int i = 0; i < SYSTEM_ORDER; ++i)
    {
        U[i].fill(0.0, init.nCell + 2);
        F[i].fill(0.0, init.nCell + 2);
        R[i].fill(0.0, init.nCell + 2);
        hlleF[i].fill(0.0, init.nCell + 1);
    }
    points.resize(init.nCell + 2);
    k.fill(0.0, init.nCell + 2);
    a.fill(0.0, init.nCell + 2);

    // Вектор для распараллеливания
    for(int i = 0; i < init.nCell + 2; i++)
    {
        parAll_v.push_back(i);
        if (i < init.nCell + 1)
        {
            parIn_v.push_back(i);
        }
    }

    // Подготовка таблиц температур
    computeT.initialize(init.t_min, init.t_max, init.t_n);

    // Расчет энергий
    EnergyDc lE;
    EnergyDc rE;
    lE.initialize();
    rE.initialize();
    lE.compute(init.lPoint);
    rE.compute(init.rPoint);

    // Исправил !Использую плотности каждого из сортов, а не их концентрации!
    for (auto i  = 0; i < solParam.nCell + 2; i++)
    {
        if (i < solParam.nCell / 2)
        {
            U[0][i] = init.lPoint.rho[0];
            U[1][i] = init.lPoint.rho[1];
            U[2][i] = (init.lPoint.rho[0] + init.lPoint.rho[1]) *
                    init.lPoint.v;
            U[3][i] = (init.lPoint.rho[0] + init.lPoint.rho[1]) *
                    (lE.fullE() + qPow(init.lPoint.v, 2.0) / 2.0);
            U[4][i] = init.lPoint.rho[0] * lE.vE12();
            U[5][i] = init.lPoint.rho[0] * lE.vE3();
            points[i] = init.lPoint;
        }
        else
        {
            U[0][i] = init.rPoint.rho[0];
            U[1][i] = init.rPoint.rho[1];
            U[2][i] = (init.rPoint.rho[0] + init.rPoint.rho[1]) *
                    init.rPoint.v;
            U[3][i] = (init.rPoint.rho[0] + init.rPoint.rho[1]) *
                    (rE.fullE() + qPow(init.rPoint.v, 2.0) / 2.0);
            U[4][i] = init.rPoint.rho[0] * rE.vE12();
            U[5][i] = init.rPoint.rho[0] * rE.vE3();
            points[i] = init.rPoint;
        }
    }
}

void MixtureCo2Ar::solve()
{
    for (int i = 0; i < solParam.nIter; ++i)
    {
        // Обновляем показатель адиабаты, скорость звука, временной шаг
        updateAK();
        updateDt();

        // Вычисляем вектор поточных и релаксационных членов + HLLE
        computeF();
        computeHlleF();
        computeR();

        // Обновляем вектор консервативных переменных
        step();

        // Возврат к макропараметрам и обновление естественного кр. условия
        updateMacroParam();
        updateBoundCond();
    }
    updateAK();
}

void MixtureCo2Ar::computeF()
{
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calc = [this](int& i)
    {
        //for (int i = 0; i < solParam.nCell + 2; ++i)
        //{
        // Временные переменные
        MacroParam p0, p1, p2;

        // Забираем известные макропараметры в (.)-ах [i-1], [i], [i+1]
        mutex.lock();
        p1 = points[i];
        if (i == 0)
        {
            p0 = p1;
            p2 = points[i + 1];
        }
        else if (i == solParam.nCell + 1)
        {
            p0 = points[i - 1];
            p2 = p1;
        }
        else
        {
            p0 = points[i - 1];
            p2 = points[i + 1];
        }
        mutex.unlock();

        // Вспомогательные величины (кончентрации и молярные доли)
        QVector<double> n0 = {p0.rho[0] / Mixture::mass(0),
                              p0.rho[1] / Mixture::mass(1)};
        QVector<double> n2 = {p2.rho[0] / Mixture::mass(0),
                              p2.rho[1] / Mixture::mass(1)};
        QVector<double> x0 = {n0[0] / (n0[0] + n0[1]),
                              n0[1] / (n0[0] + n0[1])};
        QVector<double> x2 = {n2[0] / (n2[0] + n2[1]),
                              n2[1] / (n2[0] + n2[1])};

        // Рассчитываем производные в точке i
        double dv_dx = (p2.v - p0.v) / (2.0 * solParam.dx);
        double dT_dx = (p2.t - p0.t) / (2.0 * solParam.dx);
        double dT12_dx = (p2.t12 - p0.t12) / (2.0 * solParam.dx);
        double dT3_dx = (p2.t3 - p0.t3) / (2.0 * solParam.dx);
        double dp_dx = (p2.p - p0.p) / (2.0 * solParam.dx);
        QVector<double> dx_dx = {(x2[0] - x0[0]) / (2.0 * solParam.dx),
                                 (x2[1] - x0[1]) / (2.0 * solParam.dx)};

        // Расчет поточных членов
        FlowMembersDc computer;
        computer.initialize();
        computer.compute(p1, dx_dx, dp_dx, dT_dx, dT12_dx, dT3_dx, dv_dx);

        // Обновляем значения вектора поточных членов в (.) [i]
        mutex.lock();
        for (int j = 0; j < SYSTEM_ORDER; ++j)
        {
            F[j][i] = computer.flow()[j];
        }
        mutex.unlock();
        //}
    };
    futureWatcher.setFuture(QtConcurrent::map(parAll_v, calc));
    futureWatcher.waitForFinished();
}

void MixtureCo2Ar::computeHlleF()
{
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calc = [this](int& i)
    {
        // Временные переменные
        double f_hlle, f0, f1, u0, u1, a0, a1, v0, v1, rho0, rho1, avg_k;
        double b0, b1, eta, avg_a, avg_v;

        // Забираем известные макропараметры в (.)-ах [i], [i+1]
        mutex.lock();
        a0 = qPow(a[i], 2.0);
        a1 = qPow(a[i + 1], 2.0);
        v0 = points[i].v;
        v1 = points[i + 1].v;
        rho0 = qSqrt(points[i].rho[0] + points[i].rho[1]);
        rho1 = qSqrt(points[i + 1].rho[0] + points[i + 1].rho[1]);
        avg_k = (k[i] + k[i + 1]) / 2.0;
        mutex.unlock();

        for (int j = 0; j < SYSTEM_ORDER; ++j)
        {
            // Забираем известные макропараметры в (.)-ах [i], [i+1]
            mutex.lock();
            f0 = F[j][i];
            f1 = F[j][i + 1];
            u0 = U[j][i];
            u1 = U[j][i + 1];
            mutex.unlock();

            // Расчитываем сигнальные скорости
            eta = (avg_k - 1.0) / 2.0 * rho0 * rho1 / qPow(rho0 + rho1, 2.0);
            avg_a = qSqrt((rho0 * a0 + rho1 * a1) / (rho0 + rho1) +
                          eta * qPow(v1 - v0, 2.0));
            avg_v = (rho0 * v0 + rho1 * v1) / (rho0 + rho1);

            // Расчет потока на стыке ячеек
            b0 = qMin(avg_v - avg_a, 0.0);
            b1 = qMax(avg_v + avg_a, 0.0);
            f_hlle = (b1 * f0 - b0 * f1 + b1 * b0 * (u1 - u0)) / (b1 - b0);

            // Обновляем значение [j] вектора поточных членов в (.) [i]
            mutex.lock();
            hlleF[j][i] = f_hlle;
            mutex.unlock();
        }
    };
    futureWatcher.setFuture(QtConcurrent::map(parIn_v, calc));
    futureWatcher.waitForFinished();
}

void MixtureCo2Ar::step()
{
    // Инициализация
    double dF = 0.0;
    error = 0.0;

    // Находим ошибку, обновляем вектор консервативных переменных
    for (int j = 0; j < SYSTEM_ORDER; ++j)
    {
        for (int i = 1; i < solParam.nCell + 1; ++i)
        {
            dF = hlleF[j][i] - hlleF[j][i - 1];
            U[j][i] += dt * (R[j][i] - dF / solParam.dx);
            error += qAbs(dF);
        }
    }
}

void MixtureCo2Ar::computeR()
{
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calc = [this](int& i)
    {
        // Временные переменные
        MacroParam point;

        // Берем данные извне
        mutex.lock();
        point = points[i];
        mutex.unlock();

        // Расчет времени релаксации
        double tauVT = Mixture::tauVTCO2(point.t, point.p);
        double tauVV = Mixture::tauVVCO2(point.t, point.p);

        // Расчет значений колебательных энергий
        EnergyDc e0;
        EnergyDc e1;
        e0.initialize();
        e1.initialize();
        e1.compute(point.t, point.t);
        e0.compute(point.t12, point.t3);

        // Расчет скоростей релаксации
        double r1 = (point.rho[0] + point.rho[1]) * (e1.vE12() - e0.vE12()) *
                (1.0 / tauVT + 2.0 / tauVV);
        double r2 = (point.rho[0] + point.rho[1]) * (e1.vE3() - e0.vE3()) *
                2.0 / tauVV;

        // Возвращаем новые значения
        mutex.lock();
        R[4][i] = r1;
        R[5][i] = r2;
        mutex.unlock();
    };
    futureWatcher.setFuture(QtConcurrent::map(parAll_v, calc));
    futureWatcher.waitForFinished();
}

void MixtureCo2Ar::updateMacroParam()
{
    // Распараллеливание не имеет смысла
    for (int i = 1; i < solParam.nCell + 1; ++i)
    {
        double U3 = U[3][i] - 0.5 * qPow(U[2][i], 2.0) / (U[0][i] + U[1][i]);
        points[i].rho[0] = U[0][i];
        points[i].rho[1] = U[1][i];
        points[i].v = U[2][i] / (U[0][i] + U[1][i]);
        computeT.compute(points[i], U[4][i] / U[0][i], U[5][i] / U[0][i],
                U3 / (U[0][i] + U[1][i]));
        points[i].t = computeT.T();
        points[i].t12 = computeT.T12();
        points[i].t3 = computeT.T3();
        points[i].p = K_BOLTZMANN * computeT.T() *
                (points[i].rho[0] / Mixture::mass(0) +
                points[i].rho[1] / Mixture::mass(1));
    }
}

void MixtureCo2Ar::updateAK()
{
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calc = [this](int& i)
    {
        // Временные переменные
        MacroParam point;
        double U3;

        // Берем данные извне
        mutex.lock();
        U3 = U[3][i] - 0.5 * qPow(U[2][i], 2.0) / (U[0][i] + U[1][i]);
        point = points[i];
        mutex.unlock();

        // Расчет скорости звука и показателя адиабаты
        double temp_k = 1.0 + points[i].p / U3;
        double temp_a = qSqrt(temp_k * point.p / (point.rho[0] + point.rho[1]));

        // Возвращаем новые значения
        mutex.lock();
        k[i] = temp_k;
        a[i] = temp_a;
        mutex.unlock();
    };
    futureWatcher.setFuture(QtConcurrent::map(parAll_v, calc));
    futureWatcher.waitForFinished();
}

void MixtureCo2Ar::updateDt()
{
    // Вспомогательная скорость
    double full_v = 0.0;
    double max_v  = 0.0;

    // Находим максимальную скорость распространения возмущения в потоке
    for (int i = 0; i < solParam.nCell + 2; ++i)
    {
        full_v = a[i] + points[i].v;
        if (full_v > max_v)
        {
            max_v = full_v;
        }
    }

    // Критерий Куранта-Фридрихса-Леви
    dt = solParam.cfl * solParam.dx / max_v;
}

void MixtureCo2Ar::updateBoundCond()
{
    // Продолжаем по непрерывности и сохраняем давление
    //double p = points.last().p;
    points.last() = points[solParam.nCell];
    //points.last().p = p;

    // Регулируем давление за ударной волной
    //double dRhoV = points.last().v * (points.last().rho[0] +
    //        points.last().rho[1]) - points[0].v * (points[0].rho[0] +
    //        points[0].rho[1]);
    //points.last().p += 0.5 * solParam.k * points[0].v * dRhoV;
    points.last().v = (points[0].rho[0] + points[0].rho[1]) * points[0].v /
            (points.last().rho[0] + points.last().rho[1]);
}

QVector<QVector<double>> MixtureCo2Ar::saveMacroParams()
{
    QVector<QVector<double>> table;
    table.resize(11);
    for (int i = 0; i < solParam.nCell + 2; ++i)
    {
        table[0].push_back(points[i].rho[0] / Mixture::mass(0));
        table[1].push_back(points[i].rho[1] / Mixture::mass(1));
        table[2].push_back(points[i].p);
        table[3].push_back(points[i].v);
        table[4].push_back(points[i].t);
        table[5].push_back(points[i].t12);
        table[6].push_back(points[i].t3);
        table[7].push_back(k[i]);
        table[8].push_back(a[i]);
        table[9].push_back(dt);
        table[10].push_back(error);
    }
    return table;
}

QVector<QVector<double>> MixtureCo2Ar::saveU()
{
    QVector<QVector<double>> table;
    for (int i = 0; i < SYSTEM_ORDER; ++i)
    {
        table.push_back(U[i]);
    }
    return table;
}

QVector<QVector<double>> MixtureCo2Ar::saveF()
{
    QVector<QVector<double>> table;
    for (int i = 0; i < SYSTEM_ORDER; ++i)
    {
        table.push_back(F[i]);
    }
    return table;
}

QVector<QVector<double>> MixtureCo2Ar::saveHlleF()
{
    QVector<QVector<double>> table;
    for (int i = 0; i < SYSTEM_ORDER; ++i)
    {
        table.push_back(hlleF[i]);
    }
    return table;
}

QVector<QVector<double>> MixtureCo2Ar::saveR()
{
    QVector<QVector<double>> table;
    for (int i = 0; i < SYSTEM_ORDER; ++i)
    {
        table.push_back(R[i]);
    }
    return table;
}
