#include "mixture.h"

///////////////////////////////////////////////////////////////////////////////
/// class BorderCondition
///////////////////////////////////////////////////////////////////////////////

void BorderCondition::compute(const MacroParam& lP, const int& flag)
{
    // Порядок имеет значение
    initialize(lP, flag);
    computeT();
    computeRho();
    computeP();
    computeV();
}

const MacroParam& BorderCondition::rP() const
{
    return rP_;
}

void BorderCondition::initialize(const MacroParam& lP, const int& flag)
{
    energy_.initialize();
    energy_.compute(lP);
    lRho_ = lP.rho[0] + lP.rho[1];
    y0_ = lP.rho[0] / lRho_;
    y1_ = lP.rho[1] / lRho_;
    alpha_ = y0_ / Mixture::mass(0) + y1_ / Mixture::mass(1);
    lE_ = energy_.fullE();
    lP_ = lP;
    flag_ = flag;
}

double BorderCondition::computeD(const double& t)
{
    return qPow(lP_.p + lRho_ * qPow(lP_.v, 2.0), 2.0) -
            4.0 * K_BOLTZMANN * t * alpha_ * qPow(lRho_ * lP_.v, 2.0);
}

double BorderCondition::computeRho(const double& t)
{
    return (lP_.p + lRho_ * qPow(lP_.v, 2.0) + flag_ * qSqrt(computeD(t))) /
            (2.0 * K_BOLTZMANN * t * alpha_);
}

double BorderCondition::computeF(const double& t)
{
    // Временный набор макропараметров
    MacroParam p;
    p.t = t;
    p.t12 = t;
    p.t3 = t;
    p.rho[0] = y0_;
    p.rho[1] = y1_;

    // Расчет макропараметров за УВ для данной температуры
    energy_.compute(p);
    double rRho = computeRho(t);
    double rV = lRho_ * lP_.v / rRho;
    double rP = rRho * alpha_ * K_BOLTZMANN * t;
    double lF = lRho_ * lP_.v * (lE_ + 0.5 * qPow(lP_.v, 2.0)) + lP_.p * lP_.v;
    double rF = rRho * rV * (energy_.fullE() + 0.5 * qPow(rV, 2.0)) + rP * rV;

    // Возврат отклонения от нуля
    return lF - rF;
}

void BorderCondition::computeT()
{
    // Границы области поиска
    double minT = 0.0;
    double maxT = qPow(lP_.p + lRho_ * qPow(lP_.v, 2.0), 2.0) /
            (4.0 * alpha_ * K_BOLTZMANN * qPow(lP_.v * lRho_, 2.0));

    // Текущее положение и значение нелинейной функции
    double t = 0.0;
    double f = 0.0;

    // Осуществляем метод бисекции
    while (maxT - minT > EPSILON)
    {
        t = 0.5 * (minT + maxT);
        f = computeF(t);
        if (f < 0.0)
        {
            maxT = t;
        }
        else
        {
            minT = t;
        }
    }

    // Обновляем значения температуры на правой границе T = T12 = T3
    rP_.t = 0.5 * (minT + maxT);
    rP_.t12 = rP_.t;
    rP_.t3 = rP_.t;
}

void BorderCondition::computeP()
{
    rP_.p = (rP_.rho[0] / Mixture::mass(0) + rP_.rho[1] / Mixture::mass(1)) *
            K_BOLTZMANN * rP_.t;
}

void BorderCondition::computeRho()
{
    rP_.rho[0] = y0_ * computeRho(rP_.t);
    rP_.rho[1] = y1_ * computeRho(rP_.t);
}

void BorderCondition::computeV()
{
    rP_.v = lRho_ * lP_.v / (rP_.rho[0] + rP_.rho[1]);
}

///////////////////////////////////////////////////////////////////////////////
/// class MixtureCo2Ar
///////////////////////////////////////////////////////////////////////////////

MixtureCo2Ar::MixtureCo2Ar()
{
    U.resize(SYSTEM_ORDER);
    F.resize(SYSTEM_ORDER);
    R.resize(SYSTEM_ORDER);
    hlleF.resize(SYSTEM_ORDER);
}

void MixtureCo2Ar::initialize(const MacroParam& lP, const QString& name,
                              const int& numT)
{
    // Расчет равновесных значений за УВ
    BorderCondition bc;
    bc.compute(lP);

    // Повторное использование
    dt = 0.0;
    time = 0.0;
    currIter = 0;
    shockPos = 0;
    model = numT;

    // Изменение размеров массивов
    for (int i = 0; i < SYSTEM_ORDER; ++i)
    {
        U[i].fill(0.0, N_CELL);
        F[i].fill(0.0, N_CELL);
        R[i].fill(0.0, N_CELL);
        hlleF[i].fill(0.0, N_CELL - 1);
    }
    points.resize(N_CELL);
    k.fill(0.0, N_CELL);
    a.fill(0.0, N_CELL);
    dx.fill(0.0, N_CELL);

    // Инициализация массива длин ячеек
    double alpha = (DEGREE + 1.0) * GRID_L / qPow(N_CELL, (DEGREE + 1.0)) / DX;
    for(int i = 0; i < N_CELL; i++)
    {
        dx[i] = DX * (1.0 + alpha * qPow(i, DEGREE));
    }

    // Подготовка векторов для распараллеливания
    for(int i = 0; i < N_CELL; i++)
    {
        // Векторы для распараллеливания
        parN_v.push_back(i);
        if (i < N_CELL - 1)
        {
            parN1_v.push_back(i);
        }
    }

    // Подготовка таблиц температур и энергий
    computeT.readFromFile(name, T_NUM, Y_NUM);

    // Расчет энергий
    EnergyDc lE, rE;
    lE.initialize();
    rE.initialize();
    lE.compute(lP);
    rE.compute(bc.rP());

    // Инициализация сетки макропараметров (начальное условие)
    for (int i = 0; i < N_CELL; i++)
    {
        if (i <= SHOCK_POS)
        {
            U[0][i] = lP.rho[0];
            U[1][i] = lP.rho[1];
            U[2][i] = (lP.rho[0] + lP.rho[1]) * lP.v;
            U[3][i] = (lP.rho[0] + lP.rho[1]) *
                    (lE.fullE() + qPow(lP.v, 2.0) / 2.0);
            U[4][i] = lP.rho[0] * lE.vE12();
            U[5][i] = lP.rho[0] * lE.vE3();
            points[i] = lP;
        }
        else
        {
            U[0][i] = bc.rP().rho[0];
            U[1][i] = bc.rP().rho[1];
            U[2][i] = (bc.rP().rho[0] + bc.rP().rho[1]) * bc.rP().v;
            U[3][i] = (bc.rP().rho[0] + bc.rP().rho[1]) *
                    (rE.fullE() + qPow(bc.rP().v, 2.0) / 2.0);
            U[4][i] = bc.rP().rho[0] * rE.vE12();
            U[5][i] = bc.rP().rho[0] * rE.vE3();
            points[i] = bc.rP();
        }
    }
}

void MixtureCo2Ar::solve()
{
    // Готовим progress bar
    ProgressBar bar;
    bar.initialize(FINISH_TIME);

    // Осуществляем итерационный процесс
    while (currIter < MAX_N_ITER && time < FINISH_TIME)
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

        // Возврат к основным макропараметрам
        updateMacroParam();

        // Обновляем progress bar, счетчик и таймер
        bar.update(dt);
        ++currIter;
        time += dt;
    }
    updateAK();
    findShockPos();
}

void MixtureCo2Ar::computeF()
{
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calc = [this](int& i)
    {
        // Временные переменные
        MacroParam p0, p1, p2;
        double temp_dx = 0.0;

        // Забираем известные макропараметры в (.)-ах [i-1], [i], [i+1]
        mutex.lock();
        p1 = points[i];
        switch (i)
        {
            case 0:
                p0 = p1;
                p2 = points[i + 1];
                temp_dx = 1.5 * dx[i] + 0.5 * dx[i + 1];
                break;
            case N_CELL - 1:
                p0 = points[i - 1];
                p2 = p1;
                temp_dx = 1.5 * dx[i] + 0.5 * dx[i - 1];
                break;
            default:
                p0 = points[i - 1];
                p2 = points[i + 1];
                temp_dx = dx[i] + 0.5 * (dx[i - 1] + dx[i + 1]);
                break;
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
        double dv_dx = (p2.v - p0.v) / temp_dx;
        double dT_dx = (p2.t - p0.t) / temp_dx;
        double dlnT_dx = (qLn(p2.t) - qLn(p0.t)) / temp_dx;
        double dT12_dx = (p2.t12 - p0.t12) / temp_dx;
        double dT3_dx = (p2.t3 - p0.t3) / temp_dx;
        double dlnp_dx = (qLn(p2.p) - qLn(p0.p)) / temp_dx;
        QVector<double> dx_dx = {(x2[0] - x0[0]) / temp_dx,
                                 (x2[1] - x0[1]) / temp_dx};

        // Расчет поточных членов
        FlowMembersDc computer;
        computer.initialize();
        computer.compute(p1, dx_dx, dlnp_dx, dT_dx, dlnT_dx, dT12_dx, dT3_dx,
                         dv_dx);

        // Обновляем значения вектора поточных членов в (.) [i]
        mutex.lock();
        for (int j = 0; j < SYSTEM_ORDER; ++j)
        {
            F[j][i] = computer.flow()[j];
        }
        mutex.unlock();
    };
    futureWatcher.setFuture(QtConcurrent::map(parN_v, calc));
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
        avg_v = (rho0 * v0 + rho1 * v1) / (rho0 + rho1);
        avg_k = 0.5 * (k[i] + k[i + 1]);
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

            // Обновляем значение [j] вектора поточных членов в (.) [i-1]
            mutex.lock();
            hlleF[j][i] = f_hlle;
            mutex.unlock();
        }
    };
    futureWatcher.setFuture(QtConcurrent::map(parN1_v, calc));
    futureWatcher.waitForFinished();
}

void MixtureCo2Ar::step()
{
    // Находим ошибку, обновляем вектор консервативных переменных
    for (int j = 0; j < SYSTEM_ORDER; ++j)
    {
        for (int i = 1; i < N_CELL - 1; ++i)
        {
            U[j][i] += (R[j][i] - (hlleF[j][i] - hlleF[j][i - 1]) / dx[i]) * dt;
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

        // Расчет времен релаксации VT и VV при столкновениях CO2-CO2, CO2-Ar
        double tauVTCO2 = Mixture::tauVTCO2CO2(point.t, point.rho[0]);
        double tauVVCO2 = Mixture::tauVVCO2CO2(point.t, point.rho[0]);
        double tauVTAr = Mixture::tauVTCO2Ar(point.t, point.rho[1]);
        double tauVVAr = Mixture::tauVVCO2Ar(point.t, point.rho[1]);

        // Расчет значений колебательных энергий
        EnergyDc et, et12, et3;
        et.initialize();
        et12.initialize();
        et3.initialize();
        et.compute(point.t, point.t);
        et12.compute(point.t12, point.t12);
        et3.compute(point.t3, point.t3);

        // Расчет скоростей релаксации
        double r12 = point.rho[0] *
                ((et.vE12() - et12.vE12()) * (1.0 / tauVTCO2 + 1.0 / tauVTAr) +
                 (et3.vE12() - et12.vE12()) * (2.0 / tauVVCO2 + 2.0 / tauVVAr));
        double r3 = point.rho[0] *
                (et12.vE3() - et3.vE3()) * (2.0 / tauVVCO2 + 2.0 / tauVVAr);

        // Возвращаем новые значения
        mutex.lock();
        R[4][i] = r12;
        R[5][i] = r3;
        mutex.unlock();
    };
    futureWatcher.setFuture(QtConcurrent::map(parN_v, calc));
    futureWatcher.waitForFinished();
}

void MixtureCo2Ar::updateMacroParam()
{
    // Без кинетической энергии
    double U3 = 0.0;

    // Распараллеливание не имеет смысла
    for (int i = 1; i < N_CELL - 1; ++i)
    {
        // Находим внутреннюю энергию
        U3 = U[3][i] - 0.5 * qPow(U[2][i], 2.0) / (U[0][i] + U[1][i]);
        points[i].rho[0] = U[0][i];
        points[i].rho[1] = U[1][i];
        points[i].v = U[2][i] / (U[0][i] + U[1][i]);

        // Выбираем модель: 1T, 2T, 3T
        switch (model)
        {
            case 1:
                computeT.compute(points[i], U3 / (U[0][i] + U[1][i]));
                break;
            case 2:
                computeT.compute(points[i], (U[4][i] + U[5][i]) / U[0][i],
                        U3 / (U[0][i] + U[1][i]));
                break;
            case 3:
                computeT.compute(points[i], U[4][i] / U[0][i],
                        U[5][i] / U[0][i], U3 / (U[0][i] + U[1][i]));
                break;
            default:
                computeT.compute(points[i], U[4][i] / U[0][i],
                        U[5][i] / U[0][i], U3 / (U[0][i] + U[1][i]));
        }

        // Возвращаем температуры и находим давление
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
    futureWatcher.setFuture(QtConcurrent::map(parN_v, calc));
    futureWatcher.waitForFinished();
}

void MixtureCo2Ar::updateDt()
{
    // Вспомогательная скорость
    double temp_t = 0.0;
    double min_t = 1.0;

    // Находим максимальную скорость распространения возмущения в потоке
    for (int i = 1; i < N_CELL - 1; ++i)
    {
        temp_t = dx[i] / (a[i] + points[i].v);
        if (temp_t < min_t)
        {
            min_t = temp_t;
        }
    }

    // Критерий Куранта-Фридрихса-Леви
    dt = CFL * min_t;
}

void MixtureCo2Ar::findShockPos()
{
    // Текущий и наибольший разрыв
    double maxGap = 0.0;
    double curGap = 0.0;

    // Находим положение наибольшего разрыва
    for (int i = 1; i < N_CELL - 1; ++i)
    {
        curGap = qAbs(points[i+1].p - points[i-1].p);
        if (curGap > maxGap)
        {
            maxGap = curGap;
            shockPos = i;
        }
    }
}

QVector<QVector<double>> MixtureCo2Ar::saveMacroParams()
{
    // Поиск длины свободного пробега
    double n0 = points[0].rho[0] / Mixture::mass(0);
    double n1 = points[0].rho[1] / Mixture::mass(1);
    double x0 = n0 / (n0 + n1);
    double d = x0 * Mixture::sigma(0) + (1.0 - x0) * Mixture::sigma(1);
    double L = K_BOLTZMANN * points[0].t / (qSqrt(2) * M_PI *
                                            qPow(d, 2.0) * points[0].p);

    // Подготовка таблицы
    QVector<QVector<double>> table;
    table.resize(12);

    // Заполняем таблицу
    for (int i = 0; i < N_CELL; ++i)
    {
        table[0].push_back(0.0);
        table[1].push_back(points[i].rho[0]);
        table[2].push_back(points[i].rho[1]);
        table[3].push_back(points[i].p);
        table[4].push_back(points[i].v);
        table[5].push_back(points[i].t);
        table[6].push_back(points[i].t12);
        table[7].push_back(points[i].t3);
        table[8].push_back(k[i]);
        table[9].push_back(a[i]);
        table[10].push_back(dx[i]);
        table[11].push_back(0.0);
    }
    for (int i = 1; i < N_CELL; ++i)
    {
        table[0][i] = table[0][i - 1] + 0.5 * (dx[i] + dx[i - 1]);
    }
    double ds = table[0][shockPos];
    for (int i = 0; i < N_CELL; ++i)
    {
        table[0][i] -= ds;
        table[11][i] = table[0][i] / L;
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

