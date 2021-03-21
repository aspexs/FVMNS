#include "nitrogensolver.h"


NitrogenSolver::NitrogenSolver(QObject *parent): AbstaractSolver(parent)
{

}

void NitrogenSolver::prepareSolving()
{

    U1.resize(solParam.NumCell+2);
    U2.resize(solParam.NumCell+2);
    U3.resize(solParam.NumCell+2);
    leftParam.density = leftParam.pressure /(UniversalGasConstant/molMass * leftParam.temp);
    leftParam.soundSpeed = sqrt(solParam.Gamma*leftParam.pressure/leftParam.density);
    leftParam.velocity = solParam.Ma*leftParam.soundSpeed;
    rightParam.density = ((solParam.Gamma + 1)* pow(solParam.Ma,2))/(2 + (solParam.Gamma -1)* pow(solParam.Ma,2))*leftParam.density;
    rightParam.pressure = (pow(solParam.Ma,2) * 2* solParam.Gamma - (solParam.Gamma - 1))/((solParam.Gamma +1))*leftParam.pressure;
    rightParam.temp = rightParam.pressure/(rightParam.density*UniversalGasConstant/molMass);
    auto rightParam2 = additionalSolver.BoundaryCondition[typeBC](leftParam, solParam);
    rightParam.velocity = leftParam.density*leftParam.velocity/rightParam.density;

    for(auto i  = 1; i < solParam.NumCell+1; i++)
    {
        if(i < solParam.NumCell/2+1)
        {
            U1[i] = leftParam.density;
            U2[i] = leftParam.density*leftParam.velocity;
            U3[i] = leftParam.pressure/(solParam.Gamma-1)+0.5*pow(leftParam.velocity,2)*leftParam.density;
        }

        else
        {
            U1[i] = rightParam.density;
            U2[i] = rightParam.density*rightParam.velocity;
            U3[i] = rightParam.pressure/(solParam.Gamma-1)+0.5*pow(rightParam.velocity,2)*rightParam.density;
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

    solParam.typeRightBorder = 0;
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

void NitrogenSolver::solveFlux(Matrix U1L, Matrix U2L, Matrix U3L, Matrix U1R, Matrix U2R, Matrix U3R)
{
    left_density   = U1L;
    right_density  = U1R;
    left_velocity  = U2L / left_density;
    right_velocity = U2R / right_density;
    left_pressure  = (U3L - Matrix::POW(left_velocity, 2)*0.5*left_density)*(solParam.Gamma-1);
    right_pressure = (U3R - Matrix::POW(right_velocity, 2)*0.5*right_density)*(solParam.Gamma-1);
    calcRiemanPStar();
    calcFliux();
}

void NitrogenSolver::calcRiemanPStar()
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

        rezultAfterPStart[i] = additionalSolver.ExacRiemanSolver(left,right,solParam.Gamma);
        return;
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calcPStar));
    futureWatcher.waitForFinished();
}
#include "qdebug.h"
void NitrogenSolver::calcFliux()
{
    QFutureWatcher<void> futureWatcher;
    const std::function<void(int&)> calcFlux = [this](int& i)
    {
        mutex.lock();
        auto point = rezultAfterPStart[i];
        double tempL = left_pressure[i]/(left_density[i]*UniversalGasConstant/molMass);
        double tempR = right_pressure[i]/(right_density[i]*UniversalGasConstant/molMass);
        auto du_dx = (right_velocity[i] - left_velocity[i])/delta_h;
        //auto du_dx = (right_velocity[i] - point.velocity)/delta_h*2;
        //auto du_dx = (point.velocity - left_velocity[i])/delta_h*2;
        mutex.unlock();
        double Tx = point.pressure/(point.density*UniversalGasConstant/molMass);
        double Pr = 2.0/3;
        double etta = additionalSolver.shareViscosity[typeShareVisc](leftParam.temp, Tx,0,0);
        double zetta = 0;//additionalSolver.bulkViscosity[typeBulkVisc](0,Tx,point.density, point.pressure);
        double G =  (4.0/3*etta + zetta)*du_dx;
        double k = solParam.Gamma*UniversalGasConstant/molMass*etta/(solParam.Gamma-1)/Pr;
        //qDebug() << G << Tx << du_dx;
        double dt_dx = (tempR - tempL)/delta_h;
        //double dt_dx = (tempR - Tx)/delta_h*2;
        //double dt_dx = (Tx - tempL)/delta_h*2;
        double q =   -k*dt_dx;
        F1[i] = (point.density * point.velocity);
        F2[i] = (F1.last()*point.velocity + point.pressure -G);
        F3[i] = ((point.pressure/(solParam.Gamma - 1) + point.density* pow(point.velocity,2)/2 + point.pressure) * point.velocity -G*point.velocity + q);
    };
    futureWatcher.setFuture(QtConcurrent::map(vectorForParallelSolving, calcFlux));
    futureWatcher.waitForFinished();
}

void NitrogenSolver::solve()
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
        auto pressure = (U3 - Matrix::POW(velosity,2)*0.5*U1)*(solParam.Gamma - 1);
        Matrix c = Matrix::SQRT(pressure*solParam.Gamma/U1);
        auto temp = velosity + c;
        auto max =*std::max_element(temp.begin(), temp.end());
        double dt = solParam.CFL*pow(delta_h,1)/max;
        //double dt = additionalSolver.TimeStepSolution[typeTimeStepSolution](velosity, U1, pressure, delta_h,solParam,U3);
        timeSolvind.push_back(dt);
        auto U1L = U1; U1L.removeLast();
        auto U2L = U2; U2L.removeLast();
        auto U3L = U3; U3L.removeLast();

        auto U1R = U1; U1R.removeFirst();
        auto U2R = U2; U2R.removeFirst();
        auto U3R = U3; U3R.removeFirst();
        solveFlux(U1L, U2L, U3L, U1R, U2R, U3R);

        auto res = additionalSolver.SolveEvolutionExplFirstOrder(F1, F2,F3,F4,U1, U2,U3,U4,dt,delta_h);
        U1 = res[0];
        U2 = res[1];
        U3 = res[2];
        solParam.typeRightBorder = 1;
        U1[0]=U1[1];   U1[U1.size() - 1] =rightParam.density;
        U2[0]=solParam.typeLeftBorder*U2[1];
        U2[U2.size()-1]= rightParam.density*rightParam.velocity;
        U3[0]=U3[1];   //U3[U3.size() - 1] =U3[U3.size() - 2];
        U3[U3.size() - 1] = rightParam.pressure/(solParam.Gamma-1)+0.5*pow(rightParam.velocity,2)*rightParam.density;
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

void NitrogenSolver::setTypePlot(int i)
{
    solParam.typePlot = i;
    QVector<double> values;
    switch (solParam.typePlot)
    {
        case 0: values = (U3 - Matrix::POW(U2 / U1 ,2) * 0.5 * U1) * (solParam.Gamma-1); break;
        case 1: values = U1; break;
        case 2: values = U2 / U1 ;break;
        case 3: values = ((U3 - Matrix::POW(U2 / U1 ,2) * 0.5 * U1) * (solParam.Gamma-1))/(U1 *UniversalGasConstant/molMass);break;
        default: values = U1; break;
    }
    emit updateGraph(x, values, solParam.lambda);
}
