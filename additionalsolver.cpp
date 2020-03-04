#include "additionalsolver.h"
#include "additionalsolverforco2.h"
#include <QProgressDialog>
#include <QFutureWatcher>
#include <QThread>
#include <QtConcurrent>

extern double sigma;
extern double epsilonDevK;
extern double molMass;
extern double mass;
AdditionalSolver::AdditionalSolver(QObject *parent) : QObject(parent)
{

}

void AdditionalSolver::solve()
{
    iterationVector.clear();
    rezultVector.clear();
    QVector<double> vec;
    for(int i = 0 ; i <= (stopValue - startValue)/step; i++)
    {
        iterationVector.push_back(startValue + i*step);
        vec.push_back(i);
    }
    rezultVector.resize(iterationVector.size());
    QProgressDialog *dialog = new QProgressDialog();
    dialog->setLabelText(QString("Подождите, идет расчет на %1 потоках...").arg(QThread::idealThreadCount()));
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    QFutureWatcher<void> futureWatcher;
    connect(&futureWatcher, &QFutureWatcher<void>::finished,			 dialog,		 &QProgressDialog::reset);
    connect(dialog,		    &QProgressDialog::canceled,	     		     &futureWatcher, &QFutureWatcher<void>::cancel);
    connect(&futureWatcher, &QFutureWatcher<void>::progressRangeChanged, dialog,		 &QProgressDialog::setRange);
    connect(&futureWatcher, &QFutureWatcher<void>::progressValueChanged, dialog,		 &QProgressDialog::setValue);

    const std::function<void(int)> calcCVibr = [this](int iteration)
    {
        rezultVector[iteration] = (func[typeSolve](startValue, startValue + iteration*step, density, pressure));
        return;
    };
    futureWatcher.setFuture(QtConcurrent::map(vec, calcCVibr));
    dialog->exec();

    futureWatcher.waitForFinished();

    emit completeSolution();
}

macroParam AdditionalSolver::ExacRiemanSolver(macroParam left, macroParam right, double Gamma)
{
    double maxIteration = 40; // макс число итераций
    double TOL=1e-8;
    double lambda = 0; // линия на грани КО
    macroParam ret;
    double left_soundspeed=sqrt ( Gamma*left.pressure/left.density );
    double right_soundspeed=sqrt( Gamma*right.pressure/right.density);

    double p_star= 0.5*(left.pressure+right.pressure) +
            0.125 * ( left.velocity-right.velocity ) *
            ( left.density+right.density ) *
            ( left_soundspeed+right_soundspeed );
    p_star=std::max(p_star,TOL);
    double pMin=std::min(left.pressure,right.pressure);
    double pMax=std::max(left.pressure,right.pressure);

    if ( p_star>pMax )
    {
        double temp1= sqrt ( ( 2.0/ ( Gamma+1.0 ) /left.density ) / ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *left.pressure ) );
        double temp2= sqrt ( ( 2.0/ ( Gamma+1.0 ) /right.density ) / ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *right.pressure ) );
        p_star= (temp1*left.pressure+temp2*right.pressure+ ( left.velocity-right.velocity ) ) / ( temp1+temp2 );
        p_star=std::max(p_star,TOL);
    }
    else if ( p_star<pMin )
    {
       double temp1= ( Gamma-1.0 ) / ( 2.0*Gamma );
       p_star= pow(( left_soundspeed+right_soundspeed+0.5*(Gamma-1.0 )*
                   ( left.velocity-right.velocity ) ) /
                   (left_soundspeed/pow(left.pressure,temp1) +
                   right_soundspeed/pow(right.pressure,temp1)), 1.0/temp1);
    }
    double f1 = 0, f2 = 0, f_d = 0 ;
    for(double iteration = 1;iteration < maxIteration; iteration++)
    {
        //LEFT
        double temp1 = sqrt ( ( 2.0/ ( Gamma+1.0 ) /left.density ) / ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *left.pressure ) );

        if (p_star<=left.pressure)
            f1=2.0/ ( Gamma-1.0 ) *left_soundspeed*
                    (pow(p_star/left.pressure,(Gamma-1.0 )/(2.0*Gamma))- 1.0) ;
        else
            f1= ( p_star-left.pressure ) *temp1;
        if (p_star<=left.pressure)
            f_d= pow( p_star/left.pressure,-(Gamma+1.0 )/( 2.0*Gamma ))/
                    ( left.density*left_soundspeed );
        else
            f_d=temp1* ( 1.0-0.5* ( p_star-left.pressure ) /
                         ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *left.pressure ) );
        //RIGHT
        temp1 = sqrt ( ( 2.0/ ( Gamma+1.0 ) /right.density ) /
                       ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *right.pressure ) );
        if (p_star<=right.pressure)
            f2=2.0/ ( Gamma-1.0 ) *right_soundspeed*
                    (pow(p_star/right.pressure,(Gamma-1.0 )/(2.0*Gamma))- 1.0) ;
        else
            f2= ( p_star-right.pressure ) *temp1;
        if (p_star<=right.pressure)
            f_d= f_d + pow( p_star/right.pressure,-(Gamma+1.0 )/( 2.0*Gamma ))/
                    ( right.density*right_soundspeed );
        else
            f_d=f_d + temp1* ( 1.0-0.5* ( p_star-right.pressure ) /
                         ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *right.pressure ) );
        double p_new = p_star - (f1+f2 - (left.velocity - right.velocity))/f_d;
        if(abs(p_new - p_star)/(0.5*abs(p_new + p_star)) < TOL)
            break;
        p_star = p_new;
    }
    // calculate star speed */
    double star_speed=0.5* ( left.velocity + right.velocity ) +0.5* ( f2-f1 );
    double left_star_density, left_tail_speed, left_head_speed,
            right_star_density, right_tail_speed,right_head_speed;
    //LEFT
    if ( p_star>=left.pressure ) {
            // SHOCK
        left_star_density = left.density * ( p_star / left.pressure + ( Gamma-1.0 ) / ( Gamma+1.0 ) ) /
                ( ( Gamma-1.0 ) / ( Gamma+1.0 ) * p_star / left.pressure + 1.0 );
        left_tail_speed = left.velocity -left_soundspeed * sqrt ( ( Gamma+1.0 ) / ( 2.0*Gamma ) * p_star/left.pressure +
                ( Gamma-1.0 ) / ( 2.0*Gamma ) );
        left_head_speed = left_tail_speed;
    }
    else // % left_wave_ == kRarefaction
    {
        left_star_density = left.density * pow(p_star/left.pressure,1.0/Gamma);
        left_head_speed = left.velocity - left_soundspeed;
        left_tail_speed = star_speed - sqrt ( Gamma*p_star/left_star_density );
    }
    //RIGHT
    if ( p_star>=right.pressure )
    {
        right_star_density = right.density *
                            ( p_star / right.pressure + ( Gamma-1.0 ) / ( Gamma+1.0 ) ) /
                            ( ( Gamma-1.0 ) / ( Gamma+1.0 ) * p_star / right.pressure + 1.0 );
        right_tail_speed = right.velocity +
           right_soundspeed * sqrt ( ( Gamma+1.0 ) / ( 2.0*Gamma ) * p_star/right.pressure +
           ( Gamma-1.0 ) / ( 2.0*Gamma ) );
        right_head_speed = right_tail_speed;
    }
    else // % right_wave_ == kRarefaction
    {
        right_star_density = right.density *  pow(p_star/right.pressure, 1.0/Gamma );
        right_head_speed = right.velocity + right_soundspeed;
        right_tail_speed = star_speed + sqrt ( Gamma*p_star/right_star_density );
    }

    bool is_left_of_contact = lambda  < star_speed;

    if ( is_left_of_contact )
    {// % the u is left of contact discontinuity
        if ( p_star>=left.pressure )  //the left wave is a shock
        {
            if ( lambda < left_head_speed )
            { // the u is before the shock
                ret.density  = left.density;
                ret.velocity = left.velocity;
                ret.pressure = left.pressure;
            }
            else  //% the u is behind the shock
            {
                ret.density  = left_star_density;
                ret.velocity = star_speed;
                ret.pressure = p_star;
            }
        }
        else // % the left wave is a rarefaction
        {
            if ( lambda < left_head_speed )//  % the u is before the rarefaction
            {
                ret.density  = left.density;
                ret.velocity = left.velocity;
                ret.pressure = left.pressure;
            }
            else
            {
                if ( lambda < left_tail_speed )//  % the u is inside the rarefaction
                {//% left_rarefaction (4.56)}
                    double temp1 = 2.0/ ( Gamma+1.0 ) + ( Gamma-1.0 ) / ( Gamma+1.0 )/left_soundspeed *(left.velocity - lambda);
                    ret.density = left.density *  pow(temp1, 2.0/( Gamma-1.0 ));
                    ret.pressure = left.pressure * pow(temp1, 2.0*Gamma/ ( Gamma-1.0));
                    ret.velocity = 2.0/ ( Gamma+1.0 ) * ( left_soundspeed + ( Gamma-1.0 ) /2.0*left.velocity + lambda);
                }
                else//  % the u is after the rarefaction
                {
                    ret.density  = left_star_density;
                    ret.velocity = star_speed;
                    ret.pressure = p_star;
                }
            }
        }
    }
    else// % the queried u is right of contact discontinuity
        //%------------------------------------------------------------------------
    {
        if ( p_star>=right.pressure )  //% the right wave is a shock
        {
            if ( lambda > right_head_speed )  //% the u is before the shock
            {
                ret.density  = right.density;
                ret.velocity = right.velocity;
                ret.pressure = right.pressure;
            }
            else  //% the u is behind the shock
            {
                ret.density  = right_star_density;
                ret.velocity = star_speed;
                ret.pressure = p_star;
            }
        }
        else // % the right wave is a rarefaction
        {
            if ( lambda > right_head_speed ) // % the u is before the rarefaction
            {
                ret.density  = right.density;
                ret.velocity = right.velocity;
                ret.pressure = right.pressure;
            }
            else
            {
                if ( lambda > right_tail_speed ) // % the u is inside the rarefaction
                {
                    double temp1 =2.0/ ( Gamma+1.0 ) - ( Gamma-1.0 ) / ( Gamma+1.0 ) /right_soundspeed *(right.velocity - lambda);
                    ret.density = right.density *  pow(temp1, 2.0/ ( Gamma-1.0 ) );
                    ret.pressure = right.pressure * pow(temp1, 2.0*Gamma/ ( Gamma-1.0 ) );
                    ret.velocity = 2.0/ ( Gamma+1.0 ) * ( -right_soundspeed + ( Gamma-1.0 ) /2.0*right.velocity + lambda);
                }
                else // % the u is after the rarefaction
                {
                    ret.density  = right_star_density;
                    ret.velocity = star_speed;
                    ret.pressure = p_star;
                }
            }
        }
    }
    return ret;
}

macroParam AdditionalSolver::ExacRiemanSolverCorrect(macroParam left, macroParam right, double GammaL, double GammaR)
{
    double maxIteration = 40; // макс число итераций
    double TOL=1e-8;
    double lambda = 0; // линия на грани КО
    macroParam ret;
    double left_soundspeed=sqrt ( GammaL*left.pressure/left.density );
    double right_soundspeed=sqrt( GammaR*right.pressure/right.density);

    double p_star= 0.5*(left.pressure+right.pressure) +
            0.125 * ( left.velocity-right.velocity ) *
            ( left.density+right.density ) *
            ( left_soundspeed+right_soundspeed );
    p_star=std::max(p_star,TOL);
    double pMin=std::min(left.pressure,right.pressure);
    double pMax=std::max(left.pressure,right.pressure);

    if ( p_star>pMax )
    {
        double temp1= sqrt ( ( 2.0/ ( GammaL+1.0 ) /left.density ) / ( p_star+ ( GammaL-1.0 ) / ( GammaL+1.0 ) *left.pressure ) );
        double temp2= sqrt ( ( 2.0/ ( GammaR+1.0 ) /right.density ) / ( p_star+ ( GammaR-1.0 ) / ( GammaR+1.0 ) *right.pressure ) );
        p_star= (temp1*left.pressure+temp2*right.pressure+ ( left.velocity-right.velocity ) ) / ( temp1+temp2 );
        p_star=std::max(p_star,TOL);
    }
    else if ( p_star<pMin )
    {
       double Gamma = 0.5*(GammaL + GammaR);
       double temp1= ( Gamma-1.0 ) / ( 2.0*Gamma );
       p_star= pow(( left_soundspeed+right_soundspeed+0.5*(Gamma-1.0 )*
                   ( left.velocity-right.velocity ) ) /
                   (left_soundspeed/pow(left.pressure,temp1) +
                   right_soundspeed/pow(right.pressure,temp1)), 1.0/temp1);
    }
    double f1 = 0, f2 = 0, f_d = 0 ;
    for(double iteration = 1;iteration < maxIteration; iteration++)
    {
        //LEFT
        double temp1 = sqrt ( ( 2.0/ ( GammaL+1.0 ) /left.density ) / ( p_star+ ( GammaL-1.0 ) / ( GammaL+1.0 ) *left.pressure ) );

        if (p_star<=left.pressure)
            f1=2.0/ ( GammaL-1.0 ) *left_soundspeed*
                    (pow(p_star/left.pressure,(GammaL-1.0 )/(2.0*GammaL))- 1.0) ;
        else
            f1= ( p_star-left.pressure ) *temp1;
        if (p_star<=left.pressure)
            f_d= pow( p_star/left.pressure,-(GammaL+1.0 )/( 2.0*GammaL ))/
                    ( left.density*left_soundspeed );
        else
            f_d=temp1* ( 1.0-0.5* ( p_star-left.pressure ) /
                         ( p_star+ ( GammaL-1.0 ) / ( GammaL+1.0 ) *left.pressure ) );
        //RIGHT
        temp1 = sqrt ( ( 2.0/ ( GammaR+1.0 ) /right.density ) /
                       ( p_star+ ( GammaR-1.0 ) / ( GammaR+1.0 ) *right.pressure ) );
        if (p_star<=right.pressure)
            f2=2.0/ ( GammaR-1.0 ) *right_soundspeed*
                    (pow(p_star/right.pressure,(GammaR-1.0 )/(2.0*GammaR))- 1.0) ;
        else
            f2= ( p_star-right.pressure ) *temp1;
        if (p_star<=right.pressure)
            f_d= f_d + pow( p_star/right.pressure,-(GammaR+1.0 )/( 2.0*GammaR ))/
                    ( right.density*right_soundspeed );
        else
            f_d=f_d + temp1* ( 1.0-0.5* ( p_star-right.pressure ) /
                         ( p_star+ ( GammaR-1.0 ) / ( GammaR+1.0 ) *right.pressure ) );
        double p_new = p_star - (f1+f2 - (left.velocity - right.velocity))/f_d;
        if(abs(p_new - p_star)/(0.5*abs(p_new + p_star)) < TOL)
            break;
        p_star = p_new;
    }
    // calculate star speed */
    double star_speed=0.5* ( left.velocity + right.velocity ) +0.5* ( f2-f1 );
    double left_star_density, left_tail_speed, left_head_speed,
            right_star_density, right_tail_speed,right_head_speed;
    //LEFT
    if ( p_star>=left.pressure ) {
            // SHOCK
        left_star_density = left.density * ( p_star / left.pressure + ( GammaL-1.0 ) / ( GammaL+1.0 ) ) /
                ( ( GammaL-1.0 ) / ( GammaL+1.0 ) * p_star / left.pressure + 1.0 );
        left_tail_speed = left.velocity -left_soundspeed * sqrt ( ( GammaL+1.0 ) / ( 2.0*GammaL ) * p_star/left.pressure +
                ( GammaL-1.0 ) / ( 2.0*GammaL ) );
        left_head_speed = left_tail_speed;
    }
    else // % left_wave_ == kRarefaction
    {
        left_star_density = left.density * pow(p_star/left.pressure,1.0/GammaL);
        left_head_speed = left.velocity - left_soundspeed;
        left_tail_speed = star_speed - sqrt ( GammaL*p_star/left_star_density );
    }
    //RIGHT
    if ( p_star>=right.pressure )
    {
        right_star_density = right.density *
                            ( p_star / right.pressure + ( GammaR-1.0 ) / ( GammaR+1.0 ) ) /
                            ( ( GammaR-1.0 ) / ( GammaR+1.0 ) * p_star / right.pressure + 1.0 );
        right_tail_speed = right.velocity +
           right_soundspeed * sqrt ( ( GammaR+1.0 ) / ( 2.0*GammaR ) * p_star/right.pressure +
           ( GammaR-1.0 ) / ( 2.0*GammaR ) );
        right_head_speed = right_tail_speed;
    }
    else // % right_wave_ == kRarefaction
    {
        right_star_density = right.density *  pow(p_star/right.pressure, 1.0/GammaR );
        right_head_speed = right.velocity + right_soundspeed;
        right_tail_speed = star_speed + sqrt ( GammaR*p_star/right_star_density );
    }

    bool is_left_of_contact = lambda  < star_speed;

    if ( is_left_of_contact )
    {// % the u is left of contact discontinuity
        if ( p_star>=left.pressure )  //the left wave is a shock
        {
            if ( lambda < left_head_speed )
            { // the u is before the shock
                ret.density  = left.density;
                ret.velocity = left.velocity;
                ret.pressure = left.pressure;
            }
            else  //% the u is behind the shock
            {
                ret.density  = left_star_density;
                ret.velocity = star_speed;
                ret.pressure = p_star;
            }
        }
        else // % the left wave is a rarefaction
        {
            if ( lambda < left_head_speed )//  % the u is before the rarefaction
            {
                ret.density  = left.density;
                ret.velocity = left.velocity;
                ret.pressure = left.pressure;
            }
            else
            {
                if ( lambda < left_tail_speed )//  % the u is inside the rarefaction
                {//% left_rarefaction (4.56)}
                    double temp1 = 2.0/ ( GammaL+1.0 ) + ( GammaL-1.0 ) / ( GammaL+1.0 )/left_soundspeed *(left.velocity - lambda);
                    ret.density = left.density *  pow(temp1, 2.0/( GammaL-1.0 ));
                    ret.pressure = left.pressure * pow(temp1, 2.0*GammaL/ ( GammaL-1.0));
                    ret.velocity = 2.0/ ( GammaL+1.0 ) * ( left_soundspeed + ( GammaL-1.0 ) /2.0*left.velocity + lambda);
                }
                else//  % the u is after the rarefaction
                {
                    ret.density  = left_star_density;
                    ret.velocity = star_speed;
                    ret.pressure = p_star;
                }
            }
        }
    }
    else// % the queried u is right of contact discontinuity
        //%------------------------------------------------------------------------
    {
        if ( p_star>=right.pressure )  //% the right wave is a shock
        {
            if ( lambda > right_head_speed )  //% the u is before the shock
            {
                ret.density  = right.density;
                ret.velocity = right.velocity;
                ret.pressure = right.pressure;
            }
            else  //% the u is behind the shock
            {
                ret.density  = right_star_density;
                ret.velocity = star_speed;
                ret.pressure = p_star;
            }
        }
        else // % the right wave is a rarefaction
        {
            if ( lambda > right_head_speed ) // % the u is before the rarefaction
            {
                ret.density  = right.density;
                ret.velocity = right.velocity;
                ret.pressure = right.pressure;
            }
            else
            {
                if ( lambda > right_tail_speed ) // % the u is inside the rarefaction
                {
                    double temp1 =2.0/ ( GammaR+1.0 ) - ( GammaR-1.0 ) / ( GammaR+1.0 ) /right_soundspeed *(right.velocity - lambda);
                    ret.density = right.density *  pow(temp1, 2.0/ ( GammaR-1.0 ) );
                    ret.pressure = right.pressure * pow(temp1, 2.0*GammaR/ ( GammaR-1.0 ) );
                    ret.velocity = 2.0/ ( GammaR+1.0 ) * ( -right_soundspeed + ( GammaR-1.0 ) /2.0*right.velocity + lambda);
                }
                else // % the u is after the rarefaction
                {
                    ret.density  = right_star_density;
                    ret.velocity = star_speed;
                    ret.pressure = p_star;
                }
            }
        }
    }
    return ret;
}

QVector<QVector<double> > AdditionalSolver::SolveEvolutionExplFirstOrder(Matrix F1, Matrix F2, Matrix F3,  Matrix U1old, Matrix U2old, Matrix U3old, double dt, double delta_h)
{
    int len = F1.size();
    Matrix rf;
    for(double i = 0; i < len; i++)
        rf.push_back(i);
    rf = rf * delta_h;
    Matrix U1new(len + 1) ,U2new (len + 1),U3new(len + 1), U4new(len + 1);
    auto temp_rfL = rf; temp_rfL.removeLast();
    auto temp_rfR = rf; temp_rfR.removeFirst();
    auto temp = Matrix::REVERSE(temp_rfR - temp_rfL) * dt;
    for(int i = 0 ; i < temp.size(); i++)
    {
        U1new[i+1] = U1old[i+1] - temp[i]*(F1[i+1] - F1[i]);
        U2new[i+1] = U2old[i+1] - temp[i]*(F2[i+1] - F2[i]);
        U3new[i+1] = U3old[i+1] - temp[i]*(F3[i+1] - F3[i]);
    }
    return {U1new, U2new, U3new};
}

QVector<QVector<double> > AdditionalSolver::SolveEvolutionExplFirstOrderForCO22(Matrix F1, Matrix F2, Matrix F3,Matrix F4,  Matrix U1old, Matrix U2old, Matrix U3old,Matrix U4old , double dt, double delta_h)
{
    int len = F1.size();
    Matrix rf;
    for(double i = 0; i < len; i++)
        rf.push_back(i);
    rf = rf * delta_h;
    Matrix U1new(len + 1) ,U2new (len + 1),U3new(len + 1), U4new(len + 1);
    auto temp_rfL = rf; temp_rfL.removeLast();
    auto temp_rfR = rf; temp_rfR.removeFirst();
    auto temp = Matrix::REVERSE(temp_rfR - temp_rfL) * dt;
    for(int i = 0 ; i < temp.size(); i++)
    {
        U1new[i+1] = U1old[i+1] - temp[i]*(F1[i+1] - F1[i]);
        U2new[i+1] = U2old[i+1] - temp[i]*(F2[i+1] - F2[i]);
        U3new[i+1] = U3old[i+1] - temp[i]*(F3[i+1] - F3[i]);
        U4new[i+1] = U4old[i+1] - temp[i]*(F4[i+1] - F4[i]);
    }
    return {U1new, U2new, U3new, U4new};
}

double AdditionalSolver::shareViscositySuperSimple(double startT, double currentT, double density, double pressure)
{
    double omega=0.81;
    double etha0=22.7e-6;
    return etha0*pow(currentT/startT,omega);
}

double AdditionalSolver::shareViscositySimple(double startT, double currentT, double density, double pressure)
{
    double temp1 = sqrt(kB*currentT/(M_PI*mass));
    double temp2 = 2*M_PI*pow(sigma,2);
    double omega2 = temp1*temp2;
    return (5*kB*currentT) /(8.*omega2);
}

double AdditionalSolver::shareViscosityOmega(double startT, double currentT, double density, double pressure)
{
    return (5*kB*currentT) /(8*getOmega22(currentT));
}

double AdditionalSolver::bulcViscositySimple(double startT, double currentT, double density, double pressure)
{
    double Nav = 6.02214129e23;
    double n = density/mass;
    double Crot = kB/mass;
    double Ctr = 2.0*Crot/3.0;
    double Cu = Crot + Ctr;
    double F = 1+ pow(M_PI,3/2)/2*pow(kB*currentT/2,-1/2) + (pow(M_PI,2)/4 +2)*pow(currentT/epsilonDevK,-1) + pow(M_PI,3/2)*pow(currentT/epsilonDevK,-3/2);
    double ZettaRot = ZettaInf/F;
    double tauR = (ZettaRot*M_PI*shareViscosityOmega(0, currentT)/(4*pressure));
    double Brot = (3*Crot)/(2*n*Cu*tauR);
    return (kB*currentT/Brot)* (Crot/pow(Cu,2));
}

double AdditionalSolver::bulcViscosityOld(double startT, double currentT, double density, double pressure)
{
    double zCO2Vibr = ZCO2Vibr(currentT);
    double cVibr = CVibr(currentT, zCO2Vibr);
    return (kB*currentT/Betta(currentT, cVibr))*pow((Crot() + cVibr)/(Crot() + Ctr() + cVibr),2);
}

double AdditionalSolver::bulcViscosityOld2(double CVibr, double currentT, double density, double pressure)
{
    return (kB*currentT/Betta(currentT, CVibr))*pow((Crot() + CVibr)/(Crot() + Ctr() + CVibr),2);
}

double AdditionalSolver::bulcViscosityNew(double startT, double currentT, double density, double pressure)
{
    double zCO2Vibr = ZCO2Vibr(currentT);
    double cVibr = 0;//CVibr(currentT, zCO2Vibr);
    return additionalSolverForCO2::bulcViscosity(currentT,cVibr, Ctr(), Crot(),pressure);
}

double AdditionalSolver::bulcViscosityFalse(double startT, double currentT, double density, double pressure)
{
    return 0;
}


double AdditionalSolver::vibrEnergy(double startT, double currentT, double density, double pressure)
{
    double P= 0;
    for(auto i1 = 0; i1 < 65 ; i1++)
        for(auto i2 = 0; i2 < 36; i2++)
            for(auto i3 = 0; i3 < 20; i3++)
            {
                P += ((i2+1)*(i1*e100+i2*e010+i3*e001)*exp(-((2*i1+i2)*e010+i3*e001)/(kB*currentT)));
            }
    return P/ZCO2Vibr(currentT)/mass;
}

double AdditionalSolver::fullEnergy(double startT, double currentT, double density, double pressure)
{
    return 5.0/2*kB*currentT/mass + vibrEnergy(0,currentT);
}

double AdditionalSolver::fullEnergy(double T, double Tv)
{
    return 5.0/2*kB*T/mass + vibrEnergy(0,Tv);
}

double AdditionalSolver::CVibrFunction(double startT, double currentT, double density, double pressure)
{
    return CVibr(currentT,ZCO2Vibr(currentT));
}

double AdditionalSolver::Crot()
{
    return kB/mass;
}

double AdditionalSolver::Ctr()
{
    return 3.0/2*Crot();
}

double AdditionalSolver::CVibr(double T, double ZCO2Vibr)
{
    double Q = 0, F = 0;
    for(auto i = 0; i < 65 ; i++)
        for(auto j = 0; j < 36; j++)
            for(auto l = 0; l < 20; l++)
            {
                Q += (j+1)*((i)*e100+(j)*e010+(l)*e001)*exp(-((2*i+j)*e010)/(kB*T))*exp(-l*e001/(kB*T));
                F += (j+1)*((i)*e100+(j)*e010+(l)*e001)*exp(-((2*i+j)*e010)/(kB*T))*exp(-l*e001/(kB*T))*ZCO2Vibr*((2*i+j)*e010+l*e001)/(kB*pow(T,2));
            }
    return (F-Q*DZDT(T))/(pow(ZCO2Vibr,2))/mass;
}

double AdditionalSolver::ZCO2Vibr(double T)
{
    double sum1 = 0, sum2 = 0, sum3 = 0;
       for(auto i1 = 0; i1 < 65; i1++)
       {
           sum1 += exp(-(i1*e100)/(kB*T));
           if(i1 < 36)
               sum2 += (i1+1)*exp(-(i1*e010)/(kB*T));
           if(i1 < 20)
               sum3 += exp(-(i1*e001)/(kB*T));
       }
       return sum1*sum2*sum3;
}

double AdditionalSolver::DZDT(double T)
{
    double Q = 0;
        for(auto a = 0; a < 65 ; a++)
            for(auto b = 0; b < 36; b++)
                for(auto c = 0; c < 20; c++)
                {
                    Q += (a*e100*exp(-(a*e100)/(kB*T)))/(kB*pow(T,2))*((b+1)*exp(-(b*e010)/(kB*T)))*(exp(-(c*e001)/(kB*T)))+
                            +(exp(-(a*e100)/(kB*T)))*(exp(-(c*e001)/(kB*T)))*((b+1)*b*e010*exp(-b*e010/(kB*T))/(kB*pow(T,2)))+
                            +(exp(-(a*e100)/(kB*T)))*((b+1)*exp(-(b*e010)/(kB*T)))*(c*e001*exp(-c*e001/(kB*T))/(kB*pow(T,2)));
                }
        return Q;
}

double AdditionalSolver::Betta(double T, double CVibr)
{
    double Cv = CVibr +  Crot() + Ctr();
    double F = 1+ pow(M_PI,3.0/2)/2*pow((T/epsilonDevK),-1.0/2) + (pow(M_PI,2)/4 +2)*pow(T/epsilonDevK,-1) + pow(M_PI,3.0/2)* pow(T/epsilonDevK,-3.0/2);
    double ZettaRot = ZettaInf/F;
    auto visc = shareViscosityOmega(0,T);
    double trot = ZettaRot*M_PI*visc/4;

    double a=-18.19;
    double b=40.47;
    double c=0;
    double d=0.00423;
    double tvibr = exp(a+b*pow(T,-1.0/3)+c*pow(T,-2.0/3)+d/(pow(T,-1.0/3)))*101325;
    double CIntDevTauInt = Crot()/trot + CVibr/tvibr;
    return (3.0*CIntDevTauInt)/(2.0 * Cv)* mass *UniversalGasConstant/molMass *T;
}

double AdditionalSolver::TauRot(double T, double P)
{
    double F = 1+ pow(M_PI,3.0/2)/2*pow((T/epsilonDevK),-1.0/2) + (pow(M_PI,2)/4 +2)*pow(T/epsilonDevK,-1) + pow(M_PI,3.0/2)* pow(T/epsilonDevK,-3.0/2);
    double ZettaRot = ZettaInf/F;
    return ZettaRot*M_PI*shareViscosityOmega(0,T)/(4*P);
}

double AdditionalSolver::TauVibr(double T, double P)
{
    P = P/101325;
    double a=-18.19;
    double b=40.47;
    double c=0;
    double d=0.00423;
    return exp(a+b*pow(T,-1.0/3)+c*pow(T,-2.0/3)+d/(pow(T,-1.0/3)))/P;
}

double AdditionalSolver::lambda(double T, double CVibr)
{
    double lambdaTr =(75.0*pow(kB,2)*T)/(32.0*mass*getOmega22(T));
    double lambdaVibr = (3.0*kB*T)/(8.0*getOmega11(T))*(Crot() + CVibr);
    return lambdaVibr + lambdaTr;
}

double AdditionalSolver::lambdaTr(double startT, double currentT, double density, double pressure )
{
    return (75.0*pow(kB,2)*currentT)/(32.0*mass*getOmega22(currentT));
}

double AdditionalSolver::lambdaTr_Rot(double T)
{
    return (75.0*pow(kB,2)*T)/(32.0*mass*getOmega22(T)) + (3.0*kB*T)/(8.0*getOmega11(T))*Crot()  ;
}

double AdditionalSolver::lambdaVibr2(double T, double Tv)
{
    double zCO2Vibr = ZCO2Vibr(Tv);
    double cVibr = CVibr(Tv, zCO2Vibr);
    return (3.0*kB*T)/(8.0*getOmega11(T))*cVibr;
}

double AdditionalSolver::lambdaVibr(double startT, double currentT, double density, double pressure )
{
    double zCO2Vibr = ZCO2Vibr(currentT);
    double cVibr = CVibr(currentT, zCO2Vibr);
    return (3.0*kB*currentT)/(8.0*getOmega11(currentT))*(Crot() + cVibr);
}

double AdditionalSolver::getOmega22(double T)
{
    QVector<double> f = {-0.40811, -0.05086, 0.34010, 0.70375, -0.10699, 0.00763};
    double a22= 1.5;
    double x = (log(T/epsilonDevK))+a22;

    double omegaLD = pow(f[0] + f[1]/pow(x,2) + f[2]/x + f[3]*x + f[4]*pow(x,2)+f[5]*pow(x,3),-1);
    double omegaS = sqrt(kB*T/(M_PI*mass))*2*M_PI*pow(sigma,2);
    return omegaLD*omegaS;
}

double AdditionalSolver::getOmega11(double T)
{
    QVector<double> f = {-0.16845, -0.02258, 0.19779, 0.64373, -0.09267, 0.00711};
    double a11 = 1.4;

    double x = (log(T/epsilonDevK))+a11;
    double omegaLD = pow(f[0] + f[1]/pow(x,2) + f[2]/x + f[3]*x + f[4]*pow(x,2)+f[5]*pow(x,3),-1);
    double omegaS =sqrt(kB*T/(M_PI*mass))*M_PI*pow(sigma,2);
    return omegaLD*omegaS;
}

macroParam AdditionalSolver::bondaryConditionRG(macroParam left, solverParams solParams)
{
    macroParam right;
    double zCO2Vibr = ZCO2Vibr(left.temp);
    double cVibr = CVibr(left.temp, zCO2Vibr);
    auto Cv = Crot() + Ctr() + cVibr;
    solParams.Gamma = (UniversalGasConstant/molMass + Cv)/Cv;
    right.density = ((solParams.Gamma + 1)* pow(solParams.Ma,2))/(2 + (solParams.Gamma -1)* pow(solParams.Ma,2))*left.density;
    right.pressure = (pow(solParams.Ma,2) * 2* solParams.Gamma - (solParams.Gamma - 1))/((solParams.Gamma +1))*left.pressure;
    right.temp = right.pressure/(right.density*UniversalGasConstant/molMass);
    return right;
}

macroParam AdditionalSolver::bondaryConditionRG2T(macroParam left, solverParams solParams)
{
    macroParam right;
    right.density = ((solParams.Gamma + 1)* pow(solParams.Ma,2))/(2 + (solParams.Gamma -1)* pow(solParams.Ma,2))*left.density;
    right.pressure = (pow(solParams.Ma,2) * 2* solParams.Gamma - (solParams.Gamma - 1))/((solParams.Gamma +1))*left.pressure;
    right.temp = right.pressure/(right.density*UniversalGasConstant/molMass);
    right.tempIntr = right.temp ;
    return right;
}

double AdditionalSolver::getTimeStep(Matrix velosity, Matrix density, Matrix pressure, double delta_h, solverParams solParams, Matrix gammas)
{
    Matrix c = Matrix::SQRT(pressure*solParams.Gamma/density);
    auto temp = velosity + c;
    auto max =*std::max_element(temp.begin(), temp.end());
    return solParams.CFL*pow(delta_h,2)/max*100;
}

double AdditionalSolver::getTimeStepFull(Matrix velosity, Matrix density, Matrix pressure, double delta_h, solverParams solParams, Matrix gammas)
{
    Matrix c = Matrix::SQRT(pressure*gammas/density);
    auto temp = velosity + c;
    auto max =*std::max_element(temp.begin(), temp.end());
    return solParams.CFL*pow(delta_h,2)/max*100;
}
