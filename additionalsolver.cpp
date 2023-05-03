#include "additionalsolver.h"
#include "additionalsolverforco2.h"
#include <QProgressDialog>
#include <QFutureWatcher>
#include <QThread>
#include <QtConcurrent>
#include <QProcess>
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

macroParam AdditionalSolver::ExacRiemanSolverCorrect(macroParam left, macroParam right, double GammaL, double GammaR, double lambda)
{
    double maxIteration = 40; // макс число итераций
    double TOL=1e-8;
    //double lambda = 0; // линия на грани КО
    macroParam ret;
    double left_soundspeed=sqrt ( GammaL*left.pressure/left.density );
    double right_soundspeed=sqrt( GammaR*right.pressure/right.density);
    double left_vacuum_front_speed = left.velocity + 2.0 * left_soundspeed / ( GammaL - 1.0 );
    double right_vacuum_front_speed = right.velocity - 2.0 * right_soundspeed / ( GammaR - 1.0 );
    double critical_speed =  left_vacuum_front_speed - right_vacuum_front_speed;

    if ( critical_speed < 0.0 ) //% образуется зона вукуума
    {
        double left_head_speed = left.velocity - left_soundspeed;
        double left_tail_speed = left_vacuum_front_speed;
        double right_head_speed = right.velocity + right_soundspeed;
        double right_tail_speed = right_vacuum_front_speed;
        //%-----------------------
        bool is_left_of_contact = lambda < left_tail_speed;
        if ( is_left_of_contact )// % определяем где находится искомая линия lambda слева ли от контактного разрыва
            if ( lambda < left_head_speed )
            {
                ret.density  = left.density;
                ret.velocity = left.velocity;
                ret.pressure = left.pressure;
            }
            else
            {
                //% left_rarefaction (4.56)
                double temp1 = 2.0/ ( GammaL+1.0 ) + ( GammaL-1.0 ) / ( GammaL+1.0 ) /left_soundspeed *(left.velocity - lambda);
                ret.density = left.density * pow(temp1,( 2.0/ ( GammaL-1.0 ) ));
                ret.pressure = left.pressure *  pow(temp1,( 2.0*GammaL/ ( GammaL-1.0 ) ));
                ret.velocity = 2.0/ ( GammaL+1.0 ) * ( left_soundspeed + ( GammaL-1.0 ) /2.0*left.velocity + lambda );
            }
        else
        {
            if ( lambda > right_tail_speed )
                if ( lambda > right_head_speed )
                {
                    ret.density  = right.density;
                    ret.velocity = right.velocity;
                    ret.pressure = right.pressure;
                }
                else
                    //%right_rarefaction (4.63)
                {
                    double temp1 =2.0/ ( GammaR+1.0 ) - ( GammaR-1.0 ) / ( GammaR+1.0 ) /right_soundspeed *(right.velocity - lambda);
                    ret.density = right.density * pow(temp1,( 2.0/ ( GammaR-1.0 ) ));
                    ret.pressure = right.pressure *  pow(temp1,(2.0*GammaR/ ( GammaR-1.0) ) );
                    ret.velocity = 2.0/ ( GammaR+1.0 ) * ( -right_soundspeed + ( GammaR-1.0 ) /2.0*right.velocity + lambda);
                }
            else
            {
                //% u resides inside vaccum
                ret.density=0.0;
                ret.velocity=0.0;
                ret.pressure = 0.0;
            }
        }
    }
    else
    {
    double p_star= 0.5*(left.pressure+right.pressure) +
            0.125 * ( left.velocity+right.velocity ) *
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
            }
    return ret;
}

QVector<QVector<double> > AdditionalSolver::SolveEvolutionExplFirstOrder(Matrix F1, Matrix F2, Matrix F3, Matrix F4,  Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old, double dt, double delta_h)
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

QVector<QVector<double> > AdditionalSolver::SolveEvolutionExplFirstOrderForCO22(Matrix F1, Matrix F2, Matrix F3,Matrix F4,
                                                                                Matrix U1old, Matrix U2old, Matrix U3old,Matrix U4old ,
                                                                                double dt, double delta_h, Matrix R)
{
    int len = F1.size();
    Matrix U1new(len + 1) ,U2new (len + 1),U3new(len + 1), U4new(len + 1);
    const auto temp = dt/delta_h;
    double sumDeltaF1 = 0, sumDeltaF2 = 0, sumDeltaF3 = 0, sumDeltaF4 = 0;
    for(int i = 0 ; i < len-1; i++)
    {
        const auto deltaF1 = temp*(F1[i+1]  - F1[i]);
        const auto deltaF2 = temp*(F2[i+1]  - F2[i]);
        const auto deltaF3 = temp*(F3[i+1]  - F3[i]);
        const auto deltaF4 = temp*(F4[i+1]  - F4[i]);
        sumDeltaF1 += deltaF1;
        sumDeltaF2 += deltaF2;
        sumDeltaF3 += deltaF3;
        sumDeltaF4 += deltaF4;
        U1new[i+1] = U1old[i+1] - deltaF1;
        U2new[i+1] = U2old[i+1] - deltaF2;
        U3new[i+1] = U3old[i+1] - deltaF3;
        U4new[i+1] = U4old[i+1] - deltaF4;// -  delta_h*dt*(R[i+1] - R[i]);
    }
    //U1new[0]= U1new[1];
    //U2new[0]= U2new[1];
    //U3new[0]= U3new[1];
    //U4new[0]= U4new[1];
    //U1new[len-1] = U1new[len-2];
    //U2new[len-1] = U2new[len-2];
    //U3new[len-1] = U3new[len-2];
    //U4new[len-2] = U4new[len-3];
    return {U1new, U2new, U3new, U4new};
}

QVector<QVector<double> > AdditionalSolver::SolveEvolutionExplFirstOrderForO2(Matrix F1, Matrix F2, Matrix F3,Matrix F4, Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old, double dt, double delta_h)
{
    int len = F1.size();
    Matrix U1new(len + 1) ,U2new (len + 1),U3new(len + 1), U4new(len + 1);
    auto temp = dt/delta_h;
    for(int i = 0 ; i < len-1; i++)
    {
        U1new[i+1] = U1old[i+1] - temp*(F1[i+1] - F1[i]);
        U2new[i+1] = U2old[i+1] - temp*(F2[i+1] - F2[i]);
        U3new[i+1] = U3old[i+1] - temp*(F3[i+1] - F3[i]);
        U4new[i+1] = U4old[i+1] - temp*(F4[i+1] - F4[i]);
    }

    //U1new[0] =U1new[1];
    //U2new[0] =U2new[1];
    //U3new[0] =U3new[1];
    //U4new[0] =U4new[1];
    return {U1new, U2new, U3new, U4new};
}

QVector<QVector<double> > AdditionalSolver::SEEFOForCO2(Matrix F1, Matrix F2, Matrix F3, Matrix F4, Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old, double dt, double delta_h,
                                                        Matrix F11, Matrix F22, Matrix F33, Matrix F44, Matrix R)
{
    int len = F1.size();
    Matrix U1new(len + 1) ,U2new (len + 1),U3new(len + 1), U4new(len + 1), U1(len + 1),U2(len + 1),U3(len + 1),U4(len + 1),Rnew(len + 1);
    const auto temp = dt/delta_h;
    const auto temp2 = pow(temp,2);
    const auto temp3 = dt/delta_h/delta_h;
    double sumDeltaF1 = 0, sumDeltaF2 = 0, sumDeltaF3 = 0, sumDeltaF4 = 0;
    for(int i = 1 ; i < len; i++)
    {
//        auto Q1 = U1old[i] - temp*(F1[i] - F1[i-1]);
//        auto Q1 = U1old[i] - temp*(F1[i] - F1[i-1]);
//        auto Q1 = U1old[i] - temp*(F1[i] - F1[i-1]);
//        auto Q1 = U1old[i] - temp*(F1[i] - F1[i-1]);
//        //тут вычисляем
//        auto Q = 0.5*(U1old[i] + Q1 + temp*())
        const auto deltaF1 = temp*(F1[i]  - F1[i-1]);
        const auto deltaF2 = temp*(F2[i]  - F2[i-1]);
        const auto deltaF3 = temp*(F3[i]  - F3[i-1]);
        const auto deltaF4 = temp*(F4[i]  - F4[i-1]);
        sumDeltaF1 += deltaF1;
        sumDeltaF2 += deltaF2;
        sumDeltaF3 += deltaF3;
        sumDeltaF4 += deltaF4;

        U1new[i] =  U1old[i] - deltaF1;
        U2new[i] =  U2old[i] - deltaF2;
        U3new[i] =  U3old[i] - deltaF3;
        U4new[i] =  U4old[i] - deltaF4 + dt*R[i];

    }
    return {U1new, U2new, U3new, U4new, Rnew, {sumDeltaF1+sumDeltaF2+sumDeltaF3+sumDeltaF4}};
}

double Sign1(double x)
{
    int cond  = (x>=0);
    return  cond*1 + (1-cond)*(-1);
}
double SuperBee(double r, double omega)
{
    //int cond1 = (r>0);
    //int cond2 = (r <= 0.5);
    //int cond3 = (r>1);
    //return (2*r*cond2 + (1 - cond2)*(cond3*qMin(qMin(2.0,r), 2.0/(1-omega+(1+omega)*r + 1e-12)) + (1-cond3)*1 ))*cond1;

    auto cond1 = (r>0);
    return ( qMin(2*r/(1+r + 1e-12), 2/(1-omega+(1+omega)*r + 1e-12) ) )*cond1;
}
double MUSCLlimiter(double r, double delta, double omega, int limType )
{
    auto phi = SuperBee(r,omega);
   return phi*delta;
}
double getEnergyVibrTemp(QVector<double>EnergyVibr, double energy, double energyStepTemp, double energyStartTemp)
{
    for(auto i = 0 ; i < EnergyVibr.size(); i++)
        if( energy < EnergyVibr[i])
            return i  * energyStepTemp + energyStartTemp;
    return (EnergyVibr.size()-1) * energyStepTemp + energyStartTemp;
}
macroParam point(double U1old, double U2old, double U3old, double U4old,QVector<double>EnergyVibr, double energyStepTemp, double energyStartTemp)
{
    macroParam p;
    p.density = U1old;
    p.velocity = U2old/U1old;
    auto eVibr = U4old/U1old;
    p.tempIntr = getEnergyVibrTemp(EnergyVibr,eVibr, energyStepTemp, energyStartTemp);
    auto energy_TrRot = U3old/U1old -  p.velocity*p.velocity/2 - eVibr;
    p.temp = energy_TrRot*2*mass/(5*kB);
    p.pressure = U1old* p.temp*UniversalGasConstant/molMass;
    return p;
}
QVector<Matrix> AdditionalSolver::SolveMUSCL_UL_UR(Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old, double Gamma, double dt_dx,QVector<double>EnergyVibr, double energyStepTemp, double energyStartTemp,int LimType, double omega)
{
    const auto TOL = 1e-12;
    const auto len = U1old.size();

    Matrix U1new(len + 1) ,U2new (len + 1),U3new(len + 1), U4new(len + 1);

    for(auto i = 1; i < len; i ++ )
    {
        U1new[i] = U1old[i] - U1old[i-1];
        U2new[i] = U2old[i] - U2old[i-1];
        U3new[i] = U3old[i] - U3old[i-1];
        U4new[i] = U4old[i] - U4old[i-1];
    }
    U1new[0] = U1new[1];    U1new[len] = U1new[len-1];
    U2new[0] = U2new[1];    U2new[len] = U2new[len-1];
    U3new[0] = U3new[1];    U3new[len] = U3new[len-1];
    U4new[0] = U4new[1];    U4new[len] = U4new[len-1];

    for(auto i = 0; i < U1new.size(); i ++ )
    {
        int condition = U1new[i] < TOL;
        U1new[i] = TOL*Sign1( U1new[i])*condition + (1 - condition)* U1new[i];

        condition = U2new[i] < TOL;
        U2new[i] = TOL*Sign1( U2new[i])*condition + (1 - condition)* U2new[i];

        condition = U3new[i] < TOL;
        U3new[i] = TOL*Sign1( U3new[i])*condition + (1 - condition)* U3new[i];

        condition = U4new[i] < TOL;
        U4new[i] = TOL*Sign1( U4new[i])*condition + (1 - condition)* U4new[i];
    }

    Matrix delta_i1, delta_i2, delta_i3, delta_i4, r1,r2,r3,r4;
    for(auto i = 1; i < U1new.size(); i ++ )
    {
        delta_i1.push_back(0.5*(1.0+omega)*U1new[i-1] + 0.5*(1.0-omega)*U1new[i]); r1.push_back(U1new[i-1]/U1new[i]);
         delta_i1[i-1] =  MUSCLlimiter(r1[i-1],delta_i1[i-1],omega,LimType);
        delta_i2.push_back(0.5*(1.0+omega)*U2new[i-1] + 0.5*(1.0-omega)*U2new[i]); r2.push_back(U2new[i-1]/U2new[i]);
         delta_i2[i-1] =  MUSCLlimiter(r2[i-1],delta_i2[i-1],omega,LimType);
        delta_i3.push_back(0.5*(1.0+omega)*U3new[i-1] + 0.5*(1.0-omega)*U3new[i]); r3.push_back(U3new[i-1]/U3new[i]);
         delta_i3[i-1] =  MUSCLlimiter(r3[i-1],delta_i3[i-1],omega,LimType);
        delta_i4.push_back(0.5*(1.0+omega)*U4new[i-1] + 0.5*(1.0-omega)*U4new[i]); r4.push_back(U4new[i-1]/U4new[i]);
         delta_i4[i-1] =  MUSCLlimiter(r4[i-1],delta_i4[i-1],omega,LimType);
    }
    Matrix U1BLeft = U1old - delta_i1*0.5;
    Matrix U2BLeft = U2old - delta_i2*0.5;
    Matrix U3BLeft = U3old - delta_i3*0.5;
    Matrix U4BLeft = U4old - delta_i4*0.5;

    Matrix U1BRight = U1old + delta_i1*0.5;
    Matrix U2BRight = U2old + delta_i2*0.5;
    Matrix U3BRight = U3old + delta_i3*0.5;
    Matrix U4BRight = U4old + delta_i4*0.5;

    Matrix uBLeft   = U2BLeft/U1BLeft;  //speeed
    Matrix uBRight = U2BRight/U1BRight; //speeed

    Matrix dF2, dF3,dF4;
    for(auto i = 0; i < U1new.size(); i ++ )
    {
        auto pointL = point(U1BLeft[i], U2BLeft[i],U3BLeft[i],U4BLeft[i],EnergyVibr,energyStepTemp,energyStartTemp);
        auto pointR = point(U1BRight[i], U2BRight[i],U3BRight[i],U4BRight[i],EnergyVibr,energyStepTemp,energyStartTemp);
        auto df2_i = (pointL.density * pointL.velocity*pointL.velocity + pointL.pressure) - (pointR.density * pointR.velocity*pointR.velocity + pointR.pressure);
        dF2.push_back(df2_i);
        auto HL =  pointL.pressure/pointL.density + pointL.velocity* pointL.velocity/2 + 5.0/2*kB*pointL.temp/mass + AdditionalSolver::vibrEnergy(0,pointL.tempIntr);
        auto HR =  pointR.pressure/pointR.density + pointR.velocity* pointR.velocity/2 + 5.0/2*kB*pointR.temp/mass + AdditionalSolver::vibrEnergy(0,pointR.tempIntr);
        auto df3_i = (pointL.density * pointL.velocity*(HL)) - (pointR.density * pointR.velocity*(HR));
        dF3.push_back(df3_i);
        auto df4_i = (pointL.density * pointL.velocity*(AdditionalSolver::vibrEnergy(0,pointL.tempIntr))) - (pointR.density * pointR.velocity*(AdditionalSolver::vibrEnergy(0,pointR.tempIntr)));
        dF4.push_back(df4_i);
    }

    Matrix dF1 = U2BLeft - U2BRight;

    U1BLeft = U1BLeft + dF1*0.5*dt_dx;
    U2BLeft = U2BLeft + dF2*0.5*dt_dx;
    U3BLeft = U3BLeft + dF3*0.5*dt_dx;
    U4BLeft = U4BLeft + dF4*0.5*dt_dx;
    U1BLeft.removeLast()        ;
    U2BLeft.removeLast();
    U3BLeft.removeLast();
    U4BLeft.removeLast();

    U1BRight = U1BRight + dF1*0.5*dt_dx;
    U2BRight = U2BRight + dF2*0.5*dt_dx;
    U3BRight = U3BRight + dF3*0.5*dt_dx;
    U4BRight = U4BRight + dF3*0.5*dt_dx;
    U1BRight.removeFirst()        ;
    U2BRight.removeFirst();
    U3BRight.removeFirst();
    U4BRight.removeFirst();

    return {U1BLeft,U2BLeft,U3BLeft,U4BLeft,U1BRight,U2BRight,U3BRight,U4BRight};

}

QVector<QVector<double> > AdditionalSolver::SEEFOForCO23T(Matrix F1, Matrix F2, Matrix F3, Matrix F4, Matrix F5,  Matrix F6, Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old, Matrix U5old, Matrix U6old, double dt, double delta_h, Matrix R_1, Matrix R_2)
{
    int len = F1.size();
    Matrix U1new(len + 1) ,U2new (len + 1),U3new(len + 1), U4new(len + 1), U5new(len + 1),  U6new(len + 1);
    const auto temp = dt/delta_h;
    double sumDeltaF1 = 0, sumDeltaF2 = 0, sumDeltaF3 = 0, sumDeltaF4 = 0, sumDeltaF5 = 0,  sumDeltaF6 = 0;
    for(int i = 1 ; i < len; i++)
    {
        const auto deltaF1 = temp*(F1[i]  - F1[i-1]);
        const auto deltaF2 = temp*(F2[i]  - F2[i-1]);
        const auto deltaF3 = temp*(F3[i]  - F3[i-1]);
        const auto deltaF4 = temp*(F4[i]  - F4[i-1]);
        const auto deltaF5 = temp*(F5[i]  - F5[i-1]);
        const auto deltaF6 = temp*(F6[i]  - F6[i-1]);
        sumDeltaF1 += deltaF1;
        sumDeltaF2 += deltaF2;
        sumDeltaF3 += deltaF3;
        sumDeltaF4 += deltaF4;
        sumDeltaF5 += deltaF5;
        sumDeltaF6 += deltaF6;

        U1new[i] =  U1old[i] - deltaF1;
        U2new[i] =  U2old[i] - deltaF2;
        U3new[i] =  U3old[i] - deltaF3;
        U4new[i] =  U4old[i] - deltaF4;
        U5new[i] =  U5old[i] - deltaF5 + dt*R_1[i];
        U6new[i] =  U6old[i] - deltaF6 + dt*R_2[i];

    }
    return {U1new, U2new, U3new, U4new, U5new, U6new,
        {sumDeltaF1 + sumDeltaF2 + sumDeltaF3 + sumDeltaF4 + sumDeltaF5 + sumDeltaF6}};
}


void AdditionalSolver::MackCormack(Matrix U1, Matrix U2, Matrix U3, Matrix U4,
                                   Matrix& F1, Matrix& F2, Matrix& F3, Matrix& F4,
                                   Matrix& F11, Matrix& F22, Matrix& F33, Matrix& F44,Matrix& R,
                                   QList<double>& EnergyVibr)
{
    Matrix velosity = U2/U1;
    auto EVibr = U4/U1;
    auto energyFull = U3/U1 - Matrix::POW(velosity,2)/2;
    auto E_tr_rot = energyFull - EVibr;
    Matrix T, Tvv;
    Matrix pressures;

    for (auto energy: EVibr)
    {
        auto TempTv = getEnergyVibrTemp(energy,EnergyVibr);
        if(TempTv < 300)
            TempTv = 300;
        Tvv.push_back(TempTv);
    }
    for (auto energy: E_tr_rot)
        T.push_back(energy*2*mass/(5*kB));
    F1.clear();
    F1.resize(T.size());
    F2.clear();
    F2.resize(T.size());
    F3.clear();
    F3.resize(T.size());
    F4.clear();
    F4.resize(T.size());
    F11.clear();
    F11.resize(T.size());
    F22.clear();
    F22.resize(T.size());
    F33.clear();
    F33.resize(T.size());
    F44.clear();
    F44.resize(T.size());
    R.clear();
    R.resize(T.size());
    pressures = U1*T*UniversalGasConstant/molMass;
    for(auto i = 0; i < T.size(); i++)
    {
        auto Tx = T[i];
        auto vel = velosity[i];
        auto Tv = Tvv[i];
        auto rho = U1[i];
        auto press = pressures[i];
        double omega11 = AdditionalSolver::getOmega11(Tx);
        double omega22 = AdditionalSolver::getOmega22(Tx);
        double zCO2Vibr = AdditionalSolver::ZCO2Vibr(Tv);
        double cVibr = AdditionalSolver::CVibr(Tv, zCO2Vibr);
        double etta = (5*kB*Tx) /(8*omega22);
        double Cv =5/2*kB/mass;
        double F = 1+ pow(M_PI,3.0/2)/2*pow((Tx/epsilonDevK),-1.0/2) + (pow(M_PI,2)/4 +2)*pow(Tx/epsilonDevK,-1) + pow(M_PI,3.0/2)* pow(Tx/epsilonDevK,-3.0/2);
        double ZettaRot = ZettaInf/F;
        double trot = ZettaRot*M_PI*etta/4;
        double CIntDevTauInt = kB/mass/trot;
        double betta = (3.0*CIntDevTauInt)/(2.0 * Cv)* mass *UniversalGasConstant/molMass *Tx;

        double Evibr = AdditionalSolver::vibrEnergy(0,Tx);
        double Evibr2 = AdditionalSolver::vibrEnergy(0,Tv);
        double Etr_rot = 5.0/2*kB*Tx/mass;
        double zetta = (kB*Tx/betta)*0.16;//additionalSolver.bulkViscosity[AdditionalSolver::BulkViscosity::BV_ONLY_RT_ROT](Cvibr,Tx,point.density, point.pressure);
        double P =  (4.0/3*etta + zetta);
        double qVibr =-(3.0*kB*Tx)/(8.0*omega11)*cVibr;// -additionalSolver.lambdaVibr2(Tx,Tv) * dtv_dx;
        double qTr =-((75.0*pow(kB,2)*Tx)/(32.0*mass*omega22) + (3.0*kB*Tx)/(8.0*omega11)*kB/mass); //-additionalSolver.lambdaTr_Rot(Tx)*dt_dx;
        double tauVibr = TauVibr(Tx, press);
        double entalpi = Etr_rot + Evibr2 + press/rho + pow(vel,2)/2;
        F1[i] = (rho * vel);
        F11[i] = 0;
        F2[i] = (F1[i]*vel + press);
        F22[i] = -vel*P;
        F3[i] = (F1[i]*entalpi);
        F33[i] = Tv*qVibr + Tx*qTr - P*vel * vel;
        auto deltaE = Evibr - Evibr2;
        R[i] = -rho/tauVibr * deltaE;
        F4[i] =  F1[i]*Evibr2;
        F44[i] = Tv*qVibr;
    }
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

double AdditionalSolver::bulcViscosityOld2(double CVibr, double T, double density, double pressure)
{
    //p*\tau_rot
    double F = 1+ pow(M_PI,3.0/2)/2*pow((T/epsilonDevK),-1.0/2) + (pow(M_PI,2)/4 +2)*pow(T/epsilonDevK,-1) + pow(M_PI,3.0/2)* pow(T/epsilonDevK,-3.0/2);
    double ZettaRot = ZettaInf/F;
    double visc = shareViscosityOmega(0,T);
    double t_rot = ZettaRot*M_PI*visc/4;

    //p*\tau_vibr
    double a=-18.19;
    double b=40.47;
    double c=0;
    double d=0.00423;
    double t_vibr = exp(a+b*pow(T,-1.0/3)+c*pow(T,-2.0/3)+d/(pow(T,-1.0/3)))*101325;

    double C_int = CVibr +  Crot();
    double C_v = C_int + Ctr();

    return gasConst*pow(C_int/C_v ,2) / (Crot()/t_rot + CVibr/t_vibr);
}

double AdditionalSolver::bulcViscosityNew(double startT, double currentT, double density, double pressure)
{
    //double zCO2Vibr = ZCO2Vibr(currentT);
    double cVibr = 0;//CVibr(currentT, zCO2Vibr);
    return additionalSolverForCO2::bulcViscosity(currentT,cVibr, Ctr(), Crot(),pressure);
}

double AdditionalSolver::bulcViscosityOnlyTRRot(double startT, double currentT, double density, double pressure)
{
    return (kB*currentT/Betta(currentT, 0))*pow((Crot())/(Crot() + Ctr() ),2);
}

double AdditionalSolver::bulcViscosityFalse(double startT, double currentT, double density, double pressure)
{
    return 0;
}

double E_CO2(int v1, int v2, int v3)
{
    double d[3] = { 1., 2., 1. };
   //double ECO2 = hc * (o2[0] * (v1 + d[0] / 2.0) + o2[1] * (v2 + d[1] / 2.0) + o2[2] * (v3 + d[2] / 2.0) + ox[0] * (v1 + d[0] / 2.0) * (v1 + d[0] / 2.0) + ox[1] * (v1 + d[0] / 2.0) * (v2 + d[1] / 2.0) + ox[2] * (v1 + d[0] / 2.0) * (v3 + d[2] / 2.0) + ox[3] * (v2 + d[1] / 2.0) * (v2 + d[1] / 2.0) + ox[4] * (v2 + d[1] / 2.0) * (v3 + d[2] / 2.0) + ox[5] * (v3 + d[2] / 2.0) * (v3 + d[2] / 2.0));
    double ECO2 = hc * (o2[0] * (v1 + d[0] / 2.0) + o2[1] * (v2 + d[1] / 2.0) + o2[2] * (v3 + d[2] / 2.0)); //harmonic
    //double ECO2 = v1 * e100 + v2 * e010 + v3 * e001;
    //return En_CO2_0[v1][v2]/1.e20;
    //qDebug() << v1 << v2 << v3;
    return ECO2;
}

double AdditionalSolver::vibrEnergy(double startT, double currentT, double density, double pressure)
{
    int j1, j2, j3;
        double D = 0;
        for (j1 = 0; j1 < 65; j1++)
        {
            for (j2 = 0; j2 < 36; j2++)
            {
                for (j3 = 0; j3 < 20; j3++)
                {
                    if (E_CO2(j1, j2, j3) < De)
                    {

                        D += (j2 + 1.) * (j1 * e100 + j2 * e010 + j3 * e001) * exp(-(j1 * e100 + j2 * e010 + j3 * e001) / kB / currentT);
                    }
                }
            }
        }
    //double P= 0;
    //for(auto i1 = 0; i1 < 65 ; i1++)
    //    for(auto i2 = 0; i2 < 36; i2++)
    //        for(auto i3 = 0; i3 < 20; i3++)
    //        {
    //            P += ((i2+1)*(i1*e100+i2*e010+i3*e001)*exp(-((2*i1+i2)*e010+i3*e001)/(kB*currentT)));
    //
    return D/ZCO2Vibr(currentT)/mass;
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

double AdditionalSolver::Crot(double startT, double currentT, double density, double pressure )
{
    return kB/mass;
}

double AdditionalSolver::Ctr(double startT, double currentT, double density , double pressure)
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
    static const double De = 8.83859e-19;
    int i1, i2, i3;
    double S = 0;
    for (i1 = 0; i1 < 34; i1++)
    {
        for (i2 = 0; i2 < 67; i2++)
        {
            for (i3 = 0; i3 < 20; i3++)
            {
                if ((i1 * e100 + i2 * e010 + i3 * e001) < De)
                {
                    S += (i2 + 1.) * exp(-i1 * e100 / kB / T) * exp(-i2 * e010 / kB / T) * exp(-i3 * e001 / kB / T);
                }
            }
        }
    }
    return S;
    //double sum1 = 0, sum2 = 0, sum3 = 0;
    //   for(auto i1 = 0; i1 < 65; i1++)
    //   {
    //       sum1 += exp(-(i1*e100)/(kB*T));
    //       if(i1 < 36)
    //           sum2 += (i1+1)*exp(-(i1*e010)/(kB*T));
    //       if(i1 < 20)
    //           sum3 += exp(-(i1*e001)/(kB*T));
    //   }
    //   return sum1*sum2*sum3;
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
    return (3.0* kB*T* CIntDevTauInt)/(2.0 * Cv);
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

double AdditionalSolver::lambdaForNitrogen(double gamma,double T,double density, double pressure)
{
    double Pr = 2.0/3;
    auto etta = shareViscosityOmega(0, T);
    return gamma*UniversalGasConstant/molMass*etta/(gamma-1)/Pr;
}
double AdditionalSolver::lambda(double T, double CVibr)
{
    double lambdaTr_rot = lambdaTr_Rot(T);
    double lambdaVibr = (3.0*kB*T)/(8.0*getOmega11(T))*CVibr;
    return lambdaVibr + lambdaTr_rot;
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
    double zCO2Vibr = ZCO2Vibr(left.temp);
    double cVibr = CVibr(left.temp, zCO2Vibr);
    auto Cv = Crot() + Ctr() + cVibr;
    solParams.Gamma = (UniversalGasConstant/molMass + Cv)/Cv;
    right.density = ((solParams.Gamma + 1)* pow(solParams.Ma,2))/(2 + (solParams.Gamma -1)* pow(solParams.Ma,2))*left.density;
    right.pressure = (pow(solParams.Ma,2) * 2* solParams.Gamma - (solParams.Gamma - 1))/((solParams.Gamma +1))*left.pressure;
    right.temp = right.pressure/(right.density*UniversalGasConstant/molMass);
    right.tempIntr = right.temp ;
    right.velocity = left.density*left.velocity/right.density;
    return right;
}

QStringList AdditionalSolver::runPython(macroParam left,int mlt )
{

    QStringList arguments { "Fun.py",
                            left.gas, QString::number(left.velocity), QString::number(left.density),
                            QString::number(left.temp), "Text.txt",
                            QString::number(mlt), QString::number(left.tempIntr),QString::number(left.tempIntr)};
    QProcess p;
    p.start("python", arguments);
    p.waitForFinished();
    qDebug() << p.readAllStandardError();
    return QString(p.readAll()).replace("[","").replace("]","").replace("\r\n","").split(" ");

}

double AdditionalSolver::getEnergyVibrTemp(double energy, QList<double> &e)
{
    for(auto i = 0 ; i < e.size(); i++)
        if( energy < e[i])
            return i  * 1 + 100-1;
    return (e.size()-1) * 1 + 100 - 1;
}

double AdditionalSolver::EVibr12(double sT, double T, double r, double t)
{
    int L1 = 34;
    int L2 = 67;
    int L3 = 20;
    double E = 0, D = 0;

    for (int i1 = 0; i1 <= L1; i1++)
    {
        for (int i2 = 0; i2 <= L2; i2++)
        {
            double e = (i2 + 1.) * (i1 * e100 + i2 * e010) * exp(-(i1 * e100 + i2 * e010)/ kB / T);
            D += e;
        }
    }
    double ee = ZCO2Vibr12(T);
    double w = D/ee ;
    double ww = w/mass;
    return ww;
}

double AdditionalSolver::EVibr3(double sT, double T, double r, double t)
{
    int L3 = 20;
    double E = 0, D = 0;

    for (int i3 = 0; i3 <= L3; i3++)
    {
        if ((E_CO2(0, 0, i3) < De))
        {
            D += (i3 * e001) * exp(-(i3 * e001) / kB / T);
            //E += E_CO2(0, 0, i3) * exp(-(E_CO2(0, 0, i3) - E_CO2(0, 0, 0)) / (kB * T));
        }
    }

    return D/ZCO2Vibr3(T)/massaCO2;
}

double AdditionalSolver::Lambda12(double sT, double T, double r, double t)
{
    //return tauVibrVVLosev(T, 1);
    //return (3.0*kB*T)/(8.0*getOmega11(T));
    //return (75.0*pow(kB,2)*T)/(32.0*mass*getOmega22(T)) + (3.0*kB*T)/(8.0*getOmega11(T))*Crot();
    double cVibr = CVibr12(T);
    return (3.0*kB*T)/(8.0*getOmega11(sT))*cVibr;
    //return TauVibr(T,1000);
}

double AdditionalSolver::Lambda3(double sT, double T, double r, double t)
{
    double cVibr = CVibr3(T);
    return (3.0*kB*T)/(8.0*getOmega11(sT))*cVibr;
    //return tauVibrVVLosev(T,1000);
}

double AdditionalSolver::CVibr12(double T)
{
    double Z12 = ZCO2Vibr12(T);
    long double S = 0;
    long double F = 0;
    int i1, i2;
    for (i1 = 0; i1 < 20; i1++)
    {
        for (i2 = 0; i2 < 40; i2++)
        {
            if (E_CO2(i1, i2, 0) < De)
            {

                S += (i2 + 1.) * (E_CO2(i1, 0, 0) - E_CO2(0, 0, 0) + E_CO2(0, i2, 0) - E_CO2(0, 0, 0)) * (E_CO2(i1, 0, 0) - E_CO2(0, 0, 0) + E_CO2(0, i2, 0) - E_CO2(0, 0, 0)) * exp(-(E_CO2(i1, 0, 0) - E_CO2(0, 0, 0) + E_CO2(0, i2, 0) - E_CO2(0, 0, 0)) / (kB * T)) / (kB * pow(T, 2));  // (pow(kb * T12, 2));
                //S += (i2 + 1.) * (i1 * e100 + i2 * e010) * (i1 * e100 + i2 * e010) * exp(-(i1 * e100 + i2 * e010) / (kb * T12)) / (kb * pow(T12, 2));
                F += (i2 + 1.) * (E_CO2(i1, 0, 0) - E_CO2(0, 0, 0) + E_CO2(0, i2, 0) - E_CO2(0, 0, 0)) * exp(-(E_CO2(i1, 0, 0) - E_CO2(0, 0, 0) + E_CO2(0, i2, 0) - E_CO2(0, 0, 0)) / (kB * T)) / (sqrt(kB) * T);
                //F += (i2 + 1.) * (i1 * e100 + i2 * e010) * exp(-(i1 * e100 + i2 * e010) / (kb * T12)) / (sqrt(kb) * T12);
            }
        }
    }
    return  (S / Z12 - pow(F / Z12, 2)) / massaCO2;
}

double AdditionalSolver::CVibr3(double T)
{
    int L3 = 10;
    double Z3 = ZCO2Vibr3(T);
    long double S = 0;

    for (int i3 = 0; i3 < L3; i3++)
        if (E_CO2(0, 0, i3) < De)
            S += (i3 * e001) * (i3 * e001) * exp(-(i3 * e001) / (kB * T)) / (kB * pow(T, 2));
    long double F = 0;
    for (int i3 = 0; i3 <= L3; i3++)
        if (E_CO2(0, 0, i3) < De)
        {
            F += (i3 * e001) * exp(-(i3 * e001) / (kB * T)) / (sqrt(kB) * T);
            // F += (i3 * e001) * exp(-(i3 * e001) / (kb * T3)) ;
        }
    return  (S / Z3 - pow(F / Z3, 2)) / massaCO2;
}

double AdditionalSolver::ZCO2Vibr12(double T)
{
    int L1 = 34;//20;//34;
    int L2 = 67;//40;// 67;
    //int L3 = 1;//20;
    long double Z_vibr_12 = 0;

    for (int i1 = 0; i1 < L1; i1++)
    {
        for (int i2 = 0; i2 < L2; i2++)
        {
            //for (int i3 = 0; i3 < L3; i3++) // у Алёны нет этого суммирования
            //{
                if (E_CO2(i1, i2, 0) < De)
                {
                    //Z_vibr_12 += (i2 + 1.) * exp(-(E_CO2(i1, i2, i3) - E_CO2(0, 0, 0)) / (kB * T));
                    Z_vibr_12 += (i2 + 1.) * exp(-(i1 * e100 + i2 * e010) / (kB * T));
                    // += (i2 + 1.) * exp(-i1 * e100 / kB / T) * exp(-i2 * e010 / kB / T);
                }

            //}
        }
    }
    return Z_vibr_12;
}

double AdditionalSolver::ZCO2Vibr3(double T)
{
    int L3 = 20;
    long double Z_vibr_3 = 0;


    for (int i3 = 0; i3 < L3; i3++)
    {
        if (E_CO2(0, 0, i3) < De)
        {
            //Z_vibr_3 += exp(-(E_CO2(0, 0, i3) - E_CO2(0, 0, 0)) / (kB * T));
            Z_vibr_3 += exp(-(i3 * e001) / (kB * T));
        }
    }

    return Z_vibr_3;
}

double AdditionalSolver::tauVibrVVLosev(double T, double P)
{
    P = P/101325;
    double P_tau;
    double a = -26.85;
    double b = 173.22;
    double c = -539.74;
    double d = 0.09645;
    double T13 = pow(T, -1. / 3.);
    double ptau = a + b * pow(T13, 1.) + c * pow(T13, 2.) + d / T13;
    P_tau = exp(ptau); //[atm*s] if pow(10, lgpt)*P1 [Pa*s]
    return (P_tau / P);
}

macroParam AdditionalSolver::bondaryConditionPython(macroParam left, solverParams solParams)
{
    auto mlt =8;
    macroParam right;
    while( 20 > mlt )
    {
        auto values = runPython(left, mlt);
        if(values.size() >= 3)
        {
            right.density = values[0].toDouble();
            right.velocity = values[1].toDouble();
            right.temp = values[2].toDouble();
            right.tempIntr = values.size() > 3 ? values[3].toDouble() : right.temp;
            right.pressure = right.density*UniversalGasConstant/molMass*right.temp;//  right.density/mass * kB * right.temp;
            return right;
        }
        mlt++;
    }
    return macroParam();
}

double AdditionalSolver::zVibr(double startT, double currentT, double density, double pressure)
{
    return 4*pressure*TauVibr(currentT, pressure)/(M_PI*(5*kB*currentT)/(8*getOmega22(currentT)));
}

double AdditionalSolver::getTimeStep(Matrix velosity, Matrix density, Matrix pressure, double delta_h, solverParams solParams, Matrix gammas)
{
    Matrix c = Matrix::SQRT(pressure*solParams.Gamma/density);
    auto temp = velosity + c;
    auto max =*std::max_element(temp.begin(), temp.end());
    return solParams.CFL*pow(delta_h,2)/max;
}

double AdditionalSolver::getTimeStepFull(Matrix velosity, Matrix density, Matrix pressure, double delta_h, solverParams solParams, Matrix gammas)
{
    Matrix c = Matrix::SQRT(pressure*gammas/density);
    auto temp = velosity + c;
    auto max =*std::max_element(temp.begin(), temp.end());
    return solParams.CFL*pow(delta_h,1)/max;
}

double AdditionalSolver::getTimeStep2(Matrix velosity, Matrix density, Matrix pressure, double delta_h, solverParams solParams, double x, Matrix T)
{
    Matrix c = Matrix::SQRT(pressure*solParams.Gamma/density);
    auto temp = velosity + c;
    auto E_in = temp/delta_h;
    //Matrix E_vis;
    for(auto i = 0 ; i <  T.size(); i++)
    {
        auto zco = ZCO2Vibr(T[i]);
        auto cvibr = CVibr(T[i], zco);
        auto element = 2*x/(cvibr* density[i]*delta_h*delta_h);
        E_in.push_back(element);
    }
    auto max =*std::max_element(E_in.begin(), E_in.end());
    return solParams.CFL/max;


}
