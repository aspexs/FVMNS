#include "widget.h"
#include "ui_widget.h"
#include <global.h>

extern double sigma;
extern double epsilonDevK;
extern double molMass;
extern double mass;

Widget::Widget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Widget)
{
    ui->setupUi(this);
    additionalSolver = new AdditionalSolver();
    connect(additionalSolver,SIGNAL(completeSolution()), this, SLOT(drawAdditionalSolution()));

    on_comboBox_gas_currentIndexChanged(ui->comboBox_gas->currentText());

    chart = new Chart();
    chart->setTitle("");
    chart->legend()->hide();
    chartView = new QChartView(chart);
    chart->setAnimationOptions(QChart::AllAnimations);
    chartView->setRenderHint(QPainter::Antialiasing);
    ui->widget->layout()->addWidget(chartView);
    qRegisterMetaType<QVector<double>>("QVector<double>");

}

Widget::~Widget()
{
    delete ui;
}

void Widget::on_comboBox_gas_currentIndexChanged(const QString &gas_)
{
    ui->comboBox_additionalSolvingType->clear();
    ui->comboBox_bulk->clear();
    ui->comboBox_energy->clear();
    ui->comboBox_BC->clear();
    ui->comboBox_share->clear();
    gas = gas_;
    if(gas.contains("Ar"))
    {
        molMass = 39.9e-3;
        solver = new ArgonSolver;
        ui->comboBox_additionalSolvingType->addItem("Расчет сдвиговой вязкости (Супер простая постановка) ");
        ui->comboBox_share->addItem("Супер простая постановка");
        solver->typeShareVisc = 0;
    }
    else if (gas.contains("N2"))
    {
         molMass = 14.0067e-3*2;
         sigma = 3.621e-10;
         epsilonDevK = 97.5;
         solver = new NitrogenSolver;
         ui->comboBox_additionalSolvingType->addItem("Расчет сдвиговой вязкости (Простая постановка ) ");
         ui->comboBox_additionalSolvingType->addItem("Расчет сдвиговой вязкости (Через омега интегралы) ");
         ui->comboBox_additionalSolvingType->addItem("Расчет объемной вязкости");

         ui->comboBox_share->addItem("Простая постановка");
         ui->comboBox_share->addItem("Через омега интегралы");
         solver->typeShareVisc = 1;
         ui->comboBox_bulk->addItem("Вязкость простая");
         solver->typeBulkVisc = 0;
    }
    else if (gas.contains("CO2") && !gas.contains("2T"))
    {
         molMass = 44.01e-3;
         sigma = 3.763e-10;
         epsilonDevK = 244;
         solver = new CO2Solver;
         ui->comboBox_additionalSolvingType->addItem("Расчет сдвиговой вязкости (Простая постановка ) ");
         ui->comboBox_additionalSolvingType->addItem("Расчет сдвиговой вязкости (Через омега интегралы) ");
         ui->comboBox_additionalSolvingType->addItem("Расчет объемной вязкости (old) ");
         ui->comboBox_additionalSolvingType->addItem("Расчет объемной вязкости (new) ");
         ui->comboBox_additionalSolvingType->addItem("Расчет колебательной энергии");
         ui->comboBox_additionalSolvingType->addItem("Расчет общей энергии");
         ui->comboBox_additionalSolvingType->addItem("Расчет колебательной теплоемкости");
         ui->comboBox_additionalSolvingType->addItem("Расчет поступательной теплопроводности");
         ui->comboBox_additionalSolvingType->addItem("Расчет колебательной теплопроводности");


         ui->comboBox_bulk->addItem("Old");
         ui->comboBox_bulk->addItem("New");                                 // not ready
         ui->comboBox_bulk->addItem("Without bulk viscosity");
         //ui->comboBox_energy->addItem("Из давлнеия");                       // not ready
         ui->comboBox_energy->addItem("Строгое кинетическое уравнение");
         ui->comboBox_BC->addItem("Ренкина-Гюгонио");
         //ui->comboBox_BC->addItem("Решение уравнения Оли");                 // not ready
         ui->comboBox_bulk->setCurrentIndex(0);
         ui->comboBox_energy->setCurrentIndex(0);
         ui->comboBox_BC->setCurrentIndex(0);
    }
    else
    {
        molMass = 44.01e-3;
        sigma = 3.763e-10;
        epsilonDevK = 244;
        solver = new Co22TSolver;
        ui->comboBox_bulk->addItem("Old");
        ui->comboBox_bulk->addItem("New");                                 // not ready
                                  // not ready
        //ui->comboBox_energy->addItem("Из давлнеия");                       // not ready
        ui->comboBox_energy->addItem("Строгое кинетическое уравнение");
        ui->comboBox_BC->addItem("Ренкина-Гюгонио");
        //ui->comboBox_BC->addItem("Решение уравнения Оли");                 // not ready
        ui->comboBox_bulk->setCurrentIndex(0);
        ui->comboBox_energy->setCurrentIndex(0);
        ui->comboBox_BC->setCurrentIndex(0);
    }
    ui->comboBox_additionalSolvingType->setCurrentIndex(0);

    mass = molMass/Nav;
    ui->label_sigma->setText("sigma (angstrem): " + QString::number(sigma));
    ui->label_gasMolMass->setText("Молярная масса газа: " + QString::number(molMass));
    ui->label_mass->setText("Масса молекулы газа: " + QString::number(mass));
    disconnect(solver, SIGNAL(updateGraph(QVector<double>, QVector<double>, double)), this ,SLOT(updatePlot(QVector<double>, QVector<double>, double)));
    disconnect(solver, SIGNAL(updateTime(double)), this, SLOT(updateTime(double)));
    disconnect(ui->comboBox_typePlot, SIGNAL(currentIndexChanged(int)), solver ,SLOT(setTypePlot(int )));
    disconnect(ui->pushButton_cancel,  SIGNAL(clicked()), this, SLOT(cancal()));
    disconnect(ui->pushButton_pause, SIGNAL(clicked()), this, SLOT(pause()));

    connect(solver, SIGNAL(updateGraph(QVector<double>, QVector<double>, double)), this ,SLOT(updatePlot(QVector<double>, QVector<double>, double)));
     connect(solver, SIGNAL(updateAdditionalGraph(QVector<double>, QVector<double>, double)), this ,SLOT(updateAdditionalPlot(QVector<double>, QVector<double>, double)));
    connect(solver, SIGNAL(updateTime(double)), this, SLOT(updateTime(double)));
    connect(ui->comboBox_typePlot, SIGNAL(currentIndexChanged(int)), solver ,SLOT(setTypePlot(int )));
    connect(ui->pushButton_cancel,  SIGNAL(clicked()), this, SLOT(cancal()));
    connect(ui->pushButton_pause, SIGNAL(clicked()), this, SLOT(pause()));

}

void Widget::on_comboBox_additionalSolvingType_currentIndexChanged(int index)
{
    if(index > 1)
    {
        ui->doubleSpinBox_density->setEnabled(true);
        ui->doubleSpinBox_pressure->setEnabled(true);
    }
    else
    {
        ui->doubleSpinBox_density->setEnabled(false);
        ui->doubleSpinBox_pressure->setEnabled(false);
    }
    additionalSolver->pressure = ui->doubleSpinBox_pressure->value();
    additionalSolver->density = ui->doubleSpinBox_density->value();
    if(ui->comboBox_gas->currentIndex() == 0) // Ar
        additionalSolver->typeSolve = AdditionalSolver::SHARE_VISC_SUPER_SIMPLE;
    else if (ui->comboBox_gas->currentIndex() == 1) // N2
    {
        switch (index)
        {
            case 0: additionalSolver->typeSolve = AdditionalSolver::SHARE_VISC_SIMPLE; break;
            case 1: additionalSolver->typeSolve = AdditionalSolver::SHARE_VISC_OMEGA; break;
            case 2: additionalSolver->typeSolve = AdditionalSolver::BULC_VISC_SIMPLE; break;
        }
    }
    else if(ui->comboBox_gas->currentIndex() == 2) // CO2
    {
        switch (index)
        {
            case 0: additionalSolver->typeSolve = AdditionalSolver::SHARE_VISC_SIMPLE; break;
            case 1: additionalSolver->typeSolve = AdditionalSolver::SHARE_VISC_OMEGA; break;
            case 2: additionalSolver->typeSolve = AdditionalSolver::BULC_VISC_OLD; break;
            case 3: additionalSolver->typeSolve = AdditionalSolver::BULC_VISC_NEW; break;
            case 4: additionalSolver->typeSolve = AdditionalSolver::VIBR_ENERGY; break;
            case 5: additionalSolver->typeSolve = AdditionalSolver::ALL_ENERGY; break;
            case 6: additionalSolver->typeSolve = AdditionalSolver::C_VIBR; break;
            case 7: additionalSolver->typeSolve = AdditionalSolver::LAMBDA_TR; break;
            case 8: additionalSolver->typeSolve = AdditionalSolver::LAMBDA_VIBR; break;

            default: additionalSolver->typeSolve = AdditionalSolver::SHARE_VISC_SIMPLE;
        }
    }
}

void Widget::on_pushButton_additionalSolve_clicked()
{
    isAdditionalSolve = true;
    additionalSolver->step =        ui->doubleSpinBox_step->value();
    additionalSolver->stopValue =   ui->spinBox_stop->value();
    additionalSolver->startValue = ui->spinBox_start->value();
    additionalSolver->solve();
}

void Widget::drawAdditionalSolution()
{
    chart->setData(additionalSolver->iterationVector, additionalSolver->rezultVector, additionalSolver->step);
    chart->update();
    chartView->repaint();
    repaint();
    update();
}

void Widget::on_comboBox_bulk_currentIndexChanged(int index)
{
    solver->typeBulkVisc = ui->comboBox_bulk->currentIndex();
}

void Widget::on_comboBox_energy_currentIndexChanged(int index)
{

}

void Widget::on_comboBox_BC_currentIndexChanged(int index)
{

}

void Widget::on_comboBox_typePlot_currentIndexChanged(int index)
{

}

void Widget::on_comboBox_timeStep_currentIndexChanged(int index)
{

}

void Widget::updatePlot(QVector<double> x, QVector<double> y, double lambda)
{
    _x= x;
    _y = y;
    chart->setData(x, y, lambda);
    chart->update();
    chartView->repaint();
    repaint();
    update();
}

void Widget::updateAdditionalPlot(QVector<double> x, QVector<double> y, double lambda)
{
     chart->setAdditionalData(x, y, lambda);
     chart->update();
     chartView->repaint();
     repaint();
     update();
}

void Widget::updateTime(double time)
{
     ui->label_time->setText("Время расчета :  " + QString::number(time,'g', 10));
}

void Widget::on_pushButton_start_clicked()
{
    solver->leftParam.temp = ui->doubleSpinBo_temp->value();
    solver->leftParam.pressure = ui->doubleSpinBox_pressure_2->value();
    solver->solParam.Ma = ui->doubleSpinBox_ma->value();
    solver->solParam.CFL = ui->doubleSpinBox_CFL->value();
    solver->solParam.NumCell = ui->spinBox_numCells->value();
    solver->solParam.lambdaSol = ui->spinBox_lenght->value();
    solver->solParam.lambda = ui->doubleSpinBox_lambda->value();
    solver->solParam.PlotIter = ui->spinBox_timeplot->value();
    solver->solParam.t_fin = ui->doubleSpinBox_tFin->value();
    switch (ui->comboBox_gamma->currentIndex())
    {
        case 0: solver->solParam.Gamma = 5.0/3 ; break;
        case 1: solver->solParam.Gamma = 1.4   ; break;
        case 2: solver->solParam.Gamma = 1.235 ; break;
        default:solver->solParam.Gamma = 1.4   ;
    }
    isAdditionalSolve = false;
    solver->start();
}

void Widget::on_comboBox_share_currentIndexChanged(int index)
{
    if(gas.contains("Ar"))
    {
        solver->typeShareVisc = 0;
    }
    else
    {
        solver->typeShareVisc = index + 1;
    }
}

void Widget::on_pushButton_csv_clicked()
{
    QString nameFile = "";
    if(isAdditionalSolve)
    {
        if(ui->comboBox_additionalSolvingType->currentText().contains("сдвиговой"))
            nameFile = "/shareVisc.csv";
        else if (ui->comboBox_additionalSolvingType->currentText().contains("объемной"))
            nameFile = "/bulkVisc.csv";
        else if (ui->comboBox_additionalSolvingType->currentText().contains("колебательной энергии"))
            nameFile = "/vibrEnergy.csv";
        else if (ui->comboBox_additionalSolvingType->currentText().contains("общей энергии"))
            nameFile = "/allEnergy.csv";
        else if (ui->comboBox_additionalSolvingType->currentText().contains("колебательной теплоемкости"))
            nameFile = "/CVibr.csv";
        else if (ui->comboBox_additionalSolvingType->currentText().contains("поступательной теплопроводности"))
            nameFile = "/lambdaTr.csv";
        else if (ui->comboBox_additionalSolvingType->currentText().contains("колебательной теплопроводности"))
            nameFile = "/lambdaVibr.csv";

        QFile file(QDir::currentPath() + nameFile);
        if(file.open(QFile::WriteOnly))
        {
            QTextStream out(&file);
            for(int i = 0; i <additionalSolver->iterationVector.size(); i ++)
            {
                out << additionalSolver->iterationVector[i] << ";" << additionalSolver->rezultVector[i] << "\n";
            }
        }
        file.close();
    }
    else
    {
        if(ui->comboBox_typePlot->currentIndex() == 0)
            nameFile = "/давление.csv";
        else if(ui->comboBox_typePlot->currentIndex() == 1)
            nameFile = "/плотность.csv";
        else if(ui->comboBox_typePlot->currentIndex() == 2)
            nameFile = "/скорость.csv";
        else if(ui->comboBox_typePlot->currentIndex() == 3)
            nameFile = "/температура.csv";

        QFile file(QDir::currentPath() + nameFile);
        if(file.open(QFile::WriteOnly))
        {
            QTextStream out(&file);
            for(int i = 0; i <_x.size(); i ++)
            {
                out << _x[i] << ";" << _y[i] << "\n";
            }
        }
        file.close();
    }
}

void Widget::cancal()
{
    solver->breakSolver();
}

void Widget::pause()
{
    if(solver->pauseSolve)
        ui->pushButton_pause->setText("Пауза");
    else
         ui->pushButton_pause->setText("Возоновить");
    solver->pause();
}
