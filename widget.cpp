#include "widget.h"
#include "ui_widget.h"
#include <global.h>
#include <QLineEdit>
#include <QDialogButtonBox>
#include <QFormLayout>

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
    chartView->setRenderHint(QPainter::NonCosmeticDefaultPen);
    ui->widget->layout()->addWidget(chartView);
    qRegisterMetaType<QVector<double>>("QVector<double>");
}

Widget::~Widget()
{
    delete ui;
}

void Widget::continueCalc()
{
    pause();
    QDialog dlg;
    dlg.setWindowTitle(tr("Прекратить расчет"));
    QDialogButtonBox *btn_box = new QDialogButtonBox(&dlg);
    btn_box->setStandardButtons(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);

    connect(btn_box, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    connect(btn_box, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);

    QFormLayout *layout = new QFormLayout();
    layout->addRow("Ok - Закончить расчет и сделать выгрузку", new QLabel());
    layout->addRow("Cancel - Продолжить расчет", new QLabel());
    layout->addWidget(btn_box);
    dlg.setLayout(layout);
    // В случае, если пользователь нажал "Ok".
    if(dlg.exec() == QDialog::Accepted)
    {
        solver->breakSolver();
        on_pushButton_csv_clicked();
    }
    else
        pause();
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
        ui->comboBox_additionalSolvingType->addItem("Расчет сдвиговой вязкости (Степенной закон) ");
        ui->comboBox_share->addItem("Степенной закон");
        solver->typeShareVisc = 0;
    }
    else if (gas.contains("N2"))
    {
         molMass = 14.0067e-3*2;
         sigma = 3.621e-10;
         epsilonDevK = 97.5;
         mass = molMass/Nav;
         solver = new NitrogenSolver;
         ui->comboBox_additionalSolvingType->addItem("Расчет сдвиговой вязкости (Степенной закон) ");
         ui->comboBox_additionalSolvingType->addItem("Расчет сдвиговой вязкости (По формулам КТ) ");
         ui->comboBox_additionalSolvingType->addItem("Расчет объемной вязкости");
         ui->comboBox_additionalSolvingType->addItem("C_Tr");
          ui->comboBox_additionalSolvingType->addItem("C_Rot");
          ui->comboBox_additionalSolvingType->addItem("lambda");

         ui->comboBox_share->addItem("Степенной закон");
         ui->comboBox_share->addItem("По формулам КТ");
         solver->typeShareVisc = 1;
         ui->comboBox_bulk->addItem("Только комтепательно-вращательные степени свободы");
         solver->typeBulkVisc = 0;
    }
    else if(gas.contains("CO2"))
    {
        molMass = 44.01e-3;
        sigma = 3.763e-10;
        epsilonDevK = 244;
        mass = 7.306e-26;
        ui->comboBox_additionalSolvingType->addItem("Расчет сдвиговой вязкости (Степенной закон) ");
        ui->comboBox_additionalSolvingType->addItem("Расчет сдвиговой вязкости (По формулам КТ) ");
        ui->comboBox_additionalSolvingType->addItem("Расчет объемной вязкости (Через интеграл столкновений)");
        ui->comboBox_additionalSolvingType->addItem("Расчет объемной вязкости (Через время вращательной релаксации)");
        ui->comboBox_additionalSolvingType->addItem("Расчет колебательной энергии");
        ui->comboBox_additionalSolvingType->addItem("Расчет общей энергии");
        ui->comboBox_additionalSolvingType->addItem("Расчет колебательной теплоемкости");
        ui->comboBox_additionalSolvingType->addItem("Расчет поступательной теплопроводности");
        ui->comboBox_additionalSolvingType->addItem("Расчет колебательной теплопроводности");
        ui->comboBox_additionalSolvingType->addItem("Колебательная стат.сумма");
        ui->comboBox_additionalSolvingType->addItem("Расчет колебательной энергии E12");
        ui->comboBox_additionalSolvingType->addItem("Расчет колебательной энергии E3");
        ui->comboBox_additionalSolvingType->addItem("lambda12");
        ui->comboBox_additionalSolvingType->addItem("lambda3");
         ui->comboBox_additionalSolvingType->addItem("C_Tr");
          ui->comboBox_additionalSolvingType->addItem("C_Rot");

        ui->comboBox_BC->addItem("Ренкина-Гюгонио");
        ui->comboBox_BC->addItem("Законы сохранения (Гир)");
        ui->comboBox_BC->setCurrentIndex(1);
        ui->comboBox_bulk->addItem("Простая постановка");
        ui->comboBox_bulk->addItem("Обычная постановка с учетом колебательных мод (1Т)");
        ui->comboBox_bulk->addItem("Обычная постановка без колебательных мод (2Т)");
        ui->comboBox_bulk->addItem("Без объемной вязкости");


        if(gas.contains("RS"))
             solver = new Co22TSolver;
        else if(gas.contains("HLLE"))
             solver = new Co22TSolverK;
        else if(gas.contains("3T"))
             solver = new Co23TSolver;
        else
             solver = new CO2Solver;
    }
    else if(gas.contains("O2"))
    {
        mass = molMass/Nav;
        solver = new OxygenSolver;
    }
    ui->comboBox_additionalSolvingType->setCurrentIndex(0);

    //mass = molMass/Nav;
    ui->label_sigma->setText("sigma (m): " + QString::number(sigma));
    ui->label_gasMolMass->setText("Молярная масса газа (g/mol): " + QString::number(molMass));
    ui->label_mass->setText("Масса молекулы газа (kg): " + QString::number(mass));
    disconnect(solver, SIGNAL(updateGraph(QVector<double>, QVector<double>, double)), this ,SLOT(updatePlot(QVector<double>, QVector<double>, double)));
    disconnect(solver, SIGNAL(updateTime(double, double)), this, SLOT(updateTime(double, double)));
    disconnect(ui->comboBox_typePlot, SIGNAL(currentIndexChanged(int)), solver ,SLOT(setTypePlot(int )));
    disconnect(ui->pushButton_cancel,  SIGNAL(clicked()), this, SLOT(cancal()));
    disconnect(ui->pushButton_pause, SIGNAL(clicked()), this, SLOT(pause()));
    disconnect(solver, SIGNAL(checkCalc()), this, SLOT (continueCalc()));

    connect(solver, SIGNAL(updateGraph(QVector<double>, QVector<double>, double)), this ,SLOT(updatePlot(QVector<double>, QVector<double>, double)));
    connect(solver, SIGNAL(updateAdditionalGraph(QVector<double>, QVector<double>, double)), this ,SLOT(updateAdditionalPlot(QVector<double>, QVector<double>, double)));
    connect(solver, SIGNAL(updateAdditional2Graph(QVector<double>, QVector<double>, double)), this ,SLOT(updateAdditional2Plot(QVector<double>, QVector<double>, double)));
    connect(solver, SIGNAL(updateTime(double, double)), this, SLOT(updateTime(double, double)));
    connect(ui->comboBox_typePlot, SIGNAL(currentIndexChanged(int)), solver ,SLOT(setTypePlot(int )));
    connect(ui->pushButton_cancel,  SIGNAL(clicked()), this, SLOT(cancal()));
    connect(ui->pushButton_pause, SIGNAL(clicked()), this, SLOT(pause()));
    connect(solver, SIGNAL(checkCalc()), this, SLOT (continueCalc()));
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
    additionalSolver->gamma = 5.0/4.0;
    if(ui->comboBox_gas->currentIndex() == 0) // Ar
        additionalSolver->typeSolve = AdditionalSolver::SHARE_VISC_SUPER_SIMPLE;
    else if (ui->comboBox_gas->currentIndex() == 1) // N2
    {
        switch (index)
        {
            case 0: additionalSolver->typeSolve = AdditionalSolver::SHARE_VISC_SIMPLE; break;
            case 1: additionalSolver->typeSolve = AdditionalSolver::SHARE_VISC_OMEGA; break;
            case 2: additionalSolver->typeSolve = AdditionalSolver::BULC_VISC_SIMPLE; break;
            case 3: additionalSolver->typeSolve = AdditionalSolver::C_TR; break;
            case 4: additionalSolver->typeSolve = AdditionalSolver::C_ROT; break;
            case 5: additionalSolver->typeSolve = AdditionalSolver::LAMBDA_N2; break;
        }
    }
    else if(ui->comboBox_gas->currentIndex() > 1 &&ui->comboBox_gas->currentIndex() < 6) // CO2
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
            case 9: additionalSolver->typeSolve = AdditionalSolver::Z_VIBR; break;

            case 10: additionalSolver->typeSolve = AdditionalSolver::EVIBR12; break;
            case 11: additionalSolver->typeSolve = AdditionalSolver::EVIBR3; break;
            case 12: additionalSolver->typeSolve = AdditionalSolver::LAMBDA12; break;
            case 13: additionalSolver->typeSolve = AdditionalSolver::LAMBDA3; break;
            case 14: additionalSolver->typeSolve = AdditionalSolver::C_TR; break;
            case 15: additionalSolver->typeSolve = AdditionalSolver::C_ROT; break;

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
    solver->typeBC = index;
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
    _y2.resize(y.size());
    chart->setData(x, y, lambda);
    chart->update();
    chartView->repaint();
    repaint();
    update();
}

void Widget::updateAdditionalPlot(QVector<double> x, QVector<double> y, double lambda)
{
    _y2 = y;
     chart->setAdditionalData(x, y, lambda);
     chart->update();
     chartView->repaint();
     repaint();
     update();
}

void Widget::updateAdditional2Plot(QVector<double> x, QVector<double> y, double lambda)
{
    chart->setAdditionalData2(x, y, lambda);
    chart->update();
    chartView->repaint();
    repaint();
    update();
}

void Widget::updateTime(double time, double error)
{
    _time = time;
     ui->label_time->setText("Время расчета :  " + QString::number(time,'g', 10));
     ui->label_error->setText("Суммарное отклонение :  " + QString::number(error,'g', 10));
}

void Widget::on_pushButton_start_clicked()
{
    solver->leftParam.temp = ui->doubleSpinBo_temp->value();
    solver->leftParam.pressure = ui->doubleSpinBox_pressure_2->value();
    solver->leftParam.tempIntr = ui->doubleSpinBox_temp_second->value();
    solver->solParam.Ma = ui->doubleSpinBox_ma->value();
    solver->solParam.CFL = ui->doubleSpinBox_CFL->value();
    solver->solParam.NumCell = ui->spinBox_numCells->value();
    solver->solParam.lambdaSol = ui->spinBox_lenght->value();
    solver->solParam.lambda = 0.001;
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
        else if (ui->comboBox_additionalSolvingType->currentText().contains("Расчет колебательной энергии E12"))
            nameFile = "/vibrEnergy12.csv";
        else if (ui->comboBox_additionalSolvingType->currentText().contains("Расчет колебательной энергии E3"))
            nameFile = "/vibrEnergy3.csv";
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
        else if (ui->comboBox_additionalSolvingType->currentText().contains("lambda12"))
            nameFile = "/lambda12.csv";
        else if (ui->comboBox_additionalSolvingType->currentText().contains("lambda3"))
            nameFile = "/lambda3.csv";
        else if (ui->comboBox_additionalSolvingType->currentText() == "lambda")
            nameFile = "/lambda.csv";
        QFile file(QDir::currentPath() + nameFile);
        if(file.open(QFile::WriteOnly))
        {
            QTextStream out(&file);
            for(int i = 0; i <additionalSolver->iterationVector.size(); i ++)
                out << additionalSolver->iterationVector[i] << ";" << additionalSolver->rezultVector[i] << "\n";

        }
        file.close();
    }
    else
    {
        ///if(ui->comboBox_typePlot->currentIndex() == 0)
        ///    nameFile = "/давление.csv";
        ///else if(ui->comboBox_typePlot->currentIndex() == 1)
        ///    nameFile = "/плотность.csv";
        ///else if(ui->comboBox_typePlot->currentIndex() == 2)
        ///    nameFile = "/скорость.csv";
        ///else if(ui->comboBox_typePlot->currentIndex() == 3)
        ///    nameFile = "/температура.csv";
        ///else if(ui->comboBox_typePlot->currentIndex() == 4)
        ///    nameFile = "/Напряжение.csv";
        ///else if(ui->comboBox_typePlot->currentIndex() == 5)
        ///    nameFile = "/q_tr_rot.csv";
        ///else if(ui->comboBox_typePlot->currentIndex() == 6)
        ///    nameFile = "/Q_vibr.csv";

        QFile file(QDir::currentPath() + "/AllParams_" + QString::number(ui->doubleSpinBox_ma->value()) + "_"+ ui->comboBox_gas->currentText() + ".csv");
        if(file.open(QFile::WriteOnly))
        {
            QTextStream out(&file);
            out  << "X"<< ";"<< "Press" << ";" << "Density" << ";" << "Velosity" << ";"
                << "Temp" << ";"<< "Vibt Temp" << ";"
                << "Pressure tenzor" << ";"<< "Heat vector" << ";"
                << "Vibr heat(3Т - qVibr12, 1T - C_vibr)" << ";"<< "Entalpy" << ";"
                << "vibr heat (3T - qVibr3, 1T - lambda)" << ";"<< "Temp T12 (1T - eta)" << ";" <<"Temp T3(1T - zeta)" << ";"
                << "gamma" << ";;"<< "PR" <<  "\n";
             for(int i = 0; i <_x.size(); i ++)
             {
                out << _x[i]            << ";" << solver->pres[i]   << ";" << solver->U1[i]     << ";" << solver->U2[i]/solver->U1[i] << ";"
                    << solver->T[i]     << ";" << solver->Tv[i]     << ";" << solver->P[i]      << ";" << solver->Q_t[i]              << ";"
                    << solver->Q_v[i]   << ";" << solver->Ent[i]    << ";" << solver->Q_v3[i]   << ";" << solver->T12[i]              << ";"
                    << solver->T3[i]    << ";" << solver->E_Z[i]    << ";" << solver->B_v[i]<< ";" << solver->PR[i]<<"\n" ;
             }
             out << _time;
        }
        // QFile file(QDir::currentPath() + nameFile);
        // if(file.open(QFile::WriteOnly))
        // {
        //     QTextStream out(&file);
        //     for(int i = 0; i <_x.size(); i ++)
        //     {
        //         out << _x[i] << ";" << _y[i] << ";"<< _y2[i] <<"\n";
        //     }
        //     out << _time;
        // }
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
