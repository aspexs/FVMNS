#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>
#include <abstaractsolver.h>
#include "argonsolver.h"
#include "nitrogensolver.h"
#include "co2solver.h"
#include "co22tsolver.h"
#include "co22tsolverK.h"
#include "co23tsolver.h"
#include "MixtureSolvers/mixtureco2ar.h"
#include <additionalsolver.h>
#include <chart.h>
#include <QChartView>
#include "oxygensolver.h"
namespace Ui {
class Widget;
}

class Widget : public QWidget
{
    Q_OBJECT

public:
    explicit Widget(QWidget *parent = nullptr);
    ~Widget();
public slots:
    void continueCalc();

private slots:
    void on_comboBox_gas_currentIndexChanged(const QString &gas_);
    void on_comboBox_additionalSolvingType_currentIndexChanged(int index);
    void on_pushButton_additionalSolve_clicked();
    void drawAdditionalSolution();

    void on_comboBox_bulk_currentIndexChanged(int index);

    void on_comboBox_energy_currentIndexChanged(int index);

    void on_comboBox_BC_currentIndexChanged(int index);

    void on_comboBox_typePlot_currentIndexChanged(int index);

    void on_comboBox_timeStep_currentIndexChanged(int index);
    void updatePlot(QVector<double> x, QVector<double> y, double lambda);
    void updateAdditionalPlot(QVector<double> x, QVector<double> y, double lambda);
     void updateAdditional2Plot(QVector<double> x, QVector<double> y, double lambda);
    void updateTime(double time, double error);
    void on_pushButton_start_clicked();

    void on_comboBox_share_currentIndexChanged(int index);

    void on_pushButton_csv_clicked();

    void cancal();
    void pause();

private:
    Ui::Widget *ui;
    AbstaractSolver* solver;
    AdditionalSolver* additionalSolver;
    Chart* chart;
    QChartView* chartView;
    QString gas;
    QVector<double> _x;
    QVector<double> _y;
    QVector<double> _y2;
    bool isAdditionalSolve = true;
    double _time;
};

#endif // WIDGET_H
