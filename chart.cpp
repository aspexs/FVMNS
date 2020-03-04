#include "chart.h"
#include <QtCharts/QAbstractAxis>
#include <QtCharts/QSplineSeries>
#include <QtCharts/QValueAxis>
#include <QtCore/QRandomGenerator>
#include <QtCore/QDebug>

Chart::Chart(QGraphicsItem *parent, Qt::WindowFlags wFlags):
    QChart(QChart::ChartTypeCartesian, parent, wFlags),
    m_series(nullptr),
    m_axis(new QtCharts::QValueAxis)
{

    m_series = new QLineSeries(this);
    m_series->append(0, 0);

    a_series = new QLineSeries(this);
    a_series->append(0, 0);
    QPen green(Qt::red);

    green.setWidth(3);

    QPen blue(Qt::blue);
    blue.setWidth(3);
    m_series->setPen(green);
    a_series->setPen(blue);

    addSeries(m_series);
    addSeries(a_series);
    createDefaultAxes();
    setAxisX(m_axis, m_series);
    setAxisX(m_axis, a_series);

    m_axis->setTickCount(10);

    axisX()->setRange(0, 1);
    axisY()->setRange(0, 1);
}

Chart::~Chart()
{

}

void Chart::setData(QVector<double> x, QVector<double> y, double lambda)
{
    if(x.size() != y.size())
        return;
    axisX()->setRange(*std::min_element(x.begin(), x.end())/lambda, *std::max_element(x.begin(), x.end())/lambda);
    axisY()->setRange(*std::min_element(y.begin(), y.end()), *std::max_element(y.begin(), y.end()));
    m_series->clear();
    for(auto i = 0; i < x.size(); i++)
        m_series->append(x[i]/lambda, y[i]);
}

void Chart::setAdditionalData(QVector<double> x, QVector<double> y, double lambda)
{
    if(x.size() != y.size())
        return;
    a_series->clear();
    for(auto i = 0; i < x.size(); i++)
        a_series->append(x[i]/lambda, y[i]);
}
