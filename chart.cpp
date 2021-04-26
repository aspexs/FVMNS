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

    a2_series = new QLineSeries(this);
    a2_series->append(0, 0);

    QPen green(Qt::green);
    green.setWidth(3);
    m_series->setPen(green);

    QPen blue(Qt::blue);
    blue.setWidth(3);
    a_series->setPen(blue);

    QPen magenta(Qt::magenta);
    magenta.setWidth(3);
    a2_series->setPen(magenta);

    addSeries(m_series);
    addSeries(a_series);
    addSeries(a2_series);
    createDefaultAxes();
    setAxisX(m_axis, m_series);
    setAxisX(m_axis, a_series);
    setAxisX(m_axis, a2_series);

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
    minY = *std::min_element(y.begin(), y.end());
    maxY = *std::max_element(y.begin(), y.end());
    TintMax =maxY;
    axisY()->setRange(minY, maxY);
    m_series->clear();
    for(auto i = 0; i < x.size(); i++)
        m_series->append(x[i]/lambda, y[i]);
}

void Chart::setAdditionalData(QVector<double> x, QVector<double> y, double lambda)
{
    QVector<double> e;
    e.resize(x.size());
    if(x.size() != y.size() || y == e)
        return;
    auto LocminY = *std::min_element(y.begin(), y.end());
    auto LocmaxY = *std::max_element(y.begin(), y.end());
    if(LocmaxY > maxY)
        maxY = LocmaxY;
    if(LocminY < minY)
        minY = LocminY;
    axisY()->setRange(minY, maxY);
    a_series->clear();
    if(TintMax < 1e-15)
        return;
    for(auto i = 0; i < x.size(); i++)
        a_series->append(x[i]/lambda, y[i]);
}

void Chart::setAdditionalData2(QVector<double> x, QVector<double> y, double lambda)
{
    if(x.size() != y.size())
        return;
    auto LocminY = *std::min_element(y.begin(), y.end());
    auto LocmaxY = *std::max_element(y.begin(), y.end());
    if(LocmaxY > maxY)
        maxY = LocmaxY;
    if(LocminY < minY)
        minY = LocminY;
    if(TintMax < 1e-15)
        return;
    axisY()->setRange(minY, maxY);
    //axisY()->setRange(-0.2, 1.5);
    a2_series->clear();
    for(auto i = 0; i < x.size(); i++)
        a2_series->append(x[i]/lambda, y[i]);
}
