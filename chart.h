#ifndef CHART_H
#define CHART_H

#include <QtCharts/QChart>
#include <QLineSeries>
#include <QValueAxis>
using namespace QtCharts;

class Chart: public QChart
{
    Q_OBJECT
public:
    Chart(QGraphicsItem *parent = 0, Qt::WindowFlags wFlags = 0);
    virtual ~Chart();

public slots:
    void setData(QVector<double> x, QVector<double> y, double lambda = 1);
    void setAdditionalData(QVector<double> x, QVector<double> y, double lambda = 1);
private:
    QLineSeries *m_series;
    QLineSeries *a_series;
    QValueAxis *m_axis;
    qreal m_x;

};
//![1]

#endif /* CHART_H */
