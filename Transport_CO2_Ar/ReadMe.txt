///////////////////////////////////////////////////////////////////////////////
/// Автор: Баталов Семен, 2022.
/// Инструкция для использования инструментов расчета коэффициентов переноса.
///////////////////////////////////////////////////////////////////////////////

> Основные классы для работы:
    1) tc_2::DataDc     (см. transport_m_2.h)
    2) tc_2::ComputerDc (см. transport_m_2.h)
    
> Примеры использования можно увидеть в файле main.cpp.

> Объект класса tc_2::DataDc хранит неизменяемую информацию о смеси СO2-Ar. 
    Его могут использовать несколько потоков для чтения данных из него. 

> Объект класса tc_2::ComputerDc предоставляет методы для расчета всех 
    коэффициентов переноса. Использует данные из объекта класса tc_2::DataDc. 
    Хранит все коэффициенты внутри себя, можно в любой момент к ним 
    обратиться. Предоставляет доступ не только к коэффициентам переноса, но 
    и омега-интегралам, интегральным скобкам и удельным теплоемкостям 
    для СO2.

> Пример использования классов:
    
    // Создаем объекты
    tc_2::DataDc     data;
    tc_2::ComputerDc computer;

    // Инициализируем и связываем
    dataDc.initialize();
    computer.initialize();
    computer.link(data);

    // Расчет значений в определенной точке
    computer.compute(t, t12, t3, x, p);

    // Коэффициенты переноса
    qDebug() << computer.transport().tLambda();
    qDebug() << computer.transport().rLambda();
    qDebug() << computer.transport().cLambda();
    qDebug() << computer.transport().vLambdaT12();
    qDebug() << computer.transport().vLambdaT3();
    qDebug() << computer.transport().sViscosity();
    qDebug() << computer.transport().bViscosity();
    qDebug() << computer.transport().diffusion();
    qDebug() << computer.transport().tDiffusion();
    
    // Омега-интегралы
    qDebug() << computer.omega().omega11();
    qDebug() << computer.omega().omega12();
    qDebug() << computer.omega().omega13();
    qDebug() << computer.omega().omega22();
    qDebug() << computer.omega().aa();
    qDebug() << computer.omega().bb();
    qDebug() << computer.omega().cc();
    
    // Удельные теплоемкости
    qDebug() << computer.heat().cvT12();
    qDebug() << computer.heat().cvT3();
    
    // Интегральные скобки
    qDebug() << computer.bracket().lambda();
    qDebug() << computer.bracket().lambda00();
    qDebug() << computer.bracket().lambda01();
    qDebug() << computer.bracket().lambda11();
    qDebug() << computer.bracket().eta();
    qDebug() << computer.bracket().h00();
    qDebug() << computer.bracket().beta01();
    qDebug() << computer.bracket().beta11();
    qDebug() << computer.bracket().beta0011();
    qDebug() << computer.bracket().lambdaInt();
    