#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include "splineFunc.h"
#include "tests.h"

void demonstration() {
    //Функция y = x^3 на отрезке [0, 3]
    std::vector<Real> x = { Real(0), Real(1), Real(2), Real(3), Real(4) };
    std::vector<Real> y = { Real(0), Real(1), Real(8), Real(27), Real(64) };
    std::cout << "Starting points: \n";
    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "Point " << i << " (" << x[i] << ", " << y[i] << ")\n";
    }
    std::cout << "\n";
    //Построение сплайна
    splineCoeff coef;
    SplineWorkspace ws;
    try {
        coef = calculateLoop(x, y, ws);
        std::cout << "The spline is built\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error in calculation: " << e.what() << "\n";
        return;
    }
    //Вывод полученных коэффициентов
    std::cout << "Spline coefficients\n";
    for (size_t i = 0; i < x.size(); ++i) {
        std::cout << "Section " << i << ":\n";
        std::cout << " a = " << y[i] << "\n";
        std::cout << " b = " << coef.b[i] << "\n";
        std::cout << " c = " << coef.c[i] << "\n";
        std::cout << " d = " << coef.d[i] << "\n";
    }
    std::cout << "\n";
    //Проверка в узлах интерполяции
    std::cout << "Check in the nodes\n";
    for (size_t i = 0; i < x.size(); ++i) {
        Real val = evaluateSpline(coef, x, y, x[i]);
        std::cout << "S(" << x[i] << ") = " << val
            << " | Error rate: " << std::abs(val - y[i]) << "\n";
    }
    std::cout << "\n";
    //Проверяем в промежуточных точках
    std::cout << "Check at intermediate points\n";
    std::vector<Real> test_points = { Real(0.5), Real(1.5), Real(2.5) };
    for (Real pt : test_points) {
        Real val = evaluateSpline(coef, x, y, pt);
        Real exact_val = pt * pt * pt;
        std::cout << "S(" << pt << ") = " << val
            << " | correct x^3 = " << exact_val
            << " | Error rate: " << std::abs(val - exact_val) << "\n";
    }
    std::cout << "\n";
    //Заполнение таблицы данными
    exportSplineToCSV(coef, x, y, 300);
}

int main() {
    std::cout << std::fixed << std::setprecision(PRINT_PREC);

    try {
        //Демонстрация
        demonstration();

        //Первый тест - проходит через узлы интерполяции
        std::cout << "Passes through the interpolation nodes? -> " << testInterpolationNodes() << "\n";

        //Второй тест - непрерывность производных в точках
        std::cout << "Are the derivatives continuous? -> " << testContinuity() << "\n";

        //Третий тест - вторая производная в граничных точка равна 0
        std::cout << "The second derivatives at the boundary points are equal to 0? -> " << testBoundaryPoints() << "\n";

        //Четвертый тест - восстановление линейной функции
        std::cout << "The linear function has been reconstructed? -> " << testLinearReproduction() << "\n";

        //Пятый тест - тест на время
        testTime();

    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << "\n";
    }
    return 0;
}