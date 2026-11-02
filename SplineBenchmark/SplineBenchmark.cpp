#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include "splineFunc.h"


void demonstration() {
    //Функция y = x^3 на отрезке [0, 3]
    std::vector<double> x = { 0.0, 1.0, 2.0, 3.0 };
    std::vector<double> y = { 0.0, 1.0, 8.0, 27.0 };
    std::cout << "Starting points: \n";
    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "Point " << i << " (" << x[i] << ", " << y[i] << ")\n";
    }
    std::cout << "\n";
    //Построение сплайна
    splineCoeff coef;
    try {
        coef = calculateCycle(x, y);
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
        std::cout << " a = " << coef.a[i] << "\n";
        std::cout << " b = " << coef.b[i] << "\n";
        std::cout << " c = " << coef.c[i] << "\n";
        std::cout << " d = " << coef.d[i] << "\n";
    }
    std::cout << "\n";
    //Проверка в узлах интерполяции
    std::cout << "Check in the nodes\n";
    for (size_t i = 0; i < x.size(); ++i) {
        double val = evaluateSpline(coef, x, x[i]);
        std::cout << "S(" << x[i] << ") = " << val
            << " | Error rate: " << std::abs(val - y[i]) << "\n";
    }
    std::cout << "\n";
    //Проверяем в промежуточных точках
    std::cout << "Check at intermediate points\n";
    std::vector<double> test_points = { 0.5, 1.5, 2.5 };
    for (double pt : test_points) {
        double val = evaluateSpline(coef, x, pt);
        double exact_val = pt * pt * pt;
        std::cout << "S(" << pt << ") = " << val
            << " | correct x^3 = " << exact_val
            << " | Error rate: " << std::abs(val - exact_val) << "\n";
    }
    std::cout << "\n";
    //Заполнение таблицы данными
    exportSplineToCSV(coef, x, 300);
}

//Первый тест - проверка интерполяции в узлах (над y = x^3)
bool testInterpolationNodes() {
    std::vector<double> x = { 0.0, 1.0, 2.0, 3.0 };
    std::vector<double> y = { 0.0, 1.0, 8.0, 27.0 };
    splineCoeff coef = calculateCycle(x, y);
    bool ok = true;
    for (size_t i = 0; i < x.size(); ++i) {
        double val = evaluateSpline(coef, x, x[i]);
        if (std::abs(val - y[i]) >= 1e-10) ok = false;
    }
    return ok;
}

//Второй тест - проверка непрерывности (над y = x^3)
bool testContinuity() {
    std::vector<double> x = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    std::vector<double> y = { 0.0, 1.0, 8.0, 27.0, 64.0 };
    splineCoeff coef = calculateCycle(x, y);
    bool ok = true;
    const double small = 1e-10;
    for (size_t i = 1; i + 1 < x.size(); ++i) {
        double xi = x[i];
        double leftDeriv1 = coef.b[i];
        double leftDeriv2 = coef.c[i];
        double diff = xi - x[i + 1]; 
        double rightDeriv1 = coef.b[i + 1]
            + coef.c[i + 1] * diff
            + (coef.d[i + 1] / 2.0) * diff * diff;
        double rightDeriv2 = coef.c[i + 1] + coef.d[i + 1] * diff;
        if (std::abs(leftDeriv1 - rightDeriv1) > small) ok = false;
        if (std::abs(leftDeriv2 - rightDeriv2) > small) ok = false;
    }
    return ok;
}

//Третий тест - вторые производные в граничных точках равны 0 (над y = x^3)
bool testBoundaryPoints() {
    std::vector<double> x = { 0.0, 1.0, 2.0, 3.0 };
    std::vector<double> y = { 0.0, 1.0, 8.0, 27.0 };

    splineCoeff coef = calculateCycle(x, y);
    bool ok = true;
    if (std::abs(coef.c.front()) >= 1e-10) ok = false;
    if (std::abs(coef.c.back()) >= 1e-10) ok = false;
    return ok;
}

//Четвертый тест - воспроизведение линейной функции
bool testLinearReproduction() {
    auto f = [](double x) { return 3.0 * x + 1.0; };
    std::vector<double> x = { 0.0, 1.0, 2.0, 3.0 };
    std::vector<double> y;
    for (double val : x) y.push_back(f(val));
    splineCoeff coef = calculateCycle(x, y);
    bool ok = true;
    const int intervals = 25;
    for (size_t i = 0; i < x.size() - 1; ++i) {
        for (int k = 1; k < intervals; ++k) {
            double xi = x[i] + (x[i + 1] - x[i]) * k / intervals;
            double val = evaluateSpline(coef, x, xi);
            if (std::abs(val - f(xi)) >= 1e-10) ok = false;
        }
    }
    return ok;
}

int main() {
    std::cout << std::fixed << std::setprecision(6);

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
    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << "\n";
    }
    return 0;
}