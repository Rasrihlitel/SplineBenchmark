#include <vector>
#include <chrono>
#include <algorithm>
#include <limits>
#include <iostream>
#include "splineFunc.h"

//Первый тест - проверка интерполяции в узлах (над y = x^3)
bool testInterpolationNodes() {
    std::vector<double> x = { 0.0, 1.0, 2.0, 3.0 };
    std::vector<double> y = { 0.0, 1.0, 8.0, 27.0 };
    splineCoeff coef = calculateLoop(x, y);
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
    splineCoeff coef = calculateLoop(x, y);
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

    splineCoeff coef = calculateLoop(x, y);
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
    splineCoeff coef = calculateLoop(x, y);
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

//Пятый тест - тест по времени
void testTime() {
    std::srand(42);
    int repeats = 500;
    int size = 1000000;
    std::vector<double> baseX(size);
    std::vector<double> baseY(size);
    for (int i = 0; i < size; i++) {
        baseX[i] = i * (1.0 / (size - 1));
        baseY[i] = static_cast<double>(std::rand());
    }
    long long totalTime = 0;
    for (int i = 0; i < repeats; i++) {
        std::vector<double> x(baseX), y(baseY);
        auto start = std::chrono::steady_clock::now();
        calculateLoop(x, y);
        auto end = std::chrono::steady_clock::now();
        totalTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    }
    double averageTime = static_cast<double>(totalTime) / static_cast<double>(repeats);
    std::cout << "Size: " << size << std::endl;
    std::cout << "Average time: " << averageTime << " ms" << std::endl;
}