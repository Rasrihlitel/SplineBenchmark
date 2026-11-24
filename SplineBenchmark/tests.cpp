#include <vector>
#include <chrono>
#include <algorithm>
#include <limits>
#include <iostream>
#include <omp.h>
#include "splineFunc.h"

//Первый тест - проверка интерполяции в узлах (над y = x^3)
bool testInterpolationNodes() {
    std::vector<Real> x = { Real(0), Real(1), Real(2), Real(3), Real(4) };
    std::vector<Real> y = { Real(0), Real(1), Real(8), Real(27), Real(64) };
    SplineWorkspace ws;
    splineCoeff coef = calculateLoop(x, y, ws);
    bool ok = true;
    for (size_t i = 0; i < x.size(); ++i) {
        Real val = evaluateSpline(coef, x, y, x[i]);
        if (std::abs(val - y[i]) >= TEST_EPS) ok = false;
    }
    return ok;
}

//Второй тест - проверка непрерывности (над y = x^3)
bool testContinuity() {
    std::vector<Real> x = { Real(0), Real(1), Real(2), Real(3), Real(4) };
    std::vector<Real> y = { Real(0), Real(1), Real(8), Real(27), Real(64) };
    SplineWorkspace ws;
    splineCoeff coef = calculateLoop(x, y, ws);
    bool ok = true;
    for (size_t i = 1; i + 1 < x.size(); ++i) {
        Real xi = x[i];
        Real leftDeriv1 = coef.b[i];
        Real leftDeriv2 = coef.c[i];
        Real diff = xi - x[i + 1];
        Real rightDeriv1 = coef.b[i + 1]
            + coef.c[i + 1] * diff
            + (coef.d[i + 1] / Real(2)) * diff * diff;
        Real rightDeriv2 = coef.c[i + 1] + coef.d[i + 1] * diff;
        if (std::abs(leftDeriv1 - rightDeriv1) > TEST_EPS) ok = false;
        if (std::abs(leftDeriv2 - rightDeriv2) > TEST_EPS) ok = false;
    }
    return ok;
}

//Третий тест - вторые производные в граничных точках равны 0 (над y = x^3)
bool testBoundaryPoints() {
    std::vector<Real> x = { Real(0), Real(1), Real(2), Real(3) };
    std::vector<Real> y = { Real(0), Real(1), Real(8), Real(27) };
    SplineWorkspace ws;
    splineCoeff coef = calculateLoop(x, y,ws);
    bool ok = true;
    if (std::abs(coef.c.front()) >= TEST_EPS) ok = false;
    if (std::abs(coef.c.back()) >= TEST_EPS) ok = false;
    return ok;
}

//Четвертый тест - воспроизведение линейной функции
bool testLinearReproduction() {
    auto f = [](Real x) { return Real(3) * x + Real(1); };
    std::vector<Real> x = { Real(0), Real(1), Real(2), Real(3) };
    std::vector<Real> y;
    for (Real val : x) y.push_back(f(val));
    SplineWorkspace ws;
    splineCoeff coef = calculateLoop(x, y, ws);
    bool ok = true;
    const int intervals = 25;
    for (size_t i = 0; i < x.size() - 1; ++i) {
        for (int k = 1; k < intervals; ++k) {
            Real xi = x[i] + (x[i + 1] - x[i]) * k / intervals;
            Real val = evaluateSpline(coef, x, y, xi);
            if (std::abs(val - f(xi)) >= TEST_EPS) ok = false;
        }
    }
    return ok;
}

//Пятый тест - тест по времени
void testTime() {
    omp_set_num_threads(4);
    std::srand(42);
    int repeats = 50;
    int size = 10000000;
    std::vector<Real> baseX(size);
    std::vector<Real> baseY(size);
    for (int i = 0; i < size; i++) {
        baseX[i] = Real(i) * (Real(1) / Real(size - 1));
        baseY[i] = static_cast<Real>(std::rand());
    }
    SplineWorkspace ws;

    long long timeThomas = 0;
    for (int i = 0; i < repeats; i++) {
        std::vector<Real> x(baseX), y(baseY);
        auto t0 = std::chrono::steady_clock::now();
        calculateLoop(x, y, ws, Methods::Thomas);
        auto t1 = std::chrono::steady_clock::now();
        timeThomas += std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    }
    long long timePCR = 0;
    for (int i = 0; i < repeats; i++) {
        std::vector<Real> x(baseX), y(baseY);
        auto t0 = std::chrono::steady_clock::now();
        calculateLoop(x, y, ws, Methods::PCR);
        auto t1 = std::chrono::steady_clock::now();
        timePCR += std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    }

    double avgThomas = static_cast<double>(timeThomas) / repeats;
    double avgPCR = static_cast<double>(timePCR) / repeats;

    std::cout << "Size: " << size
#ifdef USE_FLOAT
        << " (float)"
#else
        << " (double)"
#endif
        << "\n";
    std::cout << "Thomas: " << avgThomas << " ms\n";
    std::cout << "PCR: " << avgPCR << " ms\n";
}