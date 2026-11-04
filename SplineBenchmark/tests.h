#pragma once
#ifndef TESTS_H
#define TESTS_H

//Первый тест - проверка интерполяции в узлах (над y = x^3)
bool testInterpolationNodes();

//Второй тест - проверка непрерывности (над y = x^3)
bool testContinuity();

//Третий тест - вторые производные в граничных точках равны 0 (над y = x^3)
bool testBoundaryPoints();

//Четвертый тест - воспроизведение линейной функции
bool testLinearReproduction();

#endif //TESTS_H