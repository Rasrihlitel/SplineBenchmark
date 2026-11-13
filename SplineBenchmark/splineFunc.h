#pragma once
#ifndef SPLINEFUNC_H
#define SPLINEFUNC_H
#include <vector>

//Коэффициенты сплайна
struct splineCoeff {
	std::vector<double> a;
	std::vector<double> b;
	std::vector<double> c;
	std::vector<double> d;
};

// -----------------Подготовка данных-----------------//

//Подготовка входных координат
void prepareCoord(std::vector<double>& x, std::vector<double>& y);

//Подготовка входных данных
splineCoeff prepareData(const std::vector<double>& x, const std::vector<double>& y);

//-----------------Вычисление сплайна-----------------//

//Решение трехдиагональной матрицы 
void sweepMethod(
	const std::vector<double>& x,
	const std::vector<double>& y,
	std::vector<double>& cCoef
);

//Вычисление коэффициентов b и d
void computeBD(
	const std::vector<double>& x,
	const std::vector<double>& y,
	splineCoeff& data
);

//Вычисление коэффициентов сплайна
splineCoeff createSpline(
	const std::vector<double>& x,
	const std::vector<double>& y,
	splineCoeff& data
);

//Полный цикл вычисления
splineCoeff calculateLoop(std::vector<double>& x, std::vector<double>& y);

//-----------------Вычисления в точках-----------------//

//Нахождение интервала
int findInterval(const std::vector<double>& x, double xVal);

//Вычисление в точке
double evaluateSpline(
	const splineCoeff& coef,
	const std::vector<double>& x,
	double xVal
);

//-----------------Графическое представление сплайна-----------------//

//Графическое представление сплайна
void exportSplineToCSV(
	const splineCoeff& coef,
	const std::vector<double>& x,
	int numPoints
);

#endif //SPLINEFUNC_H