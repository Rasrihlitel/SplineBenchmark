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

//Подготовленные данные
struct Data {
	std::vector<double> A;
	std::vector<double> B;
	std::vector<double> C;
	std::vector<double> F;
	std::vector<double> h;
	splineCoeff coef;
};

// -----------------Проверка данных-----------------//

//Подготовка входных координат
void prepareCoord(std::vector<double>& x, std::vector<double>& y);

//Подготовка входных данных
Data prepareData(const std::vector<double>& x, const std::vector<double>& y);

//-----------------Вычисление сплайна-----------------//

//Решение трехдиагональной матрицы 
void sweepMethod(
	const std::vector<double>& A,
	const std::vector<double>& B,
	const std::vector<double>& C,
	const std::vector<double>& F,
	std::vector<double>& cCoef
);

//Вычисление коэффициентов сплайна
splineCoeff calculateSpline(
	const std::vector<double>& x,
	const std::vector<double>& y,
	Data& data
);

//Полный цикл вычисления
splineCoeff calculateCoeff(std::vector<double>& x, std::vector<double>& y);

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