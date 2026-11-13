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

// Рабочее пространство
struct SplineWorkspace {
	splineCoeff coef;
	std::vector<double> alpha;
	std::vector<double> beta;
	void resize(int size) {
		if (alpha.size() != size) {
			alpha.resize(size, 0.0);
			beta.resize(size, 0.0);
			coef.a.resize(size, 0.0);
			coef.b.resize(size, 0.0);
			coef.c.resize(size, 0.0);
			coef.d.resize(size, 0.0);
		}
	}
};

// -----------------Подготовка данных-----------------//

//Подготовка входных координат
void prepareCoord(std::vector<double>& x, std::vector<double>& y);

//-----------------Вычисление сплайна-----------------//

//Решение трехдиагональной матрицы 
void sweepMethod(
	const std::vector<double>& x,
	const std::vector<double>& y,
	SplineWorkspace& ws
);

//Вычисление коэффициентов b и d
void computeBD(
	const std::vector<double>& x,
	const std::vector<double>& y,
	SplineWorkspace& ws
);

//Вычисление коэффициентов сплайна
splineCoeff& createSpline(
	const std::vector<double>& x,
	const std::vector<double>& y,
	SplineWorkspace& ws
);

//Полный цикл вычисления
const splineCoeff& calculateLoop(
	std::vector<double>& x,
	std::vector<double>& y,
	SplineWorkspace& ws
);

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