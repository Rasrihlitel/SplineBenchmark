#pragma once
#ifndef SPLINEFUNC_H
#define SPLINEFUNC_H
#include <vector>
#include "realType.h"

//Коэффициенты сплайна
struct splineCoeff {
	std::vector<Real> b;
	std::vector<Real> c;
	std::vector<Real> d;
};

//Рабочее пространство
struct SplineWorkspace {
	splineCoeff coef;
	std::vector<Real> alpha;
	std::vector<Real> beta;
	std::vector<Real> PCR_a;
	std::vector<Real> PCR_b;
	std::vector<Real> PCR_c;
	std::vector<Real> PCR_d;
	void resize(size_t size) {
		if (alpha.size() != size) {
			alpha.assign(size, Real(0));
			beta.assign(size, Real(0));
			coef.b.assign(size, Real(0));
			coef.c.assign(size, Real(0));
			coef.d.assign(size, Real(0));
			PCR_a.assign(size, Real(0));
			PCR_b.assign(size, Real(0));
			PCR_c.assign(size, Real(0));
			PCR_d.assign(size, Real(0));
		}
	}
};

enum class Methods {
	Thomas,
	PCR
};

//-----------------Подготовка данных-----------------//

//Подготовка входных координат
void prepareCoord(std::vector<Real>& x, std::vector<Real>& y);

//-----------------Вычисление сплайна-----------------//

//Решение трехдиагональной матрицы методом прогонки Томаса
void sweepMethod(
	const std::vector<Real>& x,
	const std::vector<Real>& y,
	SplineWorkspace& ws
);

//Вычисление коэффициентов b и d
void computeBD(
	const std::vector<Real>& x,
	const std::vector<Real>& y,
	SplineWorkspace& ws
);

//Решение трехдиагональной матрицы методом циклической редукции
//outOffset - индекс в x, с которого записывается решение
void cyclicReductionParallel(
	std::vector<Real>& a,
	std::vector<Real>& b,
	std::vector<Real>& c,
	std::vector<Real>& d,
	std::vector<Real>& x
);

//Обертка-подготовка для метода циклической редукции
void cyclicReductionMethod(
	const std::vector<Real>& x,
	const std::vector<Real>& y,
	SplineWorkspace& ws
);

//Вычисление коэффициентов сплайна
const splineCoeff& createSpline(
	const std::vector<Real>& x,
	const std::vector<Real>& y,
	SplineWorkspace& ws,
	Methods method = Methods::Thomas
);

//Полный цикл вычисления
const splineCoeff& calculateLoop(
	std::vector<Real>& x,
	std::vector<Real>& y,
	SplineWorkspace& ws,
	Methods method = Methods::Thomas
);

//-----------------Вычисления в точках-----------------//

//Нахождение интервала
int findInterval(const std::vector<Real>& x, Real xVal);

//Вычисление в точке
Real evaluateSpline(
	const splineCoeff& coef,
	const std::vector<Real>& x,
	const std::vector<Real>& y,
	Real xVal
);

//-----------------Графическое представление сплайна-----------------//

//Графическое представление сплайна
void exportSplineToCSV(
	const splineCoeff& coef,
	const std::vector<Real>& x,
	const std::vector<Real>& y,
	int numPoints
);

#endif //SPLINEFUNC_H