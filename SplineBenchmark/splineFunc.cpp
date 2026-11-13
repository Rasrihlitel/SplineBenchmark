#include "splineFunc.h"
#include <algorithm>
#include <stdexcept>
#include <fstream>

//Подготовка входных координат
void prepareCoord(std::vector<double>& x, std::vector<double>& y) {
	if (x.size() != y.size())
		throw std::invalid_argument("The number of x and y coordinates does not match");
	if (x.size() < 2)
		throw std::invalid_argument("The number of coordinates cannot be less than 2");
	bool isSort = true;
	for (size_t i = 0; i < x.size() - 1; i++) {
		if (x[i + 1] < x[i]) { isSort = false; break; }
	}
	if (isSort) return;
	std::vector<std::pair<double, double>> coords(x.size());
	for (size_t i = 0; i < x.size(); i++) {
		coords[i] = { x[i], y[i] };
	}
	std::sort(coords.begin(), coords.end());
	for (size_t i = 0; i < x.size(); i++) {
		x[i] = coords[i].first;
		y[i] = coords[i].second;
	}
}

//Решение трехдиагональной матрицы 
void sweepMethod(
	const std::vector<double>& x,
	const std::vector<double>& y,
	SplineWorkspace& ws
) {
	size_t size = x.size();
	for (size_t i = 1; i < size - 1; ++i) {
		double A = x[i] - x[i - 1];
		double C = 2.0 * (x[i] - x[i - 1] + x[i + 1] - x[i]);
		double B = x[i + 1] - x[i];
		double F = 6.0 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) 
			- (y[i] - y[i - 1]) / (x[i] - x[i - 1]));

		double denom = C - ws.alpha[i - 1] * A;
		ws.alpha[i] = B / denom;
		ws.beta[i] = (F - A * ws.beta[i - 1]) / denom;
	}

	for (int i = size - 2; i >= 1; i--)
		ws.coef.c[i] = ws.beta[i] - ws.alpha[i] * ws.coef.c[i + 1];
}

void computeBD(
	const std::vector<double>& x,
	const std::vector<double>& y,
	SplineWorkspace& ws
) {
	size_t size = x.size();
	for (size_t i = 1; i < size; ++i) {
		double h = x[i] - x[i - 1];
		ws.coef.d[i] = (ws.coef.c[i] - ws.coef.c[i - 1]) / h;
		ws.coef.b[i] = (h * ws.coef.c[i] / 2.0)
			- (h * h * ws.coef.d[i] / 6.0)
			+ (y[i] - y[i - 1]) / h;
	}
}

//Вычисление коэффициентов сплайна
splineCoeff& createSpline(
	const std::vector<double>& x,
	const std::vector<double>& y,
	SplineWorkspace& ws
) {
	size_t size = x.size();
	ws.resize(size);
	ws.coef.a = y;
	sweepMethod(x, y, ws);
	computeBD(x, y, ws);
	return ws.coef;
}

//Полный цикл вычисления
const splineCoeff& calculateLoop(
	std::vector<double>& x,
	std::vector<double>& y,
	SplineWorkspace& ws
) {
	prepareCoord(x, y);
	return createSpline(x, y, ws);
}

//Нахождение интервала
int findInterval(const std::vector<double>& x, double xVal) {
	size_t size = x.size();
	int left = 0, right = size - 1;
	while (left < right - 1) {
		int mid = left + (right - left) / 2;
		if (x[mid] <= xVal) left = mid;
		else right = mid;
	}
	return right;
}

//Вычисление в точке
double evaluateSpline(
	const splineCoeff& coef,
	const std::vector<double>& x,
	double xVal
) {
	if (xVal < x.front() || xVal > x.back())
		throw std::out_of_range("Point lies outside the interpolation area");
	int interval = findInterval(x, xVal);
	double diff = xVal - x[interval];
	return coef.a[interval]
		+ coef.b[interval] * diff
		+ (coef.c[interval] / 2.0) * diff * diff
		+ (coef.d[interval] / 6.0) * diff * diff * diff;
}

//Графическое представление сплайна
void exportSplineToCSV(
	const splineCoeff& coef,
	const std::vector<double>& x,
	int numPoints
) {
	if (x.size() < 2)
		throw std::invalid_argument("At least 2 nodes required");
	if (numPoints < 2) numPoints = 2;

	double step = (x.back() - x.front()) / (numPoints - 1);
	std::ofstream file("printSpline.csv");
	if (!file.is_open()) return;

	file << "x;y\n";
	for (int i = 0; i < numPoints; ++i) {
		double xVal = x.front() + i * step;
		double yVal = evaluateSpline(coef, x, xVal);
		file << xVal << ";" << yVal << "\n";
	}
	file.close();
}