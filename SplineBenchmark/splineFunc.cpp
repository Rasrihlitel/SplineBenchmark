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
		if (x[i + 1] < x[i]) isSort = false;
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

//Подготовка входных данных
splineCoeff prepareData(const std::vector<double>& x, const std::vector<double>& y) {
	size_t size = x.size();
	splineCoeff coef;
	coef.a = y;
	coef.b.resize(size, 0.0);
	coef.c.resize(size, 0.0);
	coef.d.resize(size, 0.0);
	return coef;
}

//Решение трехдиагональной матрицы 
void sweepMethod(
	const std::vector<double>& x,
	const std::vector<double>& y,
	std::vector<double>& cCoef
) {
	size_t size = x.size();
	std::vector<double> alpha(size, 0.0);
	std::vector<double>	beta(size, 0.0);
	for (size_t i = 1; i < size - 1; ++i) {
		double A = x[i] - x[i - 1];
		double C = 2.0 * (x[i] - x[i - 1] + x[i + 1] - x[i]);
		double B = x[i + 1] - x[i];
		double F = 6.0 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) 
			- (y[i] - y[i - 1]) / (x[i] - x[i - 1]));
		double denom = C - alpha[i - 1] * A;
		alpha[i] = B / denom;
		beta[i] = (F - A * beta[i - 1]) / denom;
	}
	cCoef[size - 1] = (6.0 * ((y[size - 1] - y[size - 1]) / (x[size - 1] - x[size - 2])
		- (y[size - 2] - y[size - 3]) / (x[size - 1] - x[size - 2]))
		- (x[size - 2] - x[size - 3]) * beta[size - 2])
		/ (2.0 * (x[size - 2] - x[size - 3] + x[size - 1] - x[size - 2])
		- (x[size - 2] - x[size - 3]) * alpha[size - 2]);

	for (int i = size - 2; i >= 0; i--)
		cCoef[i] = beta[i] - alpha[i] * cCoef[i + 1];
}

void computeBD(
	const std::vector<double>& x,
	const std::vector<double>& y,
	splineCoeff& coef
) {
	size_t size = x.size();
	for (size_t i = 1; i < size; ++i) {
		coef.d[i] = (coef.c[i] - coef.c[i - 1]) / (x[i] - x[i - 1]);
		coef.b[i] = (x[i] - x[i - 1]) * coef.c[i] / 2.0
			- (x[i] - x[i - 1]) * (x[i] - x[i - 1]) * coef.d[i] / 6.0
			+ (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
	}
}

//Вычисление коэффициентов сплайна
splineCoeff createSpline(
	const std::vector<double>& x,
	const std::vector<double>& y,
	splineCoeff& coef
) {
	sweepMethod(x, y, coef.c);
	computeBD(x, y, coef);
	return coef;
}

//Полный цикл вычисления
splineCoeff calculateLoop(std::vector<double>& x, std::vector<double>& y) {
	prepareCoord(x, y);
	splineCoeff data = prepareData(x, y);
	return createSpline(x, y, data);
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