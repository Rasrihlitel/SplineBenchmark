#include "splineFunc.h"
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <fstream>
#include <omp.h>

//Подготовка входных координат
void prepareCoord(std::vector<Real>& x, std::vector<Real>& y) {
	if (x.size() != y.size())
		throw std::invalid_argument("The number of x and y coordinates does not match");
	if (x.size() < 2)
		throw std::invalid_argument("The number of coordinates cannot be less than 2");
	bool isSort = true;
	for (size_t i = 0; i < x.size() - 1; i++) {
		if (x[i + 1] < x[i]) { isSort = false; break; }
	}
	if (isSort) return;
	std::vector<std::pair<Real, Real>> coords(x.size());
	for (size_t i = 0; i < x.size(); i++) {
		coords[i] = { x[i], y[i] };
	}
	std::sort(coords.begin(), coords.end());
	for (size_t i = 0; i < x.size(); i++) {
		x[i] = coords[i].first;
		y[i] = coords[i].second;
	}
}

//Решение трехдиагональной матрицы методом прогонки Томаса
void sweepMethod(
	const std::vector<Real>& x,
	const std::vector<Real>& y,
	SplineWorkspace& ws
) {
	const size_t size = x.size();
	for (size_t i = 1; i < size - 1; ++i) {
		Real A = x[i] - x[i - 1];
		Real B = x[i + 1] - x[i];
		Real C = Real(2) * (A + B);
		Real F = Real(6) * ((y[i + 1] - y[i]) / B
			- (y[i] - y[i - 1]) / A);
		Real denom = Real(1) / (C - ws.alpha[i - 1] * A);
		ws.alpha[i] = B * denom;
		ws.beta[i] = (F - A * ws.beta[i - 1]) * denom;
	}
	for (int i = size - 2; i >= 1; i--)
		ws.coef.c[i] = ws.beta[i] - ws.alpha[i] * ws.coef.c[i + 1];
}

//Решение трехдиагональной матрицы методом циклической редукции
void cyclicReductionParallel(
    std::vector<Real>& a,
    std::vector<Real>& b,
    std::vector<Real>& c,
    std::vector<Real>& d,
    std::vector<Real>& x
) {
	size_t size = b.size() - 2;
    if (size == 0) return;
    if (size == 1) { x[1] = d[0] / b[0]; return; }
    size_t stride = 1;
    size_t num_equations = size;
    while (num_equations > 1) {
        size_t new_num = (num_equations + 1) / 2;
		if (stride < size) {
			Real f = c[0] / b[stride];
			b[0] -= f * a[stride];
			d[0] -= f * d[stride];
			c[0] = -f * c[stride];
		}

		size_t last_idx = (new_num - 1) * 2 * stride;
		size_t mid_end = new_num;

		if (last_idx > 0 && last_idx + stride >= size) {
			Real f = a[last_idx] / b[last_idx - stride];
			b[last_idx] -= f * c[last_idx - stride];
			d[last_idx] -= f * d[last_idx - stride];
			a[last_idx] = -f * a[last_idx - stride];
			mid_end = new_num - 1;
		}

#pragma omp parallel for schedule(static)
        for (int i = 1; i < mid_end; i++) {
            const size_t idx = i * 2 * stride;
            Real fl = a[idx] / b[idx - stride];
            Real fr = c[idx] / b[idx + stride];
            b[idx] -= fl * c[idx - stride] + fr * a[idx + stride];
            d[idx] -= fl * d[idx - stride] + fr * d[idx + stride];
            a[idx] = -fl * a[idx - stride];
            c[idx] = -fr * c[idx + stride];
        }
        stride *= 2;
        num_equations = new_num;
    }
    const size_t last_idx = ((size - 1) / stride) * stride;
    x[1 + last_idx] = d[last_idx] / b[last_idx];
	while (stride > 1) {
		stride /= 2;
		size_t step = 2 * stride;
#pragma omp parallel for schedule(static)
		for (size_t i = stride; i < size; i += step) {
			x[1 + i] = (d[i] - a[i] * x[1 + i - stride]) / b[i];
		}
		size_t right_edge = ((size - 1) / step) * step + stride;
		if (right_edge < size && right_edge + stride >= size) {
			x[1 + right_edge] = (d[right_edge] -
				a[right_edge] * x[1 + right_edge - stride]) / b[right_edge];
		}
	}
}

//Обертка-подготовка для метода циклической редукции
void cyclicReductionMethod(
    const std::vector<Real>& x,
    const std::vector<Real>& y,
    SplineWorkspace& ws
) {
    const size_t size = x.size();
    if (size < 3) {
        ws.coef.c[0] = ws.coef.c[size - 1] = Real(0);
        return;
    }
    const size_t N = size - 2;
    for (size_t k = 0; k < N; ++k) {
        const size_t i = k + 1;
        ws.PCR_a[k] = x[i] - x[i - 1];
		ws.PCR_c[k] = x[i + 1] - x[i];;
        ws.PCR_b[k] = Real(2) * (ws.PCR_a[k] + ws.PCR_c[k]);
        ws.PCR_d[k] = Real(6) * ((y[i + 1] - y[i]) / ws.PCR_c[k]
            - (y[i] - y[i - 1]) / ws.PCR_a[k]);
    }
    cyclicReductionParallel(ws.PCR_a, ws.PCR_b, ws.PCR_c, ws.PCR_d, ws.coef.c);
    ws.coef.c[0] = Real(0);
    ws.coef.c[size - 1] = Real(0);
}

void computeBD(
	const std::vector<Real>& x,
	const std::vector<Real>& y,
	SplineWorkspace& ws
) {
	const size_t size = x.size();
	for (size_t i = 1; i < size; ++i) {
		Real h = x[i] - x[i - 1];
		ws.coef.d[i] = (ws.coef.c[i] - ws.coef.c[i - 1]) / h;
		ws.coef.b[i] = (h * ws.coef.c[i] / Real(2)) 
			- (h * h * ws.coef.d[i] / Real(6))
			+ (y[i] - y[i - 1]) / h;
	}
}

//Вычисление коэффициентов сплайна
const splineCoeff& createSpline(
	const std::vector<Real>& x,
	const std::vector<Real>& y,
	SplineWorkspace& ws,
	Methods method
) {
	const size_t size = x.size();
	ws.resize(size);
	if (method == Methods::Thomas) sweepMethod(x, y, ws);
	else cyclicReductionMethod(x, y, ws);
	computeBD(x, y, ws);
	return ws.coef;
}

//Полный цикл вычисления
const splineCoeff& calculateLoop(
	std::vector<Real>& x,
	std::vector<Real>& y,
	SplineWorkspace& ws,
	Methods method
) {
	prepareCoord(x, y);
	return createSpline(x, y, ws, method);
}

//Нахождение интервала
int findInterval(const std::vector<Real>& x, Real xVal) {
	const size_t size = x.size();
	int left = 0, right = size - 1;
	while (left < right - 1) {
		int mid = left + (right - left) / 2;
		if (x[mid] <= xVal) left = mid;
		else right = mid;
	}
	return right;
}

//Вычисление в точке
Real evaluateSpline(
	const splineCoeff& coef,
	const std::vector<Real>& x,
	const std::vector<Real>& y,
	Real xVal
) {
	if (xVal < x.front() || xVal > x.back())
		throw std::out_of_range("Point lies outside the interpolation area");
	int interval = findInterval(x, xVal);
	Real diff = xVal - x[interval];
	return y[interval]
		+ coef.b[interval] * diff
		+ (coef.c[interval] / Real(2)) * diff * diff
		+ (coef.d[interval] / Real(6)) * diff * diff * diff;
}

//Графическое представление сплайна
void exportSplineToCSV(
	const splineCoeff& coef,
	const std::vector<Real>& x,
	const std::vector<Real>& y,
	int numPoints
) {
	if (x.size() < 2)
		throw std::invalid_argument("At least 2 nodes required");
	if (numPoints < 2) numPoints = 2;

	Real step = (x.back() - x.front()) / Real(numPoints - 1);
	std::ofstream file("printSpline.csv");
	if (!file.is_open()) return;

	file << "x;y\n";
	for (int i = 0; i < numPoints; ++i) {
		Real xVal = x.front() + Real(i) * step;
		Real yVal = evaluateSpline(coef, x, y, xVal);
		file << xVal << ";" << yVal << "\n";
	}
	file.close();
}