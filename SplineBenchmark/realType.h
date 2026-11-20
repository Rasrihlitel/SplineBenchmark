#pragma once
#ifndef REALTYPE_H
#define REALTYPE_H

#define USE_FLOAT

#ifdef USE_FLOAT
	using Real = float;
	constexpr Real TEST_EPS = 1e-3f;
	constexpr int  PRINT_PREC = 4;
#else
	using Real = double;
	constexpr Real TEST_EPS = 1e-10;
	constexpr int  PRINT_PREC = 10;
#endif

#endif //REALTYPE_H