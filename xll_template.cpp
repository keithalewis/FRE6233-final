// xll_template.cpp - Sample xll project.
#include <cassert>
#include <cmath>
#include <random>
#include "xll_template.h"

using namespace xll;

// sqrt(2 pi)
constexpr double M_SQRT2PI = 2.50662827463100050240;
// sqrt(2)
constexpr double M_SQRT2 = 1.41421356237309504880;


double n(double x)
{
	return exp(-x * x / 2) / M_SQRT2PI;
}

double N(double x)
{
	return (1 + erf(x / M_SQRT2)) / 2;
}

// Bachelier Model

double bachelier_put(double f, double sigma, double k, double t)
{
	double x = (k - f) / (f * sigma * sqrt(t));

	return (k - f)*N(x) + f*sigma*sqrt(t)*n(x);
}

int bachelier_put_test()
{
	double f = 100;
	double sigma = 0.2;
	double k = 100;
	double t = 0.25;
	double stdev = f * sigma * sqrt(t);

	//!!! Write a Monte Carlo simulation that tests the put value.
	// Use 10,000 simulations and test the error is less than 2*stdev.
	std::default_random_engine dre;
	std::normal_distribution<double> Z(0, 1);

	auto X = [f, sigma, k, t, &Z, &dre]() {
		double srt = sigma * sqrt(t);
		double F = f * (1 + srt * Z(dre));

		return std::max(k - F, 0.);
	};

	double p = fms::average(10000, X);
	double p_ = bachelier_put(f, sigma, k, t);

	assert(fabs(p - p_) < 2 * stdev);

	return 0;
}
int bachelier_put_test_ = bachelier_put_test();

AddIn xai_bachelier_put(
	Function(XLL_DOUBLE, "xll_bachelier_put", "BACHELIER.PUT")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expriation."),
		})
	.Category(CATEGORY)
	.FunctionHelp("Forward value of a Bachelier put option.")
);
double WINAPI xll_bachelier_put(double f, double sigma, double k, double t)
{
#pragma XLLEXPORT
	return bachelier_put(f, sigma, k, t);
}

double bachelier_put_delta(double f, double sigma, double k, double t)
{
	double srt = sigma * sqrt(t);
	double x = (k - f) / (f * srt);

	return -N(x) + srt*n(x); 
}

bool test_delta(double f, double sigma, double k, double t, double h)
{
	std::function<double(double)> p = [sigma, k, t](double x) {
		return bachelier_put(x, sigma, k, t);
	};
	auto dp = bachelier_put_delta(f, sigma, k, t);

	return fms::derivative_test(p, f, h, dp, 1.);
}

// !!! Test bachelier_put_delta using difference quotients for
int bachelier_put_delta_test()
{
	for (double f : {90, 100, 110})
		for (double sigma : {.1, .2, .3})
			for (double k : {90, 100, 110})
				for (double t : {.1, .2, .3})
					for (double h : {.01, .001, .0001})
						assert(test_delta(f, sigma, k, t, h));

	return 0;
}
int bachelier_put_delta_test_ = bachelier_put_delta_test();


AddIn xai_bachelier_put_delta(
	Function(XLL_DOUBLE, "xll_bachelier_put_delta", "BACHELIER.PUT.DELTA")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expriation."),
		})
		.Category(CATEGORY)
	.FunctionHelp("Forward value of a Bachelier put option.")
);
double WINAPI xll_bachelier_put_delta(double f, double sigma, double k, double t)
{
#pragma XLLEXPORT
	return bachelier_put_delta(f, sigma, k, t);
}

// !!! Create a spreadsheet with a graph of put value and put delta as a function of strike k.
// Use f = 100, sigma = 0.2, t = 0.25, and k = 80, 81, ..., 120.


