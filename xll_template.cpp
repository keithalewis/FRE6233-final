// xll_template.cpp - Sample xll project.
#include <cmath>
#include <cassert>
#include <random>
#include "xll_template.h"
#include "fms_monte_carlo.h"

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
	return ((k - f) * N(x) + f * sigma * sqrt(t) * n(x)); // !!! replace with correct formula
}

int bachelier_put_test()
{
	double f = 100;
	double sigma = 0.2;
	double k = 100;
	double t = 0.25;
	double stdev = f * sigma * sqrt(t);
	int n = 10000;

	//!!! Write a Monte Carlo simulation that tests the put value.
	// Use 10,000 simulations and test the error is less than 2*stdev.
	std::function<double(double)> payoff;
	payoff = [k](double x) { return std::max(k - x, 0.); };

	std::default_random_engine dre;
	std::normal_distribution<double> X;
	auto p = [f, sigma, t, &payoff, &X, &dre]() {
		double F = f * (1 + sigma * sqrt(t) * X(dre));
		return payoff(F);
	};

	double vn = fms::average(n, p);
	double v = bachelier_put(f, sigma, k, t);
	assert(fabs(v - vn) <= stdev * 2 / sqrt(n));
	// test assert(2 <= 1);
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
	double x = (k - f) / (f * sigma * sqrt(t));
	return (- N(x) + sigma * sqrt(t) * n(x)); // !!! replace with correct formula
}

// !!! Test bachelier_put_delta using difference quotients for
// f : {90, 100, 110}
// sigma : {.1, .2, .3}
// k : {90, 100, 110}
// t : {.1, .2, .3}
// h : {.01, .001, .0001}

int bachelier_delta_test()
{
	double fs[] = { 90, 100, 110 };
	double sigmas[] = { .1, .2, .3 };
	double ks[] = { 90, 100, 110 };
	double ts[] = { .1, .2, .3 };
	double hs[] = { .01, .001, .0001 };
	double tol = 1;

	for (double f : fs) {
		for (double sigma : sigmas) {
			for (double k : ks) {
				for (double t : ts) {
					for (double h : hs) {
						double delta_t = (bachelier_put(f+h, sigma, k, t) - bachelier_put(f-h, sigma, k, t))/(2*h);
						double delta = bachelier_put_delta(f, sigma, k, t);
						assert(fabs(delta - delta_t) <= tol);
						//assert(2 <= 1); // test 
					}
				}
			}
		}
	}

	return 0;
}
int bachelier_delta_test_ = bachelier_delta_test();

// !!! Implement BACHELIER.PUT.DELTA
AddIn xai_bachelier_put_delta(
	Function(XLL_DOUBLE, "xll_bachelier_delta", "BACHELIER.PUT.DELTA")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expriation."),
		})
		.Category(CATEGORY)
	.FunctionHelp("Delta of a Bachelier put option.")
);
double WINAPI xll_bachelier_delta(double f, double sigma, double k, double t)
{
#pragma XLLEXPORT
	return bachelier_put_delta(f, sigma, k, t);
}

// !!! Create a spreadsheet with a graph of put value and put delta as a function of strike k.
// Use f = 100, sigma = 0.2, t = 0.25, and k = 80, 81, ..., 120.


