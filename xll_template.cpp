// xll_template.cpp - Sample xll project.
#include <cmath>
#include <random>
#include <cassert>
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

namespace fms::monte_carlo {

	// Numerically stable average.
	// s_n = (x_1 + ... + x_n)/n
	// n s_n - (n-1) s_{n-1} = x_n
	// s_n = s_{n-1} + (x_n - s_{n-1})/n
	template<class X, class S = std::invoke_result<X>::type>
	inline S average(size_t n, X& x)
	{
		S s = 0;

		for (size_t m = 1; m <= n; ++m) {
			s += (x() - s) / m;
		}

		return s;
	}
}

double monte_carlo_option_value(double f, double s, double k, size_t n = 10000)
{
	std::function<double(double)> payoff;
	if (k < 0) {
		payoff = [k](double x) { return std::max(-k - x, 0.); };
	}
	else {
		payoff = [k](double x) { return std::max(x - k, 0.); };
	}

	std::default_random_engine dre;
	std::normal_distribution<double> X;
	auto p = [f, s, &payoff, &X, &dre]() {
		double F = f * (1 + s * X(dre));

		return payoff(F);
	};

	return fms::average(n, p);
}

// Bachelier Model

double bachelier_put(double f, double sigma, double k, double t)
{
	double x = (k - f) / (sigma * sqrt(t));
	return (k - f) * N(x) + f * sigma * sqrt(t) * n(x); // !!! replace with correct formula
}

int bachelier_put_test()
{
	double f = 100;
	double sigma = 0.2;
	double k = 100;
	double t = 0.25;
	double stdev = f * sigma * sqrt(t);
	int is[] = { 10000 };
	//!!! Write a Monte Carlo simulation that tests the put value.
	// Use 10,000 simulations and test the error is less than 2*stdev.
	for (size_t n : is){
		double v = bachelier_put(f, sigma, k, t);
		double vn = monte_carlo_option_value(f, sigma, k, n);
		assert(fabs(v - vn) <= 2* stdev);
	}
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
	double x = (k - f) / (sigma * sqrt(t));
	return -N(x) + sigma * sqrt(t) * n(x); // !!! replace with correct formula
}

// !!! Test bachelier_put_delta using difference quotients for
int delta_test() {
	for (int f : {90, 100, 110}) {
		for (double sigma : {.1, .2, .3}) {
			for (int k : {90, 100, 110}) {
				for (double t : {.1, .2, .3}) {
					for (double h : {.01, .001, .0001}) {
						double d = bachelier_put_delta(f, sigma, k, t);
						double dt = (bachelier_put(f + h, sigma, k, t) - bachelier_put(f - h, sigma, k, t)) / (2 * h);
						double stdev = f * sigma * sqrt(t);
						assert(fabs(d - dt) <= 2 * stdev);
					}
				}
			}
		}
		return 0;
	}
}

// !!! Implement BACHELIER.PUT.DELTA

// !!! Create a spreadsheet with a graph of put value and put delta as a function of strike k.
// Use f = 100, sigma = 0.2, t = 0.25, and k = 80, 81, ..., 120.


