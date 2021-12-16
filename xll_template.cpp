// xll_template.cpp - Sample xll project.
#include <cmath>
#include "xll_template.h"
#include <random>
#include <cassert>
#include <iostream>
#include <string.h>
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

double normal_pdf(double x)
{
	return exp(-x * x / 2) / (M_SQRT2PI);
}

double normal_cdf(double x)
{
	return (1 + erf(x / M_SQRT2)) / 2;
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
	double srt = sigma * sqrt(t);
	double d = (k - f) / (f * srt);
	return (k - f) * normal_cdf(d) + srt * f * normal_pdf(d);
}

int bachelier_put_test()
{
	//!!! Write a Monte Carlo simulation that tests the put value.
	// Use 10,000 simulations and test the error is less than 2*stdev.
	double f = 100;
	double sigma = 0.2;
	double k = 100;
	double t = 0.25;
	double stdev = f * sigma * sqrt(t);
	double sd = 2;
	size_t n = 10000;

	double v = bachelier_put(f, sigma, k, t);
	double vn = monte_carlo_option_value(f, sigma, k, n);
	//debug print
	_RPT1(0, "%f\n%f\n", v, vn);

	assert(fabs(v - vn) <= stdev * sd);

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






// bachelier_put_delta

double bachelier_put_delta(double f, double sigma, double k, double t)
{
	double srt = sigma * sqrt(t);
	double d = (k - f) / (f * srt);
	return -normal_cdf(d) + srt * normal_pdf(d);
}

AddIn xai_bachelier_put_delta(
	Function(XLL_DOUBLE, "xll_bachelier_put_delta", "BACHELIER.PUT_DELTA")
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

int bachelier_put_delta_test()
{
	// !!! Test bachelier_put_delta using difference quotients for
	// f : {90, 100, 110}
	// sigma : {.1, .2, .3}
	// k : {90, 100, 110}
	// t : {.1, .2, .3}
	// h : {.01, .001, .0001}

	double fs[] = { 90, 100, 110 };
	double sigmas[] = { .1, .2, .3 };
	double ks[] = { 90, 100, 110 };
	double ts[] = { .1, .2, .3 };
	double hs[] = { .01, .001, .0001 };
	for (double f : fs) {
		for (double k : ks) {
			for (double t : ts) {
				for (double sigma : sigmas) {
					for (double h : hs) {
						double stdev = f * sigma * sqrt(t);
						double sd = 2;
						double v = bachelier_put_delta(f, sigma, k, t);
						double vn = (monte_carlo_option_value(f + h, sigma, k) - monte_carlo_option_value(f, sigma, k)) / h;
						assert(fabs(v - vn) <= stdev * sd);
					}
				}
			}
		}
	}
	return 0;
}
int bachelier_put_delta_test_ = bachelier_put_delta_test();



// !!! Create a spreadsheet with a graph of put value and put delta as a function of strike k.
// Use f = 100, sigma = 0.2, t = 0.25, and k = 80, 81, ..., 120.


