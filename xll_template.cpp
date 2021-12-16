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

double monte_carlo_option_value(double f, double sigma, double k, int n = 10000)
{
	
	std::function<double(double)> payoff;

	payoff = [k](double x) { return std::max(k - x, 0.); };
	
	std::default_random_engine dre;
	std::normal_distribution<double> X;
	auto p = [f, sigma, &payoff, &X, &dre]() {
		double F = f * (1 + sigma * X(dre));
		return payoff(F);
	};
	return fms::average(n, p);
	
}


// Bachelier Model
// p = (k-f)*N(x) + f*sigma*sqrt(t)*n(x)
double bachelier_put(double f, double sigma, double k, double t)
{
	double x = (k - f) / (f * sigma * sqrt(t));
	return ((k-f)*N(x) + f*sigma*sqrt(t)*n(x)); // !!! replace with correct formula
}

int bachelier_put_test()
{
	double f = 100;
	double sigma = 0.2;
	double k = 100;
	double t = 0.25;
	double stdev = f * sigma * sqrt(t);
	int n = 10000;

	double p = bachelier_put(f, sigma, k, t);
	double pn = monte_carlo_option_value(f, sigma, k, n);
	assert(fabs(p-pn) <= 2 * stdev);

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
	return -N(x) + sigma*sqrt(t)*n(x); // !!! replace with correct formula
}

// !!! Test bachelier_put_delta using difference quotients for
// f : {90, 100, 110}
// sigma : {.1, .2, .3}
// k : {90, 100, 110}
// t : {.1, .2, .3}
// h : {.01, .001, .0001}
int bachelier_put_delta_test()
{
	double fs[] = { 90, 100, 110 };
	double sigmas[] = { .1, .2, .3 };
	double ks[] = { 90, 100, 110 };
	double ts[] = { .1, .2, .3 };
	double hs[] = { .01, .001, .0001 };
	int n = 10000;

	for (auto f : fs) {
		for (auto sigma : sigmas) {
			for (auto k : ks) {
				for (auto t : ts) {
					for (auto h : hs) {
						std::function<double(double)> F = [f, sigma, k, t](double x){return bachelier_put(f, sigma, k, t); };
						double df = bachelier_put_delta(f, sigma, k, t);
						double dddF = 1.0;
						double tol = 1.0;
						double x = (k - f) / (f * sigma * sqrt(t));
						fms::derivative_test(F, x, h, df, dddF = 1, tol = 1);
					}
				}
			}
		}
	}
	return 0;
}
// !!! Implement BACHELIER.PUT.DELTA
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


