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

// Bachelier Model

double bachelier_put(double f, double sigma, double k, double t)
{
	double x = (k - f) / (f * sigma * sqrt(t));
	return (k-f)*N(x)+f*sigma*sqrt(t)*n(x); // !!! replace with correct formula

}

int bachelier_put_test()
{
	double f = 100;
	double sigma = 0.2;
	double k = 100;
	double t = 0.25;
	double stdev = f * sigma * sqrt(t);
	const double TIMES = 10000;
	std::default_random_engine gen;
	std::normal_distribution<double> B(0, 1);
	double sum = 0;
	for (int i = 0; i < TIMES; i++)
	{
		double ran = B(gen);
		double F = f * (1 + sigma * sqrt(t) * ran);
		double P = 0;
		if (F < k) P = k - F;
		else P = 0;
		sum = sum + P;
	}
	double p = bachelier_put(f, sigma, k, t);
	double p_mont = sum / TIMES;
	assert(fabs(p_mont-p)<=2*stdev/sqrt(TIMES));

	//!!! Write a Monte Carlo simulation that tests the put value.
	// Use 10,000 simulations and test the error is less than 2*stdev.

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
	return sigma*sqrt(t)*n(x)-N(x); // !!! replace with correct formula
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
	.FunctionHelp("Forward delta of a Bachelier put option.")
);
double WINAPI xll_bachelier_put_delta(double f, double sigma, double k, double t)
{
#pragma XLLEXPORT
	return bachelier_put_delta(f, sigma, k, t);
}

// !!! Test bachelier_put_delta using difference quotients for
// f : {90, 100, 110}
// sigma : {.1, .2, .3}
// k : {90, 100, 110}
// t : {.1, .2, .3}
// h : {.01, .001, .0001}

// !!! Implement BACHELIER.PUT.DELTA
int bachelier_put_delta_test()
{
	double fs[3] = { 90,100,110 };
	double sigmas[3] = { 0.1,0.2,0.3 };
	double ks[3] = { 90,100,110 };
	double ts[3] = { 0.1,0.2,0.3 };
	double hs[3] = { 0.01,0.001,0.0001 };
	const double tol = 0.001;
	for (double h : hs)
	{
		for (double f : fs)
		{
			for (double k: ks)
			{
				for (double t : ts)
				{
					for (double sigma : sigmas)
					{
						double delta = bachelier_put_delta(f, sigma, k, t);
						double delta_minus = bachelier_put_delta(f-h, sigma, k, t);
						double delta_plus = bachelier_put_delta(f + h, sigma, k, t);
						double delta_test = (delta_plus - delta_minus) / (2 * h);
						assert(fabs(delta - delta_test) <= tol * fabs(delta));

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


