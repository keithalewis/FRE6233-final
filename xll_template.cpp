// xll_template.cpp - Sample xll project.
#include <cmath>
#include <random>
#include <vector>
#include "xll_template.h"

using namespace std;
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
	double x = (k - f) / (sigma * sqrt(t));
	return (k - f) * N(x) + f * sigma * sqrt(t) * n(x);
}

int bachelier_put_test()
{
	double f = 100;
	double sigma = 0.2;
	double k = 100;
	double t = 0.25;
	double stdev = f * sigma * sqrt(t);
	double p = bachelier_put(f, sigma, k, t);

	random_device rd;
	mt19937 gen(rd());
	normal_distribution<double> normal(0, t);

	//!!! Write a Monte Carlo simulation that tests the put value.
	// Use 10,000 simulations and test the error is less than 2*stdev.

	double res = 0.0, error = 0.0;
	for (int i = 0; i < 10000; i++) {
		double payoff = k - f * (1.0 + sigma * normal(gen));
		payoff = k - payoff >= 0.0 ? payoff : 0.0;
		res = (res * i + payoff) / (i + 1);
		error = (error * i + pow(payoff - p, 2)) / (i + 1);
	}
	error = sqrt(error - pow(res - p, 2));
	double difference = fabs(fabs(res - p) - 2 * stdev);
	// res is the Monte Carlo simulation result of the put value
	// error is the Monte Carlo simulation error
	// difference is the difference between Monte Carlo error and standard error.

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
	return sigma * sqrt(t) * n(x) - N(x);
}

// !!! Test bachelier_put_delta using difference quotients for
// f : {90, 100, 110}
// sigma : {.1, .2, .3}
// k : {90, 100, 110}
// t : {.1, .2, .3}
// h : {.01, .001, .0001}
int bachelier_put_delta_test() {

	vector<double> fs = { 90, 100, 110 };
	vector<double> sigmas = { .1, .2, .3 };
	vector<double> ks = { 90, 100, 110 };
	vector<double> ts = { .1, .2, .3 };
	vector<double> hs = { .01, .001, .0001 };
	vector<double> delta;
	vector<double> error;
	for (auto f : fs) {
		for (auto sigma : sigmas) {
			for (auto k : ks) {
				for (auto t : ts) {
					for (auto h : hs) {
						double calculated_delta = bachelier_put_delta(f, sigma, k, t);
						double payoff_1 = bachelier_put(f, sigma, k, t);
						double payoff_2 = bachelier_put(f + h, sigma, k, t);
						double montecarlo_delta = (payoff_2 - payoff_1) / h;
						delta.push_back(montecarlo_delta);
						error.push_back(fabs(montecarlo_delta - calculated_delta));
					}
				}
			}
		}
	}

	return 0;
}
int bachelier_put_delta_test_ = bachelier_put_delta_test();

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
	.FunctionHelp("Delta of a Bachelier put option.")
);
double WINAPI xll_bachelier_put_delta(double f, double sigma, double k, double t)
{
#pragma XLLEXPORT
	return bachelier_put_delta(f, sigma, k, t);
}

// !!! Create a spreadsheet with a graph of put value and put delta as a function of strike k.
// Use f = 100, sigma = 0.2, t = 0.25, and k = 80, 81, ..., 120.

// Sorry, professor, I cannot get my program run due to some missing files
// I am not sure whether this is caused by the last time I reinstalled the operation system or not. 
// So, I can only plot the figure by Excel. I hope you understand.
