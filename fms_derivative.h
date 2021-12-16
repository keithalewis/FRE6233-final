// fms_derivative.h - Test derivatives
#pragma once
#include <functional>
#include <limits>
namespace fms {

	inline constexpr double epsilon = std::numeric_limits<double>::epsilon();

	// symmetric difference quotient
	template<class X, class Y>
	inline auto difference_quotient(const std::function<Y(X)>& f, X dx)
	{
		return [dx, &f](X x) { 
			return (f(x + dx) - f(x - dx)) / (2 * dx); 
		};
	}

	// (f(x + h) - f(x - h))/2h = f'(x) + f'''(x) h^2/3! + O(h^3)
	template<class X, class Y>
	inline bool derivative_test(const std::function<Y(X)>& f, X x, X h, X df, X dddf, X O = 1)
	{
		dddf = std::max(fabs(dddf), 1.);
		auto Df = difference_quotient(f, h);

		double lhs = Df(x) - df;
		double rhs = dddf * h * h / 6;
		
		bool b = fabs(lhs) <= O * rhs;
		if (!b) {
			O = fabs(lhs)/fabs(rhs); // set breakpoint to check tolerance
		}

		return b;
	}
}
