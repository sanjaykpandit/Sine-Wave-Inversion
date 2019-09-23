#include "alglib/stdafx.h"
#include "alglib/interpolation.h"
#include "alglib/fasttransforms.h"
#include <iostream>

#include "CosFit.h"

using namespace alglib;

real_1d_array linspace(double start, double stop, int n)
{
	real_1d_array l;
	l.setlength(n);
	double dl = (stop - start) / (n-1);
	for (int i = 0; i < n; i++)
	{
		l[i] = start + i * dl;
	}
	return l;
}

int main()
 {
	double start = 0.0;
	double stop = 1.0;
	int n = 100;

	real_1d_array t, x;
	t = linspace(start, stop, n);

	double a = 1.0;
	double f = 3; double omega = 2.0 * pi() * f;
	double p = pi() / 4;
	double d = 0.5;
	x.setlength(n);
	for (int i = 0; i < n; i++)
	{
		x[i] = a * cos(omega * t[i] + p) + d;
	}

	for (int i = 0; i < 201; i++)
	{
		CosFit cs = CosFit(t, x);
		real_1d_array param = cs.get_parameters();

		std::cout << param.tostring(5) << std::endl;
		std::cout << a << " " << f << " " << p << " " << d << std::endl;
	}
}