#pragma once
#include<pybind11/numpy.h>
#include "alglib/interpolation.h"
#include "alglib/fasttransforms.h"
#include "alglib/statistics.h"
#include "alglib/optimization.h"
using namespace alglib;

namespace py = pybind11;

class CosFit
{
	real_1d_array t;
	real_1d_array x;

	real_1d_array get_approximate_parameters();
	double get_phase_shift(double amplitude, double frequency, double dc, int n_int);
	real_1d_array cyclicShift(real_1d_array input, int shift);

	real_1d_array linspace(double start, double stop, int n);
	
	real_1d_array abs_array(complex_1d_array);
	real_1d_array angle_array(complex_1d_array);
	real_1d_array real_array(complex_1d_array);
	real_1d_array imag_array(complex_1d_array);
	
	double max(real_1d_array x);
	double min(real_1d_array x);

	double average_sum(real_1d_array x);

	real_1d_array slice(real_1d_array x, int start, int stop);

	real_1d_array interp(real_1d_array t_int, real_1d_array x, real_1d_array y);

	void cout_array(real_1d_array x, std::string start);

public:
	CosFit(real_1d_array t, real_1d_array x);
	CosFit(py::array_t<double> t, py::array_t<double> x);
	static real_1d_array cos_array(real_1d_array t, double amplitude, double frequency, double phase, double dc);
	real_1d_array get_parameter_c();
	py::array_t<double> get_parameter();
	real_1d_array get_series_c(real_1d_array t, real_1d_array p);
	py::array_t<double> get_series(py::array_t<double> t, py::array_t<double> p);
	double get_error(py::array_t<double> p);

	void cost_function_gradient(const real_1d_array& p, double& cost_function, real_1d_array& gradient, void* ptr);
};

