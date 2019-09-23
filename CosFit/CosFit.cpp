#include "CosFit.h"

CosFit::CosFit(real_1d_array t, real_1d_array x)
{
	this->t = t;
	this->x = x;
}

CosFit::CosFit(py::array_t<double> t, py::array_t<double> x)
{
	this->t.setcontent(t.size(), &t.at(0));
	this->x.setcontent(x.size(), &x.at(0));
}

double CosFit::get_error(py::array_t<double> p)
{
	double am = p.at(0);
	double fr = p.at(1);
	double ph = p.at(2);
	double dc = p.at(3);

	int n = this->t.length();
	double omega = 2.0 * pi() * fr;

	double cost_function = 0;
	double cf, e;
	for (int i = 0; i < n; i++)
	{
		cf = am * cos(omega * t[i] + ph) + dc;
		e = cf - x[i];
		cost_function += pow(e, 2);
	}

	return cost_function;
}

void proxy_function(const real_1d_array& p, double& cost_function, real_1d_array& gradient, void* ptr)
{
	((CosFit*)ptr)->cost_function_gradient(p, cost_function, gradient, ptr);
}

real_1d_array CosFit::get_parameter_c()
{
	real_1d_array p = get_approximate_parameters();
	real_1d_array s = "[1, 1, 1, 1]";
	double epsg = 0;
	double epsf = 0;
	double epsx = 1e-9;
	ae_int_t maxits = 0;

	minlbfgsstate state;
	minlbfgscreate(1, p, state);
	minlbfgssetcond(state, epsg, epsf, epsx, maxits);
	minlbfgssetscale(state, s);

	minlbfgsreport rep;
	minlbfgsoptimize(state, &proxy_function, NULL, this);
	minlbfgsresults(state, p, rep);

	if (p[0] < 0.0)
	{
		p[0] = -p[0];
		if (p[2] > 0)
			p[2] -= pi();
		else
			p[2] += pi();
	}

	return p;
}

py::array_t<double> CosFit::get_parameter()
{
	real_1d_array p = this->get_parameter_c();
	return py::array_t<double>(p.length(), &p[0]);
}
py::array_t<double> CosFit::get_series(py::array_t<double> t, py::array_t<double> p)
{
	real_1d_array t_alglib, p_alglib;
	t_alglib.setcontent(t.size(), &t.at(0));
	p_alglib.setcontent(p.size(), &p.at(0));
	
	real_1d_array s = get_series_c(t_alglib, p_alglib);
	
	return py::array_t<double>(s.length(), &s[0]);
}

real_1d_array CosFit::get_series_c(real_1d_array t, real_1d_array p)
{
	return cos_array(t, p[0], p[1], p[2], p[3]);
}

real_1d_array CosFit::linspace(double start, double stop, int n)
{
	real_1d_array l;
	l.setlength(n);
	double dl = (stop - start) / (n - 1);
	for (int i = 0; i < n; i++)
	{
		l[i] = start + i * dl;
	}
	return l;
}

real_1d_array CosFit::get_approximate_parameters()
{
	complex_1d_array fft_x;
	fftr1d(this->x, fft_x);

	double dt = this->t[1] - this->t[0];
	int n = this->x.length();
	real_1d_array f = linspace(0, 1.0 / (2.0 * dt), int(n / 2));

	real_1d_array abs_fft_x = abs_array(fft_x);
	int max_index = -1;
	double max_value = 0;
	for (int i = 1; i < n; i++)
	{
		if (abs_fft_x[i] > max_value)
		{
			max_index = i;
			max_value = abs_fft_x[i];
		}
	}

	double frequency = f[max_index];
	double amplitude = max_value;
	double dc = average_sum(this->x);
	double phase = get_phase_shift(amplitude, frequency, dc, 100);

	real_1d_array p;
	p.setlength(4);
	p[0] = amplitude;
	p[1] = frequency;
	p[2] = phase;
	p[3] = dc;

	return p;
}

real_1d_array CosFit::abs_array(complex_1d_array c)
{
	int n = c.length();
	real_1d_array a;
	a.setlength(n);
	for (int i = 0; i < n; i++)
	{
		a[i] = abscomplex(c[i]) / n * 2.0;
	}
	return a;
}
real_1d_array CosFit::angle_array(complex_1d_array c)
{
	int n = c.length();
	real_1d_array a;
	a.setlength(n);
	for (int i = 0; i < n; i++)
	{
		a[i] = atan(c[i].y/c[i].x);
	}
	return a;
}
real_1d_array CosFit::real_array(complex_1d_array c)
{
	int n = c.length();
	real_1d_array a;
	a.setlength(n);
	for (int i = 0; i < n; i++)
	{
		a[i] = c[i].x;
	}
	return a;
}
real_1d_array CosFit::imag_array(complex_1d_array c)
{
	int n = c.length();
	real_1d_array a;
	a.setlength(n);
	for (int i = 0; i < n; i++)
	{
		a[i] = c[i].y;
	}
	return a;
}

real_1d_array CosFit::cos_array(real_1d_array t, double amplitude, double frequency, double phase, double dc)
{
	int n = t.length();
	real_1d_array x;
	x.setlength(n);
	double omega = 2.0 * pi() * frequency;
	for (int i = 0; i < n; i++)
	{
		x[i] = amplitude * cos(omega * t[i] + phase) + dc;
	}
	return x;
}

double CosFit::max(real_1d_array x)
{
	double m = 0;
	for (int i = 0; i < x.length(); i++)
	{
		if (m < x[i])
			m = x[i];
	}
	return m;
}
double CosFit::min(real_1d_array x)
{
	double m = 0;
	for (int i = 0; i < x.length(); i++)
	{
		if (m > x[i])
			m = x[i];
	}
	return m;
}

double CosFit::average_sum(real_1d_array x)
{
	double m = 0;
	int n = x.length();
	for (int i = 0; i < x.length(); i++)
	{
		m += x[i];
	}
	return m/n;
}

real_1d_array CosFit::slice(real_1d_array x, int start, int stop)
{
	real_1d_array y;
	y.setlength(stop - start + 1);
	int j = 0;
	for (int i = start; i <= stop; i++)
	{
		y[j] = x[i];
		j++;
	}
	return y;
}

real_1d_array CosFit::interp(real_1d_array t_int, real_1d_array x, real_1d_array y)
{
	int n = t_int.length();
	real_1d_array y_int;
	y_int.setlength(n);

	spline1dinterpolant s;
	spline1dbuildcubic(x, y, s);

	for (int i = 0; i < n; i++)
	{
		y_int[i] = spline1dcalc(s, t_int[i]);
	}
	return y_int;
}

void CosFit::cout_array(real_1d_array x, std::string start)
{
	std::cout << start << " = np.array(";
	std::cout << x.tostring(3) << ")" << std::endl;
}

double CosFit::get_phase_shift(double amplitude, double frequency, double dc, int n_int)
{
	double omega = 2.0 * pi() * frequency;
	int n = this->x.length();
	real_1d_array y = cos_array(this->t, amplitude, frequency, 0.0, dc);

	double dt = this->t[1] - this->t[0];
	double max_delay = pi() / omega;
	int num_steps = std::round(max_delay / dt);

	real_1d_array tsa, crr;
	tsa.setlength(2 * num_steps + 1);
	crr.setlength(2 * num_steps + 1);

	int j = 0;
	for (int ts = -num_steps; ts <= num_steps; ts++)
	{
		tsa[j] = ts;
		crr[j] = pearsoncorr2(this->x, this->cyclicShift(y, ts));
		j++;
	}

	real_1d_array tsa_int = linspace(-num_steps, num_steps, n_int);
	real_1d_array crr_int = interp(tsa_int, tsa, crr);

	int max_index = 0;
	double max_value = 0;
	for (int i = 0; i < n_int; i++)
	{
		if (max_value < crr_int[i])
		{
			max_index = i;
			max_value = crr_int[i];
		}
	}


	double time_delay = tsa_int[max_index];
	double phase_shift = -time_delay * dt * omega;
	
	return phase_shift;
}

real_1d_array CosFit::cyclicShift(real_1d_array input, int shift)
{
	real_1d_array output;
	int n = input.length();
	output.setlength(n);
	if (shift >= 0)
	{
		for (int i = 0; i < n; i++)
		{
			output[i] = input[(i + (n - shift)) % (n)];
		}
	}
	else
	{
		shift = -shift;
		for (int i = 0; i < n; i++)
		{
			output[i] = input[(i + shift) % n];
		}
	}
	return output;
}

void CosFit::cost_function_gradient(const real_1d_array& p, double& cost_function, real_1d_array& gradient, void* ptr)
{
	double am = p[0];
	double fr = p[1];
	double ph = p[2];
	double dc = p[3];

	int n = this->t.length();
	double omega = 2.0 * pi() * fr;
	
	cost_function = 0;
	
	gradient[0] = 0;
	gradient[1] = 0;
	gradient[2] = 0;
	gradient[3] = 0;

	double cf, e, sf;
	for (int i = 0; i < n; i++)
	{
		cf = am * cos(omega * t[i] + ph) + dc;
		e = cf - x[i];
		cost_function += pow(e, 2);
		gradient[0] += 2 * e * cos(omega * t[i] + ph);
		sf = am * sin(omega * t[i] + ph);
		gradient[1] += -4 * e * sf * pi() * t[i];
		gradient[2] += -2 * e * sf;
		gradient[3] += 2 * e;
	}
}