#include <pybind11/pybind11.h>
#include "CosFit.h"

namespace py = pybind11;


PYBIND11_MODULE(CosFit, m)
{
	py::class_<CosFit>(m, "CosFit")
		.def(py::init<py::array_t<double>, py::array_t<double>>())
		.def("get_parameter", &CosFit::get_parameter)
		.def("get_series", &CosFit::get_series)
		.def("get_error", &CosFit::get_error);
}