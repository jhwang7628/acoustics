#include <boost/python.hpp>

#include "../MultipoleMath.h"

BOOST_PYTHON_MODULE(_math)
{
	using namespace boost::python;
	def("spherical_harmonics", Multipole::spherical_harmonics);
	def("hankel_1st", Multipole::hankel_1st);
	def("hankel_2nd", Multipole::hankel_2nd);
	def("spherical_bessel", Multipole::spherical_bessel);
	def("regular_basis", Multipole::regular_basis);
	def("regular_basis_dir_deriv", Multipole::regular_basis_dir_deriv);
}