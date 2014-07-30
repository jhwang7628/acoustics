#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "../ModeData.h"

BOOST_PYTHON_MODULE(_mode_data)
{
	using namespace boost::python;
	typedef std::vector<REAL> RealVector;
	typedef std::vector<std::vector<REAL> > RealVectorVector;

	class_<RealVector>("RealVector")
		.def(vector_indexing_suite<RealVector>());

	class_<RealVectorVector>("RealVectorVector")
		.def(vector_indexing_suite<RealVectorVector>());

	class_<ModeData>("ModeData")
		.def_readwrite("omegaSquared", &ModeData::_omegaSquared)
		.def_readwrite("modes", &ModeData::_modes)
		.def("numModes", &ModeData::numModes)
		.def("numDOF", &ModeData::numDOF)
		.def("read", &ModeData::read)
		.def("write", &ModeData::write);
}