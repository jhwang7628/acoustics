#include <boost/python.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>

#include "../MultipoleData.h"

#include <map>

BOOST_PYTHON_MODULE(_data)
{

	typedef std::map<int, MultipoleData::MultipoleModeData> MultipoleModeDataMap;

	using namespace boost::python;
	class_<MultipoleData::MultipoleModeData>("MultipoleModeData")
		.def("mode", &MultipoleData::MultipoleModeData::mode)
		.def("frequency", &MultipoleData::MultipoleModeData::frequency)
		.def("numCoefficients", &MultipoleData::MultipoleModeData::numCoefficients)
		.def("coefficient", &MultipoleData::MultipoleModeData::coefficient);

	class_<MultipoleModeDataMap>("MultipoleModeDataMap")
		.def(map_indexing_suite<MultipoleModeDataMap>());

	class_<MultipoleData>("MultipoleData")
		.def_readwrite("c", &MultipoleData::c)
		.def_readwrite("cx", &MultipoleData::cx)
		.def_readwrite("cy", &MultipoleData::cy)
		.def_readwrite("cz", &MultipoleData::cz)
		.def_readwrite("modes", &MultipoleData::modes)
		.def("estimateModeAt", &MultipoleData::estimateModeAt);
}