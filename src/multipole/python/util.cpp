#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "../MultipoleData.h"
#include "../MultipoleUtil.h"

#include <map>

namespace gambi{
	void saveMultipoleDataToFile(const std::string & filename, const MultipoleData &  data){
		Multipole::saveToFile(filename, data);
	}

	void saveMultipoleModeDataToFile(const std::string & filename, const MultipoleData::MultipoleModeData & data){
		Multipole::saveToFile(filename, data);
	}
}

BOOST_PYTHON_MODULE(_multipole_util)
{
	using namespace boost::python;
	typedef std::vector<std::string> StringVector;

	class_<StringVector>("StringVector")
		.def(vector_indexing_suite<StringVector>());

	def("loadMultipoleDataFromFile", &Multipole::loadMultipoleDataFromFile);
	def("saveMultipoleDataToFile", &gambi::saveMultipoleDataToFile);

	def("loadMultipoleModeDataFromFile", &Multipole::loadMultipoleModeDataFromFile);
	def("saveMultipoleModeDataToFile", &gambi::saveMultipoleModeDataToFile);

	def("mergeMultipoleModeDataFiles", &Multipole::mergeMultipoleModeDataFiles);
	def("mergeMultipoleDataFiles", &Multipole::mergeMultipoleDataFile);
}