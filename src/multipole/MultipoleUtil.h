#ifndef MULTIPOLE_UTIL_H
#define MULTIPOLE_UTIL_H

#include "MultipoleData.h"
#include <string>

namespace Multipole{
	MultipoleData loadMultipoleDataFromFile(const std::string & filename);
	void saveToFile(const std::string & filename, const MultipoleData &  data);

	MultipoleData::MultipoleModeData loadMultipoleModeDataFromFile(const std::string & filename);
	void saveToFile(const std::string & filename, const MultipoleData::MultipoleModeData & data);

	//Assume that each file contains a MultipoleModeData
	MultipoleData mergeMultipoleModeDataFiles(const std::vector<std::string> & filenames, REAL c, REAL cx, REAL cy, REAL cz);
	
	//Assume that each file contains a MultipoleData
	MultipoleData mergeMultipoleDataFile(const std::vector<std::string> & filenames);
}

#endif