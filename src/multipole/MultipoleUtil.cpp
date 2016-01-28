#include "MultipoleUtil.h"

#include <iostream>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace Multipole{
	MultipoleData loadMultipoleDataFromFile(const std::string & filename){
		std::ifstream ifs(filename.c_str());
		boost::archive::binary_iarchive ia(ifs);
		MultipoleData data;
		ia >> data;
		return data;
	}

	void saveToFile(const std::string & filename, const MultipoleData &  data){
		std::ofstream ofs(filename.c_str());
		boost::archive::binary_oarchive oa(ofs);
		oa << data;
	}

	MultipoleData::MultipoleModeData loadMultipoleModeDataFromFile(const std::string & filename){
		std::ifstream ifs(filename.c_str());
		boost::archive::binary_iarchive ia(ifs);
		MultipoleData::MultipoleModeData data;
		ia >> data;
		return data;
	}

	void saveToFile(const std::string & filename, const MultipoleData::MultipoleModeData & data){
		std::ofstream ofs(filename.c_str());
		boost::archive::binary_oarchive oa(ofs);
		oa << data;		
	}

	//Assume that each file contains a MultipoleModeData
	MultipoleData mergeMultipoleModeDataFiles(const std::vector<std::string> & filenames, REAL c, REAL cx, REAL cy, REAL cz){
		MultipoleData merge;
		merge.c = c;
		merge.cx = cx;
		merge.cy = cy;
		merge.cz = cz;

		for(int i = 0; i < filenames.size(); i++){
			MultipoleData::MultipoleModeData mode = loadMultipoleModeDataFromFile(filenames[i]);
			std::cout << "Inserting mode: " << mode.mode() << std::endl;
			merge.modes.insert(std::pair<int, MultipoleData::MultipoleModeData>(mode.mode(), mode));
		}
		return merge;
	}
	
	//Assume that each file contains a MultipoleData
	MultipoleData mergeMultipoleDataFile(const std::vector<std::string> & filenames){
		MultipoleData merge;
		for(int i = 0; i < filenames.size(); i++){
			MultipoleData data = loadMultipoleDataFromFile(filenames[i]);
			if(i == 0){
				merge.c = data.c;
				merge.cx = data.cx;
				merge.cy = data.cy;
				merge.cz = data.cz;
			}
			merge.modes.insert(data.modes.begin(), data.modes.end());
		}
		return merge;
	}
}