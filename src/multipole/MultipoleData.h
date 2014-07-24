#ifndef MULTIPOLE_DATA_H
#define MULTIPOLE_DATA_H

#include <config.h>
#include <vector>
#include <map>

#include <complex>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

class MultipoleData{
public:
	class MultipoleModeData{
	public:
		MultipoleModeData(){

		}

		MultipoleModeData(int mode,
						  REAL frequency,
						  const std::vector< REAL > & coefficients)
						 :_mode(mode),
						  _frequency(frequency),
						  _coefficients(coefficients){

		}

		int mode() const{ return this->_mode; };
		REAL frequency() const {return this->_frequency; };

		const int numCoefficients() const{ return int(sqrt(_coefficients.size()/2.0)-1); }
		const std::complex<REAL> coefficient(int m, int l) const{
			return std::complex<REAL>(_coefficients[2*(l*(l+1) + m)],
									  _coefficients[2*(l*(l+1) + m)+1]);
		}

		//
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
		    ar & _mode;
		    ar & _frequency;
		    ar & _coefficients;
		}

	public://But do not use

		int _mode;
		REAL _frequency;
		std::vector< REAL > _coefficients;
	};

	REAL c;
	REAL cx, cy, cz;
	std::map<int, MultipoleModeData> modes;

	std::complex<REAL> estimateModeAt(int mode, REAL x, REAL y, REAL z);

	friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & c;
			ar & cx;
			ar & cy;
			ar & cz;
		    ar & modes;
		}

};

#endif