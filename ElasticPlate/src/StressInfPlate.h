#pragma once

#include <array>
#include <vector>
#include <complex>
#include <string>
#include "SCMap/src/types.h"

struct stress{
	double x;
	double y;
	double xy;
	double r;
	double f;
	double rf;
	double fi;
	std::complex<sc::REAL> A;
	std::complex<sc::REAL> B;
};

class StressInfPlate
{
public:
	StressInfPlate(std::vector<std::complex<sc::REAL> > const& c, double sigma_1, double sigma_2, double alpha);
	~StressInfPlate();

	stress getValue(std::complex<sc::REAL> const& p) const;


private:
	unsigned mN;
	std::complex<sc::REAL> mGamma1;
	std::complex<sc::REAL> mGamma2;

	// Series coefficents
	std::vector<std::complex<sc::REAL> > mCCoeff;
	std::vector<std::complex<sc::REAL> > mPCoeff;
	std::vector<std::complex<sc::REAL> > mQCoeff; // This starts from 0
	std::vector<std::complex<sc::REAL> > mACoeff;
	std::vector<std::complex<sc::REAL> > mKCoeff;

	void generateCoefficients();
	void setCCoeff();
	void setPCoeff();
	void setQCoeff();
	void setKCoeff();
	void solveACoeff();
	void setASystem(std::vector<std::complex<sc::REAL> > const& a, std::vector<std::complex<sc::REAL> > & b);
	sc::REAL systemError(std::vector<std::complex<sc::REAL> > const& a, std::vector<std::complex<sc::REAL> > const& b);
	std::complex<sc::REAL> getOmega(std::complex<sc::REAL> const& z ) const;
	std::complex<sc::REAL> getDOmega(std::complex<sc::REAL> const& z ) const;
	std::complex<sc::REAL> getDDOmega(std::complex<sc::REAL> const& z ) const;
	std::complex<sc::REAL> getVarphi(std::complex<sc::REAL> const& z) const;
	std::complex<sc::REAL> getDVarphi(std::complex<sc::REAL> const& z) const;
	std::complex<sc::REAL> getDDVarphi(std::complex<sc::REAL> const& z) const;
	std::complex<sc::REAL> getQ(std::complex<sc::REAL> const& z) const;
	std::complex<sc::REAL> getDQ(std::complex<sc::REAL> const& z) const;
	std::complex<sc::REAL> getK(std::complex<sc::REAL> const& z) const;
	std::complex<sc::REAL> getDK(std::complex<sc::REAL> const& z) const;

	std::vector<std::complex<sc::REAL> > zeros(unsigned N);
	void printVec(std::vector<std::complex<sc::REAL> > const& v, std::string const& name) const;
};

