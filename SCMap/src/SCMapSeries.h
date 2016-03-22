/*
 * SCMapSeries.h
 *
 *  Created on: Mar 22, 2016
 *      Author: sandro
 */

#ifndef SCMAPSERIES_H_
#define SCMAPSERIES_H_

#include <vector>
#include <complex>
#include "types.h"

class SCMapSeries
{
public:
	SCMapSeries(unsigned N, std::vector<double> const& angle, std::vector< std::complex<sc::REAL> > const& prevertex, std::complex<sc::REAL> const& C, double rot);
	SCMapSeries(double r, std::vector<double> const& angle, std::vector< std::complex<sc::REAL> > const& prevertex, std::complex<sc::REAL> const& C);
	~SCMapSeries();
	inline unsigned getSize() const { return mM; }
	std::complex<sc::REAL> getPoint(std::complex<sc::REAL> const& p) const;
	sc::REAL getRadius(std::complex<sc::REAL> const& p) const;
	std::vector< std::complex<sc::REAL> > getCoefficients() const { return mvCoefficient; }
	

private:
	unsigned mM;
	std::vector<double> mvAngle;
	std::vector< std::complex<sc::REAL> > mvPrevertex;
	std::complex<sc::REAL> mC;
	double mRotation;
	// access of mDknMatrix[n][k]: n - by prevertex entry, k - by exponent
	std::vector< std::vector< std::complex<sc::REAL> > > mDknMatrix;
	std::vector< std::complex<sc::REAL> > mvCoefficient;

	// Filling d_kn elements in cache.
	void setDknMatrix(unsigned N);

	void generate(unsigned N);
	std::complex<sc::REAL> getDkn(std::complex<sc::REAL> const& pd, double alpha, std::complex<sc::REAL> const& pv, unsigned k);
	void multiplyTwoSums(unsigned N, std::vector< std::complex<sc::REAL> > const& a, std::vector< std::complex<sc::REAL> > const& b, std::vector< std::complex<sc::REAL> >& c);
	std::complex<sc::REAL> getD1(std::complex<sc::REAL> const& z) const;
	std::complex<sc::REAL> getD2(std::complex<sc::REAL> const& z) const;

	void setSize(double r);
	double getMaxRadius();
};
#endif /* SCMAPSERIES_H_ */
