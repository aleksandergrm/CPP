// main.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <stdexcept>

#define _USE_MATH_DEFINES
#include <cmath>

#include "src/types.h"
#include "src/SCMapSeries.h"


void 
setTestPolygon(std::vector<double>& angle, std::vector< std::complex<sc::REAL> >& prevertex)
{
	angle.clear();
	prevertex.clear();
	unsigned shape = 2;

	if (shape == 1) {
		/* Square hole */
		angle.clear();
		prevertex.clear();

		sc::REAL c = 0.7071067811865475;
		for (unsigned i = 0; i < 4; i++) {
			angle.push_back(0.5);
		}
		prevertex.push_back(std::complex<sc::REAL>(c, c));
		prevertex.push_back(std::complex<sc::REAL>(-c, c));
		prevertex.push_back(std::complex<sc::REAL>(-c, -c));
		prevertex.push_back(std::complex<sc::REAL>(c, -c));
	}

	if (shape == 2) {
		/* L-Shape polygon */
		// Here we test on polygon [0, i, -1+i, -1-i, 1-i, 1]
		angle.push_back(0.5);
		angle.push_back(0.5);
		angle.push_back(0.5);
		angle.push_back(0.5);
		angle.push_back(0.5);
		angle.push_back(1.5);

		// Calculated using SC Toolbox in matlab
		prevertex.push_back(std::complex<sc::REAL>(1.0, 0.0));
		prevertex.push_back(std::complex<sc::REAL>(0.497026568112923, -0.867735322889353));
		prevertex.push_back(std::complex<sc::REAL>(-0.909257018651389, -0.416235118692774));
		prevertex.push_back(std::complex<sc::REAL>(-0.332009183382948, 0.943276153705471));
		prevertex.push_back(std::complex<sc::REAL>(0.653496652729380, 0.756929405474180));
		prevertex.push_back(std::complex<sc::REAL>(0.909257018818266, 0.416235118328236));
	}

	/* Test Milan Shape*/
	if (shape == 3) {
		angle.push_back(0.40669);
		angle.push_back(0.57213);
		angle.push_back(0.68827);
		angle.push_back(0.59756);
		angle.push_back(0.73535);

		prevertex.push_back(std::complex<sc::REAL>(1.0, 0.0));
		prevertex.push_back(std::complex<sc::REAL>(-0.06896, -0.99762));
		prevertex.push_back(std::complex<sc::REAL>(-0.93665, -0.350261));
		prevertex.push_back(std::complex<sc::REAL>(-0.73586, 0.67714));
		prevertex.push_back(std::complex<sc::REAL>(0.09187, 0.99577));
	}
}

void
checkRadius(std::vector<double> const& angle, std::vector< std::complex<sc::REAL> > const& prevertex)
{
	std::complex<sc::REAL> C(1.0, 1.0);
	//std::vector<unsigned> N = { 1000, 2000, 3000, 5000, 7000, 10000, 20000, 50000, 100000, 150000};
	std::vector<unsigned> N = { 1000, 2000, 3000, 5000};
	SCMapSeries* scMap;

	std::ofstream os("radius.txt");

	for (auto n : N)
	{
		std::cout << n << std::endl;
		scMap = new SCMapSeries(n, angle, prevertex, C);

		int i = 1;
		for (auto z : prevertex)
		{
			os << i++ << ' ' << scMap->getRadius(z) << ' ';
		}

		os << std::endl;

		delete(scMap);
	}
	os.close();
}

void
checkRadius(SCMapSeries const& scMap, std::vector< std::complex<sc::REAL> > const& prevertex)
{
	std::complex<sc::REAL> C(1.0, 1.0);

	std::ofstream os("radius.txt");

	int i = 0;
	for (auto z : prevertex)
	{
		os << i++ << ' ' << scMap.getRadius(z) << ' ';
	}

	os.close();
}

void
writeContour(unsigned N, std::vector<double> const& angle, std::vector< std::complex<sc::REAL> > const& prevertex)
{
	std::complex<sc::REAL> C(1.0, 1.0);
	SCMapSeries scMap(N, angle, prevertex, C);

	std::ofstream os("shape.txt");


	unsigned M = 5000;
	std::cout << "Start mapping points:" << std::endl;
	std::complex<sc::REAL> cp, sp;

	for (unsigned i = 0; i < M; i++) {
		cp = std::complex<sc::REAL>(std::cos(2 * i*pi / M), std::sin(2 * i*pi / M));
		sp = scMap.getPoint(cp);

		if (i % 1000 == 0)
			std::cout << "i=" << i << ": p=" << cp << ", f(p)= " << sp << std::endl;

		os << sp.real() << ' ' << sp.imag() << std::endl;
	}
	os.close();
}

void
writeContour(SCMapSeries const& scMap)
{
	std::ofstream os("shape.txt");


	unsigned M = 5000;
	std::cout << "Start mapping points:" << std::endl;
	std::complex<sc::REAL> cp, sp;

	for (unsigned i = 0; i < M; i++) {
		cp = std::complex<sc::REAL>(std::cos(2 * i*pi / M), std::sin(2 * i*pi / M));
		sp = scMap.getPoint(cp);

		if (i % 1000 == 0)
			std::cout << "i=" << i << ": p=" << cp << ", f(p)= " << sp << std::endl;

		os << sp.real() << ' ' << sp.imag() << std::endl;
	}
	os.close();
}


int 
main(int argc, char* argv[])
{
	std::vector<double> angle;
	std::vector< std::complex<sc::REAL> > prevertex, vertex;
	std::complex<sc::REAL> C(1.0, 1.0);

	setTestPolygon(angle, prevertex);

	try {
#if 0
		checkRadius(angle, prevertex);
		writeContour(1000, angle, prevertex);
#else		
		std::complex<sc::REAL> p(1.0, 1.0);
		SCMapSeries scMap(0.05, angle, prevertex, C);
		
		writeContour(scMap);
		checkRadius(scMap, prevertex);
#endif
	}
	catch (std::runtime_error& e) {
		std::cout << e.what() << std::endl;
	}

	return 0;
}

