// main.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <complex>

//#define _USE_MATH_DEFINES
#include <cmath>

#include "types.h"
#include "SCMapSeries.h"
#include "StressInfPlate.h"

void
setTestPolygon(std::vector<double>& angle, std::vector< std::complex<sc::REAL> >& prevertex)
{
	angle.clear();
	prevertex.clear();
	unsigned shape = 1;

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
	std::vector<unsigned> N = { 1000, 2000, 3000, 5000 };
	SCMapSeries* scMap;

	std::ofstream os("radius.txt");

	for each (auto n in N)
	{
		std::cout << n << std::endl;
		scMap = new SCMapSeries(n, angle, prevertex, C);

		int i = 1;
		for each (auto z in prevertex)
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
	for each (auto z in prevertex)
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

void
writeResults1D(SCMapSeries const& scMap, StressInfPlate const& sip)
{
	double R = 0.2, r;
	unsigned Nr = 100000;
	std::complex<sc::REAL> cp, mcp;
	stress ss;

	std::ofstream os("line_results.dat");
	//os << "# r sr sa t" << std::endl;

	for (unsigned i = 0; i < Nr; i++)
	{
		r = 1.0 + R*i / Nr;

		cp = std::complex<sc::REAL>(r/std::sqrt(2), r/std::sqrt(2));
		mcp = scMap.getPoint(cp);
		ss = sip.getValue(cp);

		os << r << "\t" << ss.r << "\t" << ss.f << "\t" << ss.rf << std::endl;
	}
	os.close();
}

void
writeResults2D(SCMapSeries const& scMap, StressInfPlate const& sip)
{
	double R = 6.0, r, phi, x;
	unsigned Nr = 100, Nphi = 500;
	std::complex<sc::REAL> cp, mcp;
	stress sp, sc;

	std::ofstream osCSV("results.csv");

	unsigned node = 1;
	for (unsigned i = 0; i < Nr; i++)
	{
		x = 1.0*i / Nr;
		r = 1.0 + R*std::pow(x,2);
		//std::cout << "r=" << r << std::endl;

		for (unsigned j = 0; j < Nphi; j++) 
		{
			phi = 2 * j*pi / Nphi;
			cp = std::complex<sc::REAL>(r*std::cos(phi), r*std::sin(phi));
			mcp = scMap.getPoint(cp);
			sp = sip.getValue(cp);

			// c_x, c_y, f(c)_x, f(c)_y, z, ss_r, ss_a, ss_s, thickness
			osCSV << std::real(mcp) << "," << std::imag(mcp) << "," << 1.0 - sp.x << "," << sp.x << "," << sp.y << "," << sp.x << "," << sp.y << "," << sp.xy << std::endl;
		}
	}
	osCSV.close();
	
	// matlab output
	std::ofstream osX("x.dat");
	std::ofstream osY("y.dat");
	std::ofstream osrZ("sr.dat");
	std::ofstream osfZ("sf.dat");
	std::ofstream osrfZ("srf.dat");
	std::ofstream osxZ("sx.dat");
	std::ofstream osyZ("sy.dat");
	std::ofstream osxyZ("sxy.dat");

	for (unsigned i = 0; i < Nr; i++)
	{
		x = 1.0*i / Nr;
		r = 1.0 + R*std::pow(x, 2);
		//std::cout << "r=" << r << std::endl;

		for (unsigned j = 0; j <= Nphi; j++)
		{
			phi = 2 * j*pi / Nphi;
			cp = std::complex<sc::REAL>(r*std::cos(phi), r*std::sin(phi));
			mcp = scMap.getPoint(cp);
			sp = sip.getValue(cp);

			osX << std::real(mcp) << "\t";
			osY << std::imag(mcp) << "\t";
			osrZ << sp.r << "\t";
			osfZ << sp.f << "\t";
			osrfZ << sp.rf << "\t";
			osxZ << sp.x << "\t";
			osyZ << sp.y << "\t";
			osxyZ << sp.xy << "\t";
		}
		osX << std::endl;
		osY << std::endl;
		osrZ << std::endl;
		osfZ << std::endl;
		osrfZ << std::endl;
		osxZ << std::endl;
		osyZ << std::endl;
		osxyZ << std::endl;
	}

	osX.close();
	osY.close();
	osrZ.close();
	osfZ.close();
	osrfZ.close();
	osxZ.close();
	osyZ.close();
	osxyZ.close();


	std::ofstream osTXT("hole_geo.txt");
	Nphi = 200;

	osTXT << "#Hole with N=" << scMap.getSize() << std::endl;
	osTXT << "#Group_ID Point_ID X_cord Y_cord Z_cord" << std::endl;

	for (unsigned i = 0; i < Nphi; i++)
	{
		phi = 2 * i*pi / Nphi;
		cp = std::complex<sc::REAL>(std::cos(phi), std::sin(phi));
		mcp = scMap.getPoint(cp);

		osTXT << 1 << " " << i+1 << " " << std::real(mcp) << " " << std::imag(mcp) << " " << 0.0 << std::endl;
	}
	osTXT.close();
}

void
appendPointStress(unsigned N, double sA, double sB, double sC)
{
	std::ofstream ofs;
	ofs.open("point_stress.txt", std::ofstream::out | std::ofstream::app);

	ofs << N << "\t" << sA << "\t" << sB << "\t" << sC << std::endl;

	ofs.close();
}

int
_tmain(int argc, _TCHAR* argv[])
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
		
#if 0
		double r = 0.1;
		SCMapSeries scMap(r, angle, prevertex, C);
		unsigned N = scMap.getSize() - 2;

		//writeContour(scMap);
		//checkRadius(scMap, prevertex);
#else
		//unsigned N = 3;
		//SCMapSeries scMap(N + 2, angle, prevertex, C);
#endif

#if 0
		std::vector<unsigned> nVec = {3, 7, 10, 20, 50, 100, 500, 1000, 5000, 10000, 20000, 50000, 100000};

		for (unsigned i = 0; i < nVec.size(); i++)
		{
			unsigned N = nVec[i];
			SCMapSeries scMap(N + 2, angle, prevertex, C);

			std::vector< std::complex<sc::REAL> > cc, c(scMap.getCoefficients());

			for (unsigned i = 0; i < N; i++)
				cc.push_back(c[i]);

			double s1 = 1.0;
			double s2 = 0.0;
			double a = 0;

			StressInfPlate sip(cc, s1, s2, a);

			stress ss;
			std::complex<sc::REAL> pA(1.0, 0);
			std::complex<sc::REAL> pB(1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0));
			std::complex<sc::REAL> pC(0., 1.0);

			double sA, sB, sC;
			std::cout << " -> N=" << c.size() << std::endl;

			ss = sip.getValue(pA);
			sA = ss.f;
			std::cout << " -> A: s_a=" << ss.f << ", A=" << ss.A << ", B=" << ss.B << std::endl;

			ss = sip.getValue(pB);
			sB = ss.f;
			std::cout << " -> B: s_a=" << ss.f << ", A=" << ss.A << ", B=" << ss.B << std::endl;

			ss = sip.getValue(pC);
			sC = ss.f;
			std::cout << " -> C: s_a=" << ss.f << ", A=" << ss.A << ", B=" << ss.B << std::endl;

			appendPointStress(N, sA, sB, sC);
		}
#endif

#if 1
		unsigned N = 5000;
		SCMapSeries scMap(N + 2, angle, prevertex, std::complex<sc::REAL>(1.,1.), 45.0);

		std::vector< std::complex<sc::REAL> > cc, c(scMap.getCoefficients());

		for (unsigned i = 0; i < N; i++)
			cc.push_back(c[i]);

		double s1 = 1.0;
		double s2 = 0.0;
		double a = 0.0;

		StressInfPlate sip(cc, s1, s2, a);

		std::cout << "Writing line-1D results to file... ";
		writeResults1D(scMap, sip);
		std::cout << "done!" << std::endl;
		
		std::cout << "Writing surface-2D results to file... ";
		writeResults2D(scMap, sip);
		std::cout << "done!" << std::endl;
		
#endif

#endif
	}
	catch (std::runtime_error& e) {
		std::cout << e.what() << std::endl;
	}

	std::cout << " -> Reached end of programm!" << std::endl;
	//getchar();
	return 0;
}

