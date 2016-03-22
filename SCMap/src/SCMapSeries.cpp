#include <stdexcept>
#include <iostream>
#include "SCMapSeries.h"


SCMapSeries::SCMapSeries
(
	unsigned N, 
	std::vector<double> const& angle, 
	std::vector< std::complex<sc::REAL> > const& prevertex, 
	std::complex<sc::REAL> const& C,
	double rot
)
{
	if (angle.size() != prevertex.size())
		throw std::runtime_error("Size Mismatch: angles and prevertex size does NOT match!");

	mM = angle.size();
	
	for (auto a : angle)
		mvAngle.push_back(1.0 - a);

	mvPrevertex = prevertex;
	mC = C;
	mRotation = rot;

	generate(N);
	getMaxRadius();
}

SCMapSeries::SCMapSeries
(
	double r,
	std::vector<double> const& angle,
	std::vector< std::complex<sc::REAL> > const& prevertex,
	std::complex<sc::REAL> const& C
)
{
	if (angle.size() != prevertex.size())
		throw std::runtime_error("Size Mismatch: angles and prevertex size does NOT match!");

	mM = angle.size();

	for (auto a : angle)
		mvAngle.push_back(1.0 - a);

	mvPrevertex = prevertex;
	mC = C;

	setSize(r);

	mM = mvCoefficient.size();
}


SCMapSeries::~SCMapSeries()
{
}

std::complex<sc::REAL> 
SCMapSeries::getPoint(std::complex<sc::REAL> const& p) const
{
	std::complex<sc::REAL> buf(p), d, q(1.0, 0.0);

	// We start with coefficient z^(-1) -> mvCoefficient[0] = d_1 and must set q=1/p
	q /= p;

	for (unsigned i = 0; i < mvCoefficient.size(); i++) 
	{
		d = mvCoefficient[i] * q;
		//if (std::abs(mvCoefficient[i]) < 1e-300)
		//	std::cout << "mvCoefficient[" << i << "] < 1e-300 : " << mvCoefficient[i] << std::endl;

		buf += d;
		// q is multiplied with 1/p
		q /= p;
	}

	//TODO: fix the constant A and C. 
	// A is additive constant
	//std::complex<sc::REAL> A(.0, .0);
	//return mC*(p + A/mC + buf );
	
	// rotate object for angle mRotation
	if (mRotation != 0.0) {
		std::complex<sc::REAL> A(std::cos(mRotation*pi / 180.0), std::cos(mRotation*pi / 180.0));
		buf *= A;
	}

	return buf;
}

sc::REAL
SCMapSeries::getRadius(std::complex<sc::REAL> const& p) const
{
	std::complex<sc::REAL> d1 = getD1(p);
	std::complex<sc::REAL> d2 = getD2(p);

	sc::REAL r = std::pow(std::abs(d1), 3) / std::imag(std::conj(d1)*d2);

	return r;
}

void 
SCMapSeries::generate(unsigned N)
{
	std::cout << "Start map generation ..." << std::endl;
	setDknMatrix(N);
	std::cout << " -> dnk matrix set" << std::endl;

	std::vector< std::complex<sc::REAL> > a, tmp;
	a = mDknMatrix[0];

	for (unsigned n = 0; n < (mM-1); n++) 
	{
		std::cout << " -> n=" << n + 1 << " series product started" << std::endl;
		
		multiplyTwoSums(N, a, mDknMatrix[n+1], tmp);
		a = tmp;

		std::cout << "     done." << std::endl;
	}

	std::complex<sc::REAL> d;
	sc::REAL r;

	mvCoefficient.clear();
	for (unsigned i = 2; i < tmp.size(); i++) 
	{
		r = 1.0*(i - 1);
		d = -std::complex<sc::REAL>(tmp[i].real() / r, tmp[i].imag() / r);
		mvCoefficient.push_back(d);
	}
}

void 
SCMapSeries::setDknMatrix(unsigned N)
{
	std::vector< std::complex<sc::REAL> > dn;
	std::complex<sc::REAL> buf;

	mDknMatrix.clear();

	for (unsigned n = 0; n < mM; n++)
	{
		buf = std::complex<sc::REAL>(1.0, 0.0);
		dn.push_back(buf);

		for (unsigned k = 1; k < N; k++)
		{
			
			buf = getDkn(buf, mvAngle[n], mvPrevertex[n], k);
			
			// Sign determination
			if( k % 2 == 1)
				dn.push_back(-buf);
			else
				dn.push_back(buf);	
		}
		
		mDknMatrix.push_back(dn);
		dn.clear();
	}
}

void 
SCMapSeries::multiplyTwoSums
(
	unsigned N,
	std::vector< std::complex<sc::REAL> > const& a, 
	std::vector< std::complex<sc::REAL> > const& b, 
	std::vector< std::complex<sc::REAL> >& c
)
{
	std::complex<sc::REAL> buf(0,0);
	c.clear();

	for (unsigned k = 0; k < N; k++)
	{
		buf = std::complex<sc::REAL>(0, 0);

		for (unsigned i = 0; i <= k; i++) {
			buf += a[i] * b[k - i];
		}

		if (k % 1000 == 0)
			std::cout << " -> coefficient dk=" << k << " done" << std::endl;
		
		if (std::abs(buf) < 1e-300)
			std::cout << "abs(c[" << k << "]) < 1e-300 : " << std::abs(buf) << std::endl;

		c.push_back(buf);
	}
}

std::complex<sc::REAL> 
SCMapSeries::getDkn
(
	std::complex<sc::REAL> const& dk,
	double alpha,
	std::complex<sc::REAL> const& pv,
	unsigned k
)
{
	/* We calculate term in series expansion 
	 
	 ( 1 \pm  \frac{a}{x} )^m = \sum_{i=0}^N d_i x^{-i} 
	 
	 as:

	 d_0 = 1
	 ... 
	 d_k = (-1)^k (1 - \frac{m+1}{k}) a d_{k-1}, sign is added in SCMapSeries::setDknMatrix()-"Sign determination"
	*/

	// dp: koefficient in k-1
	// pv: prevertex at position n

	sc::REAL r = (1.0 - 1.0*(alpha + 1) / k);
	return r*pv*dk;
}

std::complex<sc::REAL>
SCMapSeries::getD1(std::complex<sc::REAL> const& z) const
{
	/* This is first derivative of series function:
	f = 1 - 1*c_{-1}/z^{-2} - 2*c_{-2}/z^{-3} - 3*c_{-3}/z^{-3} - ... - n*c_{-n}/z^{-n}
	*/
	std::complex<sc::REAL> buf(1.0, 0.0), q(1.0, 0.0);

	// get pow=-2
	q /= z;
	q /= z;

	for (unsigned i = 0; i < mvCoefficient.size(); i++)
	{
		buf -= std::complex<sc::REAL>(1.0*(i + 1), 0.0) * mvCoefficient[i] * q;
		q /= z;
	}

	return buf;
}

std::complex<sc::REAL>
SCMapSeries::getD2(std::complex<sc::REAL> const& z) const
{
	/* This is second derivative of series function:
	f = 1*2*c_{-1}/z^{-3} + 2*3*c_{-2}/z^{-4} - 3*4*c_{-3}/z^{-4} + ... + n*(n+1)c_{-n}/z^{-(n+1)}
	*/
	std::complex<sc::REAL> buf(1.0, 0.0), q(1.0, 0.0);
	sc::REAL r;

	// get pow=-3
	q /= z;
	q /= z;
	q /= z;

	for (unsigned i = 0; i < mvCoefficient.size(); i++)
	{
		r = 1.0*(i + 1)*(i + 2);
		buf += std::complex<sc::REAL>(r, 0.0) * mvCoefficient[i] * q;
		q /= z;
	}

	return buf;
}

void
SCMapSeries::setSize(double r)
{
	double tol = 1e-2, q;
	unsigned N = 1000, Nmax, Nmin;
	generate(N);
	q = getMaxRadius();

	// Find N with multiplication by 2
	while ( q > r)
	{
		N *= 2;
		generate(N);
		q = getMaxRadius();
		std::cout << "First loop N=" << N << ", maxR=" << q << std::endl;
	}

	Nmax = N;
	Nmin = std::floor(N / 2);

	N = std::floor(3 * Nmax / 4);
	generate(N);
	q = getMaxRadius();

	// Start bisection serch
	while (std::abs(r - q) > tol)
	{
		if (q > r) {
			Nmin = N;
			N += (unsigned)std::floor((Nmax - N) / 2);
		}
		else {
			// Chek if difference is less than 1%
			q = (N - Nmin) / (2.0 * N);
			//std::cout << "Relative difference in N dN/N=" << q << std::endl;
			if (q < tol) {
				q = getMaxRadius();
				std::cout << "Final N=" << N << ", Rmax=" << q << std::endl;
				break;
			}

			Nmax = N;
			N -= std::floor((N - Nmin) / 2);
		}
		std::cout << "Bisection loop N=" << N << std::endl;
		generate(N);
		q = getMaxRadius();
	}
}

double
SCMapSeries::getMaxRadius()
{
	double r = 0.0, q;
	for (auto v : mvPrevertex)
	{
		q = std::abs(getRadius(v));
		if( q > r) r = q;
		std::cout << "vertex: " << v << ", r=" << q << ", maxR=" << r << std::endl;
	}

	return r;
}
