#include <iostream>
#include <cmath>
#include "StressInfPlate.h"

StressInfPlate::StressInfPlate(std::vector<std::complex<sc::REAL> > const& c, double sigma_1, double sigma_2, double alpha)
{
	mN = c.size();
	mCCoeff = zeros(mN + 1);

	double a = alpha*pi / 180;

	for (unsigned i = 1; i <= mN; i++)
		mCCoeff[i] = c[i - 1];

	// initialize vectors to be ready for filling
	mPCoeff = zeros(mN + 1);
	mQCoeff = zeros(mN + 1);
	mACoeff = zeros(mN + 1);
	mKCoeff = zeros(mN + 1);

	mGamma1 = std::complex<double>((sigma_1 + sigma_2) / 4.0, 0.0);
	mGamma2 = std::complex<double>(0.5*(sigma_2 - sigma_1)*std::cos(-2.0*a), 0.5*(sigma_2 - sigma_1)*std::sin(-2.0*a));

	//mGamma1 = std::complex<double>(0.25, 0.0);
	//mGamma2 = std::complex<double>(-0.5, 0.0);

	std::cout << " -> g=" << mGamma1 << ", g'=" << mGamma2 << std::endl;

	generateCoefficients();
}

StressInfPlate::~StressInfPlate()
{
}

stress
StressInfPlate::getValue(std::complex<sc::REAL> const& p) const
{
	std::complex<sc::REAL> omega = getOmega(p);
	std::complex<sc::REAL> Domega = getDOmega(p);
	std::complex<sc::REAL> DDomega = getDDOmega(p);

	std::complex<sc::REAL> varphi = getVarphi(p);
	std::complex<sc::REAL> Dvarphi = getDVarphi(p);
	std::complex<sc::REAL> DDvarphi = getDDVarphi(p);

	std::complex<sc::REAL> Phi(0., 0.), DPhi(0., 0.);

	if (std::abs(Domega) != 0.0)
	{
		Phi = Dvarphi / Domega;
		DPhi = (DDvarphi*Domega - Dvarphi*DDomega) / (Domega*Domega);
	}

	std::complex<sc::REAL> Q = getQ(p);
	std::complex<sc::REAL> DQ = getDQ(p);

	//std::complex<sc::REAL> psi = mGamma2*p - mGamma1/p - std::conj(Q)*Phi - getK(p);
	std::complex<sc::REAL> Dpsi = mGamma2 + mGamma1/(p*p) - DQ*Phi - Q*DPhi - getDK(p);

	std::complex<sc::REAL> A =  std::complex<sc::REAL>(2., 0.) * ( Phi + std::conj(Phi) );
	
	std::complex<sc::REAL> t_pc = std::complex<sc::REAL>(2.0, 0.0) *p*p / (std::abs(p*p) *  std::conj(Domega));
	std::complex<sc::REAL> t_cc = std::complex<sc::REAL>(2.0, 0.0) / Domega;
	std::complex<sc::REAL> q = std::conj(omega)*DPhi + Dpsi;

	std::complex<sc::REAL> B = t_pc * q;
	std::complex<sc::REAL> Bxy = t_cc * q;

	//std::cout << "z=" << p << ", A=" << A << ", B=" << B << std::endl;

	stress ss;

	ss.r = (std::real(A) - std::real(B)) / 2;
	ss.f = (std::real(A) + std::real(B)) / 2;
	ss.rf = 0.5 * std::imag(B);
	ss.A = A;
	ss.B = B;

	ss.x = (std::real(A) - std::real(Bxy)) / 2;
	ss.y = (std::real(A) + std::real(Bxy)) / 2;
	ss.xy = 0.5*std::imag(Bxy);

	ss.fi = std::arg(p*p / (std::abs(p*p) *  Domega / std::conj(Domega))) / 2.0;

	return ss;
}

void
StressInfPlate::generateCoefficients()
{
	setPCoeff();	
	setQCoeff();
	solveACoeff();
	setKCoeff();
#if 0
	printVec(mCCoeff, "C");
	printVec(mPCoeff, "P");
	printVec(mQCoeff, "Q");
	printVec(mACoeff, "A");
	printVec(mKCoeff, "K");
#endif
	std::cout << std::endl << " -> All coefficients done!" << std::endl;
}

void
StressInfPlate::setPCoeff()
{
	std::complex<sc::REAL> buf(0, 0);

	mPCoeff[mN] = mCCoeff[mN];
	mPCoeff[mN - 1] = mCCoeff[mN - 1];

	//std::cout << "N=" << mN << std::endl << std::endl;
	for (unsigned n = (mN - 2); n > 0; n--) // index runing from N-2,...,1
	{
		//std::cout << "n=" << n << std::endl;
		buf = mCCoeff[n];

		for (unsigned k = 1; k <= (mN - n - 1); k++)
		{
			//std::cout << "k=" << k << ", k+n+1=" << k+n+1 << std::endl;
			buf += std::complex<sc::REAL>(1.0*k,0.0)*(std::conj(mCCoeff[k])*mPCoeff[k + n + 1]);
		}
		//std::cout << std::endl;
		mPCoeff[n] = buf;
		//std::cout << "   p[" << n << "]: " << mPCoeff[n] << std::endl;
	}
}

void
StressInfPlate::setQCoeff()
{
	std::complex<sc::REAL> buf(0., 0.), cn1(1., 0.), z0(0., 0.);

	mQCoeff[0] += z0;
	for (unsigned k = 1; k < mN; k++)
		mQCoeff[0] += std::complex<sc::REAL>(1.0*k, 0.0)*std::conj(mCCoeff[k])*mPCoeff[k + 1];

	mQCoeff[1] += cn1;
	for (unsigned k = 1; k <= mN; k++)
		mQCoeff[1] += std::complex<sc::REAL>(1.0*k, 0.0)*std::conj(mCCoeff[k])*mPCoeff[k];

	//std::cout << "N=" << mN << std::endl << std::endl;
	for (unsigned n = 2; n <= mN; n++) // index runing from N-2,...,1
	{
		//std::cout << "n=" << n << std::endl;
		buf = std::complex<sc::REAL>(0., 0.);

		for (unsigned k = n; k <= mN; k++)
		{
			//std::cout << "k=" << k << ", k-n+1=" << k-n+1 << std::endl;
			buf += std::complex<sc::REAL>(1.0*k, 0.0)*(std::conj(mCCoeff[k])*mPCoeff[k - n + 1]);
		}
		//std::cout << std::endl;
		mQCoeff[n] = buf;
		//std::cout << "   q[" << n << "]: " << mQCoeff[n] << std::endl;
	}
}

void
StressInfPlate::solveACoeff()
{
	sc::REAL err_old = 0., err_new, EPS = 1e-50;
	std::vector<std::complex<sc::REAL> > a_old, a_new, diff;

	a_old = mCCoeff;

	setASystem(a_old, a_new);
	//printVec(a_new, "A");
	
	err_new = systemError(a_old, a_new);
	std::cout << " -> loop_err_INITIAL=" << err_new << std::endl;

	while ( std::abs(err_old - err_new) > EPS ) // iteratively solve until the differrence in error norm is less than EPS;
	{
		a_old = a_new;
		err_old = err_new;
		
		setASystem(a_old, a_new);
		//printVec(a_new, "A");
		
		err_new = systemError(a_old, a_new);

		std::cout << " -> err_diff=" << std::abs(err_old - err_new) << std::endl;
	}

	std::cout << " -> loop_err_FINAL=" << err_new << std::endl;

	mACoeff = a_new;
}

void
StressInfPlate::setASystem(
	std::vector<std::complex<sc::REAL> > const& a, 
	std::vector<std::complex<sc::REAL> > & b
	)
{
	std::complex<sc::REAL> buf(0., 0.), cn1(1., 0.), z0(0., 0.);
	b = zeros(mN + 1);

	//TODO: preveri ali tečejo vsi indeksi pravilno
	//std::cout << "N=" << mN << std::endl;

	for (unsigned k = 1; k <= (mN-2); k++)
	{
		//std::cout << "n=1, k=" << k << std::endl;
		b[1] += std::complex<sc::REAL>(1.0*k, 0.0)*mPCoeff[k + 2] * std::conj(a[k]);
	}
	b[1] -= mGamma1*mPCoeff[1] + std::conj(mGamma2);

	//getchar();

	for (unsigned n = 2; n < (mN - 1); n++)
	{
		for (unsigned k = 1; k <= (mN - n - 1); k++)
		{
			//std::cout << "n=" << n << ", k=" << k << std::endl;
			b[n] += std::complex<sc::REAL>(1.0*k, 0.0)*mPCoeff[k + n + 1] * std::conj(a[k]);
		}
		b[n] -= mGamma1*mPCoeff[n];
		//getchar();
	}

	b[mN - 1] = -mGamma1*mPCoeff[mN - 1];
	b[mN] = -mGamma1*mPCoeff[mN];
}

void
StressInfPlate::setKCoeff()
{
	//TODO: preveri ali tečejo vsi indeksi pravilno

	for (unsigned n = 1; n <= mN; n++)
	{
		for (unsigned k = n; k <= mN; k++)
		{
			//std::cout << "n=" << n << ", k=" << k << ", k-n+1=" << k-n+1 << std::endl;
			mKCoeff[n] -= std::complex<sc::REAL>(1.0*k, 0.0)*mPCoeff[k - n + 1] * std::conj(mACoeff[k]);
		}
	}
}

sc::REAL 
StressInfPlate::systemError(
	std::vector<std::complex<sc::REAL> > const& a, 
	std::vector<std::complex<sc::REAL> > const& b
	)
{
	sc::REAL err = 0.;
	unsigned N = a.size();

	for (unsigned i = 0; i < N; i++)
		err += std::abs(a[i] - b[i]);
	
	return err;
}

std::complex<sc::REAL>
StressInfPlate::getOmega(
	std::complex<sc::REAL> const& z
) const
{
	std::complex<sc::REAL> buf(z), q(1.0, 0.0);

	// We start with coefficient z^(-1) -> mvCoefficient[0] = d_1 and must set q=1/p
	q /= z;

	for (unsigned i = 1; i < mCCoeff.size(); i++)
	{
		buf += mCCoeff[i] * q;
		// q is multiplied with 1/p
		q /= z;
	}

	return buf;
}

std::complex<sc::REAL>
StressInfPlate::getDOmega(
	std::complex<sc::REAL> const& z
) const
{
	std::complex<sc::REAL> buf(1., 0.), q(1.0, 0.0);

	// We start with coefficient z^(-2)
	q /= z;
	q /= z;

	for (unsigned i = 1; i < mCCoeff.size(); i++)
	{
		buf -= std::complex<sc::REAL>(1.0*i, .0)*mCCoeff[i] * q;
		// q is multiplied with 1/p
		q /= z;
	}

	return buf;
}

std::complex<sc::REAL>
StressInfPlate::getDDOmega(
	std::complex<sc::REAL> const& z
) const
{
	std::complex<sc::REAL> buf(0., 0.), q(1.0, 0.0);

	// We start with coefficient z^(-3)
	q /= z;
	q /= z;
	q /= z;

	for (unsigned i = 1; i < mCCoeff.size(); i++)
	{
		buf += std::complex<sc::REAL>(1.0*i*(i + 1), .0)*mCCoeff[i] * q;
		// q is multiplied with 1/p
		q /= z;
	}

	return buf;
}

std::complex<sc::REAL>
StressInfPlate::getVarphi(
	std::complex<sc::REAL> const& z
) const
{
	std::complex<sc::REAL> buf(0., 0.), q(1.0, 0.0);

	buf += mGamma1*z;

	// We start with coefficient z^(-1) -> mvCoefficient[0] = d_1 and must set q=1/p
	q /= z;

	for (unsigned i = 1; i < mACoeff.size(); i++)
	{
		buf += mACoeff[i] * q;
		// q is multiplied with 1/p
		q /= z;
	}

	return buf;
}

std::complex<sc::REAL>
StressInfPlate::getDVarphi(
	std::complex<sc::REAL> const& z
) const
{
	std::complex<sc::REAL> buf(mGamma1), q(1.0, 0.0);

	// We start with coefficient z^(-2)
	q /= z;
	q /= z;

	for (unsigned i = 1; i < mACoeff.size(); i++)
	{
		buf -= std::complex<sc::REAL>(1.0*i, .0)*mACoeff[i] * q;
		// q is multiplied with 1/p
		q /= z;
	}

	return buf;
}

std::complex<sc::REAL>
StressInfPlate::getDDVarphi(
	std::complex<sc::REAL> const& z
) const
{
	std::complex<sc::REAL> buf(0., 0.), q(1.0, 0.0);

	// We start with coefficient z^(-3)
	q /= z;
	q /= z;
	q /= z;

	for (unsigned i = 1; i < mACoeff.size(); i++)
	{
		buf += std::complex<sc::REAL>(1.0*i*(i + 1), .0) * mACoeff[i] * q;
		// q is multiplied with 1/p
		q /= z;
	}

	return buf;
}

std::complex<sc::REAL>
StressInfPlate::getQ(
	std::complex<sc::REAL> const& z
) const
{
	std::complex<sc::REAL> buf(mQCoeff[0]), q(1.0, 0.0);

	q /= z;

	for (unsigned i = 1; i < mQCoeff.size(); i++)
	{
		buf += mQCoeff[i] * q;
		q /= z;
	}

	return buf;
}

std::complex<sc::REAL>
StressInfPlate::getDQ(
	std::complex<sc::REAL> const& z
) const
{
	std::complex<sc::REAL> buf(0., 0.), q(1.0, 0.0);

	q /= z;
	q /= z;

	for (unsigned i = 1; i < mQCoeff.size(); i++)
	{
		buf -= std::complex<sc::REAL>(i,0.0) * mQCoeff[i] * q;
		q /= z;
	}

	return buf;
}

std::complex<sc::REAL>
StressInfPlate::getK(
std::complex<sc::REAL> const& z
) const
{
	std::complex<sc::REAL> buf(0., 0.), q(1.0, 0.0);


	// We start with coefficient z^(-1) -> mvCoefficient[0] = d_1 and must set q=1/p
	q /= z;

	for (unsigned i = 1; i < mKCoeff.size(); i++)
	{
		buf += std::conj(mKCoeff[i]) * q;
		// q is multiplied with 1/p
		q /= z;
	}

	return buf;
}

std::complex<sc::REAL>
StressInfPlate::getDK(
std::complex<sc::REAL> const& z
) const
{
	std::complex<sc::REAL> buf(0., 0.), q(1.0, 0.0);

	// We start with coefficient z^(-2)
	q /= z;
	q /= z;

	for (unsigned i = 1; i < mKCoeff.size(); i++)
	{
		buf -= std::complex<sc::REAL>(1.0*i, .0) * std::conj(mKCoeff[i]) * q;
		// q is multiplied with 1/p
		q /= z;
	}

	return buf;
}

std::vector<std::complex<sc::REAL> > 
StressInfPlate::zeros(unsigned N)
{
	std::vector<std::complex<sc::REAL> > v;

	for (unsigned i = 0; i < N; i++)
		v.push_back(std::complex<sc::REAL>(0., 0.));

	return v;
}

void
StressInfPlate::printVec(std::vector<std::complex<sc::REAL> > const& v, std::string const& name) const
{
	unsigned N = v.size();
	std::cout << "vector size=" << N << std::endl;
	std::cout << "vector name: " << name << std::endl;
	
	for (unsigned i = 0; i < N; i++)
		std::cout << "  [" << i << "] = " << v[i] << std::endl;
}
