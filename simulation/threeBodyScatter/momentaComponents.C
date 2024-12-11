#include"crossSection.h"


double Event::Ep1() const
{
	return (m_s + m_ma*m_ma) / (2*std::sqrt(m_s));
}

double Event::Ep2() const
{
	return std::sqrt(m_s) - Ep1();
}

double Event::pk() const
{
	assert(Ep1() >= m_ma);
	return (Ep1()*Ep2() + std::sqrt(Ep1()*Ep1() - m_ma*m_ma)*Ep2()); 	
}

double Event::ppPrime(double E3, double theta, double thetaStar, double phiStar) const
{
	assert(Ep1() >= m_ma);
	return (Ep1()*q10(E3, thetaStar) -
						  q13(E3, theta, thetaStar, phiStar)*std::sqrt(Ep1()*Ep1() - m_ma*m_ma));
}

double Event::pkPrime(double E3, double theta, double thetaStar, double phiStar) const
{
	assert(Ep1() >= m_ma);
	// q3 = kaon => use q2 
	return (Ep1()*q20(E3, thetaStar) - std::sqrt(Ep1()*Ep1() - m_ma*m_ma)*q23(E3, theta, thetaStar, phiStar));
	// q3 = lepton => use q3
	//return (Ep1()*q30(E3, theta, thetaStar, phiStar) - std::sqrt(Ep1()*Ep1() - m_ma*m_ma)*q33(E3, theta, thetaStar, phiStar));
}

double Event::pPrimekPrime(double E3, double theta, double thetaStar, double phiStar) const
{
	// q3 = kaon => use q2
	return (q10(E3, thetaStar)*q20(E3, thetaStar) - 
							q11(E3, theta, thetaStar, phiStar)*q21(E3, theta, thetaStar, phiStar) -
							q12(E3, thetaStar, phiStar)*q22(E3, thetaStar, phiStar) -
							q13(E3, theta, thetaStar, phiStar)*q23(E3, theta, thetaStar, phiStar));
	//p3 = lepton => use q3
	//return (q10(E3, thetaStar)*q30(E3, theta, thetaStar, phiStar) - 
	//						q11(E3, theta, thetaStar, phiStar)*q31(E3, theta, thetaStar, phiStar) -
	//						q12(E3, thetaStar, phiStar)*q32(E3, theta, thetaStar, phiStar) -
	//						q13(E3, theta, thetaStar, phiStar)*q33(E3, theta, thetaStar, phiStar));
}

double Event::kkPrime(double E3, double theta, double thetaStar, double phiStar) const
{
	// q3 = kaon => use q2
	return (q20(E3, thetaStar)*Ep2() + q23(E3, theta, thetaStar, phiStar)*Ep2());
	// q3 = lepton => use q3
	//return (q30(E3, theta, thetaStar, phiStar)*Ep2() + q33(E3, theta, thetaStar, phiStar)*Ep2());
}

double Event::pPrimek(double E3, double theta, double thetaStar, double phiStar) const
{
	return (q10(E3, thetaStar)*Ep2() + q13(E3, theta, thetaStar, phiStar)*Ep2());
}

double Event::q12CM(double E3) const
{
	return kallen(m12Calc(E3), m_m1, m_m2) / (2*m12Calc(E3));
}

double Event::q10Star(double E3) const
{
	double m12{ m12Calc(E3) };
	return (m12*m12 + m_m1*m_m1 - m_m2*m_m2) / (2*m12);
}

double Event::q20Star(double E3) const
{
	double m12{ m12Calc(E3) };
	return (m12*m12 - m_m1*m_m1 + m_m2*m_m2) / (2*m12);
}

double Event::q11Star(double E3, double thetaStar, double phiStar) const
{
 	return q12CM(E3) * std::sin(thetaStar) * std::cos(phiStar);
}

double Event::q21Star(double E3, double thetaStar, double phiStar) const
{
	return -q11Star(E3, thetaStar, phiStar);
}

double Event::q12Star(double E3, double thetaStar, double phiStar) const
{
	return q12CM(E3) * std::sin(thetaStar) * std::sin(phiStar);
}

double Event::q22Star(double E3, double thetaStar, double phiStar) const
{
	return -q12Star(E3, thetaStar, phiStar);
}

double Event::q13Star(double E3, double thetaStar) const
{
	return q12CM(E3) * std::cos(thetaStar);
}

double Event::q23Star(double E3, double thetaStar) const
{
	return -q13Star(E3, thetaStar);
}

double Event::q10(double E3, double thetaStar) const
{
	double m12{ m12Calc(E3) };
	double gamma{ (m_s + m12*m12 - m_m3*m_m3) / (2 * std::sqrt(m_s) * m12) };
	if(fabs(1-gamma) < 0.0000001)
		    gamma = 1;
	assert(gamma >= 1);
	double beta{ std::sqrt(1 - 1/(gamma*gamma)) };
	return gamma * (q10Star(E3) + beta*q13Star(E3, thetaStar));
}

double Event::q20(double E3, double thetaStar) const
{
	double m12{ m12Calc(E3) };
	double gamma{ (m_s + m12*m12 - m_m3*m_m3) / (2 * std::sqrt(m_s) * m12) };
	if(fabs(1-gamma) < 0.0000001)
		    gamma = 1;
	assert(gamma >= 1);
	double beta{ std::sqrt(1 - 1/(gamma*gamma)) };
	return gamma * (q20Star(E3) + beta*q23Star(E3, thetaStar));
}

double Event::q11(double E3, double theta, double thetaStar, double phiStar) const
{
	double m12{ m12Calc(E3) };
	double gamma{ (m_s + m12*m12 - m_m3*m_m3) / (2 * std::sqrt(m_s) * m12) };
	if(fabs(1-gamma) < 0.0000001)
		    gamma = 1;
	assert(gamma >= 1);
	double beta{ std::sqrt(1 - 1/(gamma*gamma)) };
	double pi{ getConst(Constant::pi) };
	double firstTerm{ q11Star(E3, thetaStar, phiStar) * std::cos(theta + pi) };
	double secondTerm{ gamma * (beta*q10Star(E3) + q13Star(E3, thetaStar)) * std::sin(theta + pi) };
	return firstTerm + secondTerm;
}

double Event::q21(double E3, double theta, double thetaStar, double phiStar) const
{
	double m12{ m12Calc(E3) };
	double gamma{ (m_s + m12*m12 - m_m3*m_m3) / (2 * std::sqrt(m_s) * m12) };
	if(fabs(1-gamma) < 0.0000001)
		    gamma = 1;
	assert(gamma >= 1);
	double beta{ std::sqrt(1 - 1/(gamma*gamma)) };
	double pi{ getConst(Constant::pi) };
	double firstTerm{ q21Star(E3, thetaStar, phiStar) * std::cos(theta + pi) };
	double secondTerm{ gamma * (beta*q20Star(E3) + q23Star(E3, thetaStar)) * std::sin(theta + pi) };
	return firstTerm + secondTerm;
}

double Event::q12(double E3, double thetaStar, double phiStar) const
{
	return q12Star(E3, thetaStar, phiStar);
}

double Event::q22(double E3, double thetaStar, double phiStar) const
{
	return q22Star(E3, thetaStar, phiStar);
}

double Event::q13(double E3, double theta, double thetaStar, double phiStar) const
{
	double m12{ m12Calc(E3) };
	double gamma{ (m_s + m12*m12 - m_m3*m_m3) / (2 * std::sqrt(m_s) * m12) };
	if(fabs(1-gamma) < 0.0000001)
		    gamma = 1;
	assert(gamma >= 1);
	double beta{ std::sqrt(1 - 1/(gamma*gamma)) };
	double pi{ getConst(Constant::pi) };
	double firstTerm{ -q11Star(E3, thetaStar, phiStar) * std::sin(theta + pi) };
	double secondTerm{ gamma * (beta*q10Star(E3) + q13Star(E3, thetaStar)) * std::cos(theta + pi) };
	return firstTerm + secondTerm;
}

double Event::q23(double E3, double theta, double thetaStar, double phiStar) const
{
	double m12{ m12Calc(E3) };
	double gamma{ (m_s + m12*m12 - m_m3*m_m3) / (2 * std::sqrt(m_s) * m12) };
	if(fabs(1-gamma) < 0.0000001)
		    gamma = 1;
	assert(gamma >= 1);
	double beta{ std::sqrt(1 - 1/(gamma*gamma)) };
	double pi{ getConst(Constant::pi) };
	double firstTerm{ -q21Star(E3, thetaStar, phiStar) * std::sin(theta + pi) };
	double secondTerm{ gamma * (beta*q20Star(E3) + q23Star(E3, thetaStar)) * std::cos(theta + pi) };
	return firstTerm + secondTerm;
}

double Event::q3(double E3) const
{
	double m12{ m12Calc(E3) };
	return kallen(std::sqrt(m_s), m12, m_m3) / (2*std::sqrt(m_s));
}

double Event::q30(double E3, double theta, double thetaStar, double phiStar) const
{
	return E3;
}

double Event::q31(double E3, double theta, double thetaStar, double phiStar) const
{
	return q3(E3) * std::sin(theta);
}

double Event::q32(double E3, double theta, double thetaStar, double phiStar) const
{
	return 0;
}

double Event::q33(double E3, double theta, double thetaStar, double phiStar) const
{
	return q3(E3) * std::cos(theta);
}
