#include"crossSection.h"


double Event::getBeamEnergy() const 
{
	double gamma{ Ep1() / m_ma };
	assert(gamma > 1);
	double beta{ std::sqrt(1 - 1/(gamma*gamma)) };
	return gamma*(1 + beta)*Ep2();
	//return Ep2();
}

void Event::setS(double s)
{
	m_s = s;
	m_crossSection = std::numeric_limits<double>::min();
}

template<size_t N>
void Event::serialize(const std::array<double, N>& yVals, const std::array<double, N>& xVals,
																	 std::string& data) const
{
	for(size_t i{ 0 }; i < N; i++)
	{
		data.append(std::to_string(xVals[i]));
		data.append(" ");
		std::ostringstream stream;
		stream<<std::scientific<<std::setprecision(20)<<yVals[i];
		data.append(stream.str());
		data.append("\n");
	}
}

double Event::getConst(Constant c) 
{
	static constexpr std::array<double, static_cast<size_t>(Constant::max_constant)> constants{3.14159265359, 0.1056583755,
																																	0.51099895000e-3, 938.27208943e-3, 939.57e-3, 493.677e-3,
																																	497.611e-3, 1.16639e-5, 92.4e-3, 0.22, 299792458}; 
	return constants[static_cast<size_t>(c)];
}

double Event::kallen(double x, double y, double z) const
{
	double ret{ std::pow(x, 4) + std::pow(y, 4) + std::pow(z, 4) - 2*std::pow(x, 2)*std::pow(y, 2) 
				 - 2*std::pow(x, 2)*std::pow(z, 2) - 2*std::pow(y, 2)*std::pow(z, 2) }; 	
	if(fabs(ret) < 1e-10)
		    ret = 0;
	assert(ret >= 0);
	return std::sqrt(ret);
}

double Event::m12Calc(double E3) const
{
	double ret{ m_s + m_m3*m_m3 - 2*(std::sqrt(m_s))*E3 };
	assert(ret >= 0);
	return std::sqrt(ret);
}

double Event::m23Calc(double E3, double thetaStar) const
{
	double ret{ m_s + m_m1*m_m1 - 2*std::sqrt(m_s)*q10(E3, thetaStar) };	
	assert(ret >= 0);
	return std::sqrt(ret);
}

double Event::m13Calc(double E3, double thetaStar) const
{
	double ret{ m_s + m_m2*m_m2 - 2*std::sqrt(m_s)*q20(E3, thetaStar) };
	assert(ret >= 0);
	return std::sqrt(ret);
}

double Event::integrand(double E3, double theta, double thetaStar, double phiStar) const
{
	double m12{ m12Calc(E3) };
	double phaseSpaceFactor{ kallen(m12, m_m1, m_m2) / (m12*m12) };
	double formFactor{ dipoleFormFactor(E3, theta, thetaStar, phiStar)/*dipoleFormFactor(E3, theta, thetaStar, phiStar)*/ };
	return std::sqrt(std::pow(E3, 2) - std::pow(m_m3, 2)) * phaseSpaceFactor * ampSquared(E3, theta, thetaStar, phiStar) * 
				 std::sin(theta) * std::sin(thetaStar) * formFactor;
}

double Event::thetaPrimeIntegrand(double E3, double theta) const
{
	double gamma{ Ep1() / m_ma };
	assert(gamma > 1);
	double beta{ std::sqrt(1 - 1/(gamma*gamma)) };
	double numerator{ q3(E3) * std::sin(theta) };
	double denominator{ gamma * (q3(E3)*std::cos(theta) - beta*E3) };
	return std::atan(numerator / denominator);
}

double Event::fluxFactor() const
{
	//double p1Dotp2{ m_p1->dotProduct(*m_p2) };
	double p1Dotp2{ pk() };
	return 1/(64*std::pow(2*getConst(Constant::pi), 4) * std::sqrt(std::pow(p1Dotp2, 2) -
																																		 std::pow(m_ma, 2)*std::pow(m_mb, 2)));
}

double Event::E3Max() const
{
	double q3Max{ kallen(std::sqrt(m_s), m_m1 + m_m2, m_m3) / (2*std::sqrt(m_s)) };
	return std::sqrt(std::pow(q3Max, 2) + std::pow(m_m3, 2));
}

void Event::serializeEvent(std::string& out) const
{
	// beamEnergy crossSection
	out.append(std::to_string(getBeamEnergy()));
	out.append(" ");

	if(auto cs = getCrossSection())
	{
		std::ostringstream stream;
		double temp{ *cs * 3.89739e-28 };	// cm^2
		stream<<std::scientific<<std::setprecision(10)<<temp;
		out.append(stream.str());
	}
	
	else
		out.append("error");

	out.append("\n");
}

double Event::ampSquared(double E3, double theta, double thetaStar, double phiStar) const
{
	double a{ pk() };
	double b{ ppPrime(E3, theta, thetaStar, phiStar) };
	double c{ pkPrime(E3, theta, thetaStar, phiStar) };
	double d{ pPrimekPrime(E3, theta, thetaStar, phiStar) };
	double e{ kkPrime(E3, theta, thetaStar, phiStar) };
	double f{ pPrimek(E3, theta, thetaStar, phiStar) };
	if(m_ampSquared)
	{
		return m_ampSquared(a, b, c, d, e, f);
	}
	else
		return 0;
}

double Event::dipoleFormFactor(double E3, double theta, double thetaStar, double phiStar) const
{
	double p2Dotqi{ kkPrime(E3, theta, thetaStar, phiStar) };
	double mi{ m_m3 };
	double denominator{ 1 - (mi*mi - 2*p2Dotqi) };
	double ret{ 1 / (denominator*denominator) };
	//std::cout<<ret<<'\n';
	return ret;
	//return 1;
}

double Event::thetaPrime(double theta, double Eoverq) const
{
	double gamma{ Ep1() / m_ma };
	assert(gamma > 1);
	double beta{ std::sqrt(1 - 1/(gamma*gamma)) };
	return std::atan(std::sin(theta) / ((std::cos(theta) + beta*Eoverq)));
}

double Event::EPrime(double E) const
{
	double gamma{ Ep1() / m_ma };
	assert(gamma > 1);
	double beta{ std::sqrt(1 - 1/(gamma*gamma)) };
	return gamma * E * (1 + (q3(E) / E) * beta);
	//return 2*E;
}
