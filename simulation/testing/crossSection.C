#include"fourVector.h"
#include"crossSection/timer.h"
#include"crossSection/threadPool.h"

#include<memory>
#include<cmath>
#include<functional>
#include<array>
#include<optional>
#include<limits>
#include<fstream>
#include<iomanip>
#include<string_view>
#include<cassert>


class Event
{
public:
	enum class Constant
	{
		pi,
		muonMass,						// GeV
		electronMass,				// GeV
		protonMass,					// GeV
		neutronMass,				// GeV
		kaonChargedMass,		// GeV
		kaon0Mass,					// Gev
		GF,									// GeV^-2
		fpi,								// GeV
		Vus,
		c,									// m/s
		max_constant
	};

private:
	enum class IntegrationVariable
	{
		E3,
		theta,
		thetaStar,
		phiStar
	};

public:
	Event() = default;

	Event(int id, std::unique_ptr<FourVector> p1, std::unique_ptr<FourVector> p2, std::unique_ptr<FourVector> p3, double ma,
				double mb, double m1, double m2, double m3, double s, double (*ampSquaredPtr)(double, double, double, double, 
																																										  double, double))
		: m_id{ id },
			m_p1{ std::move(p1) },
			m_p2{ std::move(p2) },
			m_p3{ std::move(p3) },
			m_ma{ ma }, 
			m_mb{ mb },
			m_m1{ m1 },
			m_m2{ m2 },
			m_m3{ m3 },
			m_s{ s },
			m_ampSquared{ ampSquaredPtr },
			m_crossSection{ std::numeric_limits<double>::min() }           // -inf => not been calculated
			//m_beamEnergy{ getBeamEnergy() }
	{ }

	Event(int id, double ma, double mb, double m1, double m2, double m3, double s, double (*ampSquaredPtr)(double, double,
																																												double, double, double, double))
		: m_id{ id },
			m_p1{ nullptr },
			m_p2{ nullptr },
			m_p3{ nullptr },
			m_ma{ ma }, 
			m_mb{ mb },
			m_m1{ m1 },
			m_m2{ m2 },
			m_m3{ m3 },
			m_s{ s },
			m_ampSquared{ ampSquaredPtr },
			m_crossSection{ std::numeric_limits<double>::min() }           // -inf => not been calculated
			//m_beamEnergy{ getBeamEnergy() }
	{ }

	Event(const Event& e)
		:	m_id{ e.m_id },
			m_p1{ e.m_p1 ? std::make_unique<FourVector>(*(e.m_p1)) : nullptr },
			m_p2{ e.m_p2 ? std::make_unique<FourVector>(*(e.m_p2)) : nullptr },
			m_p3{ e.m_p3 ? std::make_unique<FourVector>(*(e.m_p3)) : nullptr },
			m_ma{ e.m_ma },
			m_mb{ e.m_mb },
			m_m1{ e.m_m1 },
			m_m2{ e.m_m2 },
			m_m3{ e.m_m3 },
			m_s{ e.m_s },
			m_ampSquared{ e.m_ampSquared ? e.m_ampSquared : nullptr },
			m_crossSection{ e.m_crossSection }
			//m_beamEnergy{ getBeamEnergy() }
	{ }

	~Event(){}

public:
	std::optional<double> getCrossSection() const {return m_crossSection == std::numeric_limits<double>::min() ? std::nullopt :
																								 std::optional(m_crossSection);}
	//double getBeamEnergy() const {return m_beamEnergy;}
	double getBeamEnergy() const 
	{
		double gamma{ Ep1() / m_ma };
		assert(gamma > 1);
		double beta{ std::sqrt(1 - 1/(gamma*gamma)) };
		return gamma*(1 + beta)*Ep2();
		//return Ep2();
	}
	void setS(double s)
	{
		m_s = s;
		m_crossSection = std::numeric_limits<double>::min();
	}
	static double getConst(Constant c);
	double crossSectionCalc(int n, int m);
	double E3DiffCrossSection(double value, int n, int m, int p) const;
	double thetaDiffCrossSection(double value, int n, int m, int p) const;
	double thetaStarDiffCrossSection(double value, int n, int m, int p) const;
	double phiStarDiffCrossSection(double value, int n, int m, int p) const;
	void serializeM12M23(std::string& data) const
	{
		for(const auto& entry : m_m12m23Count)
		{
			data.append(std::to_string(entry.first.first));
			data.append(" ");
			data.append(std::to_string(entry.first.second));
			data.append(" ");
			data.append(std::to_string(entry.second));
			data.append("\n");
		}
	}
	void serializeM23M13(std::string& data) const
	{
		for(auto& entry : m_m23m13Count)
		{
			data.append(std::to_string(entry.first.first));
			data.append(" ");
			data.append(std::to_string(entry.first.second));
			data.append(" ");
			data.append(std::to_string(entry.second));
			data.append("\n");
		}	
	}
	void serializeM12M13(std::string& data) const
	{
		for(auto& entry : m_m12m13Count)
		{
			data.append(std::to_string(entry.first.first));
			data.append(" ");
			data.append(std::to_string(entry.first.second));
			data.append(" ");
			data.append(std::to_string(entry.second));
			data.append("\n");
		}
	}
	void serializeEvent(std::string& out) const;
	template<size_t N>
	void serialize(const std::array<double, N>& yVals, const std::array<double, N>& xVals,
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
	double E3Max() const;
	double E3Min() const { return m_m3; }

	double thetaLab(double theta) const;
	double dThetadThetaLab(double a, double b) const;
	double dthetadthetaprime(double theta, double Eoverq) const;
	double thetaPrime(double theta, double Eoverq) const;
	double EPrime(double E) const;

private:
	int m_id;
	std::unique_ptr<FourVector> m_p1 = std::make_unique<FourVector>();
	std::unique_ptr<FourVector> m_p2 = std::make_unique<FourVector>();
	std::unique_ptr<FourVector> m_p3 = std::make_unique<FourVector>();
	double m_ma;		// nucleon
	double m_mb;		// neutrino
	double m_m1;		// nucleon
	double m_m2;		// lepton
	double m_m3;		// kaon
	double m_s;
	double (*m_ampSquared)(double, double, double, double, double, double);		// (pk, pp', pk', p'k', kk', p'k)
	double m_crossSection;
	double m_beamEnergy;
	std::map<std::pair<double,double>, int> m_m12m23Count;
	std::map<std::pair<double,double>, int> m_m23m13Count;
	std::map<std::pair<double,double>, int> m_m12m13Count;

private:
	double fluxFactor() const;
	double simpsonCalc1D(double a, double b, int n, double fixedVariable) const;
	double kallen(double x, double y, double z) const;
	double m12Calc(double E3) const;
	double m23Calc(double E3, double thetaStar) const;
	double m13Calc(double E3, double thetaStar) const;

	double integrand(double E3, double theta, double thetaStar, double phiStar) const; 
	double thetaPrimeIntegrand(double E3, double theta) const;
	void weightModify(int& weight, int i, int n) const;
	double simpsonCalc4D(double a, double b, double c, double d, double e, double f, double g, double h, int n, int m,
											 int p, int q); 
	double simpsonCalc3D(IntegrationVariable fixedParameter, double fixedParameterValue, double a, double b, double c,
											 double d, double e, double f, int n, int m, int p) const;
	double Ep1() const;
	double Ep2() const;
	double pk() const;
	double ppPrime(double E3, double theta, double thetaStar, double phiStar) const;
	double pkPrime(double E3, double theta, double thetaStar, double phiStar) const;
	double pPrimekPrime(double E3, double theta, double thetaStar, double phiStar) const;
	double kkPrime(double E3, double theta, double thetaStar, double phiStar) const;
	double pPrimek(double E3, double theta, double thetaStar, double phiStar) const;
	double q12CM(double E3) const;
	double q10Star(double E3) const;
	double q20Star(double E3) const;
	double q11Star(double E3, double thetaStar, double phiStar) const;
	double q21Star(double E3, double thetaStar, double phiStar) const;
	double q12Star(double E3, double thetaStar, double phiStar) const;
	double q22Star(double E3, double thetaStar, double phiStar) const;
	double q13Star(double E3, double thetaStar) const;
	double q23Star(double E3, double thetaStar) const;
	double q10(double E3, double thetaStar) const;
	double q20(double E3, double thetaStar) const;
	double q11(double E3, double theta, double thetaStar, double phiStar) const;
	double q21(double E3, double theta, double thetaStar, double phiStar) const;
	double q12(double E3, double thetaStar, double phiStar) const;
	double q22(double E3, double thetaStar, double phiStar) const;
	double q13(double E3, double theta, double thetaStar, double phiStar) const;
	double q23(double E3, double theta, double thetaStar, double phiStar) const;
	double q3(double E3) const;
	double q30(double E3, double theta, double thetaStar, double phiStar) const;
	double q31(double E3, double theta, double thetaStar, double phiStar) const;
	double q32(double E3, double theta, double thetaStar, double phiStar) const;
	double q33(double E3, double theta, double thetaStar, double phiStar) const;
	double ampSquared(double E3, double theta, double thetaStar, double phiStar) const;
	double dipoleFormFactor(double E3, double theta, double thetaStar, double phiStar) const;
};


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

double Event::simpsonCalc1D(double a, double b, int n, double fixedVariable) const
{
	double h = (b - a) / n; 
	double integral{ 0 };
	for(size_t i{0}; static_cast<int>(i) <= n; i++)
	{
		double x{ a + static_cast<double>(i)*h };
		double wx{ static_cast<double>((i == 0 || i == n) ? 1 : (i % 3 == 0 ? 2 : 3)) };
		integral += wx * thetaPrimeIntegrand(x, fixedVariable);
	}
	return (3.0/8.0) * h * integral;
}

void Event::weightModify(int& weight, int i, int n) const
{
	if(i == 0 || i == n)
	{
		weight *= 1;
	}
	else if(i % 3 == 0)
	{	
		weight *= 2;
	}
	else
	{
		weight *= 3;
	}
}

double Event::simpsonCalc4D(double a, double b, double c, double d, double e, double f, double g, double h, int n, int m,	
														int p, int q) 
{
	double hx{ (b - a) / (n) };
	double hy{ (d - c) / (m) };
	double hz{ (f - e) / (p) };
	double hw{ (h - g) / (q) };

	double integral{ 0 };

	for(size_t i{0}; static_cast<int>(i) <= n; i++)
	{
		double x{ a + static_cast<double>(i)*hx };
		double wx{ static_cast<double>((i == 0 || i == n) ? 1 : (i % 3 == 0 ? 2 : 3)) };
		for(size_t j{0}; static_cast<int>(j) <= m; j++)
		{
			double y{ c + static_cast<double>(j)*hy };
			double wy{ static_cast<double>((j == 0 || j == m) ? 1 : (j % 3 == 0 ? 2 : 3)) };
			for(size_t k{0}; static_cast<int>(k) <= p; k++)
			{
				double z{ e + static_cast<double>(k)*hz };
				double wz{ static_cast<double>((k == 0 || k == p) ? 1 : (k % 3 == 0 ? 2 : 3)) };
				for(size_t l{0}; static_cast<int>(l) <= q; l++)
				{
					double w{ g + static_cast<double>(l)*hw };
					double ww{ static_cast<double>((l == 0 || l == q) ? 1 : (l % 3 == 0 ? 2 : 3)) };
					
					integral += wx * wy * wz * ww * integrand(x, y, z, w);

					// dalitz plot stuff does not work concurrently
					/*
					auto round = [](double value, double precision)
					{
						return std::round(value / precision) * precision;
				  };
					double m12{ round(m12Calc(x)*m12Calc(x), 1e-3) };
					double m23{ round(m23Calc(x, z)*m23Calc(x, z), 1e-3) };
					double m13{ round(m13Calc(x, z)*m13Calc(x, z), 1e-3) };
					std::pair<double,double> m12m23 = std::make_pair(m12, m23);
					std::pair<double,double> m12m13 = std::make_pair(m12, m13);
					std::pair<double,double> m23m13 = std::make_pair(m23, m13);
					m_m12m23Count[m12m23]++;
					m_m12m13Count[m12m13]++;
					m_m23m13Count[m23m13]++;
					*/
				}
			}
		}
	}

	integral *= std::pow((3.0/8.0), 4) * hx * hy * hz * hw;
	return integral;
}

double Event::simpsonCalc3D(IntegrationVariable fixedParameter, double fixedParameterValue, double a, double b, double c,
														double d, double e, double f, int n, int m, int p) const
{
	double hx{ (b - a) / n };
	double hy{ (d - c) / m };
	double hz{ (f - e) / p };

	double integral{ 0 };

	for(size_t i{0}; static_cast<int>(i) < n; i++)
	{
		double x{ a + static_cast<double>(i)*hx };
		double wx{ static_cast<double>((i == 0 || i == n) ? 1 : (i % 3 == 0 ? 2 : 3)) };
		for(size_t j{0}; static_cast<int>(j) < m; j++)
		{
			double y{ c + static_cast<double>(j)*hy };
			double wy{ static_cast<double>((j == 0 || j == m) ? 1 : (j % 3 == 0 ? 2 : 3)) };
			for(size_t k{0}; static_cast<int>(k) < p; k++)
			{
				double z{ e + static_cast<double>(k)*hz };
				double wz{ static_cast<double>((k == 0 || k == p) ? 1 : (k % 3 == 0 ? 2 : 3)) };
				
				switch(fixedParameter)
				{
				case(IntegrationVariable::E3):
					integral += wx * wy * wz * integrand(fixedParameterValue, x, y, z);
					break;
				case(IntegrationVariable::theta):
					integral += wx * wy * wz * integrand(x, fixedParameterValue, y, z);
					break;
				case(IntegrationVariable::thetaStar):
					integral += wx * wy * wz * integrand(x, y, fixedParameterValue, z);
					break;
				case(IntegrationVariable::phiStar):
					integral += wx * wy * wz * integrand(x, y, z, fixedParameterValue);
					break;
				default:
					std::cout<<"Error: entered wrong integration variable type.\n";
					break;
				}
			}
		}
	}
	integral *= std::pow((3.0/8.0), 3) * hx * hy * hz;
	return integral;
}

double Event::crossSectionCalc(int n, int m)
{
	double E3_max{ E3Max() };//E3Max()
	double E3_min{ m_m3 };//m_m3
	double theta_max{ getConst(Constant::pi) };
	double theta_min{ 0 };
	double thetaStar_max{ getConst(Constant::pi)  };
	double thetaStar_min{ 0 };
	double phiStar_max{ 2*getConst(Constant::pi) };
	double phiStar_min{ 0 };
	double ret{ simpsonCalc4D(E3_min, E3_max, theta_min, theta_max, thetaStar_min, thetaStar_max, phiStar_min, phiStar_max, 
							n, m, m, m) };
	ret *= fluxFactor();
	m_crossSection = ret;
	return ret;
}

double Event::E3DiffCrossSection(double value, int n, int m, int p) const
{
	double theta_max{ getConst(Constant::pi) };
	double theta_min{ 0 };
	double thetaStar_max{ getConst(Constant::pi)  };
	double thetaStar_min{ 0 };
	double phiStar_max{ 2*getConst(Constant::pi) };
	double phiStar_min{ 0 };
	double ret{ simpsonCalc3D(IntegrationVariable::E3, value, theta_min, theta_max, thetaStar_min, thetaStar_max, phiStar_min,
													  phiStar_max, n, m, p) };
	ret *= fluxFactor();
	return ret;
}

double Event::thetaDiffCrossSection(double value, int n, int m, int p) const
{
	double E3_max{ E3Max() };//E3Max()
	double E3_min{ m_m3 };//m_m3
	double thetaStar_max{ getConst(Constant::pi)  };
	double thetaStar_min{ 0 };
	double phiStar_max{ 2*getConst(Constant::pi) };
	double phiStar_min{ 0 };
	double ret{ simpsonCalc3D(IntegrationVariable::theta, value, E3_min, E3_max, thetaStar_min, thetaStar_max,
													  phiStar_min, phiStar_max, n, m, p) };
	ret *= fluxFactor();
	return ret;
}

double Event::thetaStarDiffCrossSection(double value, int n, int m, int p) const
{
	double E3_max{ E3Max() };//E3Max()
	double E3_min{ m_m3 };//m_m3
	double theta_max{ getConst(Constant::pi) };
	double theta_min{ 0 };
	double phiStar_max{ 2*getConst(Constant::pi) };
	double phiStar_min{ 0 };
	double ret{ simpsonCalc3D(IntegrationVariable::thetaStar, value, E3_min, E3_max, theta_min, theta_max,
													  phiStar_min, phiStar_max, n, m, p) };
	ret *= fluxFactor();
	return ret;
}

double Event::phiStarDiffCrossSection(double value, int n, int m, int p) const
{
	double E3_max{ E3Max() };//E3Max()
	double E3_min{ m_m3 };//m_m3
	double theta_max{ getConst(Constant::pi) };
	double theta_min{ 0 };
	double thetaStar_max{ getConst(Constant::pi)  };
	double thetaStar_min{ 0 };
	double ret{ simpsonCalc3D(IntegrationVariable::phiStar, value, E3_min, E3_max, theta_min, theta_max,
													  thetaStar_min, thetaStar_max, n, m, p) };
	ret *= fluxFactor();
	return ret;
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
	//return (Ep1()*q20(E3, thetaStar) - std::sqrt(Ep1()*Ep1() - m_ma*m_ma)*q23(E3, theta, thetaStar, phiStar));
	// q3 = lepton => use q3
	return (Ep1()*q30(E3, theta, thetaStar, phiStar) - std::sqrt(Ep1()*Ep1() - m_ma*m_ma)*q33(E3, theta, thetaStar, phiStar));
}

double Event::pPrimekPrime(double E3, double theta, double thetaStar, double phiStar) const
{
	// q3 = kaon => use q2
	//return (q10(E3, thetaStar)*q20(E3, thetaStar) - 
	//						q11(E3, theta, thetaStar, phiStar)*q21(E3, theta, thetaStar, phiStar) -
	//						q12(E3, thetaStar, phiStar)*q22(E3, thetaStar, phiStar) -
	//						q13(E3, theta, thetaStar, phiStar)*q23(E3, theta, thetaStar, phiStar));
	//p3 = lepton => use q3
	return (q10(E3, thetaStar)*q30(E3, theta, thetaStar, phiStar) - 
							q11(E3, theta, thetaStar, phiStar)*q31(E3, theta, thetaStar, phiStar) -
							q12(E3, thetaStar, phiStar)*q32(E3, theta, thetaStar, phiStar) -
							q13(E3, theta, thetaStar, phiStar)*q33(E3, theta, thetaStar, phiStar));
}

double Event::kkPrime(double E3, double theta, double thetaStar, double phiStar) const
{
	// q3 = kaon => use q2
	//return (q20(E3, thetaStar)*Ep2() + q23(E3, theta, thetaStar, phiStar)*Ep2());
	// q3 = lepton => use q3
	return (q30(E3, theta, thetaStar, phiStar)*Ep2() + q33(E3, theta, thetaStar, phiStar)*Ep2());
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

namespace mat
{
	namespace con
	{
		constexpr double D{ 0.804 };
		constexpr double F{ 0.463 };
		double GF{ Event::getConst(Event::Constant::GF) };
		double fPi{ Event::getConst(Event::Constant::fpi) };
		double Vus{ Event::getConst(Event::Constant::Vus) };
	}
	
	//      a   b    c    d     e    f
	// func(pk, pp', pk', p'k', kk', p'k)

	double contactTerm(double a, double b, double c, double d, double e, double f)
	{
		constexpr double ACT{ 1 };
		constexpr double BCT{ -con::D - con::F };
		double m1{ Event::getConst(Event::Constant::neutronMass) };		// neutron -> proton
	  double m2{ Event::getConst(Event::Constant::protonMass) };
		double front{ 0.0625 * (con::GF)*(con::GF) * ACT*ACT * (con::Vus)*(con::Vus) * (1/((con::fPi)*(con::fPi))) };
		return front * ( 64*(a*d + f*c) + 64*BCT*(a*d + c*f + b*e - (3*(e*(b + m1*m2)))) );
	}
}

void write(const char* fileName, std::string_view data)
{
	std::ofstream outFile(fileName);
	if(outFile)
	{
		outFile<<std::scientific<<data;
		outFile.close();
		std::cout<<"Data written successfully\n";
	}
	else
	{
		std::cerr<<"Error opening file for writing\n";
	}
}

double Event::thetaLab(double theta) const
{
	double integral{ simpsonCalc1D(E3Min(), E3Max(), 99, theta) };
	double diff{ E3Max() - E3Min() };
	return (1/diff) * integral;
}

double Event::dThetadThetaLab(double a, double b) const
{
	double dTheta{ b - a };
	double aLab{ thetaLab(a) };
	double bLab{ thetaLab(b) };
	double dThetaLab{ bLab - aLab };
	return dTheta / dThetaLab;
}

double Event::dthetadthetaprime(double theta, double Eoverq) const
{
	double gamma{ Ep1() / m_ma };
	assert(gamma > 1);
	double beta{ std::sqrt(1 - 1/(gamma*gamma)) };
	return gamma*(1 - beta*Eoverq*std::cos(theta));
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

void crossSection()
{
	Timer t;

	// CT electron neutrino	+ neutron -> proton + muon	
	double (*func)(double, double, double, double, double, double);
	func = mat::contactTerm;	
	// q3 = kaon
	//Event event{ 0, Event::getConst(Event::Constant::neutronMass), 0, Event::getConst(Event::Constant::neutronMass),
	//				 Event::getConst(Event::Constant::muonMass), Event::getConst(Event::Constant::kaonChargedMass), 3.5, func };
	// q3 = lepton
	Event event{ 0, Event::getConst(Event::Constant::neutronMass), 0, Event::getConst(Event::Constant::protonMass),
					 Event::getConst(Event::Constant::kaon0Mass), Event::getConst(Event::Constant::muonMass), 3.5, func };
	// differential cross section 
	
	double diffMin{ 0 /*event.E3Min()*/ }; 
	double diffMax{ Event::getConst(Event::Constant::pi) /*event.E3Max()*/ };
	constexpr size_t diffN{ 300 };
	double diffStep{ (diffMax - diffMin) / (diffN - 1) };
	std::array<double, diffN> diffCrossSections;
	std::array<double, diffN> fixedParameters;
	constexpr double Eoverq{ 1.0 };
	{
		unsigned int numThreads{ std::thread::hardware_concurrency() };
		ThreadPool pool{ numThreads };
		std::cout<<"Threads: "<<numThreads<<'\n';
		for(size_t i{ 0 }; i < diffN; i++)
		{
			pool.enqueue([&, i]()
			{
				double fixedParameter{ diffMin + i*diffStep };
				fixedParameter = event.thetaPrime(fixedParameter, Eoverq);
				//fixedParameter = event.EPrime(fixedParameter);
				//fixedParameter = event.thetaLab(fixedParameter);
				//double boostedFixedParameter{ event.thetaLab(fixedParameter) }; // for dtheta
				double diffCrossSection{ std::fabs(event.thetaDiffCrossSection(fixedParameter, 36, 36, 36))/* / 2*diffMax*/ };
				//double prev{ fixedParameter - diffStep };
				//double next{ fixedParameter + diffStep };
				//double lorentzFactor{ event.dThetadThetaLab(prev, next) };  // for dtheta
				diffCrossSections[i] = diffCrossSection * 3.89739e-28 / (2*diffMax); // cm^2
				fixedParameters[i] = fixedParameter + diffMax/2;
			});
		}
	}
	std::string diffData;
	event.serialize(diffCrossSections, fixedParameters, diffData);
	std::cout<<diffData<<'\n';
	write("CTNP-l-btheta.txt", diffData);
/*
	// cross sections
	
	// elecron: 2.066953, muon: 2.38085
	constexpr double sMin{ 2.5 };  //2.25(Event::getConst(Event::neutronMass))*(Event::getConst(Event::neutronMass)) };
	constexpr double sMax{ 7.5 };	// 7.2
	constexpr size_t N{ 30 };
	constexpr double step{ (sMax - sMin) / (N-1) };
	std::array<double, N> crossSections;
	std::array<double, N> beamEnergies;
	{
		unsigned int numThreads{ std::thread::hardware_concurrency() };
		ThreadPool pool{ numThreads };
		std::cout<<"Threads: "<<numThreads<<'\n';
		for(size_t i{ 0 }; i < N; i++)
		{
			double s{ sMin + i*step };
			event.setS(s);
			event.crossSectionCalc(3, 3);
			if(auto cs = event.getCrossSection())
				crossSections[i] = *cs * 3.89739e-28; // cm^2
			else
			{
				crossSections[i] = NAN;
				std::cout<<"no cs\n";
			}
			beamEnergies[i] = event.getBeamEnergy();
		}
			
			double s{ sMin + i*step };
			pool.enqueue([&, i, s]()
			{
				Event e{ event };
				e.setS(s);
				e.crossSectionCalc(36, 36);
				if(auto cs = e.getCrossSection())
					crossSections[i] = *cs * 3.89739e-28; // cm^2
				else
				{
					crossSections[i] = NAN;
					std::cout<<"no cs\n";
				}
				beamEnergies[i] = e.getBeamEnergy();
			});
		}
		
	}
	std::string data;
	event.serialize(crossSections, beamEnergies, data);
	std::cout<<data<<'\n';
	write("CTNP-f.txt", data);	
	
*/	
	// dalitz
/*	
	event.crossSectionCalc(99, 99);
	std::string data1;
	event.serializeM12M23(data1);
	write("CTNPm12m23-long.txt", data1);

	std::cout<<"Time elapsed: "<<t.elapsed()<<'\n';
	*/
//---------------------------------------------------------------------------------------------------------------------------	
}
