#include"fourVector.h"
#include"timer.h"
#include"threadPool.h"

#include<memory>
#include<cmath>
#include<functional>
#include<array>
#include<optional>
#include<limits>
#include<fstream>
#include<iomanip>


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
			m_crossSection{ std::numeric_limits<double>::min() },           // -inf => not been calculated
			m_beamEnergy{ std::sqrt(s) - ma }
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
			m_crossSection{ std::numeric_limits<double>::min() },           // -inf => not been calculated
			m_beamEnergy{ std::sqrt(s) - Ep() }
	{ }

	Event(const Event& e)
		:	m_id{ e.m_id },
			m_p1{ std::make_unique<FourVector>(*(e.m_p1)) },
			m_p2{ std::make_unique<FourVector>(*(e.m_p2)) },
			m_p3{ std::make_unique<FourVector>(*(e.m_p3)) },
			m_ma{ e.m_ma },
			m_mb{ e.m_mb },
			m_m1{ e.m_m1 },
			m_m2{ e.m_m2 },
			m_m3{ e.m_m3 },
			m_s{ e.m_s },
			m_ampSquared{ e.m_ampSquared },
			m_crossSection{ e.m_crossSection },
			m_beamEnergy{ e.m_beamEnergy }
	{ }

	~Event(){}

public:
	std::optional<double> getCrossSection() const {return m_crossSection == std::numeric_limits<double>::min() ? std::nullopt :
																								 std::optional(m_crossSection);}
	double getBeamEnergy() const {return m_beamEnergy;}
	static double getConst(Constant c);
	double crossSectionCalc(int n, int m);
	void serializeEvent(std::string& out) const;

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

public:
	double fluxFactor() const;
	double simpsonCalc1D(double a, double b, int n, double (*func)(double)) const;
	double kallen(double x, double y, double z) const;
	double m12Calc(double E3) const;
	double integrand(double E3, double theta, double thetaStar, double phiStar) const; 
	double E3Max() const;
	void weightModify(int& weight, int i, int n) const;
	double simpsonCalc4D(double a, double b, double c, double d, double e, double f, double g, double h, int n, int m,
											 int p, int q) const; 
	double Ep() const;
	double cosAlpha(double theta, double thetaStar, double phiStar) const;
	double qSquared(double theta, double thetaStar, double phiStar) const;
	double pk() const;
	double ppPrime(double theta, double thetaStar, double phiStar) const;
	double pkPrime(double theta, double thetaStar, double phiStar) const;
	double pPrimekPrime() const;
	double kkPrime(double theta, double thetaStar, double phiStar) const;
	double pPrimek(double theta, double thetaStar, double phiStar) const;
	double ampSquared(double theta, double thetaStar, double phiStar) const;
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
	return std::sqrt(std::pow(x, 4) + std::pow(y, 4) + std::pow(z, 4) - 2*std::pow(x, 2)*std::pow(y, 2) 
				 - 2*std::pow(x, 2)*std::pow(z, 2) - 2*std::pow(y, 2)*std::pow(z, 2)); 	
}

double Event::m12Calc(double E3) const
{
	return std::sqrt(m_s + std::pow(m_m3, 2) - 2*(std::sqrt(m_s))*E3); 	
}

double Event::integrand(double E3, double theta, double thetaStar, double phiStar) const
{
	//std::cout<<"E3: "<<E3<<'\n';
	double m12{ m12Calc(E3) };
	double phaseSpaceFactor{ kallen(m12, m_m1, m_m2) / (m12*m12) };
	//std::cout<<"m12: "<<m12<<", phase space factor: "<<phaseSpaceFactor<<'\n';
	//std::cout<<"matrix element: "<<ampSquared(theta, thetaStar, phiStar)<<'\n';
	return std::sqrt(std::pow(E3, 2) - std::pow(m_m3, 2)) * phaseSpaceFactor * ampSquared(theta, thetaStar, phiStar) * 
				 std::sin(theta) * std::sin(thetaStar);
}

double Event::fluxFactor() const
{
	//double p1Dotp2{ m_p1->dotProduct(*m_p2) };
	double p1Dotp2{ pk() };
	return 1/(64*std::pow(2*getConst(Constant::pi), 4) * std::sqrt(std::pow(p1Dotp2, 2) - std::pow(m_ma, 2)*std::pow(m_mb, 2)));
}

double Event::E3Max() const
{
	double q3Max{ kallen(std::sqrt(m_s), m_m1 + m_m2, m_m3) / (2*std::sqrt(m_s)) };
	//std::cout<<"q3 max: "<<q3Max<<'\n';
	return std::sqrt(std::pow(q3Max, 2) + std::pow(m_m3, 2));
	//return std::sqrt(m_s) - Ep() + m_m3;
}

double Event::simpsonCalc1D(double a, double b, int n, double (*func)(double)) const
{
	double ret = func(b) + func(a);
	double h = (b - a) / n; 
	for(size_t i{1}; static_cast<int>(i) < n - 1; i++)
	{
		double k = a + static_cast<double>(i)*h;
		if(i%3==0)
		{
			ret += 2*func(k);
		}
		else
		{
			ret += 3*func(k);
		}
	}
	return (3.0/8.0) * h * ret;
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
														int p, int q) const
{
	double hx{ (b - a) / n };
	double hy{ (d - c) / m };
	double hz{ (f - e) / p };
	double hw{ (h - g) / q };

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
				for(size_t l{0}; static_cast<int>(l) < q; l++)
				{
					double w{ g + static_cast<double>(l)*hw };
					double ww{ static_cast<double>((l == 0 || l == q) ? 1 : (l % 3 == 0 ? 2 : 3)) };
					
					integral += wx * wy * wz * ww * integrand(x, y, z, w);
				}
			}
		}
	}

	integral *= std::pow((3.0/8.0), 4) * hx * hy * hz * hw;
	//std::cout<<"Integral: "<<integral<<'\n';
	return integral;
}

double Event::crossSectionCalc(int n, int m)
{
	double E3_max{ E3Max() };
	double E3_min{ m_m3 };
	double theta_max{ getConst(Constant::pi) };
	double theta_min{ 0 };
	double thetaStar_max{ getConst(Constant::pi)  };
	double thetaStar_min{ 0 };
	double phiStar_max{ 2*getConst(Constant::pi) };
	double phiStar_min{ 0 };

	//std::cout<<"E3 max: "<<E3_max<<", E3 min: "<<E3_min<<'\n';	

	double ret{ simpsonCalc4D(E3_min, E3_max, theta_min, theta_max, thetaStar_min, thetaStar_max, phiStar_min, phiStar_max, 
							n, m, m, m) };

	ret *= fluxFactor();

	m_crossSection = ret;
	return ret;
}

void Event::serializeEvent(std::string& out) const
{
	// beamEnergy crossSection
	out.append(std::to_string(m_beamEnergy));
	out.append(" ");

	if(auto cs = getCrossSection())
	{
		std::ostringstream stream;
		double temp{ *cs * 3.894e-27 };	// cm^2
		stream<<std::scientific<<std::setprecision(10)<<temp;
		out.append(stream.str());
	}
	
	else
		out.append("error");

	out.append("\n");
}

double Event::Ep() const
{
	//std::cout<<"ma: "<<m_ma<<'\n';
	//std::cout<<"Ep: "<< (m_s + m_ma*m_ma) / (2*std::sqrt(m_s))<<'\n';
	return (m_s + m_ma*m_ma) / (2*std::sqrt(m_s));
}

double Event::cosAlpha(double theta, double thetaStar, double phiStar) const
{
	//return 0.5*(std::cos(phiStar) + (std::cos(phiStar) + 1)*std::cos(theta + thetaStar + getConst(Constant::pi)) - 1);
	return std::cos(phiStar); 
	//return std::cos(getConst(Constant::pi) - theta - thetaStar);
}

double Event::qSquared(double theta, double thetaStar, double phiStar) const
{
	double EpSqrd{ Ep()*Ep() };
	double mNSqrd{ m_ma*m_ma };
	double mNPrimeSqrd{ m_m1*m_m1 };
	double cosA{ cosAlpha(theta, thetaStar, phiStar) };
	double ret{ mNSqrd + mNPrimeSqrd - 2*(-std::sqrt(EpSqrd - mNSqrd) * std::sqrt(EpSqrd - mNPrimeSqrd)*cosA + EpSqrd) };
	//std::cout<<"q sqrd: "<<ret<<'\n';
	return ret;
}

double Event::pk() const
{
	return 0.5*(m_s - m_mb*m_mb - m_ma*m_ma);
}

double Event::ppPrime(double theta, double thetaStar, double phiStar) const
{
	return 0.5*(-qSquared(theta, thetaStar, phiStar) + 2*m_ma*m_ma);
}

double Event::pkPrime(double theta, double thetaStar, double phiStar) const
{
	return 0.5*(m_s + qSquared(theta, thetaStar, phiStar) - m_ma*m_ma - m_mb*m_mb);
}

double Event::pPrimekPrime() const
{
	return 0.5*(m_s - m_ma*m_ma - m_mb*m_mb);
}

double Event::kkPrime(double theta, double thetaStar, double phiStar) const
{
	return 0.5*(-qSquared(theta, thetaStar, phiStar) + 2*m_mb*m_mb);
}

double Event::pPrimek(double theta, double thetaStar, double phiStar) const
{
	return 0.5*(m_s + qSquared(theta, thetaStar, phiStar) - m_ma*m_ma - m_mb*m_mb);
}

double Event::ampSquared(double theta, double thetaStar, double phiStar) const
{
	double a{ pk() };
	double b{ ppPrime(theta, thetaStar, phiStar) };
	double c{ pkPrime(theta, thetaStar, phiStar) };
	double d{ pPrimekPrime() };
	double e{ kkPrime(theta, thetaStar, phiStar) };
	double f{ pPrimek(theta, thetaStar, phiStar) };
	//std::cout<<"p'k: "<<f<<'\n';
	//std::cout<<a<<", "<<b<<", "<<c<<", "<<d<<", "<<e<<", "<<f<<'\n';
	if(m_ampSquared)
	{
		//std::cout<<m_ampSquared(a, b, c, d, e, f)<<'\n';
		return m_ampSquared(a, b, c, d, e, f);
	}
	else
		return 0;
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
		constexpr double BCT{ con::D - con::F };
		double m1{ Event::getConst(Event::Constant::neutronMass) };		// neutron -> neutron
	  double m2{ Event::getConst(Event::Constant::protonMass) };
		double front{ 0.0625 * (con::GF)*(con::GF) * ACT*ACT * (con::Vus)*(con::Vus) * (1/((con::fPi)*(con::fPi))) };
		return front * ( 64*(a*d + f*c) + 64*BCT*(a*d + c*f + b*e - 3*(e*(b + m1*m2))) );
	}
}

void crossSection()
{
	Timer t;

	int id{ 0 };
	auto makeEvent2 = [&id](double ma, double mb, double m1, double m2, double m3, double s, 
													double (*ampSquared)(double, double, double, double, double, double))
	{
		Event e{ id, ma, mb, m1, m2, m3, s, ampSquared };
		id++;
		return e;
	};

//--------------------------------------------------------------------------------------------------------------------------//		
	// CT electron neutrino	+ neutron -> proton + electron, s = 3.25GeV^2	
	double (*func)(double, double, double, double, double, double);
	func = mat::contactTerm;	
	/*
	// single event
	Event e = makeEvent2(0.9395654, 0, 0.938272, 
											 0.10566, 0.4947648, 10.7, func);
	e.crossSectionCalc(1000, 20);
	std::string data;
	e.serializeEvent(data);
	std::cout<<data<<'\n';
	*/
	// elecron: 2.066953
	constexpr double sMin{ 2.38085 };  //(Event::getConst(Event::neutronMass))*(Event::getConst(Event::neutronMass)) };
	constexpr double sMax{ 16 };
	constexpr size_t N{ 50 };
	constexpr double step{ (sMax - sMin) / N };
	std::array<double, N> crossSections;
	std::array<double, N> beamEnergies;
	std::string data;
	{
		ThreadPool pool{ 50 };
		std::mutex updateDataMutex;
		for(size_t i{ 0 }; i < N; i++)
		{
			double s{ sMin + i*step };
			pool.enqueue([&, i, s]()
			{
				Event e = makeEvent2(0.9395654, 0, 0.938272, 
												 		 0.10566, 0.4947648, s, func);
				//std::cout<<"made event\n";
				e.crossSectionCalc(100, 100);
				if(auto cs = e.getCrossSection())
				{
					std::lock_guard<std::mutex> lock(updateDataMutex);
					crossSections[i] = *cs;
					beamEnergies[i] = e.getBeamEnergy();
					e.serializeEvent(data);
				}
			});
		}
	}
	
	std::cout<<data<<'\n';
	std::cout<<"here\n";
	//for(size_t i{ 0 }; i < N; i++)
	//{
	//	std::cout<<beamEnergies[i]<<" "<<crossSections[i]<<'\n';
	//}
	
	std::ofstream outFile("CT.txt");
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

	std::cout<<"Time elapsed: "<<t.elapsed()<<'\n';
	
//---------------------------------------------------------------------------------------------------------------------------	
}
