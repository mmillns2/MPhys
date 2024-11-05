#include"fourVector.h"

#include<memory>
#include<cmath>
#include<functional>
#include<array>
#include<optional>


class Event
{
public:
	enum Constant
	{
		pi,
		muonMass,						// MeV
		electronMass,				// MeV
		protonMass,					// MeV
		neutronMass,				// MeV
		kaonChargedMass,		// MeV
		kaon0MassMass,			// Mev
		max_constant
	};

public:
	Event() = default;

	Event(int id, std::unique_ptr<FourVector> p1, std::unique_ptr<FourVector> p2, std::unique_ptr<FourVector> p3, double ma,
				double mb, double m1, double m2, double m3, double s, double (*ampSquared)(double, double, double, double))
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
			m_ampSquared{ ampSquared },
			m_crossSection{ -1 },           // -1 => not been calculated
			m_beamEnergy{ std::sqrt(s - (ma*ma)) }
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
	std::optional<double> getCrossSection() const {return m_crossSection == -1 ? std::nullopt : std::optional(m_crossSection);}
	static double getConst(Constant c);
	void crossSectionCalc(int n, int m);
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
	double (*m_ampSquared)(double, double, double, double);
	double m_crossSection;
	double m_beamEnergy;

private:
	double fluxFactor() const;
	double simpsonCalc1D(double a, double b, int n, double (*func)(double)) const;
	double kallen(double x, double y, double z) const;
	double m12Calc(double s, double m3, double E3) const;
	double integrand(double E3, double theta, double thetaStar, double phiStar) const; 
	double E3Max() const;
	void weightModify(int& weight, int i, int n) const;
	double simpsonCalc4D(double a, double b, double c, double d, double e, double f, double g, double h, int n, int m,
											 int p, int q) const; 
};


double Event::getConst(Constant c)
{
	static constexpr std::array<double, static_cast<size_t>(Constant::max_constant)> constants{3.14159265359, 0.1134289259,
																																	0.51099895000, 938.27208943, 939.57, 493.677, 497.611}; 
	return constants[static_cast<size_t>(c)];
}

double Event::kallen(double x, double y, double z) const
{
	return std::sqrt(std::pow(x, 4) + std::pow(y, 4) + std::pow(z, 4) - 2*std::pow(x, 2)*std::pow(y, 2) 
				 - 2*std::pow(x, 2)*std::pow(z, 2) - 2*std::pow(y, 2)*std::pow(z, 2)); 	
}

double Event::m12Calc(double s, double m3, double E3) const
{
	return s + std::pow(m3, 2) - 2*(std::sqrt(s))*E3; 	
}

double Event::integrand(double E3, double theta, double thetaStar, double phiStar) const
{
	double m12{ m12Calc(m_s, m_m3, E3) };
	double phaseSpaceFactor{ kallen(m12, m_m1, m_m2) / (m_m1 * m_m2) };
	return std::sqrt(std::pow(E3, 2) - std::pow(m_m3, 2)) * phaseSpaceFactor * m_ampSquared(1.0,1.0,1.0,1.0);
}

double Event::fluxFactor() const
{
	double p1Dotp2{ m_p1->dotProduct(*m_p2) };
	return 1/(64*std::pow(2*getConst(Constant::pi), 4) * std::sqrt(std::pow(p1Dotp2, 2) - std::pow(m_ma, 2)*std::pow(m_mb, 2)));
}

double Event::E3Max() const
{
	double q3Max{ kallen(std::sqrt(m_s), m_m1 + m_m2, m_m3) / (2*std::sqrt(m_s)) };
	return std::sqrt(std::pow(q3Max, 2) + std::pow(m_m3, 2));
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
		for(size_t j{0}; static_cast<int>(j) < m; j++)
		{
			for(size_t k{0}; static_cast<int>(k) < p; k++)
			{
				for(size_t l{0}; static_cast<int>(l) < q; l++)
				{
					double x{ a + static_cast<double>(i)*hx };
					double y{ c + static_cast<double>(j)*hy };
					double z{ e + static_cast<double>(k)*hz };
					double w{ g + static_cast<double>(l)*hw };

					int weight{ 1 };
					weightModify(weight, static_cast<int>(i), n);
					weightModify(weight, static_cast<int>(j), m);
					weightModify(weight, static_cast<int>(k), p);
					weightModify(weight, static_cast<int>(l), q);

					integral += weight * integrand(x, y, z, w);
				}
			}
		}
	}

	integral *= (3/8) * hx * hy * hz * hw;
	return integral;
}

void Event::crossSectionCalc(int n, int m)
{
	double E3_max{ E3Max() };
	double E3_min{ m_m3 };
	double theta_max{ getConst(Constant::pi) };
	double theta_min{ 0 };
	double thetaStar_max{ getConst(Constant::pi)  };
	double thetaStar_min{ 0 };
	double phiStar_max{ 2*getConst(Constant::pi) };
	double phiStar_min{ 0 };

	double ret{ simpsonCalc4D(E3_min, E3_max, theta_min, theta_max, thetaStar_min, thetaStar_max, phiStar_min, phiStar_max, 
							n, m, m, m) };

	m_crossSection = ret;
}

void Event::serializeEvent(std::string& out) const
{
	// beamEnergy crossSection
	out.append(std::to_string(m_beamEnergy));
	out.append(" ");

	if(getCrossSection())
		out.append(std::to_string(m_crossSection));
	
	else
		out.append("error");

	out.append("\n");
}


void crossSection()
{
	auto makeVector = [](double E, double px, double py, double pz)
	{
		FourVector v{E, px, py, pz};
		return make_unique<FourVector>(v);
	};

	int id{ 0 };

	auto makeEvent = [&makeVector](int& id,
											double E1, double px1, double py1, double pz1, 
											double E2, double px2, double py2, double pz2,
											double E3, double px3, double py3, double pz3,
											double ma, double mb, double m1, double m2, double m3, double s,
											double (*ampSquared)(double, double, double, double))
	{
		Event e{id,
						makeVector(E1, px1, py1, pz1),
						makeVector(E2, px2, py2, pz2),
						makeVector(E3, px3, py3, pz3),
						ma, mb, m1, m2, m3, s, ampSquared};

		id++;
		return e;
	};

	
	
	
}
