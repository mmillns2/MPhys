#include"fourVector.h"

#include<memory>
#include<cmath>
#include<functional>


namespace Constants
{
	double pi{3.14159265359};
}


class Event
{
public:
	Event() = default;

	Event(double ma, double mb, double m1, double m2, double m3, double s, double (*ampSquared)(double, double, double, double))
		: m_ma{ma}, m_mb{mb}, m_m1{m1}, m_m2{m2}, m_m3{m3}, m_s{s}, m_ampSquared{ampSquared}
	{ }

	~Event(){}

	double crossSection(int n, int m) const;

private:
	std::unique_ptr<FourVector> m_p1;
	std::unique_ptr<FourVector> m_p2;
	std::unique_ptr<FourVector> m_p3;
	double m_ma;
	double m_mb;
	double m_m1;
	double m_m2;
	double m_m3;
	double m_s;
	double (*m_ampSquared)(double, double, double, double);

	double fluxFactor() const;
	double simpsonCalc1D(double a, double b, int n, double (*func)(double)) const;
	double kallen(double x, double y, double z) const;
	double m12Calc(double s, double m3, double E3) const;
	double integrand(double E3, double theta, double thetaStar, double phiStar) const; 
	double E3Max() const;
	void weightModify(int& weight, int i, int n) const;
	double simpsonCalc4D(double a, double b, double c, double d, double e, double f, double g, double h, int n, int m,
											 int p, int q, std::function<double(double, double, double, double)> integrandPtr) const; 
};


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
	return 1/(64*std::pow(2*Constants::pi, 4) * std::sqrt(std::pow(p1Dotp2, 2) - std::pow(m_ma, 2)*std::pow(m_mb, 2)));
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
	for(int i{1}; i < n - 1; i++)
	{
		double k = a + i*h;
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
														int p, int q, std::function<double(double, double, double, double)> integrandPtr) const
{
	double hx{ (b - a) / n };
	double hy{ (d - c) / m };
	double hz{ (f - e) / p };
	double hw{ (h - g) / q };

	double integral{ 0 };

	for(int i{0}; i < n; i++)
	{
		for(int j{0}; j < m; j++)
		{
			for(int k{0}; k < p; k++)
			{
				for(int l{0}; l < q; l++)
				{
					double x{ a + i*hx };
					double y{ c + j*hy };
					double z{ e + k*hz };
					double w{ g + l*hw };

					int weight{ 1 };
					weightModify(weight, i, n);
					weightModify(weight, j, m);
					weightModify(weight, k, p);
					weightModify(weight, l, q);

					integral += weight * integrandPtr(x, y, z, w);
				}
			}
		}
	}

	integral *= (3/8) * hx * hy * hz * hw;
	return integral;
}

double Event::crossSection(int n, int m) const
{
	double E3_max{ E3Max() };
	double E3_min{ m_m3 };
	double theta_max{ Constants::pi };
	double theta_min{ 0 };
	double thetaStar_max{ Constants::pi  };
	double thetaStar_min{ 0 };
	double phiStar_max{ 2*Constants::pi };
	double phiStar_min{ 0 };

	std::function<double(double, double, double, double)> integrandPtr;
	integrandPtr = [this](double a, double b, double c, double d){ return this->integrand(a, b, c, d); };

	double ret{ simpsonCalc4D(E3_min, E3_max, theta_min, theta_max, thetaStar_min, thetaStar_max, phiStar_min, phiStar_max, 
							n, m, m, m, integrandPtr) };
	
	return ret;
}


void crossSection()
{

}
