#include"crossSection.h"


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
