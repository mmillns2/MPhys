#include"crossSection.h"


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
