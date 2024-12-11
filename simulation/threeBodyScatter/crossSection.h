#pragma once

#include"threadPool.h"
#include"timer.h"

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

		Event(int id, double ma, double mb, double m1, double m2, double m3, double s, double (*ampSquaredPtr)(double, double,
																																												double, double, double, double))
		: m_id{ id },
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
	double getBeamEnergy() const; 
	void setS(double s);
	static double getConst(Constant c);
	double crossSectionCalc(int n, int m);
	double E3DiffCrossSection(double value, int n, int m, int p) const;
	double thetaDiffCrossSection(double value, int n, int m, int p) const;
	double thetaStarDiffCrossSection(double value, int n, int m, int p) const;
	double phiStarDiffCrossSection(double value, int n, int m, int p) const;
	void serializeEvent(std::string& out) const;
	template<size_t N>
	void serialize(const std::array<double, N>& yVals, const std::array<double, N>& xVals,
																	 std::string& data) const;
	double E3Max() const;
	double E3Min() const { return m_m3; }

	double thetaPrime(double theta, double Eoverq) const;
	double EPrime(double E) const;

private:
	int m_id;
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
