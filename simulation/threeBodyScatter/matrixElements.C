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
		double m1{ Event::getConst(Event::Constant::neutronMass) };		// neutron -> proton
	  double m2{ Event::getConst(Event::Constant::neutronMass) };
		double front{ 0.0625 * (con::GF)*(con::GF) * ACT*ACT * (con::Vus)*(con::Vus) * (1/((con::fPi)*(con::fPi))) };
		return front * ( 64*(a*d + f*c) + 64*BCT*(a*d + c*f + b*e - (3*(e*(b + m1*m2)))) );
	}
}
