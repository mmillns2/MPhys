#pragma once

#include<array>
#include<iostream>

class FourVector
{
public:
	FourVector() = default;

	FourVector(double E, double p1, double p2, double p3)
	{
		setE(E);
		setP1(p1);
		setP2(p2);
		setP3(p3);
	}

	FourVector(std::array<double, 4> vec) : m_vec{vec} {}

	~FourVector(){}

	double getE() const {return m_vec[0];}
	double getP1() const {return m_vec[1];}	
	double getP2() const {return m_vec[2];}
	double getP3() const {return m_vec[3];}

	void setE(const double E) {m_vec[0] = E;}
	void setP1(const double p1) {m_vec[1] = p1;}
	void setP2(const double p2) {m_vec[2] = p2;}
	void setP3(const double p3) {m_vec[3] = p3;}

	double dotProduct(const FourVector& vec) const
	{
		double ret{0};
		ret += getE() * vec.getE();
		ret -= getP1() * vec.getP1();
		ret -= getP2() * vec.getP2();
		ret -= getP3() * vec.getP3();
		return ret;
	}

	double length() const
	{
		double ret{0};
		ret += (getE() * getE());
		ret -= (getP1() * getP1());
		ret -= (getP2() * getP2());
		ret -= (getP3() * getP3());
		return ret;
	}

	FourVector& operator+=(const FourVector& rhs)
	{
		this->setE(this->getE() + rhs.getE());
		this->setP1(this->getP1() + rhs.getP1());
		this->setP2(this->getP2() + rhs.getP2());
		this->setP3(this->getP3() + rhs.getP3());
		return *this;
	}

	FourVector operator+(const FourVector& rhs)
	{
		FourVector ret;
		ret.setE(this->getE() + rhs.getE());
		ret.setP1(this->getP1() + rhs.getP1());
		ret.setP2(this->getP2() + rhs.getP2());
		ret.setP3(this->getP3() + rhs.getP3());
		return ret;
	}

	friend std::ostream& operator<<(std::ostream& out, const FourVector& vec)
	{
		out<<"4-momentum: ("<<vec.getE()<<", "<<vec.getP1()<<", "<<vec.getP2()<<", "<<vec.getP3()<<")\n";
		return out;
	}

private:
	std::array<double, 4> m_vec;
};
