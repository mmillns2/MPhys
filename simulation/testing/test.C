#include"fourVector.h"
#include<iostream>

void test()
{
	FourVector v1{1, 1, 1, 1};
	FourVector v2{{2, 2, 2, 2}};
	std::cout<<"v1: "<<v1<<"v2: "<<v2<<'\n';

	FourVector v3 = v1 + v2;
	std::cout<<"v3: "<<v3<<'\n';

	std::cout<<"v1: "<<v1<<"v2: "<<v2<<'\n';

	std::cout<<"length v2: "<<v2.length()<<'\n';
	
	std::cout<<"v1 dot v2: "<<v1.dotProduct(v2)<<'\n';

	double length2{v2.length()};
	std::cout<<"length: "<<length2<<'\n';
}
