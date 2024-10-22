#include<iostream>
#include<vector>


double simpsonCalc(std::vector<double> xData, double a, double b, int n, double (*func)(double))
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

void simpson()
{
	int n = 1000000;
	std::vector<double> xVals;
	for(int i = 0; i <= 100; i++)
		xVals.push_back(i);

	double integral = simpsonCalc(xVals, 0, 100, n, [](double x){return x*x;});
	
	std::cout<<integral<<'\n';

	return;
}
