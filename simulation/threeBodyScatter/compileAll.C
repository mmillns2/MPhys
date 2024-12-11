#include"crossSection.h"

#include"crossSection.C"
#include"crossSectionCalculations.C"
#include"integration.C"
#include"momentaComponents.C"
#include"matrixElements.C"


void write(const char* fileName, std::string_view data)
{
	std::ofstream outFile(fileName);
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
}

void compileAll()
{
	Timer t;

	double (*func)(double, double, double, double, double, double);
	func = mat::contactTerm;	
	// q3 = kaon
	//Event event{ 0, Event::getConst(Event::Constant::neutronMass), 0, Event::getConst(Event::Constant::neutronMass),
	//				 Event::getConst(Event::Constant::muonMass), Event::getConst(Event::Constant::kaonChargedMass), 3.5, func };
	// q3 = lepton
	Event event{ 0, Event::getConst(Event::Constant::neutronMass), 0, Event::getConst(Event::Constant::protonMass),
					 Event::getConst(Event::Constant::kaon0Mass), Event::getConst(Event::Constant::muonMass), 3.5, func };

	// cross sections
	
	constexpr double sMin{ 2.5 };  
	constexpr double sMax{ 7.5 };
	constexpr size_t N{ 30 };
	constexpr double step{ (sMax - sMin) / (N-1) };
	std::array<double, N> crossSections;
	std::array<double, N> beamEnergies;
	{
		unsigned int numThreads{ std::thread::hardware_concurrency() };
		ThreadPool pool{ numThreads };
		std::cout<<"Threads: "<<numThreads<<'\n';
		for(size_t i{ 0 }; i < N; i++)
		{
			double s{ sMin + i*step };
			event.setS(s);
			event.crossSectionCalc(3, 3);
			if(auto cs = event.getCrossSection())
				crossSections[i] = *cs * 3.89739e-28; // cm^2
			else
			{
				crossSections[i] = NAN;
				std::cout<<"no cs\n";
			}
			beamEnergies[i] = event.getBeamEnergy();
		}
	}
	std::string data;
	event.serialize(crossSections, beamEnergies, data);
	std::cout<<data<<'\n';
	write("test.txt", data);	
	
	std::cout<<"Time elapsed: "<<t.elapsed()<<'\n';
}
