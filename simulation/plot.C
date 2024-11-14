#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TApplication.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>



//  Function to read whitespace-delimited data from a text file
std::vector<std::vector<double>> readDataFromTxt(const std::string& filename) 
{
  std::vector<std::vector<double>> ret;
  std::ifstream file(filename);

  // Check if the file was opened successfully
  if(!file.is_open()) 
  {                         
    std::cerr<<"Error: Could not open file "<<filename<<std::endl;
    return ret;                                     
  }
                                         
  std::string line;
  while(std::getline(file, line)) 
  {
    std::istringstream iss(line);
		std::string energyStr;
		std::string csStr;
		iss>>energyStr;
		iss>>csStr;
		double energy{ std::stod(energyStr) };
		double cs{ std::stod(csStr) };
		std::vector<double> row{ energy, cs };

    // Add the row to the data vector if its not empty
    if (!row.empty())
    {
      ret.push_back(row);
    }
  }
  file.close();
  return ret;
} 


void plot()
{
  // Import data from txt file (whitespace delimiter)
  std::string filename{ "CT.txt" };
  std::vector<std::vector<double>> data{ readDataFromTxt(filename) };
	const size_t n{ data.size() };

  // Separate each column into individual vectors
	double energyVals[n];
	double csVals[n];
	for(size_t i{ 0 }; i < n; i++)
	{
		energyVals[i] = data[i][0];
		csVals[i] = data[i][1];
		std::cout<<energyVals[i]<<" "<<csVals[i]<<'\n';
	}

  // Create a canvas
  TCanvas *c1 = new TCanvas("c1","Our Plot", 2048, 1536);
 
  c1->SetGrid();
  c1->GetFrame()->SetBorderSize(12);
 
  //  Instantiate a graph with errors
  TGraphErrors *graph = new TGraphErrors(n, energyVals, csVals, nullptr, nullptr);

  //  Change how the plot looks
  graph->SetTitle("CT neutron -> proton");
  graph->SetMarkerColor(kBlue);
  graph->SetMarkerStyle(kOpenCircle);
	graph->SetLineColor(kBlue);
	graph->GetXaxis()->SetRangeUser(0.5, 2.5);
 
  //  Plot the graph
  graph->DrawClone("APE");
 
  c1->Update();
  c1->Print("CTPlot.png");
}
