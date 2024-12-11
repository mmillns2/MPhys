#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TH2D.h>
#include <TStyle.h>
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
		if(csStr == "nan" || csStr == "-nan")
			continue;
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

std::vector<std::vector<double>> readDataFromTxt2(const std::string& filename) 
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
		std::string mij1Str;
		std::string mij2Str;
		std::string countStr;
		iss>>mij1Str;
		iss>>mij2Str;
		iss>>countStr;
		double mij1{ std::stod(mij1Str) };
		double mij2{ std::stod(mij2Str) };
		double count{ std::stod(countStr) };
		std::vector<double> row{ mij1, mij2, count };

    // Add the row to the data vector if its not empty
    if (!row.empty())
    {
      ret.push_back(row);
    }
  }
  file.close();
  return ret;
}
void crossSectionPlot()
{
	// Import data from txt file (whitespace delimiter)
  std::string filename{ "CTPP-k-btheta.txt" };
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
  graph->SetTitle("CT proton -> proton q3=kaon -b diff theta");
  graph->SetMarkerColor(kBlue);
  graph->SetMarkerStyle(kOpenCircle);
	graph->SetLineColor(kBlue);
	//graph->GetXaxis()->SetRangeUser(0.5, 2.5);
	//graph->GetYaxis()->SetRangeUser(0, 2.5e-41);
 
  //  Plot the graph
  graph->DrawClone("APE");
 
  c1->Update();
  c1->Print("CTPP-k-bthetaPlot.png");
}

void dalitzPlot()
{
	// Import data from txt file (whitespace delimiter)
  std::string filename{ "CTNPm12m23-long.txt" };
  std::vector<std::vector<double>> data{ readDataFromTxt2(filename) };
	const size_t n{ data.size() };

  // Separate each column into individual vectors
	double mij1Vals[n];
	double mij2Vals[n];
	double countVals[n];
	for(size_t i{ 0 }; i < n; i++)
	{
		mij1Vals[i] = data[i][0];
		mij2Vals[i] = data[i][1];
		countVals[i] = data[i][2];
		//std::cout<<mij1Vals[i]<<" "<<mij2Vals[i]<<" "<<countVals[i]<<'\n';
	}

  // Create a canvas
  TCanvas *c1 = new TCanvas("c1","Our Plot", 2048, 1536);
 
  //  Instantiate a graph with errors
	int bins{ static_cast<int>(std::sqrt(n)) - 10 };
	TH2F* hDalitz = new TH2F("hDalitz", "Dalitz Plot;M_{12}^{2} (GeV^{2});M_{23}^{2} (GeV^{2})",
									bins,
									*(std::min_element(mij1Vals, mij1Vals + n)), *(std::max_element(mij1Vals, mij1Vals + n)), // x axis
									bins,
									*(std::min_element(mij2Vals, mij2Vals + n)), *(std::max_element(mij2Vals, mij2Vals + n)));// y axis
	TH2F* h2Dalitz = new TH2F("hDalitz", "Dalitz Plot;M_{12}^{2} (GeV^{2});M_{23}^{2} (GeV^{2})",
									bins,
									*(std::min_element(mij1Vals, mij1Vals + n)), *(std::max_element(mij1Vals, mij1Vals + n)), // x axis
									bins,
									*(std::min_element(mij2Vals, mij2Vals + n)), *(std::max_element(mij2Vals, mij2Vals + n)));// y axis


	//Fill the histogram with your data
  for(size_t i{ 0 }; i < n; i++)
	{
  	hDalitz->Fill(mij1Vals[i], mij2Vals[i], countVals[i]);
  	h2Dalitz->Fill(mij1Vals[i], mij2Vals[i], countVals[i]);
	}

	// set range
	hDalitz->GetXaxis()->SetRangeUser(*(std::min_element(mij1Vals, mij1Vals + n))-0.2,
																		*(std::max_element(mij1Vals, mij1Vals + n))+0.2);
	hDalitz->GetYaxis()->SetRangeUser(*(std::min_element(mij2Vals, mij2Vals + n))-0.1,
																		*(std::max_element(mij2Vals, mij2Vals + n))+0.1);
	h2Dalitz->GetXaxis()->SetRangeUser(*(std::min_element(mij1Vals, mij1Vals + n))-0.2,
																		*(std::max_element(mij1Vals, mij1Vals + n))+0.2);
	h2Dalitz->GetYaxis()->SetRangeUser(*(std::min_element(mij2Vals, mij2Vals + n))-0.1,
																		*(std::max_element(mij2Vals, mij2Vals + n))+0.1);


  // Plot the graph
  hDalitz->Draw("SURF2Z");
	h2Dalitz->Draw("SURF2Z same");
	gStyle->SetOptStat(0); // get rid of stats box
 
  c1->Update();
  c1->Print("CTNPm12m23Plot-long.png");
}

void plot()
{
	//dalitzPlot(); 
	crossSectionPlot();
}
