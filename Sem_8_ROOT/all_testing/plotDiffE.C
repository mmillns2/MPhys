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


int plotDiffE()
{
  //  Import data from txt file (whitespace delimiter)
  std::string filename3{ "dataFiles/CTNN-l-b3E3.txt" };
  std::vector<std::vector<double>> data3{ readDataFromTxt(filename3) };
	const size_t n3{ data3.size() };
  //  Import data from txt file (whitespace delimiter)
  std::string filename2{ "dataFiles/CTNP-l-b2E3.txt" };
  std::vector<std::vector<double>> data2{ readDataFromTxt(filename2) };
	const size_t n2{ data2.size() };
  //  Import data from txt file (whitespace delimiter)
  std::string filename1{ "dataFiles/CTNN-l-b1E3.txt" };
  std::vector<std::vector<double>> data1{ readDataFromTxt(filename1) };
	const size_t n1{ data1.size() };



  //  Separate each column into individual vectors
	double energyVals1[n1];
	double csVals1[n1];
	for(size_t i{ 0 }; i < n1; i++)
	{
		energyVals1[i] = data1[i][0];
		csVals1[i] = data1[i][1];
		std::cout<<energyVals1[i]<<" "<<csVals1[i]<<'\n';
	}

//  Separate each column into individual vectors
	double energyVals2[n2];
	double csVals2[n2];
	for(size_t i{ 0 }; i < n2; i++)
	{
		energyVals2[i] = data2[i][0];
		csVals2[i] = data2[i][1];
		std::cout<<energyVals2[i]<<" "<<csVals2[i]<<'\n';
	}

//  Separate each column into individual vectors
	double energyVals3[n3];
	double csVals3[n3];
	for(size_t i{ 0 }; i < n3; i++)
	{
		energyVals3[i] = data3[i][0];
		csVals3[i] = data3[i][1];
		std::cout<<energyVals3[i]<<" "<<csVals3[i]<<'\n';
	}

  //  Create a canvas
  TCanvas *c1 = new TCanvas("c1","Our Plot", 2048, 1536);
 
  c1->SetGrid();
  c1->GetFrame()->SetBorderSize(25);
 
  //  Instantiate graph1
  TGraphErrors *graph1 = new TGraphErrors(n1, energyVals1, csVals1, nullptr, nullptr);
  //  Instantiate graph2
  TGraphErrors *graph2 = new TGraphErrors(n2, energyVals2, csVals2, nullptr, nullptr);
  //  Instantiate graph3
  TGraphErrors *graph3 = new TGraphErrors(n3, energyVals3, csVals3, nullptr, nullptr);

  //  Change how the plot looks
  graph2->SetTitle("#nu_{#mu} + N #rightarrow P + #mu^{-} + K^{0}");
  gStyle->SetTitleFont(32, "t");

  //  Set axis titles
  graph2->GetXaxis()->SetTitle("lepton Scattering Energy E_{3} (GeV)");
  graph2->GetYaxis()->SetTitle("Differential Cross-Section #frac{d#sigma}{dE_{3}} (cm^{2}GeV^{-1})");

  graph2->GetXaxis()->SetTitleFont(42);
  graph2->GetYaxis()->SetTitleFont(42);

  //  Customize axis labels
  graph2->GetXaxis()->SetLabelFont(42);
  graph2->GetXaxis()->SetLabelSize(0.025); // Set label size
  graph2->GetYaxis()->SetLabelFont(42);
  graph2->GetYaxis()->SetLabelSize(0.025);  

//  Set the number of divisions (ticks) on the X and Y axes
  graph2->GetXaxis()->SetNdivisions(505);  // More ticks on X-axis
  graph2->GetYaxis()->SetNdivisions(505);  // Fewer ticks on Y-axis
//  Set tick font for X and Y axes
//  gStyle->SetTickFont(132); 
//  gStyle->SetTickFont(132);  // Bold Times Roman for axis ticks



  graph1->SetMarkerColor(kRed);
  graph1->SetMarkerStyle(20);
  graph1->SetLineColor(kBlack);

  graph2->SetMarkerColor(kBlue);
  graph2->SetMarkerStyle(22);

  graph3->SetMarkerColor(kBlue);
  graph3->SetMarkerStyle(22);

/*
//  Customize the X-axis to show ticks in terms of pi 
  TGaxis *axis = new TGaxis(graph->GetXaxis()->GetXmin(), 0, 
			    graph->GetXaxis()->GetXmax(), 0, 
			    graph->GetXaxis()->GetXmin()/TMath::Pi(), 
			    graph->GetXaxis()->GetXmax()/TMath::Pi(), 505, "+L"); // Show labels in terms of pi
*/
//  graph->GetXaxis()->SetTickLength(0.02);

//  graph2->GetXaxis()->SetRangeUser(0, 1);
//  graph->GetYaxis()->SetRangeUser(0, 5e-42);  //  PP scale
//  graph->GetYaxis()->SetRangeUser(0,120e-42);  // NP or  NN scale
 
  //  Draw the customised axis
//  axis->Draw();

  //  Plot the graph
  graph2->DrawClone("APE");
//  graph1->DrawClone("P SAME");
//  graph3->DrawClone("P SAME");

// Build and Draw a legend
  TLatex text;
  text.SetTextSize(0.04);
//  leg.SetFillColor(0);
//  graph2->SetFillColor(0);
//  leg.AddEntry(graph3,"s = 7.5GeV^{2}");
  text.SetTextAlign(13);
  text.SetNDC();
  text.DrawLatex(.15,.85,"s = 3.5GeV^{2} ");

 
  c1->Update();
  c1->Print("plotting/testPlotCTNP-l-bE3.png");
return 0;
}
