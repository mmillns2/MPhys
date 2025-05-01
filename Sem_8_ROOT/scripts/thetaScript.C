R__LOAD_LIBRARY(event1.so);


#include <limits>
#include <cmath>


#include"event1.h"


void wModelsScript()
{
// Input files, add as many as necessary
  std::vector<std::string> inputFiles = {
    "/gluster/data/theory/ogregory/Argon_Mono_14_2000.root",
    "/gluster/data/theory/ogregory/Oxygen_Mono_14_2000.root",
    "/gluster/data/theory/ogregory/Hydrogen_Mono_14_2000.root",
    "/gluster/data/theory/ogregory/Argon_Mono_anti14_2000.root",
    "/gluster/data/theory/ogregory/Oxygen_Mono_anti14_2000.root",
    "/gluster/data/theory/ogregory/Hydrogen_Mono_anti14_2000.root"
    "/gluster/data/theory/ogregory/Argon_DUNE.root"
  };

// Outer vector takes beam energy, Inner vector is simply: variable "x" against weight differential Xsec/"x"
  // { 3000, 2500, 2000, 1500 }
  std::vector<std::vector<std::pair<double,double>>> W_vals(7);

// Set boundaries for each variable
  double minTheta = std::numeric_limits<double>::max();
  double maxTheta = std::numeric_limits<double>::min();

// Initialise variables, need different bins for each plotting variable 'x'
	int bins{ 100 };
	int N{ };

// Loops through each file
  for(size_t i{ 0 }; i < inputFiles.size(); ++i) 
  {
    // open file and get the tree
    TFile* f = TFile::Open(inputFiles[i].c_str());
    TTree* t = static_cast<TTree*>(f->Get("treeout"));

    event* e = new event();
    t->SetBranchAddress("e", &e);
    N = t->GetEntries();
    
    // First Loop through real events to determine min and max values
    for(Long64_t eventIndex = 0; eventIndex < N; eventIndex++)
    {
      t->GetEntry(eventIndex);

      double theta = e->costheta();

      minTheta[i] = std::min(minTheta[i], theta);
      maxTheta[i] = std::max(maxTheta[i], theta);
    }
    // Populates the result vectors defined above, can be merged with the above potentially
    for(Long64_t eventIndex = 0; eventIndex < N; eventIndex++)
    {
      t->GetEntry(eventIndex);

      double theta = e->costheta();
      double weight = e->weight / N;

      W_vals[i].push_back(std::make_pair(theta, bins*weight / (maxTheta[i] - minTheta[i]))); 
    }
    // File cleanup
    delete e;
    f->Close();
    delete f;
  }

  // Canvas for our plot
  TCanvas* cW_14 = new TCanvas("cW_Ar_14", "Energy Plot", 2048, 1536);
  TCanvas* cW_anti_14 = new TCanvas("cW_Ar_14anti", "Energy Plot", 2048, 1536);
  TCanvas* cW_DUNE = new TCanvas("cW_DUNE", "Energy Plot", 2048, 1536);
  // Legend for our plot
  auto legendW = new TLegend(0.75, 0.75, 0.95, 0.95);
	
  // Generate histograms
  TH1D* W_1 = new TH1D("W1", "Leptonic Angular Distribution for Mono-energetic Neutrinos Incident on Nuclei;cos(#theta);d#sigma/cos(#theta) [cm^{2}]", bins, minTheta[0], maxTheta[0]);
  TH1D* W_2 = new TH1D("W2", "Leptonic Angular Distribution for Mono-energetic Neutrinos Incident on Nuclei;cos(#theta);d#sigma/cos(#theta) [cm^{2}]", bins, minTheta[1], maxTheta[1]);
  TH1D* W_3 = new TH1D("W3", "Leptonic Angular Distribution for Mono-energetic Neutrinos Incident on Nuclei;cos(#theta);d#sigma/cos(#theta) [cm^{2}]", bins, minTheta[2], maxTheta[2]);
  TH1D* W_4 = new TH1D("W4", "Leptonic Angular Distribution for Mono-energetic Antineutrinos Incident on Nuclei;cos(#theta);d#sigma/cos(#theta) [cm^{2}]", bins, minTheta[3], maxTheta[3]);
  TH1D* W_5 = new TH1D("W5", "Leptonic Angular Distribution for Mono-energetic Antineutrinos Incident on Nuclei;cos(#theta);d#sigma/cos(#theta) [cm^{2}]", bins, minTheta[4], maxTheta[4]);
  TH1D* W_6 = new TH1D("W6", "Leptonic Angular Distribution for Mono-energetic Antineutrinos Incident on Nuclei;cos(#theta);d#sigma/cos(#theta) [cm^{2}]", bins, minTheta[5], maxTheta[5]);
  TH1D* W_7 = new TH1D("W7", "Leptonic Angular Distribution for the DUNE FHC_FD_Flux Beam Incident on Argon;cos(#theta);d#sigma/cos(#theta) [cm^{2}]", bins, minTheta[6], maxTheta[6]);
//  TH1D* hW_model_7 = new TH1D("W 7", "Hadronic Invariant mass Distribution;cos(#theta);d#sigma/cos(#theta) [cm^{2}]", bins, minTheta, maxTheta);

  
  // Fill histograms with event data .first is x .second is y
  for(Long64_t eventIndex = 0; eventIndex < N; eventIndex++)
  {
    W_1->Fill(W_vals[0][eventIndex].first, W_vals[0][eventIndex].second);
    W_3->Fill(W_vals[1][eventIndex].first, W_vals[1][eventIndex].second);
    W_2->Fill(W_vals[2][eventIndex].first, W_vals[2][eventIndex].second);
    W_4->Fill(W_vals[3][eventIndex].first, W_vals[3][eventIndex].second);
    W_5->Fill(W_vals[4][eventIndex].first, W_vals[4][eventIndex].second);
    W_6->Fill(W_vals[5][eventIndex].first, W_vals[5][eventIndex].second);
    W_7->Fill(W_vals[6][eventIndex].first, W_vals[6][eventIndex].second);
  }

  // Global pointer needs to point to the canvas we're adding to

  // mu neu (PDG: 14)
  cW_14->cd();
  // Legend for our plot
  auto legendNeu = new TLegend(0.75, 0.75, 0.95, 0.95);

  // Removes the deafult ROOT stats
  gStyle->SetOptStat(0);
  // Aesthetics
  W_1->SetLineColor(kMagenta);
  W_1->Draw("HIST E1");

  W_2->SetLineColor(kOrange);
  W_2->Draw("HIST E1");

  W_3->SetLineColor(kCyan);
  W_3->Draw("HIST E1");

  legendNeu->SetHeader("Mono-energetic {#nu_{#mu}} beam at {#sqrt{s}=2} GeV");
  legendNeu->AddEntry(W_1, "{Ar_18}^{40} Target ", "l");
  legendNeu->AddEntry(W_2, "{O_8}^{16} Target", "l");
  legendNeu->AddEntry(W_3, "{H_1}^1 Target", "l");

  // Output
  cW_14->Update();
  cW_14->Print("ThetaL_14.png");


  // anti mu neu (PDG: -14)
  cW_14anti->cd();
  // Legend for our plot
  auto legendAnti = new TLegend(0.75, 0.75, 0.95, 0.95);
  // Removes the deafult ROOT stats
  gStyle->SetOptStat(0);
  // Aesthetics
  W_4->SetLineColor(kMagenta);
  W_4->Draw("HIST E1");

  W_5->SetLineColor(kOrange);
  W_5->Draw("HIST E1");

  W_6->SetLineColor(kCyan);
  W_6->Draw("HIST E1");

  legendAnti->SetHeader("Mono-energetic {#nu_{#mu}} beam at {#sqrt{s}=2} GeV");
  legendAnti->AddEntry(W_3, "{Ar_18}^{40} Target ", "l");
  legendAnti->AddEntry(W_4, "{O_8}^{16} Target", "l");
  legendAnti->AddEntry(W_6, "{H_1}^1 Target", "l");

  // Output
  cW_14anti->Update();
  cW_14anti->Print("ThetaL_anti_14.png");



  // DUNE beam
  cW_DUNE->cd();
  // Legend for our plot
  auto legendDUNE = new TLegend(0.75, 0.75, 0.95, 0.95);
  // Removes the deafult ROOT stats
  gStyle->SetOptStat(0);
  // Aesthetics
  W_7->SetLineColor(kMagenta);
  W_7->Draw("HIST E1");


  legendAnti->SetHeader("DUNE FHC_FD_Flux Beam");
  legendAnti->AddEntry(W_7, "{Ar_18}^{40} Target ", "l");

  // Output
  cW_DUNE->Update();
  cW_DUNE->Print("ThetaL_plot_DUNE.png");



  delete W_1;
  delete W_2;
  delete W_3;
  delete W_4;
  delete W_5;
  delete W_6;
  delete W_7;
  delete cW_14;
  delete cW_anti_14;
  delete cW_DUNE;
}
