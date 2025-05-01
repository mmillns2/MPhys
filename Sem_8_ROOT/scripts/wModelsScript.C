R__LOAD_LIBRARY(event1.so);


#include <limits>
#include <cmath>


#include"event1.h"


void wModelsScript()
{
// Input files, add as many as necessary
  std::vector<std::string> inputFiles = {
    "/gluster/data/theory/ogregory/w_model_1.root",
    "/gluster/data/theory/ogregory/w_model_2.root",
    "/gluster/data/theory/ogregory/w_model_3.root",
    "/gluster/data/theory/ogregory/w_model_4.root",
    "/gluster/data/theory/ogregory/w_model_5.root",
    "/gluster/data/theory/ogregory/w_model_6.root"
  };

// Set boundaries for each variable
std::vector<double> minW (7, std::numeric_limits<double>::max());
std::vector<double> maxW (7, std::numeric_limits<double>::min());

// Set boundaries for each variable
  double minW = std::numeric_limits<double>::max();
  double maxW = std::numeric_limits<double>::min();

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
      particle N0 = e->W();

      double W = std::sqrt((vect(N0)+vect(K))*(vect(N0)+vect(K)));

      minW[i] = std::min(minW[i], W);
      maxW[i] = std::max(maxW[i], W);
    }
    // Populates the result vectors defined above, can be merged with the above potentially
    for(Long64_t eventIndex = 0; eventIndex < N; eventIndex++)
    {
      t->GetEntry(eventIndex);

      double W = e->W();

      double weight = e->weight / N;
  
      W_vals[i].push_back(std::make_pair(W, weight / (maxW[i] - minW[i]))); 
    }
    // File cleanup
    delete e;
    f->Close();
    delete f;
  }

  // Canvas for our plot
  TCanvas* cW = new TCanvas("cEnergy", "Energy Plot", 2048, 1536);
  // Legend for our plot
  auto legendW = new TLegend(0.75, 0.75, 0.95, 0.95);
	
  // Generate histograms
  TH1D* hW_model_1 = new TH1D("hW model 1", "Hadronic Invariant mass Distribution;W [MeV^{2}];d#sigma/dW [cm^{2}MeV^{-2}]", bins, minW[0], maxW[0]);
  TH1D* hW_model_2 = new TH1D("hW model 2", "Hadronic Invariant mass Distribution;W [MeV^{2}];d#sigma/dW [cm^{2}MeV^{-2}]", bins, minW[1], maxW[1]);
  TH1D* hW_model_3 = new TH1D("hW model 3", "Hadronic Invariant mass Distribution;W [MeV^{2}];d#sigma/dW [cm^{2}MeV^{-2}]", bins, minW[2], maxW[2]);
  TH1D* hW_model_4 = new TH1D("hW model 4", "Hadronic Invariant mass Distribution;W [MeV^{2}];d#sigma/dW [cm^{2}MeV^{-2}]", bins, minW[3], maxW[3]);
  TH1D* hW_model_5 = new TH1D("hW model 5", "Hadronic Invariant mass Distribution;W [MeV^{2}];d#sigma/dW [cm^{2}MeV^{-2}]", bins, minW[4], maxW[4]);
  TH1D* hW_model_6 = new TH1D("hW model 6", "Hadronic Invariant mass Distribution;W [MeV^{2}];d#sigma/dW [cm^{2}MeV^{-2}]", bins, minW[5], maxW[5]);
//  TH1D* hW_model_7 = new TH1D("hW model 7", "Hadronic Invariant mass Distribution;W [MeV^{2}];d#sigma/dW [cm^{2}MeV^{-2}]", bins, minW, maxW);

  
  // Fill histograms with event data .first is x .second is y
  for(Long64_t eventIndex = 0; eventIndex < N; eventIndex++)
  {
    hW_model_1->Fill(W_vals[0][eventIndex].first, W_vals[0][eventIndex].second);
    hW_model_3->Fill(W_vals[1][eventIndex].first, W_vals[1][eventIndex].second);
    hW_model_2->Fill(W_vals[2][eventIndex].first, W_vals[2][eventIndex].second);
    hW_model_4->Fill(W_vals[3][eventIndex].first, W_vals[3][eventIndex].second);
    hW_model_5->Fill(W_vals[4][eventIndex].first, W_vals[4][eventIndex].second);
    hW_model_6->Fill(W_vals[5][eventIndex].first, W_vals[5][eventIndex].second);
  //  hW_model_7->Fill(W_vals[6][eventIndex].first, W_vals[6][eventIndex].second);
  }

  // Global pointer needs to point to the canvas we're adding to
  cW->cd();
  // Removes the deafult ROOT stats
  gStyle->SetOptStat(0);

  // Aesthetics
  hW_model_1->SetLineColor(kRed);
  hW_model_1->Draw("HIST E1");

  hW_model_2->SetLineColor(kBlue);
  hW_model_2->Draw("HIST E1 SAME");
  
  hW_model_3->SetLineColor(kGreen);
  hW_model_3->Draw("HIST E1 SAME");

  hW_model_4->SetLineColor(kCyan);
  hW_model_4->Draw("HIST E1 SAME");

  hW_model_5->SetLineColor(kMagenta);
  hW_model_5->Draw("HIST E1 SAME");

  hW_model_6->SetLineColor(kOrange);
  hW_model_6->Draw("HIST E1 SAME");
  
//  hW_model_7->SetLineColor(kBlack);
//  hW_model_7->Draw("SAME");

  legendW->SetHeader("Mono-energetic {#nu_{#mu}} beam at {#sqrt{s}=2} GeV on a {C_6}^{12} Target");
  legendW->AddEntry(hW_model_1, "Free Target", "l");
  legendW->AddEntry(hW_model_2, "Fermi Gas", "l");
  legendW->AddEntry(hW_model_3, "Local Fermi Gas", "l");
  legendW->AddEntry(hW_model_4, "Bodek-Ritchie Fermi Gas", "l");
  legendW->AddEntry(hW_model_5, "Effective Spectral Function", "l");
  legendW->AddEntry(hW_model_6, "Effective Nuclear Potential", "l");
//  legendW->AddEntry(hW_model_7, "Deuterium", "l");
  legendW->Draw();

  cW->Update();
  cW->Print("W-models-plot.png");

  delete hW_model_1;
  delete hW_model_2;
  delete hW_model_3;
  delete hW_model_4;
  delete hW_model_5;
  delete hW_model_6;
//  delete hW_model_7;
  delete cW;
}
