R__LOAD_LIBRARY(event1.so);


#include <limits>
#include <cmath>


#include"event1.h"


void template() 
{
// Input files, add as many as necessary
  std::vector<std::string> inputFiles = {
    "/gluster/data/theory/ogregory/Argon_Mono_14_2000.root",
    "/gluster/data/theory/ogregory/Argon_Mono_anti14_2000.root",
    "/gluster/data/theory/ogregory/Argon_DUNE.root"
  };

  std::vector<std::vector<std::pair<double,double>>> m12_vals(4);
  std::vector<std::vector<std::pair<double,double>>> m23_vals(4);

// Set boundaries for each variable
  double minM12 = std::numeric_limits<double>::max();
  double maxM12 = std::numeric_limits<double>::min();
  double minM23 = std::numeric_limits<double>::max();
  double maxM23 = std::numeric_limits<double>::min();

// Initialise variables, need different bins for each plotting variable 'x'
	int bins{ };
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
    
    // Set bins here
    bins = std::max(10, static_cast<int>((N)));
    
    // First Loop through real events to determine min and max values
    for(Long64_t eventIndex = 0; eventIndex < N; eventIndex++)
    {
      t->GetEntry(eventIndex);

      vect lepton4Mom= e->out.at(0);
      vect kaon4Mom= e->out.at(2);
      vect N14Mom= e->out.at(1);

      double m12 = (N14Mom + kaon4Mom) * (N14Mom + kaon4Mom);
      double m23 = (kaon4Mom + lepton4Mom) * (kaon4Mom + lepton4Mom);

      minM12 = std::min(minM12, m12);
      maxM12 = std::max(maxM12, m12);
      minM23 = std::min(minM23, m23);
      maxM23 = std::max(maxM23, m23);
    }

    // Populates the result vectors defined above, can be merged with the above potentially
    for(Long64_t eventIndex = 0; eventIndex < N; eventIndex++)
    {
      t->GetEntry(eventIndex);

      vect lepton4Mom= e->out.at(0);
      vect kaon4Mom= e->out.at(2);
      vect N14Mom= e->out.at(1);

      double m12 = (N14Mom + kaon4Mom) * (N14Mom + kaon4Mom);
      double m23 = (kaon4Mom + lepton4Mom) * (kaon4Mom + lepton4Mom);
      double weight = e->weight / N;
 
      m12_vals[i].push_back(std::make_pair(m12, weight)); 
      m23_vals[i].push_back(std::make_pair(m23, weight)); 
    }
    
    // File cleanup
    delete e;
    f->Close();
    delete f;
  }

  // Canvas for each plot
  TCanvas* cDalitzNeu = new TCanvas("cDalitzNeu", "Dalitz Plot", 2048, 1536);
  TCanvas* cDalitzAnti = new TCanvas("cDalitzAnti", "Dalitz Plot", 2048, 1536);
  TCanvas* cDalitzDUNE = new TCanvas("cDalitzDUNE", "Dalitz Plot", 2048, 1536);

  // Generate Dalitz
  TH2D* hDalitzNeu = new TH2D("hDalitzNeu", "Dalitz Plot for Neutrinos;W^{2} [MeV^{2}];m_{Kl}^{2} [MeV^{2}];Events", bins, minM12, maxM12, bins, minM23, maxM23);
  TH2D* hDalitzAnti = new TH2D("hDalitzAnti", "Dalitz Plot for AntiNeutrinos;W^{2} [MeV^{2}];m_{Kl}^{2} [MeV^{2}];Events", bins, minM12, maxM12, bins, minM23, maxM23);
  TH2D* hDalitzDUNE = new TH2D("hDalitzDUNE", "Dalitz Plot for DUNE;W^{2} [MeV^{2}];m_{Kl}^{2} [MeV^{2}];Events", bins, minM12, maxM12, bins, minM23, maxM23);
  
  // Fill histograms with event data .first is x .second is y
  for(Long64_t eventIndex = 0; eventIndex < N; eventIndex++)
  {
    hDalitzNeu->Fill(m12_vals[0][eventIndex].first, m23_vals[0][eventIndex].first, m12_vals[0][eventIndex].second);
    hDalitzAnti->Fill(m12_vals[1][eventIndex].first, m23_vals[1][eventIndex].first, m12_vals[1][eventIndex].second);
    hDalitzNeu->Fill(m12_vals[2][eventIndex].first, m23_vals[2][eventIndex].first, m12_vals[2][eventIndex].second);
  }


   
  cDalitzNeu->cd();
  hDalitzNeu->Draw("COLZ");
  cDalitzNeu->Update();
  cDalitzNeu->Print("dalitz-neutrinos.png");

  cDalitzNeu->cd();
  hDalitzNeu->Draw("COLZ");
  cDalitzNeu->Update();
  cDalitzNeu->Print("dalitz-antineutrinos.png");

  cDalitzNeu->cd();
  hDalitzNeu->Draw("COLZ");
  cDalitzNeu->Update();
  cDalitzNeu->Print("dalitz-DUNE.png");
    
  delete hDalitzNeu;
  delete hDalitzAnti;
  delete hDalitzDUNE;
  delete cDalitzNeu;
  delete cDalitzAnti;
  delete cDalitzDUNE;
}
