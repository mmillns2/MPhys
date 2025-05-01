R__LOAD_LIBRARY(event1.so);


#include <limits>
#include <cmath>


#include"event1.h"


void wModelsScript() 
{
// Input files, add as many as necessary
  std::vector<std::string> inputFiles = {
    "/gluster/data/theory/mmillns/myeventfile_3000.root",
    "/gluster/data/theory/mmillns/myeventfile_2500.root",
    "/gluster/data/theory/mmillns/myeventfile_2000.root",
    "/gluster/data/theory/mmillns/myeventfile_1500.root",
  };

// Outer vector takes beam energy, Inner vector is simply: variable "x" against weight differential Xsec/"x"
  // { 3000, 2500, 2000, 1500 }
  std::vector<std::vector<std::pair<double,double>>> q2_vals(4);
  std::vector<std::vector<std::pair<double,double>>> theta_vals(4);
  std::vector<std::vector<std::pair<double,double>>> W_vals(4);
  std::vector<std::vector<std::pair<double,double>>> m12_vals(4);
  std::vector<std::vector<std::pair<double,double>>> m23_vals(4);

// Set boundaries for each variable
  double minQ2 = std::numeric_limits<double>::max();
  double maxQ2 = std::numeric_limits<double>::min();
	double minTheta = std::numeric_limits<double>::max();
	double maxTheta = std::numeric_limits<double>::min();
  double minW = std::numeric_limits<double>::max();
  double maxW = std::numeric_limits<double>::min();
  double minM12 = std::numeric_limits<double>::max();
  double maxM12 = std::numeric_limits<double>::min();
  double minM23 = std::numeric_limits<double>::max();
  double maxM23 = std::numeric_limits<double>::min();

// Initialise variables, need different bins for each plotting variable 'x'
	int bins{ };
	int N{ };

// Loops through each file
  for(size_t i{ 0 }; i < inputFiles.size(); ++i) {
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
      particle N0 = e->in.at(1);

      vect neutrino4Mom = e->in.at(0);
      vect lepton4Mom= e->out.at(0);
      vect kaon4Mom= e->out.at(1);
      vect N14Mom= e->out.at(2);
      vect q = neutrino4Mom - lepton4Mom;

      double q2 = -q*q;
      double theta = acos(e->costheta()) * 180/M_PI;
      double W = e->W();
      double m12 = (N14Mom + kaon4Mom) * (N14Mom + kaon4Mom);
      double m23 = (kaon4Mom + lepton4Mom) * (kaon4Mom + lepton4Mom);

      minQ2 = std::min(minQ2, q2);
      maxQ2 = std::max(maxQ2, q2);
      minTheta = std::min(minTheta, theta);
      maxTheta = std::max(maxTheta, theta);
      minW = std::min(minW, W);
      maxW = std::max(maxW, W);
      minM12 = std::min(minM12, m12);
      maxM12 = std::max(maxM12, m12);
      minM23 = std::min(minM23, m23);
      maxM23 = std::max(maxM23, m23);
    }

    // Populates the result vectors defined above, can be merged with the above potentially
    for(Long64_t eventIndex = 0; eventIndex < N; eventIndex++)
    {
      t->GetEntry(eventIndex);
      particle N0 = e->in.at(1);

      vect neutrino4Mom = e->in.at(0);
      vect lepton4Mom= e->out.at(0);
      vect kaon4Mom= e->out.at(1);
      vect N14Mom= e->out.at(2);
      vect q = neutrino4Mom - lepton4Mom;

      double q2 = -q*q;
      double theta = acos(e->costheta()) * 180/M_PI;
      double W = e->W();
      double m12 = (N14Mom + kaon4Mom) * (N14Mom + kaon4Mom);
      double m23 = (kaon4Mom + lepton4Mom) * (kaon4Mom + lepton4Mom);
      double weight = e->weight / N;
 
      q2_vals[i].push_back(std::make_pair(q2, weight / (maxQ2 - minQ2))); 
      W_vals[i].push_back(std::make_pair(W, weight / (maxW - minW))); 
      theta_vals[i].push_back(std::make_pair(theta, weight / (maxTheta - minTheta))); 
      m12_vals[i].push_back(std::make_pair(m12, weight)); 
      m23_vals[i].push_back(std::make_pair(m23, weight)); 
    }
    
    // File cleanup
    delete e;
    f->Close();
    delete f;
  }

  // Canvas for each plot
  TCanvas* cQ2 = new TCanvas("cQ2", "Q^{2} Plot", 2048, 1536);
  TCanvas* cTheta = new TCanvas("cTheta", "Angular Plot", 2048, 1536);
  TCanvas* cWtheta = new TCanvas("cQ2theta", "Energy-Angular Plot", 2048, 1536);
  TCanvas* cW = new TCanvas("cEnergy", "Energy Plot", 2048, 1536);
  TCanvas* cDalitz = new TCanvas("cDalitz", "Dalitz Plot", 2048, 1536);
  
  // Legend for each plot
  auto legendQ2 = new TLegend(0.75, 0.75, 0.95, 0.95);
  auto legendTheta = new TLegend(0.75, 0.75, 0.95, 0.95);
  auto legendWtheta = new TLegend(0.75, 0.75, 0.95, 0.95);
  auto legendW = new TLegend(0.75, 0.75, 0.95, 0.95);
	
  // Generate histograms at 3000
  TH1D* hQ2_3000 = new TH1D("hQ2 3000", "Q^{2} Distribution;Q^{2} [MeV^{2}];d#sigma/dQ^{2} [cm^{2}MeV^{-2}]", bins, minQ2, maxQ2);
  TH1D* hTheta_3000 =  new TH1D("hTheta 3000", "Angular Distribution;Theta [degrees];d#sigma/#theta [cm^{2}]", bins, minTheta, maxTheta);
  TH1D* hW_3000 = new TH1D("hW 3000", "Hadronic Invariant mass Distribution;W [MeV^{2}];d#sigma/dW [cm^{2}MeV^{-2}]", bins, minW, maxW);
  TH2D* hWtheta_3000 = new TH2D("hWTheta 3000", "Double Differential (W, #theta);W [MeV^{2}];#theta [degrees];d#sigma/dWd#theta [cm^{2}MeV^{-2}]", bins, minW, maxW, bins, minTheta, maxTheta);

  // Generate histograms at 2500
  TH1D* hQ2_2500 = new TH1D("hQ2 2500", "", bins, minQ2, maxQ2);
  TH1D* hTheta_2500 =  new TH1D("hTheta 2500", "", bins, minTheta, maxTheta);
  TH1D* hW_2500 = new TH1D("hW 2500", "", bins, minW, maxW);
  //TH2D* hWtheta_2500 = new TH2D("hWTheta 2500", "", bins, minW, maxW, bins, minTheta, maxTheta);

  // Generate histograms at 2000
  TH1D* hQ2_2000 = new TH1D("hQ2 2000", "", bins, minQ2, maxQ2);
  TH1D* hTheta_2000 =  new TH1D("hTheta 2000", "", bins, minTheta, maxTheta);
  TH1D* hW_2000 = new TH1D("hW 2000", "", bins, minW, maxW);
  //TH2D* hWtheta_2000 = new TH2D("hWTheta 2000", "", bins, minW, maxW, bins, minTheta, maxTheta);

  // Generate histograms at 1500
  TH1D* hQ2_1500 = new TH1D("hQ2 1500", "", bins, minQ2, maxQ2);
  TH1D* hTheta_1500 =  new TH1D("hTheta 1500", "", bins, minTheta, maxTheta);
  TH1D* hW_1500 = new TH1D("hW 1500", "", bins, minW, maxW);
  //TH2D* hWtheta_1500 = new TH2D("hWTheta 1500", "", bins, minW, maxW, bins, minTheta, maxTheta);

  // Generate Dalitz
  TH2D* hDalitz = new TH2D("hDalitz", "Dalitz Plot;W^{2} [MeV^{2}];m_{Kl}^{2} [MeV^{2}];Events", bins, minM12, maxM12, bins, minM23, maxM23);
  
  // Fill histograms with event data .first is x .second is y
  for(Long64_t eventIndex = 0; eventIndex < N; eventIndex++)
  {
    hQ2_3000->Fill(q2_vals[0][eventIndex].first, q2_vals[0][eventIndex].second);
    hW_3000->Fill(W_vals[0][eventIndex].first, W_vals[0][eventIndex].second);
    hTheta_3000->Fill(theta_vals[0][eventIndex].first, theta_vals[0][eventIndex].second);
    
    hQ2_2500->Fill(q2_vals[1][eventIndex].first, q2_vals[1][eventIndex].second);
    hW_2500->Fill(W_vals[1][eventIndex].first, W_vals[1][eventIndex].second);
    hTheta_2500->Fill(theta_vals[1][eventIndex].first, theta_vals[1][eventIndex].second);

    hQ2_2000->Fill(q2_vals[2][eventIndex].first, q2_vals[2][eventIndex].second);
    hW_2000->Fill(W_vals[2][eventIndex].first, W_vals[2][eventIndex].second);
    hTheta_2000->Fill(theta_vals[2][eventIndex].first, theta_vals[2][eventIndex].second);

    hQ2_1500->Fill(q2_vals[3][eventIndex].first, q2_vals[3][eventIndex].second);
    hW_1500->Fill(W_vals[3][eventIndex].first, W_vals[3][eventIndex].second);
    hTheta_1500->Fill(theta_vals[3][eventIndex].first, theta_vals[3][eventIndex].second);

    // 2D plots only filled at 3000 MeV
    assert(fabs(theta_vals[0][eventIndex].second - W_vals[0][eventIndex].second) < 1e-5);
    hWtheta_3000->Fill(W_vals[0][eventIndex].first, theta_vals[0][eventIndex].first, theta_vals[0][eventIndex].second);

    assert(fabs(m12_vals[0][eventIndex].second - m23_vals[0][eventIndex].second) < 1e-5);
    hDalitz->Fill(m12_vals[0][eventIndex].first, m23_vals[0][eventIndex].first, m12_vals[0][eventIndex].second);
  }

  // Global pointer needs to point to the canvas we're adding to
  cQ2->cd();
  // Removes the deafult ROOT stats
  gStyle->SetOptStat(0);

  // Aesthetics
  hQ2_3000->SetLineColor(kRed);
  hQ2_3000->GetXaxis()->SetRangeUser(0.0, 2500e3);
  hQ2_3000->GetYaxis()->SetRangeUser(0.0, 0.2e-48);
  hQ2_3000->Draw();

  hQ2_2500->SetLineColor(kBlue);
  hQ2_2500->GetXaxis()->SetRangeUser(0.0, 2500e3);
  hQ2_2500->GetYaxis()->SetRangeUser(0.0, 0.2e-48);
  hQ2_2500->Draw("SAME");
  
  hQ2_2000->SetLineColor(kGreen);
  hQ2_2000->GetXaxis()->SetRangeUser(0.0, 2500e3);
  hQ2_2000->GetYaxis()->SetRangeUser(0.0, 0.2e-48);
  hQ2_2000->Draw("SAME");

  hQ2_1500->SetLineColor(kOrange);
  hQ2_1500->GetXaxis()->SetRangeUser(0.0, 2500e3);
  hQ2_1500->GetYaxis()->SetRangeUser(0.0, 0.2e-48);
  hQ2_1500->Draw("SAME");

  legendQ2->AddEntry(hQ2_3000, "3000 MeV", "l");
  legendQ2->AddEntry(hQ2_2500, "2500 MeV", "l");
  legendQ2->AddEntry(hQ2_2000, "2000 MeV", "l");
  legendQ2->AddEntry(hQ2_1500, "1500 MeV", "l");
  legendQ2->Draw();

  cQ2->Update();
  cQ2->Print("q2-nuwro-test.png");
  

  cW->cd();
  gStyle->SetOptStat(0);

  hW_3000->SetLineColor(kRed);
  hW_3000->Draw();

  hW_2500->SetLineColor(kBlue);
  hW_2500->Draw("SAME");
  
  hW_2000->SetLineColor(kGreen);
  hW_2000->Draw("SAME");

  hW_1500->SetLineColor(kOrange);
  hW_1500->Draw("SAME");

  legendW->AddEntry(hW_3000, "3000 MeV", "l");
  legendW->AddEntry(hW_2500, "2500 MeV", "l");
  legendW->AddEntry(hW_2000, "2000 MeV", "l");
  legendW->AddEntry(hW_1500, "1500 MeV", "l");
  legendW->Draw();

  cW->Update();
  cW->Print("W-nuwro-test.png");

  cTheta->cd();
  gStyle->SetOptStat(0);

  hTheta_3000->SetLineColor(kRed);
  hTheta_3000->Draw();

  hTheta_2500->SetLineColor(kBlue);
  hTheta_2500->Draw("SAME");
  
  hTheta_2000->SetLineColor(kGreen);
  hTheta_2000->Draw("SAME");

  hTheta_1500->SetLineColor(kOrange);
  hTheta_1500->Draw("SAME");

  legendTheta->AddEntry(hTheta_3000, "3000 MeV", "l");
  legendTheta->AddEntry(hTheta_2500, "2500 MeV", "l");
  legendTheta->AddEntry(hTheta_2000, "2000 MeV", "l");
  legendTheta->AddEntry(hTheta_1500, "1500 MeV", "l");
  legendTheta->Draw();

  cTheta->Update();
  cTheta->Print("theta-nuwro-test.png");

  cWtheta->cd();
  hWtheta_3000->Draw("CONT1");
  cWtheta->Update();
  cWtheta->Print("W-theta-nuwro-test.png");
   
  cDalitz->cd();
  hDalitz->Draw("COLZ");
  cDalitz->Update();
  cDalitz->Print("dalitz-nuwro-test.png");
    
  delete hQ2_3000;
  delete hTheta_3000;
  delete hW_3000;
  delete hQ2_2500;
  delete hTheta_2500;
  delete hW_2500;
  //delete hWtheta_2500;
  delete hDalitz;
  delete cQ2;
  delete cTheta;
  delete cWtheta;
  delete cDalitz;
  delete cW;

}
