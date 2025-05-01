R__LOAD_LIBRARY(event1.so);


#include <limits>
#include <cmath>


#include"event1.h"


void normalize(TH1D* hist) {
  int nBins = hist->GetNbinsX();
  for(int i = 1; i <= nBins; ++i) {
    double binContent = hist->GetBinContent(i);
    double binWidth = hist->GetBinWidth(i);
    if(binWidth > 0) {
      // normalize
      hist->SetBinContent(i, binContent / binWidth); 
      // adjust errors
      hist->SetBinError(i, hist->GetBinError(i) / binWidth);
    }
  }
}

void normalize(TH2D* hist) {
  int nBinsX = hist->GetNbinsX();
  int nBinsY = hist->GetNbinsY();
  for(int i = 1; i <= nBinsX; ++i) {
    for(int j = 1; i <= nBinsY; ++j) {
      double binContent = hist->GetBinContent(i, j);
      double binWidthX = hist->GetBinWidth(i);
      double binWidthY = hist->GetBinWidth(j);
      if(binWidthX > 0 && binWidthY > 0) {
        hist->SetBinContent(i, j, binContent / (binWidthX * binWidthY);
        hist->SetBinError(i, j, hist->GetBinError(i, j) / (binWidthX * binWidthY));
      }
    }
  }
}

void plots() {
  std::vector<std::string> inputFiles = {
    "/gluster/data/theory/mmillns/myeventfile.root",
  };
  for(std::string file : inputFiles) {
    // open file and get the tree
    TFile* f = TFile::Open(file.c_str());
    TTree* t = static_cast<TTree*>(f->Get("treeout"));

    event* e = new event();
    t->SetBranchAddress("e", &e);
    int N = t->GetEntries();

    TCanvas* cQ2 = new TCanvas("cQ2", "Q^{2} Plot", 2048, 1536);
    TCanvas* cTheta = new TCanvas("cTheta", "Angular Plot", 2048, 1536);
    TCanvas* cQ2theta = new TCanvas("cQ2theta", "Energy-Angular Plot", 2048, 1536);
    TCanvas* cEnergy = new TCanvas("cEnergy", "Energy Plot", 2048, 1536);
    TCanvas* cDalitz = new TCanvas("cDalitz", "Dalitz Plot", 2048, 1536);

    int bins{ std::max(10, static_cast<int>(std::sqrt(N) - 10)) };

    double minQ2 = std::numeric_limits<double>::max();
		double maxQ2 = std::numeric_limits<double>::min();
		double minTheta = std::numeric_limits<double>::max();
		double maxTheta = std::numeric_limits<double>::min();
    double minEbeam = std::numeric_limits<double>::max();
    double maxEbeam = std::numeric_limits<double>::min();
    double minM12 = std::numeric_limits<double>::max();
    double maxM12 = std::numeric_limits<double>::min();
    double minM23 = std::numeric_limits<double>::max();
    double maxM23 = std::numeric_limits<double>::min();

    for(Long64_t eventIndex = 0; eventIndex < N; eventIndex++) {
			t->GetEntry(eventIndex);
      particle N0 = e->in.at(1);
      vect neutrino4Mom = e->in.at(0);
      vect lepton4Mom= e->out.at(0);
      vect kaon4Mom= e->out.at(1);
      vect N14Mom= e->out.at(2);
      vect q = neutrino4Mom - lepton4Mom;
      double q2 = -q*q;
		  double theta = acos(e->costheta()) * 180/M_PI;
      double Ebeam = std::sqrt(e->s()) - N0.mass();
      double m12 = (N14Mom + kaon4Mom) * (N14Mom + kaon4Mom);
      double m23 = (kaon4Mom + lepton4Mom) * (kaon4Mom + lepton4Mom);

			minQ2 = std::min(minQ2, q2);
			maxQ2 = std::max(maxQ2, q2);
	    minTheta = std::min(minTheta, theta);
			maxTheta = std::max(maxTheta, theta);
      minEbeam = std::min(minEbeam, Ebeam);
      maxEbeam = std::max(maxEbeam, Ebeam);
      minM12 = std::min(minM12, m12);
      minM12 = std::max(maxM12, m12);
      minM23 = std::min(minM23, m23);
      minM23 = std::max(maxM23, m23);
    }

    TH1D* hQ2 = new TH1D("hQ2", "Q^{2} Distribution;Q^{2} [MeV^{2}];Events", bins, minQ2, maxQ2);
    TH1D* hTheta =  new TH1D("hTheta", "Theta Distribution;Theta [degrees];Events", bins, minTheta, maxTheta);
    TH1D* hEnergy = new TH1D("hEnergy", "Energy Distribution;E [MeV^{2}];Events", bins, minQ2, maxQ2);
    TH2D* hQ2theta = new TH2D("hQ2Theta", "Differential Q^{2} vs #theta;Q^{2} [MeV^{2}];#theta [degrees];dN/dQ^{2}d#theta", bins, minQ2, maxQ2, bins, minTheta, maxTheta);
    TH2D* hDalitz = new TH2D("hDalitz", "Dalitz Plot;m_{lK}^{2} [MeV^{2}];m_{K#pi}^{2} [MeV^{2}];Events", bins, minMlk2, maxMlk2, bins, minMkpi2, maxMkpi2);
		
    for(Long64_t eventIndex = 0; eventIndex < N; eventIndex++) {
			t->GetEntry(eventIndex);
      particle N0 = e->in.at(1);
      vect neutrino4Mom = e->in.at(0);
      vect lepton4Mom= e->out.at(0);
      vect kaon4Mom= e->out.at(1);
      vect N14Mom= e->out.at(2);
      vect q = neutrino4Mom - lepton4Mom;
      double q2 = -q*q;
			double theta = acos(e->costheta()) * 180/3.14;
      double Ebeam = std::sqrt(e->s()) - N0.mass();
      double m12 = (N14Mom + kaon4Mom) * (N14Mom + kaon4Mom);
      double m23 = (kaon4Mom + lepton4Mom) * (kaon4Mom + lepton4Mom);
      double weight = e->weight;

      hQ2->Fill(q2, weight);
      hEnergy->Fill(Ebeam, weight);
      hTheta->Fill(theta, weight);
      hQ2theta->Fill(q2, theta, weight);
      hDalitz->Fill(m12, m23, weight);
    }
    
    cQ2->cd();
    normalize(hQ2);
    hQ2->Draw();
    cQ2->Update();
    cQ2->Print("q2-nuwro-test.png");
    
    cEnergy->cd();
    normalize(hEnergy);
    hEnergy->Draw();
    cEnergy->Update();
    cEnergy->Print("energy-nuwro-test.png");

    cTheta->cd();
    normalize(hTheta);
    hTheta->Draw();
    cTheta->Update();
    cTheta->Print("theta-nuwro-test.png");

    cQ2theta->cd();
    normalize(hQ2theta);
    hQ2theta->Draw("COLZ");
    cQ2theta->Update();
    cQ2theta->Print("q2-theta-nuwro-test.png");
     
    cDalitz->cd();
    normalize(hDalitz);
    hDalitz->Draw("COLZ");
    cDalitz->Update();
    cDalitz->Print("dalitz-nuwro-test.png");
      
    delete hQ2;
    delete hTheta;
    delete hEnergy;
    delete hQ2theta;
    delete hDalitz;
    delete cQ2;
    delete cTheta;
    delete cQ2theta;
    delete cDalitz;
    delete cEnergy;
    delete e;
    f->Close();
    delete f;
  }
}
