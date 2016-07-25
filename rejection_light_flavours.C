#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TMatrixD.h"
#include "TMath.h"
#include <iostream>

#include "functions.C"
using namespace std;

void rejection_light_flavours() {

  TCanvas *canvas = new TCanvas("Histocanvas","Graph2D example",0,0,800,800);
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gStyle->SetPalette(1);

  Float_t minPt = 350.;
  Float_t maxPt = 1000.;
  Float_t minTrkPt = 0.9;
  Float_t maxEta = 2.1;
  Float_t bEff = 0.70;
  Float_t jetIP3DProb_light;

  /////////////////////////////
  // HISTOs
  TH1F * IP3D_llr_b = new TH1F("llr", "Log Likelihood Ratio for b with IP3D", 1000, -20, 50);
    	IP3D_llr_b->SetLineColor(2);
  TH1F * IP3D_llr_l = new TH1F("llr", "Log Likelihood Ratio for l with IP3D", 1000, -20, 50);

  ////////////////////////////
  // TREE READER
  
  TFile *myStdFile = TFile::Open("flav_Akt4EMTo_tt_Std.root");
  TTreeReader myStdReader("bTag_AntiKt4EMTopoJets", myStdFile);
  
  TTreeReaderValue<vector<float> > jetPt(myStdReader, "jet_pt");
  TTreeReaderValue<vector<float> > jetEta(myStdReader, "jet_eta");
  TTreeReaderValue<vector<float> > jetPhi(myStdReader, "jet_phi");
  TTreeReaderValue<vector<int> >   jetFlav(myStdReader, "jet_truthflav");
  TTreeReaderValue<vector<float> > secVertB(myStdReader, "bH_Lxy");
  TTreeReaderValue<vector<float> > secVertC(myStdReader, "cH_Lxy");

  //SV values
  TTreeReaderValue<vector<float> > jetIP3DProbQ(myStdReader, "jet_ip3d_pu");
  TTreeReaderValue<vector<float> > jetIP3DProbB(myStdReader, "jet_ip3d_pb");
  TTreeReaderValue<vector<float> > jetIP3DProbC(myStdReader, "jet_ip3d_pc");
  TTreeReaderValue<vector<float> > jetIP3DProb(myStdReader, "jet_ip3d_llr");

  ///////////////////////////////
  // FILL HISTOs
  // Loop over TTree entries
  while (myStdReader.Next()) {
  
    //loop over jets
    for(int ij = 0; ij < jetPt->size(); ij++){
    
      jetIP3DProb_light = (jetIP3DProbQ->at(ij) / (jetIP3DProbB->at(ij) + jetIP3DProbC->at(ij)));
    
    	if(jetFlav->at(ij) == 0) {
      		IP3D_llr_l -> Fill(jetIP3DProb_light);
    	}
    	if(jetFlav->at(ij) == 5) {
      		IP3D_llr_b -> Fill(jetIP3DProb->at(ij));
    	}
    }
  }
  cout << IP3D_llr_b -> GetBinContent(500) << endl;
  cout << FindEff(IP3D_llr_b, IP3D_llr_l) << endl;
  IP3D_llr_b -> Draw();
  IP3D_llr_l -> Draw("SAME");
}
