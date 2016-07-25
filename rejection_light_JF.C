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

void rejection_light_JF() {
  
  TCanvas *canvas = new TCanvas("Histocanvas","Graph2D example",0,0,800,800);
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gStyle->SetPalette(1);


  Float_t minPt = 350.;
  Float_t maxPt = 1000.;
  Float_t minTrkPt = 0.9;
  Float_t maxEta = 2.1;
  Float_t bEff = 0.70;
  Double_t jetJFProb_light;
  
  /////////////////////////////
  // HISTOs
  TH1F * JF_llr_b = new TH1F("llr", "Log Likelihood Ratio for b with JF", 1000, -20, 50);
  TH1F * JF_llr_l = new TH1F("llr", "Log Likelihood Ratio for light with JF", 1000, -20, 50);

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

  //JF values
  TTreeReaderValue<vector<float> > jetJFProbQ(myStdReader, "jet_jf_pu");
  TTreeReaderValue<vector<float> > jetJFProbB(myStdReader, "jet_jf_pb");
  TTreeReaderValue<vector<float> > jetJFProbC(myStdReader, "jet_jf_pc");
  TTreeReaderValue<vector<float> > jetJFProb(myStdReader, "jet_jf_llr");

  ///////////////////////////////
  // FILL HISTOs
  // Loop over TTree entries
  while (myStdReader.Next()) {
  
    //loop over jets
    for(int ij = 0; ij < jetPt->size(); ij++){
    
      jetJFProb_light = (jetJFProbQ->at(ij) / (jetJFProbB->at(ij) + jetJFProbC->at(ij)));
    
    	if(jetFlav->at(ij) == 0) {
      		JF_llr_l -> Fill(jetJFProb_light);
    	}
    	if(jetFlav->at(ij) == 5) {
      		JF_llr_b -> Fill(jetJFProb->at(ij));
    	}
    }
  }
  //////////////////////////////////////
  // Working Point evaluation
  
  float working_point_JF = IntegrateEff(JF_llr_b, bEff);
  TLine * line = new TLine(working_point_JF, 0, working_point_JF, 250);
  //JF_llr->Draw();
  float rej = EvaluateEff(JF_llr_l, bEff, working_point_JF);
  
  cout << rej << endl;
  
  JF_llr_l -> Draw();
}


