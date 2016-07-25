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

#ifndef FUNCTIONS
#include "functions.C"
#endif
/*
#ifndef CANVAS
TCanvas *canvas = new TCanvas("EffCanvas","Graph2D example",0,0,800,800);
#define CANVAS
#endif
*/
using namespace std;

void EffVsLxySV_b() {
  
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gStyle->SetPalette(1);


  Float_t minPt = 350.;
  Float_t maxPt = 1000.;
  Float_t minTrkPt = 0.9;
  Float_t maxEta = 2.5;
  Float_t bEff = 0.70;
  
  /////////////////////////////
  // HISTOs
  TH1F * SV_llr = new TH1F("llr", "Log Likelihood Ratio for b with SV", 1000, -20, 50);
  TH1F * SV_llr_beampipe_03 = new TH1F("", "llr for b with SV in [0,3[ mm range", 1000, -20, 50);
  TH1F * SV_llr_beampipe_06 = new TH1F("", "llr for b with SV in [0,3[ mm range", 1000, -20, 50);
  TH1F * SV_llr_beampipe_10 = new TH1F("", "llr for b with SV in [3,10[ mm range", 1000, -20, 50);
  TH1F * SV_llr_beampipe_18 = new TH1F("", "llr for b with SV in [0,3[ mm range", 1000, -20, 50);
  TH1F * SV_llr_beampipe_25 = new TH1F("", "llr for b with SV in [10,25[ cm range", 1000, -20, 50);
  TH1F * SV_llr_IBL = new TH1F("", "llr for b with SV in [25,35[ cm range (around IBL)", 1000, -20, 50);
  TH1F * SV_llr_B = new TH1F("", "llr for b with SV after 35 mm (IBL to B layer)", 1000, -20, 50);
  TH1F * SV_llr_afterB = new TH1F("", "llr for b with SV after 35 mm (IBL to B layer)", 1000, -20, 50);

  TGraph * SV_efficiency_Lxy = new TGraph(8);
  	SV_efficiency_Lxy->SetLineColor(3);
  	SV_efficiency_Lxy->SetLineWidth(4);
  	SV_efficiency_Lxy->SetMarkerColor(4);
  	SV_efficiency_Lxy->SetMarkerSize(1.5);
  	SV_efficiency_Lxy->SetMarkerStyle(21);
  	SV_efficiency_Lxy->SetTitle("Efficiency degradation for SV reconstruction");
  	SV_efficiency_Lxy->GetXaxis()->SetTitle("Distance from IP [mm]");
  	SV_efficiency_Lxy->GetYaxis()->SetTitle("eff(Lxy)/eff(TOT)");

  ////////////////////////////
  // TREE READER
  
  TFile *myStdFile = TFile::Open("ntuples/flav_Akt4EMTo_tt_Std.root");
  TTreeReader myStdReader("bTag_AntiKt4EMTopoJets", myStdFile);
  
  TTreeReaderValue<vector<float> > jetPt(myStdReader, "jet_pt");
  TTreeReaderValue<vector<float> > jetEta(myStdReader, "jet_eta");
  TTreeReaderValue<vector<float> > jetPhi(myStdReader, "jet_phi");
  TTreeReaderValue<vector<int> >   jetFlav(myStdReader, "jet_truthflav");
  TTreeReaderValue<vector<float> > secVertB(myStdReader, "bH_Lxy");
  TTreeReaderValue<vector<float> > secVertC(myStdReader, "cH_Lxy");

  //SV values
  TTreeReaderValue<vector<float> > jetSVProbQ(myStdReader, "jet_sv1_pu");
  TTreeReaderValue<vector<float> > jetSVProbB(myStdReader, "jet_sv1_pb");
  TTreeReaderValue<vector<float> > jetSVProbC(myStdReader, "jet_sv1_pc");
  TTreeReaderValue<vector<float> > jetSVProb(myStdReader, "jet_sv1_llr");

  ///////////////////////////////
  // FILL HISTOs
  // Loop over TTree entries
  while (myStdReader.Next()) {

    //loop over jets
    for(int ij = 0; ij < jetPt->size(); ij++){
    
    if(jetFlav->at(ij) == 5 && jetEta->at(ij) < maxEta){	//bottom flavour jet
    
    	SV_llr->Fill(jetSVProb->at(ij));
    	
    	if(secVertB->at(ij) < 3) {
    	  SV_llr_beampipe_03->Fill(jetSVProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 6) {
    	  SV_llr_beampipe_06->Fill(jetSVProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 10) {
    	  SV_llr_beampipe_10->Fill(jetSVProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 18) {
    	  SV_llr_beampipe_18->Fill(jetSVProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 25) {
    	  SV_llr_beampipe_25->Fill(jetSVProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 35) {
    	  SV_llr_IBL->Fill(jetSVProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 50){
    	  SV_llr_B->Fill(jetSVProb->at(ij));
    	}
    	else {
    	  SV_llr_afterB->Fill(jetSVProb->at(ij));
    	}
      }
    }
  }
  //////////////////////////////////////
  // Working Point evaluation
  
  float working_point_SV = IntegrateEff(SV_llr, bEff);
  TLine * line = new TLine(working_point_SV, 0, working_point_SV, 250);
  //SV_llr->Draw();
  float eff_03 = EvaluateEff(SV_llr_beampipe_03, bEff, working_point_SV);
  float eff_06 = EvaluateEff(SV_llr_beampipe_06, bEff, working_point_SV);
  float eff_10 = EvaluateEff(SV_llr_beampipe_10, bEff, working_point_SV);
  float eff_18 = EvaluateEff(SV_llr_beampipe_18, bEff, working_point_SV);
  float eff_25 = EvaluateEff(SV_llr_beampipe_25, bEff, working_point_SV);
  float eff_IBL = EvaluateEff(SV_llr_IBL, bEff, working_point_SV);
  float eff_B = EvaluateEff(SV_llr_B, bEff, working_point_SV);
  float eff_afterB = EvaluateEff(SV_llr_afterB, bEff, working_point_SV);
  
  SV_efficiency_Lxy -> SetPoint(0, 3, eff_03);
  SV_efficiency_Lxy -> SetPoint(1, 6, eff_06);
  SV_efficiency_Lxy -> SetPoint(2, 10, eff_10);
  SV_efficiency_Lxy -> SetPoint(3, 18, eff_18);
  SV_efficiency_Lxy -> SetPoint(4, 25, eff_25);
  SV_efficiency_Lxy -> SetPoint(5, 35, eff_IBL);
  SV_efficiency_Lxy -> SetPoint(6, 50, eff_B);
  SV_efficiency_Lxy -> SetPoint(7, 60, eff_afterB);
  SV_efficiency_Lxy -> Draw("same lp");
  
  TText * lineName = new TText(60, eff_afterB, "SV");
  lineName->SetTextColor(3);
  lineName -> Draw("SAME");
}


