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

void EffVsLxyMV2_b() {
  
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gStyle->SetPalette(1);


  Float_t minPt = 350.;
  Float_t maxPt = 1000.;
  Float_t minTrkPt = 0.9;
  Float_t maxEta = 2.5;
  Float_t bEff = 0.70;
  Double_t MV2_llr_b = 0;
  
  /////////////////////////////
  // HISTOs
  TH1F * MV2_llr = new TH1F("llr", "Log Likelihood Ratio for b with MV2", 1000, -20, 50);
  TH1F * MV2_llr_beampipe_03 = new TH1F("", "llr for b with MV2 in [0,3[ mm range", 1000, -20, 50);
  TH1F * MV2_llr_beampipe_06 = new TH1F("", "llr for b with MV2 in [0,3[ mm range", 1000, -20, 50);
  TH1F * MV2_llr_beampipe_10 = new TH1F("", "llr for b with MV2 in [3,10[ mm range", 1000, -20, 50);
  TH1F * MV2_llr_beampipe_18 = new TH1F("", "llr for b with MV2 in [0,3[ mm range", 1000, -20, 50);
  TH1F * MV2_llr_beampipe_25 = new TH1F("", "llr for b with MV2 in [10,25[ cm range", 1000, -20, 50);
  TH1F * MV2_llr_IBL = new TH1F("", "llr for b with MV2 in [25,35[ cm range (around IBL)", 1000, -20, 50);
  TH1F * MV2_llr_B = new TH1F("", "llr for b with MV2 after 35 mm (IBL to B layer)", 1000, -20, 50);
  TH1F * MV2_llr_afterB = new TH1F("", "llr for b with MV2 after 35 mm (IBL to B layer)", 1000, -20, 50);

  TGraph * MV2_efficiency_Lxy = new TGraph(8);

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
  TTreeReaderValue<vector<double> > jetMV2ProbQ(myStdReader, "jet_mv2m_pu");
  TTreeReaderValue<vector<double> > jetMV2ProbB(myStdReader, "jet_mv2m_pb");
  TTreeReaderValue<vector<double> > jetMV2ProbC(myStdReader, "jet_mv2m_pc");

  ///////////////////////////////
  // FILL HISTOs
  // Loop over TTree entries
  while (myStdReader.Next()) {

    //loop over jets
    for(int ij = 0; ij < jetPt->size(); ij++){
    
    if(jetFlav->at(ij) == 5 && jetEta->at(ij) < maxEta){	//bottom flavour jet
    
        MV2_llr_b = (jetMV2ProbB->at(ij) / (jetMV2ProbQ->at(ij) + jetMV2ProbC->at(ij)));
    
    	MV2_llr->Fill(MV2_llr_b);
    	
    	if(secVertB->at(ij) < 3) {
    	  MV2_llr_beampipe_03->Fill(MV2_llr_b);
    	}
    	else if(secVertB->at(ij) < 6) {
    	  MV2_llr_beampipe_06->Fill(MV2_llr_b);
    	}
    	else if(secVertB->at(ij) < 10) {
    	  MV2_llr_beampipe_10->Fill(MV2_llr_b);
    	}
    	else if(secVertB->at(ij) < 18) {
    	  MV2_llr_beampipe_18->Fill(MV2_llr_b);
    	}
    	else if(secVertB->at(ij) < 25) {
    	  MV2_llr_beampipe_25->Fill(MV2_llr_b);
    	}
    	else if(secVertB->at(ij) < 35) {
    	  MV2_llr_IBL->Fill(MV2_llr_b);
    	}
    	else if(secVertB->at(ij) < 50){
    	  MV2_llr_B->Fill(MV2_llr_b);
    	}
    	else {
    	  MV2_llr_afterB->Fill(MV2_llr_b);
    	}
      }
    }
  }
  //////////////////////////////////////
  // Working Point evaluation
  
  float working_point_MV2 = IntegrateEff(MV2_llr, bEff);
  TLine * line = new TLine(working_point_MV2, 0, working_point_MV2, 250);
  //SV_llr->Draw();
  float eff_03 = EvaluateEff(MV2_llr_beampipe_03, bEff, working_point_MV2);
  float eff_06 = EvaluateEff(MV2_llr_beampipe_06, bEff, working_point_MV2);
  float eff_10 = EvaluateEff(MV2_llr_beampipe_10, bEff, working_point_MV2);
  float eff_18 = EvaluateEff(MV2_llr_beampipe_18, bEff, working_point_MV2);
  float eff_25 = EvaluateEff(MV2_llr_beampipe_25, bEff, working_point_MV2);
  float eff_IBL = EvaluateEff(MV2_llr_IBL, bEff, working_point_MV2);
  float eff_B = EvaluateEff(MV2_llr_B, bEff, working_point_MV2);
  float eff_afterB = EvaluateEff(MV2_llr_afterB, bEff, working_point_MV2);
  
  MV2_efficiency_Lxy -> SetPoint(0, 3, eff_03);
  MV2_efficiency_Lxy -> SetPoint(1, 6, eff_06);
  MV2_efficiency_Lxy -> SetPoint(2, 10, eff_10);
  MV2_efficiency_Lxy -> SetPoint(3, 18, eff_18);
  MV2_efficiency_Lxy -> SetPoint(4, 25, eff_25);
  MV2_efficiency_Lxy -> SetPoint(5, 35, eff_IBL);
  MV2_efficiency_Lxy -> SetPoint(6, 50, eff_B);
  MV2_efficiency_Lxy -> SetPoint(7, 60, eff_afterB);
  ////////////////////////////////////
  MV2_efficiency_Lxy->SetLineColor(4);
  MV2_efficiency_Lxy->SetLineWidth(4);
  MV2_efficiency_Lxy->SetMarkerColor(4);
  MV2_efficiency_Lxy->SetMarkerSize(1.5);
  MV2_efficiency_Lxy->SetMarkerStyle(21);
  MV2_efficiency_Lxy->SetTitle("Efficiency degradation for MV2 reconstruction");
  MV2_efficiency_Lxy->GetXaxis()->SetTitle("Distance from IP [mm]");
  MV2_efficiency_Lxy->GetYaxis()->SetTitle("eff(Lxy)/eff(TOT)");
  ///////////////////////////////////
  TLine * IBLline = new TLine(33.5, 0.15, 33.5, 1.3);
  
  
  // Draw calls
  MV2_efficiency_Lxy -> Draw("same lp");
	TText * lineName = new TText(60, eff_afterB, "MV2");
	lineName->SetTextColor(4);
  lineName -> Draw("SAME");

}


