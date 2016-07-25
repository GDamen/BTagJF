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

using namespace std;

void EffVsLxyJF_b() {

	/*
	#ifndef CANVAS
	TCanvas *canvas = new TCanvas("JFEffCanvas","Graph2D example",0,0,800,800);
	#define CANVAS
	#endif
	*/
  
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gStyle->SetPalette(1);


  Float_t minPt = 350.;
  Float_t maxPt = 1000.;
  Float_t minTrkPt = 0.9;
  Float_t maxEta = 2.5;
  Float_t bEff = 0.70;
  
  /////////////////////////////
  // HISTOs
  TH1F * JF_llr = new TH1F("llr", "Log Likelihood Ratio for b with JF", 1000, -20, 50);
  TH1F * JF_llr_beampipe_03 = new TH1F("", "llr for b with JF in [0,3[ mm range", 1000, -20, 50);
  TH1F * JF_llr_beampipe_06 = new TH1F("", "llr for b with JF in [0,3[ mm range", 1000, -20, 50);
  TH1F * JF_llr_beampipe_10 = new TH1F("", "llr for b with JF in [3,10[ mm range", 1000, -20, 50);
  TH1F * JF_llr_beampipe_18 = new TH1F("", "llr for b with JF in [0,3[ mm range", 1000, -20, 50);
  TH1F * JF_llr_beampipe_25 = new TH1F("", "llr for b with JF in [10,25[ cm range", 1000, -20, 50);
  TH1F * JF_llr_IBL = new TH1F("", "llr for b with JF in [25,35[ cm range (around IBL)", 1000, -20, 50);
  TH1F * JF_llr_B = new TH1F("", "llr for b with JF after 35 mm (IBL to B layer)", 1000, -20, 50);
  TH1F * JF_llr_afterB = new TH1F("", "llr for b with JF after 35 mm (IBL to B layer)", 1000, -20, 50);

  TGraph * JF_efficiency_Lxy = new TGraph(8);
  	

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
    
    if(jetFlav->at(ij) == 5 && jetEta->at(ij) < maxEta){	//bottom flavour jet
    
    	JF_llr->Fill(jetJFProb->at(ij));
    	
    	if(secVertB->at(ij) < 3) {
    	  JF_llr_beampipe_03->Fill(jetJFProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 6) {
    	  JF_llr_beampipe_06->Fill(jetJFProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 10) {
    	  JF_llr_beampipe_10->Fill(jetJFProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 18) {
    	  JF_llr_beampipe_18->Fill(jetJFProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 25) {
    	  JF_llr_beampipe_25->Fill(jetJFProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 35) {
    	  JF_llr_IBL->Fill(jetJFProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 50){
    	  JF_llr_B->Fill(jetJFProb->at(ij));
    	}
    	else {
    	  JF_llr_afterB->Fill(jetJFProb->at(ij));
    	}
      }
    }
  }
  //////////////////////////////////////
  // Working Point evaluation
  
  float working_point_JF = IntegrateEff(JF_llr, bEff);
  TLine * line = new TLine(working_point_JF, 0, working_point_JF, 250);
  //JF_llr->Draw();
  float eff_03 = EvaluateEff(JF_llr_beampipe_03, bEff, working_point_JF);
  float eff_06 = EvaluateEff(JF_llr_beampipe_06, bEff, working_point_JF);
  float eff_10 = EvaluateEff(JF_llr_beampipe_10, bEff, working_point_JF);
  float eff_18 = EvaluateEff(JF_llr_beampipe_18, bEff, working_point_JF);
  float eff_25 = EvaluateEff(JF_llr_beampipe_25, bEff, working_point_JF);
  float eff_IBL = EvaluateEff(JF_llr_IBL, bEff, working_point_JF);
  float eff_B = EvaluateEff(JF_llr_B, bEff, working_point_JF);
  float eff_afterB = EvaluateEff(JF_llr_afterB, bEff, working_point_JF);
  
  JF_efficiency_Lxy -> SetPoint(0, 3, eff_03);
  JF_efficiency_Lxy -> SetPoint(1, 6, eff_06);
  JF_efficiency_Lxy -> SetPoint(2, 10, eff_10);
  JF_efficiency_Lxy -> SetPoint(3, 18, eff_18);
  JF_efficiency_Lxy -> SetPoint(4, 25, eff_25);
  JF_efficiency_Lxy -> SetPoint(5, 35, eff_IBL);
  JF_efficiency_Lxy -> SetPoint(6, 50, eff_B);
  JF_efficiency_Lxy -> SetPoint(7, 60, eff_afterB);
  
  JF_efficiency_Lxy->SetLineColor(2);
  	JF_efficiency_Lxy->SetLineWidth(4);
  	JF_efficiency_Lxy->SetMarkerColor(4);
  	JF_efficiency_Lxy->SetMarkerSize(1.5);
  	JF_efficiency_Lxy->SetMarkerStyle(21);
  	JF_efficiency_Lxy->SetTitle("Efficiency degradation for b hadrons reconstruction");
  	JF_efficiency_Lxy->GetXaxis()->SetTitle("Distance from IP [mm]");
  	JF_efficiency_Lxy->GetYaxis()->SetTitle("eff(Lxy)/eff(TOT)");
  JF_efficiency_Lxy -> Draw("same");
  
  TText * lineName = new TText(60, eff_afterB, "JF");
     lineName->SetTextColor(2);
  lineName -> Draw("same");
}


