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

void EffVsLxyIP3D_b() {

	/*
	#ifndef CANVAS
	TCanvas *canvas = new TCanvas("IP3DEffCanvas","Graph2D example",0,0,800,800);
	#define CANVAS
	#endif*/
  
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gStyle->SetPalette(1);


  Float_t minPt = 350.;
  Float_t maxPt = 1000.;
  Float_t minTrkPt = 0.9;
  Float_t maxEta = 2.5;
  Float_t bEff = 0.70;
  
  /////////////////////////////
  // HISTOs
  TH1F * IP3D_llr = new TH1F("llr", "Log Likelihood Ratio for b with IP3D", 1000, -20, 50);
  TH1F * IP3D_llr_beampipe_03 = new TH1F("", "llr for b with IP3D in [0,3[ mm range", 1000, -20, 50);
  TH1F * IP3D_llr_beampipe_06 = new TH1F("", "llr for b with IP3D in [0,3[ mm range", 1000, -20, 50);
  TH1F * IP3D_llr_beampipe_10 = new TH1F("", "llr for b with IP3D in [3,10[ mm range", 1000, -20, 50);
  TH1F * IP3D_llr_beampipe_18 = new TH1F("", "llr for b with IP3D in [0,3[ mm range", 1000, -20, 50);
  TH1F * IP3D_llr_beampipe_25 = new TH1F("", "llr for b with IP3D in [10,25[ cm range", 1000, -20, 50);
  TH1F * IP3D_llr_IBL = new TH1F("", "llr for b with IP3D in [25,35[ cm range (around IBL)", 1000, -20, 50);
  TH1F * IP3D_llr_B = new TH1F("", "llr for b with IP3D after 35 mm (IBL to B layer)", 1000, -20, 50);
  TH1F * IP3D_llr_afterB = new TH1F("", "llr for b with IP3D after 35 mm (IBL to B layer)", 1000, -20, 50);

  TGraph * IP3D_efficiency_Lxy = new TGraph(8);

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

  //IP3D values
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
    
    if(jetFlav->at(ij) == 5 && jetEta->at(ij) < maxEta){	//bottom flavour jet
    
    	IP3D_llr->Fill(jetIP3DProb->at(ij));
    	
    	if(secVertB->at(ij) < 3) {
    	  IP3D_llr_beampipe_03->Fill(jetIP3DProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 6) {
    	  IP3D_llr_beampipe_06->Fill(jetIP3DProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 10) {
    	  IP3D_llr_beampipe_10->Fill(jetIP3DProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 18) {
    	  IP3D_llr_beampipe_18->Fill(jetIP3DProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 25) {
    	  IP3D_llr_beampipe_25->Fill(jetIP3DProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 35) {
    	  IP3D_llr_IBL->Fill(jetIP3DProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 50){
    	  IP3D_llr_B->Fill(jetIP3DProb->at(ij));
    	}
    	else {
    	  IP3D_llr_afterB->Fill(jetIP3DProb->at(ij));
    	}
      }
    }
  }
  //////////////////////////////////////
  // Working Point evaluation
  
  float working_point_IP3D = IntegrateEff(IP3D_llr, bEff);
  //IP3D_llr->Draw();
  float eff_03 = EvaluateEff(IP3D_llr_beampipe_03, bEff, working_point_IP3D);
  float eff_06 = EvaluateEff(IP3D_llr_beampipe_06, bEff, working_point_IP3D);
  float eff_10 = EvaluateEff(IP3D_llr_beampipe_10, bEff, working_point_IP3D);
  float eff_18 = EvaluateEff(IP3D_llr_beampipe_18, bEff, working_point_IP3D);
  float eff_25 = EvaluateEff(IP3D_llr_beampipe_25, bEff, working_point_IP3D);
  float eff_IBL = EvaluateEff(IP3D_llr_IBL, bEff, working_point_IP3D);
  float eff_B = EvaluateEff(IP3D_llr_B, bEff, working_point_IP3D);
  float eff_afterB = EvaluateEff(IP3D_llr_afterB, bEff, working_point_IP3D);
  
  IP3D_efficiency_Lxy -> SetPoint(0, 3, eff_03);
  IP3D_efficiency_Lxy -> SetPoint(1, 6, eff_06);
  IP3D_efficiency_Lxy -> SetPoint(2, 10, eff_10);
  IP3D_efficiency_Lxy -> SetPoint(3, 18, eff_18);
  IP3D_efficiency_Lxy -> SetPoint(4, 25, eff_25);
  IP3D_efficiency_Lxy -> SetPoint(5, 35, eff_IBL);
  IP3D_efficiency_Lxy -> SetPoint(6, 50, eff_B);
  IP3D_efficiency_Lxy -> SetPoint(7, 60, eff_afterB);
  
  /////////////////////////////////////////
  // Efficiency plot draw
  IP3D_efficiency_Lxy->SetTitle("B Tagging efficiency degradation");
  IP3D_efficiency_Lxy->SetLineColor(1);
  IP3D_efficiency_Lxy->SetLineWidth(4);
  IP3D_efficiency_Lxy->SetMarkerColor(4);
  IP3D_efficiency_Lxy->SetMarkerSize(1.5);
  IP3D_efficiency_Lxy->SetMarkerStyle(21);
  IP3D_efficiency_Lxy->GetXaxis()->SetTitle("Distance from IP [mm]");
  IP3D_efficiency_Lxy->GetYaxis()->SetTitle("eff(Lxy)/eff(TOT)");
  IP3D_efficiency_Lxy -> Draw("alp");
  
  TLine * line = new TLine(33.5, 0.15, 33.5, 1.3);
  line -> Draw("SAME");
  
  TText * lineName = new TText(60, eff_afterB, "IP3D");
  lineName->SetTextColor(1);
  lineName -> Draw("SAME");
}




//CALCOLO ERRORI?
